
#include "controller.h"
#include <unistd.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <thread>
#include "absl/strings/str_cat.h"
#include "src/agd/errors.h"
#include "src/common/all_all_executor.h"
#include "src/common/cluster_set.h"
#include "src/dist/checkpoint.h"
#include "src/dist/all_all_dist.h"

using std::cout;
using std::string;
using std::thread;
using namespace std::chrono_literals;

long int timestamp() {
  time_t t = std::time(0);
  long int now = static_cast<long int>(t);
  return now;
}

// https://stackoverflow.com/questions/2209135/safely-prompt-for-yes-no-with-cin
bool PromptForChar(const string& prompt, char& readch) {
  std::string tmp;
  std::cout << prompt << std::endl;
  if (std::getline(std::cin, tmp)) {
    // Only accept single character input
    if (tmp.length() == 1) {
      readch = tmp[0];
    } else {
      // For most input, char zero is an appropriate sentinel
      readch = '\0';
    }
    return true;
  }
  return false;
}

agd::Status Controller::Run(const Params& params,
                            const Parameters& aligner_params,
                            std::vector<std::unique_ptr<Dataset>>& datasets) {
  checkpoint_timer_ = timestamp();
  std::atomic_int_fast32_t outstanding_requests{0};
  cout << "Num seqs threshold: " << params.nseqs_threshold << std::endl;
  // index all sequences
  agd::Status s = Status::OK();
  const char* data;
  size_t length;
  size_t id = 0;  // absolute ID
  for (auto& dataset : datasets) {
    size_t dataset_index = 0;
    s = dataset->GetNextRecord(&data, &length);
    while (s.ok()) {
      auto seq = Sequence(absl::string_view(data, length), dataset->Name(),
                          dataset->Size(), dataset_index, id);
      sequences_.push_back(std::move(seq));
      id++;
      dataset_index++;
      s = dataset->GetNextRecord(&data, &length);
    }
  }

  if (params.checkpoint_interval) {
    cout << "controller will checkpoint with interval of "
         << params.checkpoint_interval << std::endl;
  }

  auto total_merges = sequences_.size() - 1;
  outstanding_merges_ = total_merges;
  if (params.dataset_limit > 0) {
    outstanding_merges_ = params.dataset_limit - 1;
  }
  cout << "outstanding merges to complete: " << outstanding_merges_ << std::endl;
  cout << "dup removal thresh is " << params.dup_removal_thresh << std::endl;
  cout << "Using " << params.num_threads << " threads to merge partials\n";

  // create envs, params
  string logpam_json_file = absl::StrCat(params.data_dir_path, "logPAM1.json");
  string all_matrices_json_file =
      absl::StrCat(params.data_dir_path, "all_matrices.json");

  std::ifstream logpam_stream(logpam_json_file);
  std::ifstream allmat_stream(all_matrices_json_file);

  if (!logpam_stream.good()) {
    return agd::errors::Internal("File ", logpam_json_file, " not found.");
  }

  if (!allmat_stream.good()) {
    return agd::errors::Internal("File ", all_matrices_json_file,
                                 " not found.");
  }

  json logpam_json;
  logpam_stream >> logpam_json;

  json all_matrices_json;
  allmat_stream >> all_matrices_json;

  AlignmentEnvironments envs;

  // initializing envs is expensive, so don't copy this
  cout << "Worker initializing environments from " << params.data_dir_path
       << std::endl;
  envs.InitFromJSON(logpam_json, all_matrices_json, aligner_params.min_score);
  cout << "Done.\n";

  // done init envs

  // connect to zmq queues
  auto address = std::string("tcp://*:");
  auto response_queue_address =
      absl::StrCat(address, params.response_queue_port);
  auto request_queue_address = absl::StrCat(address, params.request_queue_port);
  auto incomplete_request_queue_address =
      absl::StrCat(address, params.incomplete_request_queue_port);
  auto set_request_channel_address =
      absl::StrCat(address, params.set_request_port);

  context_ = zmq::context_t(1);
  context_sink_ = zmq::context_t(2);
  context_ir_sink_ = zmq::context_t(1);
  context_set_request_ = zmq::context_t(1);
  try {
    zmq_recv_socket_.reset(new zmq::socket_t(context_sink_, ZMQ_PULL));
  } catch (...) {
    return agd::errors::Internal("Could not create zmq PULL socket ");
  }

  try {
    zmq_send_socket_.reset(new zmq::socket_t(context_, ZMQ_REP));
  } catch (...) {
    return agd::errors::Internal("Could not create zmq REP socket");
  }

  try {
    zmq_set_request_socket_.reset(
        new zmq::socket_t(context_set_request_, ZMQ_REP));
  } catch (...) {
    return agd::errors::Internal(
        "Could not create a zmp REP socket -- large_pm");
  }

  // zmq_recv_socket_->setsockopt(ZMQ_RCVHWM, 2);
  // TODO make sendhwm a param
  // zmq_send_socket_->setsockopt(ZMQ_SNDHWM, 3);
  // int val = zmq_send_socket_->getsockopt<int>(ZMQ_SNDHWM);
  // cout << "snd hwm value is " << val << " \n";

  // waits for timeout time and then returns with EAGAIN
  zmq_set_request_socket_->setsockopt(ZMQ_RCVTIMEO, 1000);
  int val = zmq_set_request_socket_->getsockopt<int>(ZMQ_RCVTIMEO);
  cout << "recv timeout set to " << val << " ms \n";

  try {
    zmq_incomplete_request_socket_.reset(
        new zmq::socket_t(context_ir_sink_, ZMQ_PULL));
  } catch (...) {
    return agd::errors::Internal("Could not create zmq INCOMP REQ socket ");
  }

  try {
    zmq_recv_socket_->bind(response_queue_address.c_str());
  } catch (...) {
    return agd::errors::Internal("Could not connect to zmq at ",
                                 response_queue_address);
  }

  try {
    zmq_send_socket_->bind(request_queue_address.c_str());
  } catch (...) {
    return agd::errors::Internal("Could not connect to zmq at ",
                                 request_queue_address);
  }

  try {
    zmq_set_request_socket_->bind(set_request_channel_address.c_str());
  } catch (...) {
    return agd::errors::Internal("Could not connect to zmq at ",
                                 set_request_channel_address);
  }

  try {
    zmq_incomplete_request_socket_->bind(
        incomplete_request_queue_address.c_str());
  } catch (...) {
    return agd::errors::Internal("Could not connect to zmq at ",
                                 incomplete_request_queue_address);
  }

  request_queue_.reset(
      new ConcurrentQueue<MarshalledRequest>(params.queue_depth));
  response_queue_.reset(
      new ConcurrentQueue<MarshalledResponse>(params.queue_depth));
  incomplete_request_queue_.reset(
      new ConcurrentQueue<MarshalledRequest>(params.queue_depth));
  sets_to_merge_queue_.reset(
      new ConcurrentQueue<MarshalledClusterSet>(sequences_.size()));

  int total_sent = 0;
  request_queue_thread_ = thread([this, &total_sent]() {
    // get msg from work queue (output)
    // encode
    // put in work queue zmq
    // repeat
    MarshalledRequest merge_request;

    while (run_) {
      if (!request_queue_->pop(merge_request)) {
        continue;
      }

      zmq::message_t message;
      zmq_send_socket_->recv(&message);

      auto size = merge_request.buf.size();
      zmq::message_t msg(merge_request.buf.release_raw(), size, free_func,
                         NULL);
      /*cout << "pushing request of size " << size << " of type "
           << (merge_request.has_batch() ? "batch " : "partial ") << std::endl;*/

      bool success = zmq_send_socket_->send(std::move(msg));
      if (!success) {
        cout << "Thread failed to send request over zmq!\n";
      }
      total_sent++;
    }

    cout << "Work queue thread ending. Total sent: " << total_sent << std::endl;
  });

  // std::thread queue_measure_thread = std::thread([this](){
  //     // get timestamp, queue size
  //     //cout << "queue measure thread starting ...\n";
  //     while(run_) {
  //       time_t result = std::time(nullptr);
  //       //timestamps_.push_back(static_cast<long int>(result));
  //       //queue_sizes_.push_back(work_queue_->size());
  //       std::cout << static_cast<long int>(result) << ": " <<
  //       response_queue_->size() << std::endl;
  //       std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  //       // if (queue_sizes_.size() >= 1000000) {
  //       //  break;  // dont run forever ...
  //       // }
  //     }
  //     cout << "queue measure thread finished\n";
  //   }
  // );

  set_request_thread_ = thread([this]() {
    // receive requests from worker and send them cluster sets

    while (run_) {
      zmq::message_t message;
      bool success = zmq_set_request_socket_->recv(&message);
      // returns 0 on timeout
      if (!success) {
        continue;
      }

      // extract the id
      int id = *reinterpret_cast<int*>(message.data());
      // cout << "Got set request with id [" << id << "]\n";

      // search in partial merge map
      PartialMergeItem* partial_item;
      {
        absl::MutexLock l(&mu_);
        auto partial_it = partial_merge_map_.find(id);
        if (partial_it == partial_merge_map_.end()) {
          cout << "the set with " << id << " is not in the map \n";
          exit(0);
        }
        partial_item = &partial_it->second;
        // cout << "Some data --> " << partial_item->buf.size() << std::endl;

        // creating a copy since zmq takes ownership without creating a copy
        agd::Buffer buf;
        buf.AppendBuffer(partial_item->marshalled_set_buf.data(),
                         partial_item->marshalled_set_buf.size());
        zmq::message_t msg(buf.release_raw(), buf.size(), free_func, NULL);
        // cout << "Size of message = " << msg.size() << std::endl;
        bool success = zmq_set_request_socket_->send(std::move(msg));
        if (!success) {
          cout << "Thread failed to send cluster set over zmq!\n";
        }
        //cout << "Set sent with id: [" << id << "] \n";
      }
    }
  });

  int total_received = 0;
  response_queue_thread_ = thread([this, &total_received]() {
    // get msg from work queue (output)
    // encode
    // put in work queue zmq
    // repeat
    while (run_) {
      MarshalledResponse response;
      bool msg_received = zmq_recv_socket_->recv(&response.msg, ZMQ_NOBLOCK);
      // basically implements polling to avoid blocking recv calls
      // not ideal but whatever, this is a research project! :-D
      if (!msg_received) {
        // cout << "res qt -- Stuck here.\n";
        continue;
      }
      total_received++;
      // cout << "response set size "  << response.Set().NumClusters() << "
      // clusters\n";

      response_queue_->push(std::move(response));
    }

    cout << "Work queue thread ending. Total received: " << total_received
         << std::endl;
  });

  incomplete_request_queue_thread_ = thread([this, &total_received]() {
    // get msg from incomplete request queue
    // put into work queue
    while (run_) {
      zmq::message_t msg;
      bool msg_received =
          zmq_incomplete_request_socket_->recv(&msg, ZMQ_NOBLOCK);
      if (!msg_received) {
        // cout << "irqt -- Stuck here.\n";
        continue;
      }
      total_received++;
      // read as Marshalled request and push into request_queue_
      MarshalledRequest request;
      request.buf.AppendBuffer(reinterpret_cast<const char*>(msg.data()),
                               msg.size());
      PartialRequestHeader* h =
          reinterpret_cast<PartialRequestHeader*>(msg.data());
      cout << "Received incomplete request with ID: " << h->id << std::endl;
      request_queue_->push(std::move(request));
    }
  });

  // partial mergers in this thread may need to be more fully parallelized to
  // prevent bottlenecks. May be required to use a different structure for
  // tracking partial mergers rather than the current map, which needs to be
  // locked
    
  AllAllDist allalldist(request_queue_.get(), sequences_, std::string(params.output_dir));

  auto worker_func = [this, &outstanding_requests, &allalldist]() {
    // read from result queue
    // if is a batch result and is small enough, push to WorkManager
    // if is partial result (ID will be
    // in merge map), lookup and merge with partial_set, if partial set
    // complete,
    //    push to WorkManager
    MarshalledResponse response;
    while (run_) {
      if (!response_queue_->pop(response)) {
        // cout << "worker_func -- no response\n";
        continue;
      }

      auto id = response.ID();
      auto type = response.Type();
      PartialMergeItem* partial_item;

      if (type == RequestType::Batch) {
        MarshalledClusterSet new_set(response);
        // no need to sort for batches
        sets_to_merge_queue_->push(std::move(new_set));
        outstanding_requests--;

      } else if (type == RequestType::Partial) {
        {
          absl::MutexLock l(&mu_);
          auto partial_it = partial_merge_map_.find(id);
          if (partial_it == partial_merge_map_.end()) {
            cout << "error thing in map was not -1, was " << id << std::endl;
            exit(0);
          }
          partial_item = &partial_it->second;
        }

        auto indexes = response.Indexes();
        int start_index = std::get<0>(indexes);
        int end_index = std::get<1>(indexes);
        int cluster_index = std::get<2>(indexes);

        partial_item->partial_set.MergeClusterSet(response.Set(), start_index,
                                                  end_index, cluster_index);
        auto val = partial_item->num_received++;

        // check if all have been received
        if (partial_item->num_expected - 1 == val) {
          MarshalledClusterSet set;
          if (outstanding_merges_ == 0) {
            cout << "building final marshalled set, time: "
                 << static_cast<long int>(std::time(0)) << std::endl;
          }
          partial_item->partial_set.BuildMarshalledSet(&set);
          set.SortSet();
          if (outstanding_merges_ == 0) {
            cout << "done, time: " << static_cast<long int>(std::time(0))
                 << std::endl;
          }
          {
            if (outstanding_merges_ == 1) {
              cout << "last request complete, 1 merge left, time: "
                   << static_cast<long int>(std::time(0)) << std::endl;
            }
            absl::MutexLock l(&mu_);
            partial_merge_map_.erase(id);
          }
          sets_to_merge_queue_->push(std::move(set));
          outstanding_requests--;
        }
      } else if (type == RequestType::Alignment) {
       //cout << "received matches\n";
        const char* matches_buf = reinterpret_cast<const char*>(response.msg.data());
        allalldist.ProcessResult(matches_buf);

      } else {
        cout << "Response was not of any type!!!!\n";
        exit(0);
      }

    }  // end while
  };

  worker_threads_.reserve(params.num_threads);
  for (size_t i = 0; i < params.num_threads; i++) {
    worker_threads_.push_back(std::thread(worker_func));
  }

  cout << "loading to marshalled sets\n";
  // dump all sequences in single cluster sets into the queue
  // for now assumes all of this will fit in memory
  // even a million sequences would just be a few MB

  // here is where we load in checkpointed state
  // if checkpoint dir detected, ask if user wants to load
  // else, is there a set of clusters to add more data to,
  // load here (as last thing in sets to merge queue)
  char response = 'n';
  if (params.checkpoint_interval &&
      CheckpointFileExists(params.checkpoint_dir)) {
    if (params.load_checkpoint_auto) {
       response = 'y';
    } else {
       auto prompt = absl::StrCat("Checkpoint found at ", params.checkpoint_dir,
                                  ". Do you want to load it? (y/n):");
       while (PromptForChar(prompt, response)) {
         if ((response == 'y') | (response == 'n')) {
           break;
         }
       }
    }
  }

  if (response == 'n') {
    int i = 0;
    for (const auto& s : sequences_) {
      if (params.dataset_limit > 0) {
        if (i == params.dataset_limit) break;
      }
      MarshalledClusterSet set(s.ID());
      sets_to_merge_queue_->push(std::move(set));
      i++;
    }
  } else {
    // load the checkpoint
    agd::Status stat =
        LoadCheckpointFile(params.checkpoint_dir, sets_to_merge_queue_);
    if (!stat.ok()) {
      return stat;
    }
    outstanding_merges_ = sets_to_merge_queue_->size() - 1;
    cout << "Loaded checkpoint, outstanding merges: " << outstanding_merges_
         << std::endl;
  }
  cout << "done\n";

  auto t0 = std::chrono::high_resolution_clock::now();
  // just use 'this' thread to schedule outgoing work
  // take stuff from sets_to_merge_ and schedule the work in
  // batches or split partial merges for larger sets
  std::vector<MarshalledClusterSet> sets;
  sets.resize(2);
  while (outstanding_merges_ > 0) {
    // check checkpoint timer
    // if time to checkpoint, wait till all in flight are complete, and
    // the sets to merge queue is stable
    // then dump sets to merge state to checkpoint file
    if (params.checkpoint_interval > 0 &&
        timestamp() - checkpoint_timer_ > params.checkpoint_interval) {
      cout << "Checkpointing, waiting for outstanding requests...\n";
      auto start_wait_checkpoint_time = std::chrono::high_resolution_clock::now();

      while (outstanding_requests.load() > 0) {
        cout << "Waiting to checkpoint, " << outstanding_requests.load()
             << " requests outstanding ...\n";
        std::this_thread::sleep_for(500ms);
      }
      auto end_wait_checkpoint_time = std::chrono::high_resolution_clock::now();
      auto wait_duration = end_wait_checkpoint_time - start_wait_checkpoint_time;
      cout << "Writing checkpoint after " 
           << std::chrono::duration_cast<std::chrono::seconds>(wait_duration).count()
           << "sec waiting for outstanding requests ..." << std::endl;
      // write sets to merge queue
      agd::Status stat =
          WriteCheckpointFile(params.checkpoint_dir, sets_to_merge_queue_);
      if (!stat.ok()) {
        return stat;
      }
      checkpoint_timer_ = timestamp();
      cout << "Checkpoint complete" << std::endl;
    }

    if (!sets_to_merge_queue_->pop(sets[0])) {
      continue;
    }
    // get next so we have two to merge
    if (!sets_to_merge_queue_->pop(sets[1])) {
      return agd::errors::Internal(
          "error: did not get set for second to merge with.");
    }

    // form request, push to queue
    if (sets[0].NumClusters() < params.batch_size ||
        sets[1].NumClusters() < params.batch_size) {
      // create a batch, they are small
      // cout << "two sets are small, batching ...\n";

      uint32_t total_clusters = sets[0].NumClusters() + sets[1].NumClusters();
      // cout << "total clusters: " << total_clusters << std::endl;
      MarshalledRequest request;
      request.CreateBatchRequest(-1);
      request.AddSetToBatch(sets[0]);
      request.AddSetToBatch(sets[1]);

      outstanding_merges_--;
      while (total_clusters < params.batch_size && outstanding_merges_ > 0) {
        // add more cluster sets to batch
        if (!sets_to_merge_queue_->pop(sets[0])) {
          return agd::errors::Internal(
              "ERROR: did not get set for second to merge with.");
        }

        total_clusters += sets[0].NumClusters();
        request.AddSetToBatch(sets[0]);
        outstanding_merges_--;
      }
      // cout << "Sending a batch request -1\n";
      request_queue_->push(std::move(request));
      outstanding_requests++;

    } else {
      outstanding_merges_--;
      PartialMergeItem item;
      item.num_received = 0;

      //cout << "Swapping\n";
      // use outstanding merges as id
      if (sets[0].NumClusters() < sets[1].NumClusters()) {
        std::swap(sets[0], sets[1]);
      }

      //cout << "starting parse\n";
      // iterate through the second cluster set to determine num_expected
      MarshalledClusterView cluster, cluster2;

      uint32_t num_chunks = 0, num_seqs = 0;
      uint32_t pi = -1, ci = 0;  // ci: current index, pi: previous index
      sets[1].Reset();
      sets[0].Reset();

      while (sets[0].NextCluster(&cluster)) {
        sets[1].Reset();
        num_seqs = 0;
        pi = -1, ci = 0;  // ci: current index, pi: previous index
        while (sets[1].NextCluster(&cluster2)) {
          num_seqs += cluster2.NumSeqs();
          if (num_seqs > params.nseqs_threshold - cluster.NumSeqs()) {
            if (ci == pi + 1) {
              pi = ci;
              num_seqs = 0;
            } else {
              pi = ci - 1;
              num_seqs = cluster2.NumSeqs();
            }
            num_chunks++;
          }
          ci++;
        }
        num_chunks += 1;
      }

      item.num_expected = num_chunks;
      /*std::cout << "Num expected: " << item.num_expected << " "
                << " set1 clusters: " << sets[0].NumClusters() << std::endl;*/
      // Reset calls done in function
      item.partial_set.Init(sets[0], sets[1]);
      item.marshalled_set_buf.AppendBuffer(sets[1].buf.data(),
                                           sets[1].buf.size());

      {
        absl::MutexLock l(&mu_);
        partial_merge_map_.insert_or_assign(outstanding_merges_,
                                            std::move(item));
      }

      uint32_t total_cluster = sets[0].NumClusters();
      uint32_t cluster_index = 0;

      sets[0].Reset();
      while (sets[0].NextCluster(&cluster)) {
        sets[1].Reset();
        num_seqs = 0;
        pi = -1, ci = 0;  // ci: current index, pi: previous index

        while (sets[1].NextCluster(&cluster2)) {
          num_seqs += cluster2.NumSeqs();
          if (num_seqs > params.nseqs_threshold - cluster.NumSeqs()) {
            MarshalledRequest request;
            if (ci == pi + 1) {
              request.CreatePartialRequest(outstanding_merges_, cluster,
                                                   pi + 1, ci, cluster_index);
              pi = ci;
              num_seqs = 0;
            } else {
              request.CreatePartialRequest(
                  outstanding_merges_, cluster, pi + 1, ci - 1, cluster_index);
              pi = ci - 1;
              num_seqs = cluster2.NumSeqs();
            }
            request_queue_->push(std::move(request));
          }
          ci++;
        }

        // form the request for the last chunk
        MarshalledRequest request;
        request.CreatePartialRequest(outstanding_merges_, cluster,
                                             pi + 1, ci - 1, cluster_index);
        request_queue_->push(std::move(request));

        cluster_index++;
      }
      assert(cluster_index == total_cluster);

      outstanding_requests++;
      if (outstanding_merges_ == 1) {
        cout << "last request sent, 1 merge left, time: "
             << static_cast<long int>(std::time(0)) << std::endl;
      }
    }
    time_t now_time = std::time(0);
    cout << "[" << std::put_time(std::localtime(&now_time), "%F %T") << "] " 
         << "outstanding merges: " << outstanding_merges_ << std::endl;
  }

  // we must now wait for the last results to come in
  // wait for worker thread to push last merged set
  // TODO add a timeout or something?
  cout << "done and waiting for final result...\n";
  // while (sets_to_merge_queue_->size() != 1);;

  MarshalledClusterSet final_set;
  sets_to_merge_queue_->peek(final_set);
  auto t1 = std::chrono::high_resolution_clock::now();
  cout << "final set size is " << final_set.NumClusters() << " clusters" << std::endl;
  cout << "partial merge map size is " << partial_merge_map_.size() << std::endl;
  cout << "sets to merge size is " << sets_to_merge_queue_->size() << std::endl;
  cout << "request queue size " << request_queue_->size() << std::endl;
  cout << "response queue size " << response_queue_->size() << std::endl;

  auto duration = t1 - t0;
  auto sec = std::chrono::duration_cast<std::chrono::seconds>(duration);
  cout << "Clustering execution time: " << sec.count() << " seconds." << std::endl;

  std::ofstream timing_file("dist_timing.txt", std::ofstream::out);
  timing_file << sec.count() << std::endl;
  timing_file.close();
  

  ClusterSet set(final_set, sequences_);
  std::vector<string> placeholder = {"dist_placeholder"};
  string json_output_file = absl::StrCat(params.output_dir, "dist_clusters.json");
  set.DumpJson(json_output_file, placeholder);

  if (!params.exclude_allall) {
    cout << "scheduling all-all alignments on workers..." << std::endl;
    /*AllAllExecutor executor(std::thread::hardware_concurrency(), 500, &envs,
                            &aligner_params);*/

    set.ScheduleAlignments(&allalldist, sequences_);
    allalldist.Finish();
    /*executor.Initialize();
    set.ScheduleAlignments(&executor);
    executor.FinishAndOutput("dist_output_dir");*/
  } else {
    cout << "Skipping all-all alignments ..." << std::endl;
  }

  cout << "clustering complete!! Joining threads ..." << std::endl;

  run_ = false;
  response_queue_->unblock();
  request_queue_->unblock();
  sets_to_merge_queue_->unblock();
  incomplete_request_queue_->unblock();
  for (auto& t : worker_threads_) {
    t.join();
  }
  request_queue_thread_.join();
  response_queue_thread_.join();
  incomplete_request_queue_thread_.join();
  set_request_thread_.join();

  cout << "All threads joined." << std::endl;

  return agd::Status::OK();
}
