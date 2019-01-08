
#include "controller.h"
#include <iostream>
#include <fstream>
#include "src/agd/errors.h"
#include "src/common/cluster_set.h"
#include "src/common/all_all_executor.h"

using std::cout;
using std::thread;
using std::string;

// merge other into
void MergePartials(cmproto::ClusterSet& set, const cmproto::ClusterSet& other,
                   uint32_t original_size) {
  // relies on clusters in the sets being in the same order
  // where any new cluster is the last element
  // we do not "fully merge" any of the clusters in `set`

  cout << "merging partial clusters ...\n";
  for (uint32_t i = 0; i < original_size; i++) {
    auto* mut_cluster = set.mutable_clusters(i);
    auto cluster = other.clusters(i);

    for (auto index : cluster.indexes()) {
      if (std::find(mut_cluster->indexes().begin(),
                    mut_cluster->indexes().end(),
                    index) == mut_cluster->indexes().end()) {
        // doesnt exist, add
        mut_cluster->add_indexes(index);
      }
    }
  }
  // if the partial merge generated a new cluster, add it to the set
  if (other.clusters_size() > original_size) {
    assert(original_size == other.clusters_size() - 1);
    auto* c = set.add_clusters();
    c->CopyFrom(other.clusters(other.clusters_size() - 1));
  }
}

agd::Status Controller::Run(size_t num_threads, size_t queue_depth,
                            const std::string& controller_ip,
                            int request_queue_port, int response_queue_port,
                            const std::string& data_dir_path,
                            std::vector<std::unique_ptr<Dataset>>& datasets) {
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

  auto total_merges = sequences_.size() - 1;
  outstanding_merges_ = total_merges;
  cout << "outstanding merges to complete: " << outstanding_merges_ << "\n";

  // create envs, params
  string logpam_json_file = data_dir_path + "logPAM1.json";
  string all_matrices_json_file = data_dir_path + "all_matrices.json";

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
  cout << "Worker initializing environments from " << data_dir_path << "\n";
  envs.InitFromJSON(logpam_json, all_matrices_json);
  cout << "Done.\n";

  // done init envs

  Parameters params;  // using default params for now




  // connect to zmq queues
  auto address = std::string("tcp://*:");
  auto response_queue_address = absl::StrCat(address, response_queue_port);
  auto request_queue_address = absl::StrCat(address, request_queue_port);

  context_ = zmq::context_t(1);
  try {
    zmq_recv_socket_.reset(new zmq::socket_t(context_, ZMQ_PULL));
  } catch (...) {
    return agd::errors::Internal("Could not create zmq PULL socket ");
  }

  try {
    zmq_send_socket_.reset(new zmq::socket_t(context_, ZMQ_PUSH));
  } catch (...) {
    return agd::errors::Internal("Could not create zmq PUSH socket ");
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

  request_queue_.reset(new ConcurrentQueue<cmproto::MergeRequest>(queue_depth));
  response_queue_.reset(new ConcurrentQueue<cmproto::Response>(queue_depth));
  sets_to_merge_queue_.reset(
      new ConcurrentQueue<cmproto::ClusterSet>(sequences_.size()));

  request_queue_thread_ = thread([this]() {
    // get msg from work queue (output)
    // encode
    // put in work queue zmq
    // repeat
    cmproto::MergeRequest merge_request;
    while (run_) {
      if (!request_queue_->pop(merge_request)) {
        continue;
      }
      auto size = merge_request.ByteSizeLong();
      zmq::message_t msg(size);
      cout << "pushing request of size " << size << " of type "
           << (merge_request.has_batch() ? "batch " : "partial ") << "\n";
      auto success = merge_request.SerializeToArray(msg.data(), size);
      if (!success) {
        cout << "Thread failed to serialize request protobuf!\n";
      }

      success = zmq_send_socket_->send(std::move(msg));
      if (!success) {
        cout << "Thread failed to send request over zmq!\n";
      }
    }

    cout << "Work queue thread ending.\n";
  });

  response_queue_thread_ = thread([this]() {
    // get msg from work queue (output)
    // encode
    // put in work queue zmq
    // repeat
    cmproto::Response response;
    zmq::message_t msg;
    while (run_) {
      zmq_recv_socket_->recv(&msg);

      if (!response.ParseFromArray(msg.data(), msg.size())) {
        cout << "Failed to parse merge request protobuf!!\n";
        return;
      }
      cout << "parsed a zmq response with " << msg.size() << " bytes and "
           << response.set().clusters_size() << " clusters\n";

      response_queue_->push(response);
    }

    cout << "Work queue thread ending.\n";
  });

  worker_thread_ = thread([this]() {
    // read from result queue
    // if is a batch result and is small enough, push to WorkManager
    // if is partial result (ID will be
    // in merge map), lookup and merge with partial_set, if partial set
    // complete,
    //    push to WorkManager
    cmproto::Response response;
    while (run_) {
      if (!response_queue_->pop(response)) {
        continue;
      }

      auto id = response.id();
      cout << "repsonse id is " << id << "\n";
      // will need to lock this section
      auto partial_it = partial_merge_map_.find(id);
      if (partial_it != partial_merge_map_.end()) {
        cout << "preparing to merge partial result... \n";
        // its part of a larger merger
        auto& partial_item = partial_it->second;
        // merge response_set into partial_item.partial_set
        MergePartials(partial_item.partial_set, response.set(),
                      partial_item.original_size);
        // check if this was the last one
        partial_item.num_received++;
        if (partial_item.num_expected == partial_item.num_received) {
          sets_to_merge_queue_->push(std::move(partial_item.partial_set));
          // remove partial it, its done now
          partial_merge_map_.erase(partial_it);
        }
      } else {
        // it's a full result
        cout << "pushing full result \n";
        sets_to_merge_queue_->push(std::move(response.set()));
      }
    }
  });

  // dump all sequences in single cluster sets into the queue
  // for now assumes all of this will fit in memory
  // even a million sequences would just be a few MB
  for (const auto& s : sequences_) {
    cmproto::ClusterSet set;
    auto* c = set.add_clusters();
    c->add_indexes(s.ID());
    sets_to_merge_queue_->push(std::move(set));
  }

  cout << "sets to merge queue size is " << sets_to_merge_queue_->size() << "\n";
  cout << "outstanding merges: " << outstanding_merges_ << "\n";
  cout << "sets to merge queue is loaded and ready, make sure workers are ready, press key to compute...\n";
  std::cin.get();

  // just use 'this' thread to schedule outgoing work
  // take stuff from sets_to_merge_ and schedule the work in
  // batches or split partial merges for larger sets
  cmproto::MergeRequest request;
  std::vector<cmproto::ClusterSet> sets;
  sets.resize(2);
  while (outstanding_merges_ > 1) {
    sets[0].Clear();
    sets[1].Clear();
    if (!sets_to_merge_queue_->pop(sets[0])) {
      continue;
    }
    if (outstanding_merges_ == 1) {
      // we are done
      cout << "last set is completed.\n";
    } else if (outstanding_merges_ > 1) {
      // get next so we have two to merge
      if (!sets_to_merge_queue_->pop(sets[1])) {
        return agd::errors::Internal(
            "error: did not get set for second to merge with.");
      }
      cout << "processing two sets ...\n";
      cout << "set one size: " << sets[0].clusters_size()
           << ", set two size: " << sets[1].clusters_size() << "\n";

      // form request, push to queue
      if (sets[0].clusters_size() < 10 || sets[1].clusters_size() < 10) {
        // create a batch, they are small
        cout << "two sets are small, batching ...\n";
        cmproto::MergeBatch* batch = request.mutable_batch();
        uint32_t total_clusters =
            sets[0].clusters_size() + sets[1].clusters_size();

        auto* c = batch->add_sets();
        c->CopyFrom(sets[0]);
        c = batch->add_sets();
        c->CopyFrom(sets[1]);
        while (total_clusters < batch_size_ && outstanding_merges_ > 1) {
          // add more cluster sets to batch
          if (!sets_to_merge_queue_->pop(sets[0])) {
            return agd::errors::Internal(
                "ERROR: did not get set for second to merge with.");
          }

          c = batch->add_sets();
          c->CopyFrom(sets[0]);
          outstanding_merges_--;
        }
        request.set_id(0);
        // if the queue uses copy semantics im not sure how protobufs
        // with submessages will behave
        request_queue_->push(std::move(request));
        request.Clear();

      } else {
        // either set is large enough, split the computation into multiple
        // requests
        cout << "splitting merger of two large sets into partial mergers\n";
        outstanding_merges_--;
        // make a map entry for this multi-part request
        PartialMergeItem item;
        item.num_received = 0;
        // use the outstanding merges as id
        request.set_id(outstanding_merges_);
        if (sets[0].clusters_size() < sets[1].clusters_size()) {
          sets[0].Swap(&sets[1]);
        }
        // each work item does a partial merge of one cluster in sets[0] to
        // all clusters in sets[1]
        item.num_expected = sets[0].clusters_size();
        item.partial_set.CopyFrom(sets[1]);
        item.original_size = sets[1].clusters_size();
        partial_merge_map_.insert_or_assign(outstanding_merges_, item);

        for (const auto& c : sets[0].clusters()) {
          cmproto::MergePartial* partial_request = request.mutable_partial();
          auto* cluster = partial_request->mutable_cluster();
          cluster->CopyFrom(c);
          auto* cluster_set = partial_request->mutable_set();
          cluster_set->CopyFrom(sets[1]);
          cout << "pushing partial request with " << partial_request->set().clusters_size() << " clusters in set\n";
          request_queue_->push(std::move(request));
          request.Clear();
        }
      }
    }
  }

  // we must now wait for the last results to come in
  // wait for worker thread to push last merged set
  // TODO add a timeout or something?
  while (sets_to_merge_queue_->size() != 1);;

  if (sets_to_merge_queue_->size() != 1) {
    cout << "where did the last set go??\n";
  } else {
    cout << "scheduling final alignments on controller...\n";
    cmproto::ClusterSet final_set;
    sets_to_merge_queue_->peek(final_set);

    ClusterSet set(final_set, sequences_);
    AllAllExecutor executor(std::thread::hardware_concurrency(), 500, &envs, &params);
    set.ScheduleAlignments(&executor);
    executor.FinishAndOutput("output_dir");
  }

  cout << "clustering complete!!\n";

  return agd::Status::OK();
}
