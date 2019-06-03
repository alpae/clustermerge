
#include "bottom_up_merge.h"
#include <iostream>

using std::cout;
using std::string;

BottomUpMerge::BottomUpMerge(std::vector<std::unique_ptr<Dataset>>& datasets,
                             ProteinAligner* aligner) {
  aligner_ = aligner;

  const char* data;
  size_t size;
  uint32_t id = 0;
  for (auto& dataset : datasets) {
    cout << "Parsing dataset " << dataset->Name() << " ...\n";

    auto s = dataset->GetNextRecord(&data, &size);
    uint32_t genome_index = 0;
    while (s.ok()) {
      // coverages_.push_back(string());
      // coverages_.back().resize(size);

      // cout << "Adding sequence id " << id << "\n";
      if (size > 60000) {
        cout << "over size " << size << "\n";
        exit(0);
      }
      Sequence seq(absl::string_view(data, size),  // coverages_.back(),
                   dataset->Name(), dataset->Size(), genome_index++, id++);

      ClusterSet cs(seq);

      sets_.push_back(std::move(cs));
      s = dataset->GetNextRecord(&data, &size);
    }
  }
}

// Add by akash
BottomUpMerge::BottomUpMerge(
    nlohmann::json dataset_json_obj,
    std::vector<std::unique_ptr<Dataset>>& datasets_old,
    std::vector<std::unique_ptr<Dataset>>& datasets, ProteinAligner* aligner) {
  // Add by akash : create a map from the of dataset and use that to store into
  // old cluster set data structure
  const char* data_old;
  size_t size_old;
  uint32_t id_old = 0;

  // Create a map
  //std::map<std::string, std::vector<Sequence>> dataset_old_map;
  std::vector<Sequence> old_sequences;

  for (auto& dataset_old : datasets_old) {
    cout << "Parsing dataset " << dataset_old->Name() << " ...\n";

    auto s_old = dataset_old->GetNextRecord(&data_old, &size_old);
    uint32_t genome_index_old = 0;
    while (s_old.ok()) {
      // coverages_.push_back(string());
      // coverages_.back().resize(size);

      // cout << "Adding sequence id " << id << "\n";
      if (size_old > 60000) {
        cout << "over size " << size_old << "\n";
        exit(0);
      }

      Sequence seq(absl::string_view(data_old, size_old),  // coverages_.back(),
                   dataset_old->Name(), dataset_old->Size(), genome_index_old++,
                   id_old++);

      // ClusterSet cs(seq);

      // sets_.push_back(std::move(cs));
      //dataset_old_map[dataset_old->Name()].push_back(seq);
      old_sequences.push_back(std::move(seq));
      s_old = dataset_old->GetNextRecord(&data_old, &size_old);
    }
  }

  //std::map<std::string, std::vector<Sequence>>::iterator it;

  // debugging
  // for(it=dataset_old_map.begin();it!=dataset_old_map.end();it++)
  // std::cout<<it->first<<std::endl;

  // Create old cluster set, by traversing through json obj
  long long int counter = 0;
  long long int clusters_size = dataset_json_obj["clusters"].size();

  for (const auto& cluster : dataset_json_obj["clusters"]) {
    Cluster c;
    for (const auto& seq : cluster) {
      int abs_index = seq["AbsoluteIndex"];
      c.AddSequence(old_sequences[abs_index]);
    }
    old_set_.AddCluster(c);
  }

  /*for (; counter < clusters_size; counter++) {
    Cluster c;
    long long int seq_num = dataset_json_obj["clusters"][counter].size();
    for (long long int seq_ind = 0; seq_ind < seq_num; seq_ind++) {
      int abs_index = dataset_json_obj["clusters"][counter][seq_ind]["AbsoluteIndex"];
      //std::string genome_name =
      //    dataset_json_obj["clusters"][counter][seq_ind]["Genome"];
      //long long int genome_index =
       //   dataset_json_obj["clusters"][counter][seq_ind]["Index"];
      //c.AddSequence(dataset_old_map[genome_name][genome_index]);
      c.AddSequence(old_sequences[abs_index]);
      // std::cout<<"Add genome: "<<genome_name<<" with index:
      // "<<genome_index<<std::endl;
    }

    old_set_.AddCluster(c);
  }*/

  aligner_ = aligner;

  const char* data;
  size_t size;
  uint32_t id = id_old;
  for (auto& dataset : datasets) {
    cout << "Parsing dataset " << dataset->Name() << " ...\n";

    auto s = dataset->GetNextRecord(&data, &size);
    uint32_t genome_index = 0;
    while (s.ok()) {
      // coverages_.push_back(string());
      // coverages_.back().resize(size);

      // cout << "Adding sequence id " << id << "\n";
      if (size > 60000) {
        cout << "over size " << size << "\n";
        exit(0);
      }
      Sequence seq(absl::string_view(data, size),  // coverages_.back(),
                   dataset->Name(), dataset->Size(), genome_index++, id++);

      ClusterSet cs(seq);

      sets_.push_back(std::move(cs));
      s = dataset->GetNextRecord(&data, &size);
    }
  }
}

void BottomUpMerge::DebugDump() {
  cout << "Dumping merger ... \n";
  for (const auto& cs : sets_) {
    cs.DebugDump();
  }
}

agd::Status BottomUpMerge::RunMulti(
    size_t num_threads, size_t dup_removal_threshold, AllAllExecutor* executor,
    MergeExecutor* merge_executor, bool do_allall,
    std::vector<std::string>& datasetsFileName) {
  cluster_sets_left_ = sets_.size();

  // launch threads, join threads
  auto cluster_worker = [this, &merge_executor, &dup_removal_threshold]() {
    // cout << "cluster worker starting\n";
    // need own aligner per thread
    ProteinAligner aligner(aligner_->Envs(), aligner_->Params());

    // is atomic actually needed here?
    while (cluster_sets_left_.load() > 1) {
      queue_mu_.Lock();
      if (sets_.size() > 1 && cluster_sets_left_.load() > 1) {  // enough to pop
        auto s1 = std::move(sets_.front());
        sets_.pop_front();
        auto s2 = std::move(sets_.front());
        sets_.pop_front();
        cluster_sets_left_.fetch_sub(1);
        // cout << "cluster sets left: " << cluster_sets_left_.load() << "\n";
        queue_mu_.Unlock();

        // this part takes a while for larger sets
        auto t0 = std::chrono::high_resolution_clock::now();

        // eventually we may want to call single thread mergeClusters for
        // small cluster sets as it may be more efficient

        // swap so we have the larger set first, this results
        // in a larger number of smaller work items
        if (s1.Size() < s2.Size()) {
          s1.Swap(&s2);
          // cout << "swapped89\n";
        }

        auto merged_set = s1.MergeClustersParallel(s2, merge_executor);
        auto t1 = std::chrono::high_resolution_clock::now();

        auto duration = t1 - t0;
        auto sec = std::chrono::duration_cast<std::chrono::seconds>(duration);
        // cout << "merge took: " << sec.count() << " seconds.\n";

        if (merged_set.Size() > dup_removal_threshold) {
          merged_set.RemoveDuplicates();
        }

        queue_mu_.Lock();
        sets_.push_back(std::move(merged_set));
        queue_pop_cv_.SignalAll();
        queue_mu_.Unlock();

      } else if (sets_.size() <= 1) {  // wait until enough
        while (sets_.size() <= 1 && cluster_sets_left_.load() > 1) {
          /*cout << "waiting, set size: " << sets_.size()
               << " sets left: " << cluster_sets_left_.load() << "\n";*/
          queue_pop_cv_.Wait(&queue_mu_);
        }
        if (cluster_sets_left_.load() <= 1) {
          queue_mu_.Unlock();
          break;
        }

        auto s1 = std::move(sets_.front());
        sets_.pop_front();
        auto s2 = std::move(sets_.front());
        sets_.pop_front();
        cluster_sets_left_.fetch_sub(1);
        // cout << "cluster sets left: " << cluster_sets_left_.load() << "\n";
        queue_mu_.Unlock();

        // swap so we have the larger set first, this results
        // in a larger number of smaller work items
        if (s1.Size() < s2.Size()) {
          s1.Swap(&s2);
          // cout << "swapped127\n";
        }

        // this part takes a while for larger sets
        auto t0 = std::chrono::high_resolution_clock::now();
        auto merged_set = s1.MergeClustersParallel(s2, merge_executor);
        auto t1 = std::chrono::high_resolution_clock::now();

        auto duration = t1 - t0;
        auto sec = std::chrono::duration_cast<std::chrono::seconds>(duration);
        // cout << "merge took: " << sec.count() << " seconds.\n";

        if (merged_set.Size() > dup_removal_threshold) {
          merged_set.RemoveDuplicates();
        }

        queue_mu_.Lock();
        sets_.push_back(std::move(merged_set));
        queue_pop_cv_.SignalAll();
        queue_mu_.Unlock();

      } else {  // only one cluster set left, done.
        queue_mu_.Unlock();
        break;
      }
    }

    // cout << "Cluster eval thread finishing...\n";
  };

  cout << "Clustering " << cluster_sets_left_.load() << " sequences...\n";

  auto t0 = std::chrono::high_resolution_clock::now();
  threads_.reserve(num_threads);
  for (size_t i = 0; i < num_threads; i++) {
    threads_.push_back(std::thread(cluster_worker));
  }

  for (auto& t : threads_) {
    t.join();
  }

  if (old_set_.Size() >= 1) {
    std::cout << "Merging data of older set with the new result\n";
    auto s1 = std::move(sets_.front());
    sets_.pop_front();
    //auto s2 = std::move(old_set_);
    auto merged_set = s1.MergeClustersParallel(old_set_, merge_executor);
    sets_.push_back(std::move(merged_set));
  }

  auto t1 = std::chrono::high_resolution_clock::now();

  auto duration = t1 - t0;
  auto sec = std::chrono::duration_cast<std::chrono::seconds>(duration);
  cout << "Clustering execution time (s): " << sec.count() << "\n";

  assert(sets_.size() == 1);
  // now we are all finished clustering
  auto& final_set = sets_[0];

  // Add by akash
  final_set.DumpJson("clusters.json", datasetsFileName);
  cout << "Total clusters: " << final_set.Size()
       << ", dumped to clusters.json.\n";

  // for all clusters in final set, schedule all-all alignments with executor
  if (do_allall) {
    cout << "Scheduling all-all alignments ...\n";
    final_set.ScheduleAlignments(executor);
    cout << "Finished alignment scheduling. \n";
  } else {
    cout << "Skipping all all computation ...\n";
  }

  return agd::Status::OK();
}

// Add by akash
agd::Status BottomUpMerge::Run(AllAllExecutor* executor,
                               size_t dup_removal_threshold, bool do_allall,
                               std::vector<std::string>& datasetsFileName) {
  auto t0 = std::chrono::high_resolution_clock::now();
  while (sets_.size() > 1) {
    // dequeue 2 sets
    // merge the sets into one
    // push onto queue

    auto& s1 = sets_[0];
    auto& s2 = sets_[1];

    /*cout << "Merging cluster sets of size " << s1.Size()
      << " and " << s2.Size() << "\n";
    cout << "Cluster sets remaining: " << sets_.size() - 2 << "\n";*/
    // s1.DebugDump();
    // cout << "\nand\n";
    // s2.DebugDump();
    auto merged_set = s1.MergeClusters(s2, aligner_);

    if (merged_set.Size() > dup_removal_threshold) {
      merged_set.RemoveDuplicates();
    }
    sets_.pop_front();
    sets_.pop_front();
    sets_.push_back(std::move(merged_set));
  };

  auto t1 = std::chrono::high_resolution_clock::now();

  auto duration = t1 - t0;
  auto sec = std::chrono::duration_cast<std::chrono::seconds>(duration);
  cout << "Clustering execution time: " << sec.count() << " seconds.\n";

  auto& final_set = sets_[0];

  // Add by akash
  final_set.DumpJson("clusters.json", datasetsFileName);

  // for all clusters in final set, schedule all-all alignments with executor
  if (do_allall) {
    final_set.ScheduleAlignments(executor);
  } else {
    cout << "Skipping all all computation ...\n";
  }

  cout << "Total clusters: " << final_set.Size() << "\n";

  return agd::Status::OK();
}
