
#include <fstream>
#include <iostream>
#include "args.h"
#include "dataset/agd_protein_dataset.h"
#include "dataset/fasta_dataset.h"
#include "dataset/load_dataset.h"
#include "src/common/aligner.h"
#include "src/common/all_all_executor.h"
#include "src/common/bottom_up_merge.h"
#include "src/common/debug.h"

using std::cout;
using std::string;
using std::unique_ptr;

int main(int argc, char** argv) {
  args::ArgumentParser parser("ClusterMerge",
                              "Bottom up protein cluster merge.");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::ValueFlag<unsigned int> threads_arg(
      parser, "threads",
      absl::StrCat("Number of threads to use for all-all [",
                   std::thread::hardware_concurrency(), "]"),
      {'t', "threads"});
  args::ValueFlag<unsigned int> cluster_threads_arg(
      parser, "cluster_threads",
      absl::StrCat("Number of threads to use for clustering [",
                   std::thread::hardware_concurrency(), "]"),
      {'c', "cluster-threads"});
  args::ValueFlag<unsigned int> merge_threads_arg(
      parser, "merge_threads",
      absl::StrCat("Number of threads to use for merging [",
                   std::thread::hardware_concurrency(), "]"),
      {'m', "merge-threads"});
  args::ValueFlag<unsigned int> dup_removal_threshold_arg(
      parser, "duplicate removal threshold",
      "How big a set of clusters should be before duplicates are filtered out "
      "[MAX_INT]",
      {'r', "dup_removal_thresh"});
  args::ValueFlag<std::string> json_data_dir(
      parser, "data_dir",
      "Directory containing alignment environment matrices in JSON "
      "(logPAM1.json, all_matrices.json) [data/matrices/json]",
      {'d', "data_dir"});
  args::ValueFlag<std::string> output_dir(
      parser, "output_dir",
      "Output directory. Will be overwritten if exists."
      "[./output_matches]",
      {'o', "output_dir"});
  args::ValueFlag<std::string> aligner_params_arg(
      parser, "aligner parameters",
      "JSON containing alignment and clustering parameters "
      "[data/default_aligner_params.json].",
      {'a', "aligner_params"});
  args::ValueFlag<std::string> input_file_list(
      parser, "file_list", "JSON containing list of input AGD/FASTA datasets.",
      {'i', "input_list"});
  args::PositionalList<std::string> datasets_opts(
      parser, "datasets",
      "AGD/FASTA Protein datasets to cluster. If present, will override "
      "`input_list` "
      "argument.");
  args::Flag exclude_allall(
      parser, "exclude_allall",
      "Don't perform intra-cluster all-all alignment, just do the clustering.",
      {'x', "exclude_allall"});

  args::ValueFlag<std::string> file_name(
      parser, "file",
      "Adds clustering data from an already clustered json file and merges "
      "with clusters from current dataset",
      {'f', "file"});

  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  } catch (args::ValidationError e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  std::string input_file_name_temp = args::get(input_file_list);
  std::vector<std::string> dataset_file_names;
  std::vector<unique_ptr<Dataset>> datasets_old;
  json dataset_json_obj;

  if (file_name) {
    agd::Status s_old;

    std::ifstream dataset_stream(args::get(file_name));

    if (!dataset_stream.good()) {
      s_old = agd::errors::NotFound("No such file: ", args::get(file_name));
    }

    if (!s_old.ok()) {
      cout << s_old.ToString() << "\n";
      return 0;
    }

    dataset_stream >> dataset_json_obj;

    for (const auto& old_dataset : dataset_json_obj["datasets"]) {
      s_old = LoadDatasetsJSON(old_dataset, &datasets_old);
      dataset_file_names.push_back(old_dataset);

      if (!s_old.ok()) {
        cout << s_old.ToString() << "\n";
        return 0;
      }
    }
  }

  dataset_file_names.push_back(input_file_name_temp);

  unsigned int threads = std::thread::hardware_concurrency();
  if (threads_arg) {
    threads = args::get(threads_arg);
    // do not allow more than hardware threads
    if (threads > std::thread::hardware_concurrency()) {
      threads = std::thread::hardware_concurrency();
    }
  }
  // cout << "Using " << threads << " hardware threads for alignment.\n";

  unsigned int cluster_threads = std::thread::hardware_concurrency();
  if (cluster_threads_arg) {
    cluster_threads = args::get(cluster_threads_arg);
    // do not allow more than hardware threads
    if (cluster_threads > std::thread::hardware_concurrency()) {
      cluster_threads = std::thread::hardware_concurrency();
    }
  }
  // cout << "Using " << cluster_threads << " hardware threads for
  // clustering.\n";

  unsigned int merge_threads = std::thread::hardware_concurrency();
  if (merge_threads_arg) {
    merge_threads = args::get(merge_threads_arg);
    // do not allow more than hardware threads
    if (merge_threads > std::thread::hardware_concurrency()) {
      merge_threads = std::thread::hardware_concurrency();
    }
  }
  // cout << "Using " << merge_threads << " hardware threads for merging.\n";

  uint32_t dup_removal_threshold = UINT_MAX;
  if (dup_removal_threshold_arg) {
    dup_removal_threshold = args::get(dup_removal_threshold_arg);
  }

  // get output dir to use
  string dir("output_matches");
  if (output_dir) {
    dir = args::get(output_dir);
  }
  cout << "Using " << dir << " for output.\n";

  json aligner_params_json;
  if (aligner_params_arg) {
    string aligner_params_file = args::get(aligner_params_arg);

    std::ifstream aligner_params_stream(aligner_params_file);

    if (!aligner_params_stream.good()) {
      std::cerr << "Provided file " << aligner_params_file << " not found.\n";
      return 1;
    }

    aligner_params_stream >> aligner_params_json;
  } else {
    string aligner_params_file("data/default_aligner_params.json");
    std::ifstream aligner_params_stream(aligner_params_file);

    if (!aligner_params_stream.good()) {
      std::cerr << "Default file " << aligner_params_file << " not found.\n";
      return 1;
    }

    aligner_params_stream >> aligner_params_json;
  }

  Parameters aligner_params;

  auto min_score_it = aligner_params_json.find("min_score");
  if (min_score_it != aligner_params_json.end()) {
    aligner_params.min_score = *min_score_it;
  }

  auto max_aa_uncovered_it = aligner_params_json.find("max_aa_uncovered");
  if (max_aa_uncovered_it != aligner_params_json.end()) {
    aligner_params.max_n_aa_not_covered = *max_aa_uncovered_it;
  }

  auto min_full_merge_score_it =
      aligner_params_json.find("min_full_merge_score");
  if (min_full_merge_score_it != aligner_params_json.end()) {
    aligner_params.min_full_merge_score = *min_full_merge_score_it;
  }

  auto blosum_it = aligner_params_json.find("blosum");
  if (blosum_it != aligner_params_json.end()) {
    aligner_params.use_blosum = *blosum_it;
  }

  // load alignment envs and initialize (this is for SWPS3)
  string json_dir_path = "data/matrices/json/";
  if (json_data_dir) {
    json_dir_path = args::get(json_data_dir);
    // OSX  STILL  doesn't have experimental/filesystem so this isn't portable
    if (json_dir_path.back() != '/') {
      json_dir_path += '/';
    }
  }

  string logpam_json_file = json_dir_path + "logPAM1.json";
  string blosum_json_file = json_dir_path + "BLOSUM62.json";
  string all_matrices_json_file = json_dir_path + "all_matrices.json";

  std::ifstream logpam_stream(logpam_json_file);
  std::ifstream blosum_stream(blosum_json_file);
  std::ifstream allmat_stream(all_matrices_json_file);

  if (!logpam_stream.good()) {
    std::cerr << "File " << logpam_json_file << " not found.\n";
    return 1;
  }

  if (!blosum_stream.good()) {
    std::cerr << "File " << blosum_json_file << " not found.\n";
    return 1;
  }

  if (!allmat_stream.good()) {
    std::cerr << "File " << all_matrices_json_file << " not found.\n";
    return 1;
  }

  json logpam_json;
  logpam_stream >> logpam_json;

  json blosum_json;
  blosum_stream >> blosum_json;

  json all_matrices_json;
  allmat_stream >> all_matrices_json;

  AlignmentEnvironments envs;

  // initializing envs is expensive, so don't copy this

  cout << "Initializing alignment environments from " << json_dir_path << "\n";
  envs.InitFromJSON(logpam_json, all_matrices_json, aligner_params.min_score);
  if (aligner_params.use_blosum) {
    cout << "using blosum62 matrix ...\n";
    envs.UseBlosum(blosum_json, aligner_params.min_score);
  }
  cout << "Done.\n";

  // done init envs

  std::vector<unique_ptr<Dataset>> datasets;
  agd::Status s;

  if (datasets_opts) {
    // load and parse protein datasets
    // cluster merge sequences are simply string_views over this data
    // so these structures must live until computations complete
    if (input_file_list) {
      cout << "WARNING: ignoring input file list and using positionals!\n";
    }
    s = LoadDatasetsPositional(datasets_opts, &datasets);
  } else if (input_file_list) {
    s = LoadDatasetsJSON(args::get(input_file_list), &datasets);
  } else {
    cout << "No datasets provided. See usage: \n";
    std::cerr << parser;
    exit(0);
  }

  if (!s.ok()) {
    cout << s.ToString() << "\n";
    return 0;
  }

  // init aligner object
  ProteinAligner aligner(&envs, &aligner_params);

  // build initial clustersets
  // one sequence, in one cluster, in one set
  // then, we bottom-up merge them

  cout << "Datasets loaded ...\n";

  // Add by akash
  // TODO deduplicate this code
  if (file_name) {
    BottomUpMerge merger(dataset_json_obj, datasets_old, datasets, &aligner);

    AllAllExecutor executor(threads, 1000, &envs, &aligner_params);
    executor.Initialize();

    auto t0 = std::chrono::high_resolution_clock::now();
    MergeExecutor merge_executor(merge_threads, 200, &envs, &aligner_params);

    // Add by akash
    merger.RunMulti(cluster_threads, dup_removal_threshold, &executor,
                    &merge_executor, !exclude_allall, dataset_file_names);

    // Call recombine and execute RunMulti again

    // merger.DebugDump();
    // wait and finish call on executor
    // which dumps final results to disk
    executor.FinishAndOutput(dir);
    auto t1 = std::chrono::high_resolution_clock::now();

    auto duration = t1 - t0;
    auto sec = std::chrono::duration_cast<std::chrono::seconds>(duration);

    cout << "Execution time: " << sec.count() << " seconds.\n";

  } else {
    BottomUpMerge merger(datasets, &aligner);

    AllAllExecutor executor(threads, 1000, &envs, &aligner_params);
    executor.Initialize();

    auto t0 = std::chrono::high_resolution_clock::now();
    MergeExecutor merge_executor(merge_threads, 200, &envs, &aligner_params);

    // Add by akash
    merger.RunMulti(cluster_threads, dup_removal_threshold, &executor,
                    &merge_executor, !exclude_allall, dataset_file_names);

    // merger.DebugDump();
    // wait and finish call on executor
    // which dumps final results to disk
    executor.FinishAndOutput(dir);
    auto t1 = std::chrono::high_resolution_clock::now();

    auto duration = t1 - t0;
    auto sec = std::chrono::duration_cast<std::chrono::seconds>(duration);

    cout << "Execution time: " << sec.count() << " seconds.\n";
  }
  return (0);
}
