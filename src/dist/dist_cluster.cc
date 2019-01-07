
// Entry point for the master server for the distributed
// clustering implementation

#include <fstream>
#include <iostream>
#include <thread>
#include "args.h"
#include "controller.h"
#include "src/common/alignment_environment.h"
#include "src/common/params.h"
#include "src/dataset/load_dataset.h"
#include "worker.h"
/*#include "src/dist/proto/requests.pb.h"
#include "zmq.hpp"*/

using namespace std;

constexpr char cluster_format[] =
    "{\n"
    " \"controller\": \"<ip/addr>\",\n"
    " \"request_queue_port\": <port num>,\n"
    " \"response_queue_port\": <port num>,\n"
    "}\n";

/*
Server cluster format example
{
  "controller": "<ip/addr>",
  "request_queue_port": <port num>,
  "pull_port": <port num>
}
*/

#define DEFAULT_QUEUE_DEPTH 5
#define DEFAULT_RESPONSE_QUEUE_PORT 5556
#define DEFAULT_REQUEST_QUEUE_PORT 5555

int main(int argc, char* argv[]) {
  args::ArgumentParser parser("ClusterMerge",
                              "Bottom up protein cluster merge.");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::ValueFlag<std::string> json_data_dir(
      parser, "data_dir",
      "Directory containing alignment environment matrices in JSON "
      "(logPAM1.json, all_matrices.json) [data/matrices/json]",
      {'d', "data_dir"});
  args::ValueFlag<unsigned int> threads_arg(
      parser, "threads",
      absl::StrCat("Number of threads to use for all-all [",
                   std::thread::hardware_concurrency(), "]"),
      {'t', "threads"});
  args::ValueFlag<unsigned int> queue_depth_arg(
      parser, "queue_depth", "Depth of the local work and response queue",
      {'q', "queue_depth"});
  args::ValueFlag<std::string> input_file_list(
      parser, "file_list", "JSON containing list of input AGD datasets.",
      {'i', "input_list"});
  args::ValueFlag<std::string> output_dir(
      parser, "output_dir",
      "Output directory. Will be overwritten if exists."
      "[./output_matches]",
      {'o', "output_dir"});
  args::ValueFlag<std::string> server_config_file(
      parser, "server_config",
      absl::StrCat(
          "JSON file containing the cluster server configuration, example: \n",
          cluster_format),
      {'s', "server_config"});
  args::PositionalList<std::string> datasets_opts(
      parser, "datasets",
      "AGD Protein datasets to cluster. If present, will override `input_list` "
      "argument.");
  args::Flag controller(parser, "controller",
                        "Designate this process as the cluster controller.",
                        {'C', "controller"});

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

  // parse the server cluster config file to see if we are a worker or the
  // controller

  bool is_controller = false;
  if (controller) {
    is_controller = args::get(controller);
  }

  json server_config_json;
  string controller_ip;
  int request_queue_port;
  int response_queue_port;
  if (server_config_file) {
    string server_config_path = args::get(server_config_file);
    std::ifstream server_config_stream(server_config_path);

    if (!server_config_stream.good()) {
      std::cerr << "File " << server_config_path << " not found.\n";
      return 1;
    }

    server_config_stream >> server_config_json;

  } else {
    cout << "The server config file (--server_config) is required.\n";
    return 1;
  }

  auto controller_it = server_config_json.find("controller");
  if (controller_it == server_config_json.end() && !is_controller) {
    cout << "This process is not controller and server config lacks controller "
            "addr\n";
    return 1;
  } else {
    controller_ip = *controller_it;
  }

  auto request_queue_port_it = server_config_json.find("request_queue_port");
  if (request_queue_port_it == server_config_json.end()) {
    request_queue_port = DEFAULT_REQUEST_QUEUE_PORT;  // default
  } else {
    request_queue_port = *request_queue_port_it;
  }

  auto response_queue_port_it = server_config_json.find("response_queue_port");
  if (response_queue_port_it == server_config_json.end()) {
    response_queue_port = DEFAULT_RESPONSE_QUEUE_PORT;  // default
  } else {
    response_queue_port = *response_queue_port_it;
  }

  unsigned int threads = std::thread::hardware_concurrency();
  if (threads_arg) {
    threads = args::get(threads_arg);
    // do not allow more than hardware threads
    if (threads > std::thread::hardware_concurrency()) {
      threads = std::thread::hardware_concurrency();
    }
  }
  cout << "Using " << threads << " hardware threads for alignment.\n";

  unsigned int queue_depth = DEFAULT_QUEUE_DEPTH;
  if (queue_depth_arg) {
    queue_depth = args::get(queue_depth_arg);
  }

  // get output dir to use, only needed if controller
  string dir("output_matches");
  if (output_dir) {
    dir = args::get(output_dir);
  }
  cout << "Using " << dir << " for output.\n";

  // load alignment envs and initialize (this is for SWPS3)
  string json_dir_path = "data/matrices/json/";
  if (json_data_dir) {
    json_dir_path = args::get(json_data_dir);
    // OSX  STILL  doesn't have experimental/filesystem so this isn't portable
    if (json_dir_path.back() != '/') {
      json_dir_path += '/';
    }
  }

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

  if (is_controller) {
    // launch controller(push_port, pull_port)
    Controller controller;
    Status stat = controller.Run(threads, queue_depth, controller_ip, request_queue_port,
                   response_queue_port, json_dir_path, datasets);
    if (!stat.ok()) {
      cout << "Error: " << stat.error_message() << "\n";
    }
  } else {
    // load datasets, launch worker
    // launch worker(args, controller_ip, push_port, pull_port)
    Worker worker;
    Status stat = worker.Run(threads, queue_depth, json_dir_path, controller_ip,
               request_queue_port, response_queue_port, datasets);
    if (!stat.ok()) {
      cout << "Error: " << stat.error_message() << "\n";
    }
  }

  return 0;
}