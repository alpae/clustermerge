
cc_library(
    name = "agd",
    srcs = glob(["*.cc"]),
    hdrs = glob(["*.h", "*.hpp"]),
    linkopts = ["-lz"],
    copts = ["-std=c++14"],
    visibility = ["//visibility:public"],
    deps = ["@com_google_absl//absl/strings",
        "@com_google_absl//absl/synchronization",
        "@json//:json-cpp"
    ]
)