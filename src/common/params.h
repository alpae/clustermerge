
#pragma once

// TODO rename to something less generic e.g AlignerParams
struct Parameters {
  int min_score = 181;
  //bool subsequence_homology = true;
  //size_t max_representatives = 1;
  size_t max_n_aa_not_covered = 15;
  // min score for full merge
  float min_full_merge_score = 250.0f;
  bool use_blosum = false;
};
