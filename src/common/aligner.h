#pragma once

#include "alignment_environment.h"
#include "params.h"
#include "swps3/extras.h"

// copied aligner class from Persona, so we are using AGD Status here
// quite a bit

class ProteinAligner {
 public:
  const static double pam_list_[];

  struct Alignment {
    int seq1_min;
    int seq1_max;
    int seq2_min;
    int seq2_max;
    const AlignmentEnvironment* env;
    double score;
    double pam_distance;
    double pam_variance;
    int seq1_length;  // len of aligned string minus '_'
    int seq2_length;  // len of aligned string minus '_'
    std::string ToString() {
      std::stringstream ss;
      ss << "score: " << score << " s1min " << seq1_min << " s1max " << seq1_max
         << " s2min " << seq2_min << " s2max " << seq2_max;
      return ss.str();
    }
  };

  ProteinAligner(const AlignmentEnvironments* envs, const Parameters* params)
      : envs_(envs), params_(params) {
    bt_data_ = new BTData();
  }

  ~ProteinAligner() {
    if (buf1_) {
      delete[] buf1_;
    }
    delete bt_data_;
  }

  agd::Status AlignLocal(const char* seq1, const char* seq2, int seq1_len,
                         int seq2_len, Alignment& result);

  bool PassesThreshold(const char* seq1, const char* seq2, int seq1_len,
                       int seq2_len);

  bool LogPamPassesThreshold(const char* seq1, const char* seq2, int seq1_len,
                       int seq2_len);

  // with full range calc
  agd::Status AlignDouble(const char* seq1, const char* seq2, int seq1_len,
                          int seq2_len, bool stop_at_threshold,
                          Alignment& result, const AlignmentEnvironment& env);

  // const AlignmentEnvironments* Envs() { return envs_; }

  agd::Status AlignSingle(const char* seq1, const char* seq2, int seq1_len,
                          int seq2_len, Alignment& result);

  const Parameters* Params() { return params_; }
  const AlignmentEnvironments* Envs() { return envs_; }

  size_t NumAlignments() { return num_alignments_; }

 private:
  const AlignmentEnvironments* envs_;
  const Parameters* params_;
  size_t num_alignments_ = 0;

  struct StartPoint {
    Alignment alignment;
    double estimated_pam;
    char* seq1;
    char* seq2;
    int seq1_len;
    int seq2_len;
  };

  void FindStartingPoint(const char* seq1, const char* seq2, int seq1_len,
                         int seq2_len, StartPoint& point);

  int AlignStrings(double* matrix, char* s1, int len1, char* s2, int len2,
                   double escore, char* o1, char* o2, double maxerr,
                   double gap_open, double gap_ext, BTData* data);

  double c_align_double_global(double* matrix, const char* s1, int ls1,
                               const char* s2, int ls2, double gap_open,
                               double gap_ext, BTData* data);

  // these are for reuse over align double/local methods
  char* buf1_ = nullptr;      // MAXSEQLEN
  char* buf2_ = nullptr;      // MAXSEQLEN
  char* savebuf1_ = nullptr;  // MAXSEQLEN
  char* savebuf2_ = nullptr;  // MAXSEQLEN

  BTData* bt_data_;
};
