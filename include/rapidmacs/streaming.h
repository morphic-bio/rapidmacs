#ifndef RAPIDMACS_STREAMING_H_
#define RAPIDMACS_STREAMING_H_

#include <cstdint>
#include <string>
#include <vector>

#include "rapidmacs/call_peaks.h"
#include "rapidmacs/fragments.h"

namespace macs3 {

struct StreamingParams {
  int32_t bin_width = 50;
  int32_t ext_size = 150;
  int32_t local_window_bins = 200;
  double fdr = 0.05;
  double p_value = 0.01;
  int32_t merge_bin_gap = 0;
};

// Per-run statistics filled in by StreamingPeakCaller::CallPeaks.
struct StreamingStats {
  int64_t fragments_consumed = 0;
  int64_t saturated_bins = 0;  // bins that hit UINT16_MAX
  int64_t total_bins = 0;
};

// Single-pass peak caller. Drains a FragmentIterator into per-chrom uint16
// bin arrays (saturating at 65535/bin), then runs CallPeaksOnBins per chrom.
// Fragments are dropped immediately after binning — peak RSS is dominated by
// the bin arrays themselves (~124 MB for hg38 at 50 bp bins).
class StreamingPeakCaller {
 public:
  explicit StreamingPeakCaller(const StreamingParams& params);

  // Drains `it`. On success returns true and fills *out_peaks; on iterator
  // parse / I/O error returns false, sets *err to the iterator's error
  // message, and leaves *out_peaks containing peaks derived from any
  // partial stream consumed before the error (useful for debugging).
  // *out_stats (if non-null) is filled in regardless of success.
  bool CallPeaks(FragmentIterator& it,
                 std::vector<chromap::peaks::Peak>* out_peaks,
                 std::string* err,
                 StreamingStats* out_stats = nullptr);

 private:
  StreamingParams params_;
};

}  // namespace macs3

#endif  // RAPIDMACS_STREAMING_H_
