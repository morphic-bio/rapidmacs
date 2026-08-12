#include "rapidmacs/binned_signal.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

namespace chromap {
namespace peaks {

namespace {

int32_t CutLeft(int32_t start) { return start; }
int32_t CutRight(int32_t end) { return end - 1; }

}  // namespace

bool BuildBinnedCutSignal(const ChromFragments& ch, int32_t bin_width_bp,
                          int32_t ext_size_bp, std::vector<uint16_t>* obs,
                          int32_t* padded_max_coord) {
  obs->clear();
  *padded_max_coord = 0;
  if (ch.frags.empty() || bin_width_bp < 1 || ext_size_bp < 1) {
    return false;
  }
  const int32_t half = ext_size_bp / 2;
  int32_t max_pos = 0;
  for (const Fragment& f : ch.frags) {
    for (int32_t c : {CutLeft(f.start), CutRight(f.end)}) {
      const int32_t hi = c + half;
      if (hi > max_pos) {
        max_pos = hi;
      }
    }
  }
  if (max_pos < 1) {
    return false;
  }
  *padded_max_coord = max_pos;
  const int64_t n_bins =
      (static_cast<int64_t>(max_pos) + bin_width_bp - 1) / bin_width_bp;
  if (n_bins < 1 || n_bins > (1 << 30)) {
    return false;
  }
  std::vector<int64_t> diff(static_cast<size_t>(n_bins) + 1, 0);
  for (const Fragment& f : ch.frags) {
    for (int32_t c : {CutLeft(f.start), CutRight(f.end)}) {
      int32_t L = c - half;
      const int32_t R = c + half;
      if (R <= 0) {
        continue;
      }
      if (L < 0) {
        L = 0;
      }
      int32_t b0 = L / bin_width_bp;
      int32_t b1 = (R - 1) / bin_width_bp;
      if (b0 > b1) {
        continue;
      }
      b1 = std::min(b1, static_cast<int32_t>(n_bins) - 1);
      if (b0 < 0) {
        b0 = 0;
      }
      const int64_t w = static_cast<int64_t>(f.count);
      diff[static_cast<size_t>(b0)] += w;
      if (static_cast<int64_t>(b1) + 1 < static_cast<int64_t>(diff.size())) {
        diff[static_cast<size_t>(b1) + 1] -= w;
      }
    }
  }
  obs->resize(static_cast<size_t>(n_bins), 0);
  int64_t acc = 0;
  constexpr int64_t kSatMax = std::numeric_limits<uint16_t>::max();
  for (int64_t i = 0; i < n_bins; ++i) {
    acc += diff[static_cast<size_t>(i)];
    int64_t v = acc;
    if (v < 0) v = 0;          // belt-and-suspenders; diffs balance to non-negative
    if (v > kSatMax) v = kSatMax;
    (*obs)[static_cast<size_t>(i)] = static_cast<uint16_t>(v);
  }
  return n_bins > 0;
}

}  // namespace peaks
}  // namespace chromap
