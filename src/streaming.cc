#include "rapidmacs/streaming.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>

namespace macs3 {

namespace {

constexpr uint32_t kBinSatMax = std::numeric_limits<uint16_t>::max();

inline void SaturatingAdd(uint16_t& bin, uint32_t delta) {
  uint32_t v = static_cast<uint32_t>(bin) + delta;
  if (v > kBinSatMax) v = kBinSatMax;
  bin = static_cast<uint16_t>(v);
}

void AddCutToBins(int32_t cut, int32_t half, int32_t bin_width,
                  uint32_t weight, std::vector<uint16_t>* bins,
                  int32_t* max_bin_inclusive) {
  int32_t L = cut - half;
  const int32_t R = cut + half;
  if (R <= 0) return;
  if (L < 0) L = 0;
  int32_t b0 = L / bin_width;
  int32_t b1 = (R - 1) / bin_width;
  if (b0 > b1) return;
  if (b0 < 0) b0 = 0;
  if (b1 >= static_cast<int32_t>(bins->size())) {
    bins->resize(static_cast<size_t>(b1) + 1, 0);
  }
  for (int32_t b = b0; b <= b1; ++b) {
    SaturatingAdd((*bins)[static_cast<size_t>(b)], weight);
  }
  if (b1 > *max_bin_inclusive) *max_bin_inclusive = b1;
}

}  // namespace

StreamingPeakCaller::StreamingPeakCaller(const StreamingParams& params)
    : params_(params) {}

bool StreamingPeakCaller::CallPeaks(FragmentIterator& it,
                                    std::vector<chromap::peaks::Peak>* out_peaks,
                                    std::string* err,
                                    StreamingStats* out_stats) {
  if (out_peaks == nullptr) return false;
  out_peaks->clear();
  StreamingStats stats;

  const int32_t half = params_.ext_size / 2;
  std::unordered_map<int32_t, std::vector<uint16_t>> bins_by_chrom;
  std::unordered_map<int32_t, int32_t> max_bin;

  FragmentRecord rec;
  while (it.Next(&rec)) {
    if (rec.chrom_id < 0) continue;
    if (rec.end <= rec.start) continue;
    auto& bins = bins_by_chrom[rec.chrom_id];
    auto mb_ins = max_bin.emplace(rec.chrom_id, -1);
    int32_t& mb = mb_ins.first->second;
    AddCutToBins(rec.start, half, params_.bin_width, rec.count, &bins, &mb);
    AddCutToBins(rec.end - 1, half, params_.bin_width, rec.count, &bins, &mb);
    ++stats.fragments_consumed;
  }
  const bool stream_ok = it.Error().empty();
  if (err != nullptr) *err = it.Error();

  const auto& names = it.ChromNames();
  for (auto& kv : bins_by_chrom) {
    const int32_t cid = kv.first;
    auto& obs = kv.second;
    auto mb_it = max_bin.find(cid);
    if (mb_it == max_bin.end() || mb_it->second < 0) continue;
    obs.resize(static_cast<size_t>(mb_it->second) + 1, 0);
    for (uint16_t v : obs) {
      ++stats.total_bins;
      if (v == kBinSatMax) ++stats.saturated_bins;
    }
    if (cid < 0 || cid >= static_cast<int32_t>(names.size())) continue;
    std::vector<chromap::peaks::Peak> out;
    chromap::peaks::CallPeaksOnBins(names[static_cast<size_t>(cid)],
                                    params_.bin_width, obs,
                                    params_.local_window_bins, params_.fdr,
                                    params_.p_value, params_.merge_bin_gap,
                                    &out);
    for (const auto& p : out) out_peaks->push_back(p);
  }

  if (out_stats != nullptr) *out_stats = stats;
  return stream_ok;
}

}  // namespace macs3
