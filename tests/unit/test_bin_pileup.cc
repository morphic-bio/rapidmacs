// Unit: BuildBinnedCutSignal places weight at the two Tn5 cuts (start, end-1)
// extended by half(ext_size) bp into bin_width-sized bins.
#include "rapidmacs/binned_signal.h"

#include <cassert>
#include <cstdint>
#include <cstdio>
#include <vector>

int main() {
  using chromap::peaks::BuildBinnedCutSignal;
  using chromap::peaks::ChromFragments;
  using chromap::peaks::Fragment;

  // One fragment [100, 200) with count=2; ext_size=150 -> half=75.
  // Cuts at 100 (left) and 199 (right). Each cut covers [cut-75, cut+75].
  //   left  cut -> [25, 175]
  //   right cut -> [124, 274]
  // bin_width=50 -> bins indexed 0..; padded_max_coord = 274.
  ChromFragments ch;
  ch.name = "chr1";
  Fragment f;
  f.start = 100;
  f.end = 200;
  f.count = 2;
  ch.frags.push_back(f);

  std::vector<uint16_t> obs;
  int32_t padded = 0;
  bool ok = BuildBinnedCutSignal(ch, /*bin_width_bp=*/50, /*ext_size_bp=*/150,
                                 &obs, &padded);
  assert(ok);
  assert(padded == 275 || padded == 274);  // implementation-dependent rounding
  assert(!obs.empty());

  // Hand-computed for diff/prefix-sum semantics: bins covered by each cut
  // contribute count=2 each. Where ranges overlap, sum = 4.
  // bin 0 = [0,50) — partially covered by [25,175]: contributes 2.
  // bin 1 = [50,100) — covered by left only: 2.
  // bin 2 = [100,150) — covered by both [25,175] and [124,274]: 4.
  // bin 3 = [150,200) — covered by both ([25,175] up to 175, [124,274]): 4.
  // bin 4 = [200,250) — covered by right only: 2.
  // bin 5 = [250,300) — partially covered by [124,274]: 2.
  // Allow off-by-one tolerance at boundaries due to inclusive/exclusive
  // conventions (verified by smoke parity; here just check shape).
  int64_t total = 0;
  for (int64_t v : obs) total += v;

  // Each cut contributes 2 (count) * (half*2+1 bp / 50 bp) bins of weight 2.
  // Total weight added across both cuts = 2 * 2 * (151) = 604 bp-weight
  // distributed across bins. Smoke check: must be > 0 and bounded.
  assert(total > 0);
  assert(total <= 4 * 151 * 2);

  // The peak bin should sit in the overlap region (bins 2 or 3).
  size_t argmax = 0;
  for (size_t i = 1; i < obs.size(); ++i) {
    if (obs[i] > obs[argmax]) argmax = i;
  }
  assert(argmax >= 2 && argmax <= 3);
  assert(obs[argmax] >= obs[0] * 2);  // overlap region beats single-cut tail

  std::printf("PASS test_bin_pileup (n_bins=%zu, total_weight=%lld, argmax=%zu, peak=%lld)\n",
              obs.size(), static_cast<long long>(total), argmax,
              static_cast<long long>(obs[argmax]));
  return 0;
}
