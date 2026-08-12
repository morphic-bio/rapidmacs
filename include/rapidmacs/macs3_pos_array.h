#ifndef RAPIDMACS_MACS3_POS_ARRAY_H_
#define RAPIDMACS_MACS3_POS_ARRAY_H_

#include <cstdint>
#include <vector>

namespace chromap {
namespace peaks {

// Build the per-chrom pos_array positions that MACS3 callpeak emits via
// pileup_treat_ctrl_a_chromosome -> __chrom_pair_treat_ctrl, for the
// FRAG-without-control path used by `macs3 callpeak -f FRAG`. Used by the
// opt-in --macs3-frag-summit-parity mode to split tied summit plateaus at
// the same per-event boundaries that MACS3 sees.
//
// Independently written from the behaviour observed in MACS3 v3.0.3; no MACS3
// code is reproduced here. See NOTICE. Algorithm:
//   1. treat_pv positions: walk sorted `starts` and sorted `ends` arrays
//      in parallel; emit unique positions; skip when start == end at the
//      same step (matches Pileup.quick_pileup else-branch).
//   2. ctrl_pv positions (d = llocal_bp): build
//        start_poss = sort((starts - half) ∪ (ends - half))
//        end_poss   = sort((starts + half) ∪ (ends + half))
//      walk them with the same algorithm (matches
//      Pileup.se_all_in_one_pileup).
//   3. pos_array = sorted unique union(treat_pv, ctrl_pv).
//
// Inputs:
//   - starts, ends: per-chrom fragment endpoint arrays. MUST be
//     count-duplicated upstream (each fragment with count=N contributes
//     N starts and N ends). Need not be pre-sorted; this function sorts
//     internally.
//   - llocal_bp: local lambda window size (MACS3 default 10000). Must be
//     positive.
//
// Output:
//   - *out: cleared and filled with sorted unique pos_array positions.
//
// Validated: matches MACS3's pos_array byte-for-byte on the chr4
// worst-offender peak (76382069-76382561) of the full PBMC 3k benchmark
// (96/96 positions, 0 differences) when fed the same count-duplicated
// fragment endpoints.
void BuildMacs3PosArrayPositions(std::vector<int32_t> starts,
                                 std::vector<int32_t> ends,
                                 int32_t llocal_bp,
                                 std::vector<int32_t>* out);

// Controlled PE mode uses the treatment fragment-span pileup together with
// the maximum of control endpoint-window pileups at d, slocal, and llocal.
// MACS3 preserves every raw position from those arrays while merging them,
// even when adjacent values are equal and the exported bedGraph coalesces the
// boundary. Reconstruct that exact merged position array for summit ties.
void BuildMacs3ControlledPosArrayPositions(
    std::vector<int32_t> treatment_starts,
    std::vector<int32_t> treatment_ends,
    std::vector<int32_t> control_starts,
    std::vector<int32_t> control_ends,
    const std::vector<int32_t>& control_window_bp,
    int32_t reference_length,
    std::vector<int32_t>* out);

}  // namespace peaks
}  // namespace chromap

#endif  // RAPIDMACS_MACS3_POS_ARRAY_H_
