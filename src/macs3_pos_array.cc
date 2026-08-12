#include "rapidmacs/macs3_pos_array.h"

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <vector>

namespace chromap {
namespace peaks {

namespace {

// Mirror MACS3.Signal.Pileup.quick_pileup / se_all_in_one_pileup main loop:
// walk sorted_starts and sorted_ends in parallel, emit unique positions,
// skip when start[i_s] == end[i_e] at the same step (the else-branch that
// drops "no-op" coincidences where a +1 start and a -1 end at the same
// position cancel out).
void EmitParallelWalkPositions(const std::vector<int32_t>& sorted_starts,
                               const std::vector<int32_t>& sorted_ends,
                               std::vector<int32_t>* out) {
  out->clear();
  const size_t ls = sorted_starts.size();
  const size_t le = sorted_ends.size();
  if (ls == 0 || le == 0) {
    return;
  }
  out->reserve(ls + le);

  // MACS3 emits an initial "first chunk of 0" position at min(starts[0],
  // ends[0]) when that position > 0. Mirror this even though it lands at
  // chrom-leftmost and rarely affects peak summits.
  int32_t pre_p = std::min(sorted_starts[0], sorted_ends[0]);
  if (pre_p > 0) {
    out->push_back(pre_p);
  } else {
    pre_p = 0;
  }

  size_t i_s = 0;
  size_t i_e = 0;
  while (i_s < ls && i_e < le) {
    const int32_t s = sorted_starts[i_s];
    const int32_t e = sorted_ends[i_e];
    if (s < e) {
      if (s != pre_p) {
        out->push_back(s);
        pre_p = s;
      }
      ++i_s;
    } else if (s > e) {
      if (e != pre_p) {
        out->push_back(e);
        pre_p = e;
      }
      ++i_e;
    } else {
      // s == e: drop the no-op pair (start +1 and end -1 cancel; MACS3
      // does not emit a pos_array entry here).
      ++i_s;
      ++i_e;
    }
  }
  while (i_e < le) {
    const int32_t e = sorted_ends[i_e];
    if (e != pre_p) {
      out->push_back(e);
      pre_p = e;
    }
    ++i_e;
  }
  while (i_s < ls) {
    const int32_t s = sorted_starts[i_s];
    if (s != pre_p) {
      out->push_back(s);
      pre_p = s;
    }
    ++i_s;
  }
}

}  // namespace

void BuildMacs3PosArrayPositions(std::vector<int32_t> starts,
                                 std::vector<int32_t> ends,
                                 int32_t llocal_bp,
                                 std::vector<int32_t>* out) {
  if (out == nullptr) {
    return;
  }
  out->clear();
  if (starts.empty() || ends.empty() ||
      starts.size() != ends.size() || llocal_bp <= 0) {
    return;
  }

  std::sort(starts.begin(), starts.end());
  std::sort(ends.begin(), ends.end());

  // treat_pv positions (PE-mode quick_pileup): walk starts vs ends.
  std::vector<int32_t> treat_pos;
  EmitParallelWalkPositions(starts, ends, &treat_pos);

  // ctrl_pv positions (se_all_in_one_pileup with d = llocal_bp):
  //   start_poss = sort((starts - half) ∪ (ends - half))
  //   end_poss   = sort((starts + half) ∪ (ends + half))
  // The shifted-down positions get clamped to 0 (MACS3 fix_coordinates).
  // We don't clamp the upper bound to chromosome length; downstream
  // callers consume only positions that fall inside merged ppois regions
  // (well within chrom bounds), so clamping at the right edge has no
  // effect on summit-splitting.
  const int32_t half = llocal_bp / 2;
  const size_t n = starts.size();
  std::vector<int32_t> ctrl_starts;
  std::vector<int32_t> ctrl_ends;
  ctrl_starts.reserve(2 * n);
  ctrl_ends.reserve(2 * n);
  for (size_t i = 0; i < n; ++i) {
    const int32_t s = starts[i];
    const int32_t e = ends[i];
    ctrl_starts.push_back(s - half < 0 ? 0 : s - half);
    ctrl_starts.push_back(e - half < 0 ? 0 : e - half);
    ctrl_ends.push_back(s + half);
    ctrl_ends.push_back(e + half);
  }
  std::sort(ctrl_starts.begin(), ctrl_starts.end());
  std::sort(ctrl_ends.begin(), ctrl_ends.end());

  std::vector<int32_t> ctrl_pos;
  EmitParallelWalkPositions(ctrl_starts, ctrl_ends, &ctrl_pos);

  // Free intermediates before the union.
  std::vector<int32_t>().swap(ctrl_starts);
  std::vector<int32_t>().swap(ctrl_ends);

  // pos_array = sorted union of treat_pv and ctrl_pv positions. Both
  // inputs are already sorted+unique (the parallel-walk emit does the
  // dedup), so std::set_union produces the merged sorted+unique result.
  out->reserve(treat_pos.size() + ctrl_pos.size());
  std::set_union(treat_pos.begin(), treat_pos.end(),
                 ctrl_pos.begin(), ctrl_pos.end(),
                 std::back_inserter(*out));
}

void BuildMacs3ControlledPosArrayPositions(
    std::vector<int32_t> treatment_starts,
    std::vector<int32_t> treatment_ends,
    std::vector<int32_t> control_starts,
    std::vector<int32_t> control_ends,
    const std::vector<int32_t>& control_window_bp,
    int32_t reference_length,
    std::vector<int32_t>* out) {
  if (out == nullptr) {
    return;
  }
  out->clear();
  if (treatment_starts.empty() || treatment_ends.empty() ||
      treatment_starts.size() != treatment_ends.size() ||
      control_starts.empty() || control_ends.empty() ||
      control_starts.size() != control_ends.size() ||
      control_window_bp.empty() || reference_length <= 0) {
    return;
  }

  std::sort(treatment_starts.begin(), treatment_starts.end());
  std::sort(treatment_ends.begin(), treatment_ends.end());
  std::vector<int32_t> treatment_pos;
  EmitParallelWalkPositions(treatment_starts, treatment_ends, &treatment_pos);

  std::vector<int32_t> control_pos;
  for (int32_t window_bp : control_window_bp) {
    if (window_bp <= 0) {
      continue;
    }
    const int32_t half = window_bp / 2;
    std::vector<int32_t> shifted_starts;
    std::vector<int32_t> shifted_ends;
    shifted_starts.reserve(control_starts.size() * 2);
    shifted_ends.reserve(control_ends.size() * 2);
    auto clamp = [&](int64_t value) {
      if (value < 0) return static_cast<int32_t>(0);
      if (value > reference_length) return reference_length;
      return static_cast<int32_t>(value);
    };
    for (size_t i = 0; i < control_starts.size(); ++i) {
      shifted_starts.push_back(
          clamp(static_cast<int64_t>(control_starts[i]) - half));
      shifted_starts.push_back(
          clamp(static_cast<int64_t>(control_ends[i]) - half));
      shifted_ends.push_back(
          clamp(static_cast<int64_t>(control_starts[i]) + half));
      shifted_ends.push_back(
          clamp(static_cast<int64_t>(control_ends[i]) + half));
    }
    std::sort(shifted_starts.begin(), shifted_starts.end());
    std::sort(shifted_ends.begin(), shifted_ends.end());
    std::vector<int32_t> window_pos;
    EmitParallelWalkPositions(shifted_starts, shifted_ends, &window_pos);
    if (control_pos.empty()) {
      control_pos.swap(window_pos);
    } else {
      std::vector<int32_t> merged;
      merged.reserve(control_pos.size() + window_pos.size());
      std::set_union(control_pos.begin(), control_pos.end(),
                     window_pos.begin(), window_pos.end(),
                     std::back_inserter(merged));
      control_pos.swap(merged);
    }
  }

  out->reserve(treatment_pos.size() + control_pos.size());
  std::set_union(treatment_pos.begin(), treatment_pos.end(),
                 control_pos.begin(), control_pos.end(),
                 std::back_inserter(*out));
}

}  // namespace peaks
}  // namespace chromap
