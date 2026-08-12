#include "rapidmacs/macs3_frag_workspace.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <iterator>
#include <limits>
#include <map>
#include <queue>
#include <tuple>

#include <htslib/khash.h>
#include <omp.h>

#include "rapidmacs/macs3_pos_array.h"
#include "rapidmacs/peak_io.h"

namespace chromap {
namespace peaks {

namespace {

float Quantize5Float(double x) {
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%.5f", x);
  return std::strtof(buf, nullptr);
}

double Round5Score(double x) {
  return std::round(x * 1e5) / 1e5;
}

int32_t EffectiveFragWeight(int32_t raw_count, int32_t max_count,
                            bool macs3_uint8_counts) {
  int64_t w = raw_count;
  if (max_count > 0) {
    w = std::min<int64_t>(w, max_count);
  }
  if (w <= 0) {
    return 0;
  }
  if (!macs3_uint8_counts) {
    if (w > std::numeric_limits<int32_t>::max()) {
      return std::numeric_limits<int32_t>::max();
    }
    return static_cast<int32_t>(w);
  }
  const auto u = static_cast<uint8_t>(static_cast<unsigned int>(w));
  return static_cast<int32_t>(u);
}

int32_t RowCountForMacsTotals(int32_t raw_count, int32_t max_count,
                              bool macs3_uint8_counts) {
  int64_t t = raw_count;
  if (max_count > 0) {
    t = std::min<int64_t>(t, max_count);
  }
  if (t <= 0) {
    return 0;
  }
  if (macs3_uint8_counts) {
    return static_cast<int32_t>(static_cast<uint8_t>(static_cast<unsigned int>(t)));
  }
  if (t > std::numeric_limits<int32_t>::max()) {
    return std::numeric_limits<int32_t>::max();
  }
  return static_cast<int32_t>(t);
}

void ReleaseEvents(std::vector<Macs3FragEvent>* ev) {
  std::vector<Macs3FragEvent>().swap(*ev);
}

int NormalizeThreads(int requested, size_t n_tasks) {
  if (requested < 1 || n_tasks <= 1) {
    return 1;
  }
  const size_t rt = static_cast<size_t>(requested);
  return static_cast<int>(std::min(rt, n_tasks));
}

bool CoalesceEventPositions(std::vector<Macs3FragEvent>* ev) {
  if (ev->empty()) {
    return true;
  }
  std::sort(ev->begin(), ev->end(),
            [](const Macs3FragEvent& a, const Macs3FragEvent& b) {
              if (a.pos != b.pos) {
                return a.pos < b.pos;
              }
              return a.delta < b.delta;
            });
  size_t wri = 0;
  for (size_t i = 0; i < ev->size(); ++i) {
    if (wri > 0 && ev->at(i).pos == ev->at(wri - 1).pos) {
      const int64_t delta = static_cast<int64_t>(ev->at(wri - 1).delta) +
                            static_cast<int64_t>(ev->at(i).delta);
      if (delta < std::numeric_limits<int32_t>::min() ||
          delta > std::numeric_limits<int32_t>::max()) {
        return false;
      }
      ev->at(wri - 1).delta = static_cast<int32_t>(delta);
    } else {
      ev->at(wri++) = ev->at(i);
    }
  }
  ev->resize(wri);
  return true;
}

void TrackFromEvents(const std::vector<Macs3FragEvent>& ev, Macs3BdgTrack* out) {
  out->p.clear();
  out->v.clear();
  if (ev.empty()) {
    return;
  }
  double z = 0.0;
  double pre_z = -1.0e4;
  int32_t s = 0;
  out->p.reserve(ev.size());
  out->v.reserve(ev.size());
  for (const auto& e : ev) {
    const int32_t pos = e.pos;
    const double v = static_cast<double>(e.delta);
    if (pos != s) {
      if (z == pre_z && !out->p.empty()) {
        out->p.back() = pos;
      } else {
        out->p.push_back(pos);
        out->v.push_back(Quantize5Float(z));
        pre_z = z;
      }
    }
    z += v;
    s = pos;
  }
}

struct WeightedEndpoint {
  WeightedEndpoint() : pos(0), weight(0) {}
  WeightedEndpoint(int32_t p, int32_t w) : pos(p), weight(w) {}

  int32_t pos;
  int32_t weight;
};

struct LocalBdgRun {
  LocalBdgRun() : end(0), value(0.f) {}
  LocalBdgRun(int32_t e, float v) : end(e), value(v) {}

  int32_t end;
  float value;
};

static_assert(sizeof(LocalBdgRun) == 8,
              "controlled-lambda run must remain one 64-bit word");

struct LocalBdgTrack {
  std::vector<LocalBdgRun> runs;
};

bool EnvValueIsOne(const char* name) {
  const char* value = std::getenv(name);
  return value != nullptr && std::strcmp(value, "1") == 0;
}

void BuildSortedControlEndpoints(const ChromFragments& ch, int32_t max_count,
                                 bool macs3_uint8_counts,
                                 std::vector<WeightedEndpoint>* starts,
                                 std::vector<WeightedEndpoint>* ends,
                                 bool sort_ends = true,
                                 double* end_sort_seconds = nullptr) {
  starts->clear();
  ends->clear();
  starts->reserve(ch.frags.size());
  ends->reserve(ch.frags.size());
  for (const Fragment& f : ch.frags) {
    const int32_t w =
        EffectiveFragWeight(f.count, max_count, macs3_uint8_counts);
    if (w <= 0) continue;
    starts->push_back(WeightedEndpoint{f.start, w});
    ends->push_back(WeightedEndpoint{f.end, w});
  }
  const auto by_pos = [](const WeightedEndpoint& a,
                         const WeightedEndpoint& b) {
    return a.pos < b.pos;
  };
  // Fragment inputs are contractually sorted by (start, end). Retain a
  // defensive fallback for programmatic callers that construct the vectors
  // themselves, while avoiding an unnecessary N log N start-coordinate sort
  // for the BAMPE/BEDPE readers.
  if (!std::is_sorted(starts->begin(), starts->end(), by_pos)) {
    if (sort_ends) {
      std::sort(starts->begin(), starts->end(), by_pos);
    } else {
      // The lazy sweep needs each start to remain paired with its fragment's
      // unsorted end. This branch is only a defensive fallback: file-backed
      // inputs are contractually start-sorted.
      std::vector<std::pair<WeightedEndpoint, WeightedEndpoint>> paired;
      paired.reserve(starts->size());
      for (size_t i = 0; i < starts->size(); ++i) {
        paired.push_back(std::make_pair(starts->at(i), ends->at(i)));
      }
      std::sort(
          paired.begin(), paired.end(),
          [&](const std::pair<WeightedEndpoint, WeightedEndpoint>& a,
              const std::pair<WeightedEndpoint, WeightedEndpoint>& b) {
            return by_pos(a.first, b.first);
          });
      for (size_t i = 0; i < paired.size(); ++i) {
        starts->at(i) = paired[i].first;
        ends->at(i) = paired[i].second;
      }
    }
  }
  if (sort_ends) {
    if (end_sort_seconds == nullptr) {
      std::sort(ends->begin(), ends->end(), by_pos);
    } else {
      const auto t0 = std::chrono::steady_clock::now();
      std::sort(ends->begin(), ends->end(), by_pos);
      const auto t1 = std::chrono::steady_clock::now();
      *end_sort_seconds = std::chrono::duration<double>(t1 - t0).count();
    }
  } else if (end_sort_seconds != nullptr) {
    *end_sort_seconds = 0.0;
  }
}

void TrackFromSortedSpanEndpoints(
    const std::vector<WeightedEndpoint>& starts,
    const std::vector<WeightedEndpoint>& ends, Macs3BdgTrack* out) {
  out->p.clear();
  out->v.clear();
  if (starts.empty() || starts.size() != ends.size()) return;

  size_t is = 0;
  size_t ie = 0;
  const size_t n = starts.size();
  out->p.reserve(n * 2);
  double z = 0.0;
  double pre_z = -1.0e4;
  int32_t previous_pos = 0;
  while (is < n || ie < n) {
    int32_t pos = std::numeric_limits<int32_t>::max();
    if (is < n) pos = std::min(pos, starts[is].pos);
    if (ie < n) pos = std::min(pos, ends[ie].pos);
    if (pos != previous_pos) {
      if (z == pre_z && !out->p.empty()) {
        out->p.back() = pos;
      } else {
        out->p.push_back(pos);
        out->v.push_back(Quantize5Float(z));
        pre_z = z;
      }
    }
    int64_t delta = 0;
    while (is < n && starts[is].pos == pos) delta += starts[is++].weight;
    while (ie < n && ends[ie].pos == pos) delta -= ends[ie++].weight;
    z += static_cast<double>(delta);
    previous_pos = pos;
  }
}

// Treatment-only MACS3 local lambda: merge the four monotone endpoint-event
// streams directly.  The streaming heap path used the same unclamped event
// coordinates, so keeping negative left-window coordinates here preserves
// its exact track alignment semantics.
void TrackFromSortedEndpointWindowsUnclamped(
    const std::vector<WeightedEndpoint>& starts,
    const std::vector<WeightedEndpoint>& ends, int32_t window_bp,
    Macs3BdgTrack* out) {
  out->p.clear();
  out->v.clear();
  if (window_bp < 1 || starts.empty() || starts.size() != ends.size()) return;

  const int32_t half = window_bp / 2;
  const int32_t right_shift = window_bp - half;
  size_t ss = 0, se = 0, es = 0, ee = 0;
  const size_t n = starts.size();
  auto source_pos = [&](int source, size_t index) {
    const std::vector<WeightedEndpoint>& src = source < 2 ? starts : ends;
    const int32_t shift = (source == 0 || source == 2) ? -half : right_shift;
    return static_cast<int32_t>(static_cast<int64_t>(src[index].pos) + shift);
  };
  out->p.reserve(n * 4);
  double z = 0.0;
  double pre_z = -1.0e4;
  int32_t previous_pos = 0;
  while (ss < n || se < n || es < n || ee < n) {
    int32_t pos = std::numeric_limits<int32_t>::max();
    if (ss < n) pos = std::min(pos, source_pos(0, ss));
    if (se < n) pos = std::min(pos, source_pos(1, se));
    if (es < n) pos = std::min(pos, source_pos(2, es));
    if (ee < n) pos = std::min(pos, source_pos(3, ee));
    if (pos != previous_pos) {
      if (z == pre_z && !out->p.empty()) {
        out->p.back() = pos;
      } else {
        out->p.push_back(pos);
        out->v.push_back(Quantize5Float(z));
        pre_z = z;
      }
    }
    int64_t delta = 0;
    while (ss < n && source_pos(0, ss) == pos) delta += starts[ss++].weight;
    while (se < n && source_pos(1, se) == pos) delta -= starts[se++].weight;
    while (es < n && source_pos(2, es) == pos) delta += ends[es++].weight;
    while (ee < n && source_pos(3, ee) == pos) delta -= ends[ee++].weight;
    z += static_cast<double>(delta);
    previous_pos = pos;
  }
}

bool TrackFromSortedEndpointWindows(
    const std::vector<WeightedEndpoint>& starts,
    const std::vector<WeightedEndpoint>& ends, int32_t reference_length,
    int32_t window_bp, float scale, float baseline, LocalBdgTrack* out) {
  out->runs.clear();
  if (window_bp < 1 || starts.empty() || starts.size() != ends.size()) {
    return false;
  }
  const int32_t half = window_bp / 2;
  auto clamp = [&](int64_t p) {
    if (p < 0) return static_cast<int32_t>(0);
    if (reference_length > 0 && p > reference_length) {
      return reference_length;
    }
    if (p > std::numeric_limits<int32_t>::max()) {
      return std::numeric_limits<int32_t>::max();
    }
    return static_cast<int32_t>(p);
  };

  // Each source is already sorted, and adding a constant plus clamping is
  // monotone. Merge the four event streams directly instead of materializing
  // and comparison-sorting 4*N (position, delta) records for every window.
  size_t ss = 0, se = 0, es = 0, ee = 0;
  const size_t n = starts.size();
  auto source_pos = [&](int source, size_t index) {
    const std::vector<WeightedEndpoint>& src =
        source < 2 ? starts : ends;
    const int64_t shift = (source == 0 || source == 2) ? -half : half;
    return clamp(static_cast<int64_t>(src[index].pos) + shift);
  };
  out->runs.reserve(n * 4);
  double z = 0.0;
  double pre_z = -1.0e4;
  int32_t previous_pos = 0;
  while (ss < n || se < n || es < n || ee < n) {
    int32_t pos = std::numeric_limits<int32_t>::max();
    if (ss < n) pos = std::min(pos, source_pos(0, ss));
    if (se < n) pos = std::min(pos, source_pos(1, se));
    if (es < n) pos = std::min(pos, source_pos(2, es));
    if (ee < n) pos = std::min(pos, source_pos(3, ee));

    if (pos != previous_pos) {
      if (z == pre_z && !out->runs.empty()) {
        out->runs.back().end = pos;
      } else {
        out->runs.push_back(LocalBdgRun(pos, Quantize5Float(z)));
        pre_z = z;
      }
    }

    int64_t delta = 0;
    while (ss < n && source_pos(0, ss) == pos) {
      delta += starts[ss++].weight;
    }
    while (se < n && source_pos(1, se) == pos) {
      delta -= starts[se++].weight;
    }
    while (es < n && source_pos(2, es) == pos) {
      delta += ends[es++].weight;
    }
    while (ee < n && source_pos(3, ee) == pos) {
      delta -= ends[ee++].weight;
    }
    if (delta < std::numeric_limits<int32_t>::min() ||
        delta > std::numeric_limits<int32_t>::max()) {
      return false;
    }
    z += static_cast<double>(delta);
    previous_pos = pos;
  }

  for (LocalBdgRun& run : out->runs) {
    const float scaled = static_cast<float>(run.value * scale);
    run.value = std::max(scaled, baseline);
  }
  return true;
}

struct LazyWindowEnd {
  LazyWindowEnd() : open(0), close(0), weight(0) {}
  LazyWindowEnd(int32_t o, int32_t c, int32_t w)
      : open(o), close(c), weight(w) {}

  int32_t open;
  int32_t close;
  int32_t weight;
};

struct LazyWindowEndGreater {
  bool operator()(const LazyWindowEnd& a, const LazyWindowEnd& b) const {
    if (a.open != b.open) return a.open > b.open;
    return a.close > b.close;
  }
};

// One controlled local-window sweep. Fragment starts arrive monotonically;
// only their not-yet-opened end windows require a heap. Start-window closes
// are monotone in input order, and end-window closes become monotone as their
// opens leave the heap, so flat cursor queues suffice for those two streams.
class LazyControlWindowSweep {
 public:
  LazyControlWindowSweep()
      : reference_length_(0),
        half_(0),
        scale_(0.f),
        baseline_(0.f),
        out_(nullptr),
        z_(0.0),
        pre_z_(-1.0e4),
        previous_pos_(0),
        have_pending_start_open_(false),
        pending_start_open_pos_(0),
        pending_start_open_delta_(0),
        max_heap_depth_(0),
        ok_(true) {}

  bool Init(int32_t reference_length, int32_t window_bp, float scale,
            float baseline, size_t n, LocalBdgTrack* out) {
    if (window_bp < 1 || out == nullptr) return false;
    reference_length_ = reference_length;
    half_ = window_bp / 2;
    scale_ = scale;
    baseline_ = baseline;
    out_ = out;
    out_->runs.clear();
    out_->runs.reserve(n * 4);
    start_closes_.clear();
    end_closes_.clear();
    z_ = 0.0;
    pre_z_ = -1.0e4;
    previous_pos_ = 0;
    have_pending_start_open_ = false;
    pending_start_open_pos_ = 0;
    pending_start_open_delta_ = 0;
    max_heap_depth_ = 0;
    ok_ = true;
    EndHeap empty;
    end_heap_.swap(empty);
    return true;
  }

  void AddFragment(const WeightedEndpoint& start,
                   const WeightedEndpoint& end) {
    if (!ok_) return;
    const int32_t open = Clamp(static_cast<int64_t>(start.pos) - half_);
    if (!have_pending_start_open_ || open != pending_start_open_pos_) {
      FlushPendingStartOpen();
      DrainBefore(open);
      have_pending_start_open_ = true;
      pending_start_open_pos_ = open;
      pending_start_open_delta_ = 0;
    }
    pending_start_open_delta_ += start.weight;

    const int32_t start_close =
        Clamp(static_cast<int64_t>(start.pos) + half_);
    start_closes_.push_back(WeightedEndpoint(start_close, start.weight));

    const int32_t end_open = Clamp(static_cast<int64_t>(end.pos) - half_);
    const int32_t end_close = Clamp(static_cast<int64_t>(end.pos) + half_);
    end_heap_.push(LazyWindowEnd(end_open, end_close, end.weight));
    max_heap_depth_ = std::max(max_heap_depth_, end_heap_.size());
  }

  bool Finish() {
    FlushPendingStartOpen();
    while (ok_ && HasPending()) {
      const int32_t pos = NextPendingPosition();
      ProcessPosition(pos, 0);
    }
    if (!ok_) return false;
    for (LocalBdgRun& run : out_->runs) {
      const float scaled = static_cast<float>(run.value * scale_);
      run.value = std::max(scaled, baseline_);
    }
    return true;
  }

  size_t max_heap_depth() const { return max_heap_depth_; }

 private:
  using EndHeap = std::priority_queue<LazyWindowEnd,
                                      std::vector<LazyWindowEnd>,
                                      LazyWindowEndGreater>;

  int32_t Clamp(int64_t p) const {
    if (p < 0) return 0;
    if (reference_length_ > 0 && p > reference_length_) {
      return reference_length_;
    }
    if (p > std::numeric_limits<int32_t>::max()) {
      return std::numeric_limits<int32_t>::max();
    }
    return static_cast<int32_t>(p);
  }

  int32_t NextPendingPosition() const {
    int32_t pos = std::numeric_limits<int32_t>::max();
    if (!start_closes_.empty()) pos = std::min(pos, start_closes_.front().pos);
    if (!end_heap_.empty()) pos = std::min(pos, end_heap_.top().open);
    if (!end_closes_.empty()) pos = std::min(pos, end_closes_.front().pos);
    return pos;
  }

  bool HasPending() const {
    return !start_closes_.empty() || !end_heap_.empty() ||
           !end_closes_.empty();
  }

  void EmitClose(int32_t pos) {
    if (z_ == pre_z_ && !out_->runs.empty()) {
      out_->runs.back().end = pos;
    } else {
      out_->runs.push_back(LocalBdgRun(pos, Quantize5Float(z_)));
      pre_z_ = z_;
    }
  }

  void ProcessPosition(int32_t pos, int64_t direct_delta) {
    if (!ok_) return;
    if (pos != previous_pos_) EmitClose(pos);

    int64_t delta = direct_delta;
    while (!start_closes_.empty() && start_closes_.front().pos == pos) {
      delta -= start_closes_.front().weight;
      start_closes_.pop_front();
    }
    while (!end_heap_.empty() && end_heap_.top().open == pos) {
      const LazyWindowEnd e = end_heap_.top();
      end_heap_.pop();
      delta += e.weight;
      end_closes_.push_back(WeightedEndpoint(e.close, e.weight));
    }
    while (!end_closes_.empty() && end_closes_.front().pos == pos) {
      delta -= end_closes_.front().weight;
      end_closes_.pop_front();
    }
    if (delta < std::numeric_limits<int32_t>::min() ||
        delta > std::numeric_limits<int32_t>::max()) {
      ok_ = false;
      return;
    }
    z_ += static_cast<double>(delta);
    previous_pos_ = pos;
  }

  void FlushPendingStartOpen() {
    if (!have_pending_start_open_ || !ok_) return;
    ProcessPosition(pending_start_open_pos_, pending_start_open_delta_);
    have_pending_start_open_ = false;
    pending_start_open_delta_ = 0;
  }

  void DrainBefore(int32_t target) {
    while (ok_ && HasPending()) {
      const int32_t pos = NextPendingPosition();
      if (pos >= target) break;
      ProcessPosition(pos, 0);
    }
  }

  int32_t reference_length_;
  int32_t half_;
  float scale_;
  float baseline_;
  LocalBdgTrack* out_;
  std::deque<WeightedEndpoint> start_closes_;
  std::deque<WeightedEndpoint> end_closes_;
  EndHeap end_heap_;
  double z_;
  double pre_z_;
  int32_t previous_pos_;
  bool have_pending_start_open_;
  int32_t pending_start_open_pos_;
  int64_t pending_start_open_delta_;
  size_t max_heap_depth_;
  bool ok_;
};

// Advance up to three controlled windows during one pass over the aligned
// start/end arrays. Each window retains its own pending-end heap.
bool TrackFromLazyEndpointWindowsShared(
    const std::vector<WeightedEndpoint>& starts,
    const std::vector<WeightedEndpoint>& ends, int32_t reference_length,
    const int32_t* window_bp, const float* scale, float baseline,
    LocalBdgTrack* const* outputs, size_t window_count,
    size_t* max_heap_depths) {
  if (starts.empty() || starts.size() != ends.size() || window_count < 1 ||
      window_count > 3) {
    return false;
  }
  std::array<LazyControlWindowSweep, 3> sweeps;
  for (size_t w = 0; w < window_count; ++w) {
    if (!sweeps[w].Init(reference_length, window_bp[w], scale[w], baseline,
                        starts.size(), outputs[w])) {
      return false;
    }
  }
  for (size_t i = 0; i < starts.size(); ++i) {
    for (size_t w = 0; w < window_count; ++w) {
      sweeps[w].AddFragment(starts[i], ends[i]);
    }
  }
  bool ok = true;
  for (size_t w = 0; w < window_count; ++w) {
    ok = sweeps[w].Finish() && ok;
    if (max_heap_depths != nullptr) {
      max_heap_depths[w] = sweeps[w].max_heap_depth();
    }
  }
  return ok;
}

// MACS3.Signal.Pileup.over_two_pv_array(..., func="max"). It deliberately
// stops when either track ends and retains every input breakpoint; the latter
// is required for byte-identical tied-summit splitting.
void MaxOverTwoPvArrays(const LocalBdgTrack& a, const LocalBdgTrack& b,
                        LocalBdgTrack* out) {
  out->runs.clear();
  size_t ia = 0;
  size_t ib = 0;
  out->runs.reserve(a.runs.size() + b.runs.size());
  while (ia < a.runs.size() && ib < b.runs.size()) {
    const float value = std::max(a.runs[ia].value, b.runs[ib].value);
    if (a.runs[ia].end < b.runs[ib].end) {
      out->runs.push_back(LocalBdgRun(a.runs[ia++].end, value));
    } else if (a.runs[ia].end > b.runs[ib].end) {
      out->runs.push_back(LocalBdgRun(b.runs[ib++].end, value));
    } else {
      out->runs.push_back(LocalBdgRun(a.runs[ia].end, value));
      ++ia;
      ++ib;
    }
  }
}

void MaxOverThreePvArrays(const LocalBdgTrack& a, const LocalBdgTrack& b,
                          const LocalBdgTrack& c, LocalBdgTrack* out) {
  out->runs.clear();
  size_t ia = 0, ib = 0, ic = 0;
  out->runs.reserve(a.runs.size() + b.runs.size() + c.runs.size());
  while (ia < a.runs.size() && ib < b.runs.size() && ic < c.runs.size()) {
    const float value = std::max(
        a.runs[ia].value,
        std::max(b.runs[ib].value, c.runs[ic].value));
    const int32_t next = std::min(
        a.runs[ia].end, std::min(b.runs[ib].end, c.runs[ic].end));
    out->runs.push_back(LocalBdgRun(next, value));
    if (a.runs[ia].end == next) ++ia;
    if (b.runs[ib].end == next) ++ib;
    if (c.runs[ic].end == next) ++ic;
  }
}

// Controlled calling keeps the large, short-lived control tracks as a single
// array of 8-byte runs. Convert only the final max-control track to the public
// split-vector lambda representation, while building the summit breakpoint
// union in the same coordinate walk.
void BuildLambdaTrackFromTreatAndLocalControl(
    const std::vector<int32_t>& t_p, const LocalBdgTrack& ctrl,
    Macs3BdgTrack* out, std::vector<int32_t>* merged_positions) {
  out->p.clear();
  out->v.clear();
  if (merged_positions != nullptr) {
    merged_positions->clear();
    merged_positions->reserve(t_p.size() + ctrl.runs.size());
  }
  if (t_p.empty() || ctrl.runs.empty()) {
    return;
  }

  size_t it = 0;
  size_t ic = 0;
  bool have_prev = false;
  int32_t prev_p = 0;
  float prev_v = 0.f;
  out->p.reserve(t_p.size() + ctrl.runs.size());
  out->v.reserve(t_p.size() + ctrl.runs.size());
  while (it < t_p.size() && ic < ctrl.runs.size()) {
    int32_t cur_p = 0;
    const float cur_v = ctrl.runs[ic].value;
    if (t_p[it] < ctrl.runs[ic].end) {
      cur_p = t_p[it++];
    } else if (t_p[it] > ctrl.runs[ic].end) {
      cur_p = ctrl.runs[ic++].end;
    } else {
      cur_p = t_p[it++];
      ++ic;
    }
    if (merged_positions != nullptr &&
        (merged_positions->empty() || merged_positions->back() != cur_p)) {
      merged_positions->push_back(cur_p);
    }
    if (!have_prev) {
      prev_v = cur_v;
      prev_p = cur_p;
      have_prev = true;
      continue;
    }
    if (prev_v != cur_v) {
      out->p.push_back(prev_p);
      out->v.push_back(prev_v);
      prev_v = cur_v;
    }
    prev_p = cur_p;
  }
  if (have_prev) {
    out->p.push_back(prev_p);
    out->v.push_back(prev_v);
  }
  if (merged_positions != nullptr) {
    while (it < t_p.size()) {
      const int32_t p = t_p[it++];
      if (merged_positions->empty() || merged_positions->back() != p) {
        merged_positions->push_back(p);
      }
    }
    while (ic < ctrl.runs.size()) {
      const int32_t p = ctrl.runs[ic++].end;
      if (merged_positions->empty() || merged_positions->back() != p) {
        merged_positions->push_back(p);
      }
    }
  }
}

void ScoreRunPush(Macs3BdgTrack* out, int32_t seg_lo, int32_t seg_hi,
                  double val, bool quantize_values, bool* has, int32_t* lo,
                  int32_t* hi, double* last) {
  if (seg_hi <= seg_lo) {
    return;
  }
  if (!*has) {
    *has = true;
    *lo = seg_lo;
    *hi = seg_hi;
    *last = val;
    return;
  }
  const bool same_value = quantize_values
                              ? std::fabs(val - *last) <= 1e-5
                              : static_cast<float>(val) == static_cast<float>(*last);
  if (same_value && seg_lo == *hi) {
    *hi = seg_hi;
    return;
  }
  out->p.push_back(*hi);
  out->v.push_back(quantize_values ? Quantize5Float(*last)
                                   : static_cast<float>(*last));
  *lo = seg_lo;
  *hi = seg_hi;
  *last = val;
}

void ScoreRunFlush(Macs3BdgTrack* out, bool* has, int32_t* lo, int32_t* hi,
                   double* last, bool quantize_values) {
  (void)lo;
  if (!*has) {
    return;
  }
  out->p.push_back(*hi);
  out->v.push_back(quantize_values ? Quantize5Float(*last)
                                   : static_cast<float>(*last));
  *has = false;
}

uint32_t FloatBits(float v) {
  uint32_t out = 0;
  std::memcpy(&out, &v, sizeof(out));
  return out;
}

KHASH_MAP_INIT_INT64(ppois_score_lut, double)

class KhashPpoisScoreLut {
 public:
  explicit KhashPpoisScoreLut(size_t expected_buckets)
      : table_(kh_init(ppois_score_lut)) {
    if (table_ != nullptr && expected_buckets > 0) {
      const size_t capped = std::min<size_t>(
          expected_buckets, std::numeric_limits<khint_t>::max());
      (void)kh_resize(ppois_score_lut, table_, static_cast<khint_t>(capped));
    }
  }

  ~KhashPpoisScoreLut() {
    if (table_ != nullptr) kh_destroy(ppois_score_lut, table_);
  }

  bool Find(uint64_t key, double* value) const {
    if (table_ == nullptr) return false;
    const khint_t it = kh_get(ppois_score_lut, table_, key);
    if (it == kh_end(table_)) return false;
    *value = kh_value(table_, it);
    return true;
  }

  void Insert(uint64_t key, double value) {
    if (table_ == nullptr) return;
    int inserted = 0;
    const khint_t it = kh_put(ppois_score_lut, table_, key, &inserted);
    if (inserted >= 0 && it != kh_end(table_)) kh_value(table_, it) = value;
  }

 private:
  khash_t(ppois_score_lut)* table_ = nullptr;
};

struct PpoisCacheStats {
  size_t segments = 0;
  size_t hits = 0;
  size_t misses = 0;
};

template <bool kCollectStats, typename ScoreCache>
bool BuildMacs3FragPpoisTrackCached(
    const Macs3BdgTrack& treat, const Macs3BdgTrack& lambda,
    double pseudocount, Macs3BdgTrack* out, bool quantize_values,
    ScoreCache* score_cache, PpoisCacheStats* stats) {
  if (out == nullptr || score_cache == nullptr ||
      treat.p.size() != treat.v.size() || lambda.p.size() != lambda.v.size()) {
    return false;
  }
  out->p.clear();
  out->v.clear();
  if (treat.p.empty() || lambda.p.empty()) return true;
  const float pseudocount_f = static_cast<float>(pseudocount);
  size_t i1 = 0;
  size_t i2 = 0;
  int32_t prev = 0;
  bool has = false;
  int32_t lo = 0;
  int32_t hi = 0;
  double last = 0.0;
  while (i1 < treat.p.size() && i2 < lambda.p.size()) {
    const int32_t p1 = treat.p[i1];
    const int32_t p2 = lambda.p[i2];
    const int32_t nx = std::min(p1, p2);
    if (nx > prev) {
      const int32_t seg_lo = std::max(prev, static_cast<int32_t>(0));
      const int32_t seg_hi = std::max(nx, static_cast<int32_t>(0));
      if (seg_hi > seg_lo) {
        const uint64_t key =
            (static_cast<uint64_t>(FloatBits(treat.v[i1])) << 32) |
            static_cast<uint64_t>(FloatBits(lambda.v[i2]));
        double sp = 0.0;
        if (!score_cache->Find(key, &sp)) {
          const double raw = Macs3BdgcmpPpoisNegLog10P(
              treat.v[i1], lambda.v[i2], static_cast<double>(pseudocount_f));
          sp = quantize_values ? Round5Score(raw)
                               : static_cast<double>(static_cast<float>(raw));
          score_cache->Insert(key, sp);
          if (kCollectStats) ++stats->misses;
        } else if (kCollectStats) {
          ++stats->hits;
        }
        if (kCollectStats) ++stats->segments;
        ScoreRunPush(out, seg_lo, seg_hi, sp, quantize_values, &has, &lo, &hi,
                     &last);
      }
    }
    if (p1 <= p2) ++i1;
    if (p2 <= p1) ++i2;
    prev = nx;
  }
  ScoreRunFlush(out, &has, &lo, &hi, &last, quantize_values);
  return true;
}

}  // namespace

void BuildLambdaTrackFromTreatAndControl(const std::vector<int32_t>& t_p,
                                         const Macs3BdgTrack& ctrl,
                                         Macs3BdgTrack* out,
                                         bool quantize_values,
                                         std::vector<int32_t>* merged_positions) {
  out->p.clear();
  out->v.clear();
  if (merged_positions != nullptr) {
    merged_positions->clear();
    merged_positions->reserve(t_p.size() + ctrl.p.size());
  }
  if (t_p.empty() || ctrl.p.empty() || ctrl.p.size() != ctrl.v.size()) {
    return;
  }
  size_t it = 0;
  size_t ic = 0;
  bool have_prev = false;
  int32_t prev_p = 0;
  float pre_v = 0.f;
  out->p.reserve(t_p.size() + ctrl.p.size());
  out->v.reserve(t_p.size() + ctrl.p.size());
  while (it < t_p.size() && ic < ctrl.p.size()) {
    int32_t cur_p = 0;
    const float cur_v = ctrl.v[ic];
    if (t_p[it] < ctrl.p[ic]) {
      cur_p = t_p[it];
      ++it;
    } else if (t_p[it] > ctrl.p[ic]) {
      cur_p = ctrl.p[ic];
      ++ic;
    } else {
      cur_p = t_p[it];
      ++it;
      ++ic;
    }
    if (merged_positions != nullptr &&
        (merged_positions->empty() || merged_positions->back() != cur_p)) {
      merged_positions->push_back(cur_p);
    }
    if (!have_prev) {
      pre_v = cur_v;
      prev_p = cur_p;
      have_prev = true;
      continue;
    }
    const bool value_changed = quantize_values
                                   ? std::fabs(static_cast<double>(pre_v) -
                                               static_cast<double>(cur_v)) > 1e-5
                                   : pre_v != cur_v;
    if (value_changed) {
      out->p.push_back(prev_p);
      out->v.push_back(quantize_values ? Quantize5Float(pre_v) : pre_v);
      pre_v = cur_v;
    }
    prev_p = cur_p;
  }
  if (have_prev) {
    out->p.push_back(prev_p);
    out->v.push_back(quantize_values ? Quantize5Float(pre_v) : pre_v);
  }
  if (merged_positions != nullptr) {
    while (it < t_p.size()) {
      const int32_t p = t_p[it++];
      if (merged_positions->empty() || merged_positions->back() != p) {
        merged_positions->push_back(p);
      }
    }
    while (ic < ctrl.p.size()) {
      const int32_t p = ctrl.p[ic++];
      if (merged_positions->empty() || merged_positions->back() != p) {
        merged_positions->push_back(p);
      }
    }
  }
}

bool BuildMacs3FragWorkspaceFromFragments(
    const std::vector<ChromFragments>& by_chrom,
    const Macs3FragWorkspaceParams& params,
    Macs3FragPeakWorkspace* workspace) {
  if (workspace == nullptr || params.frag_pileup_max_count < 0 ||
      params.effective_genome_size < 1 || params.llocal_bp < 1) {
    return false;
  }
  std::vector<std::string> chrom_names;
  chrom_names.reserve(by_chrom.size());
  for (const ChromFragments& ch : by_chrom) {
    chrom_names.push_back(ch.name);
  }
  if (!InitMacs3FragWorkspace(chrom_names, params, workspace)) {
    return false;
  }
  for (size_t i = 0; i < by_chrom.size(); ++i) {
    const ChromFragments& ch = by_chrom[i];
    workspace->events_by_chrom[i].treat_span_events.reserve(ch.frags.size() * 2);
    if (params.build_treatment_local_lambda_events) {
      workspace->events_by_chrom[i].lambda_end_window_events.reserve(
          ch.frags.size() * 4);
    }
    if (params.byte_parity_summits &&
        params.capture_raw_pos_array_endpoints) {
      workspace->pos_array_starts_capture[i].reserve(ch.frags.size());
      workspace->pos_array_ends_capture[i].reserve(ch.frags.size());
    }
    for (const Fragment& f : ch.frags) {
      if (!AddMacs3FragWorkspaceFragment(
              workspace, static_cast<uint32_t>(i), f.start, f.end, f.count)) {
        return false;
      }
    }
  }
  // Keep captured endpoints until no-control pos-array finalization. The
  // controlled path disables capture and reuses finalized track breakpoints.
  return true;
}

bool InitMacs3FragWorkspace(const std::vector<std::string>& chrom_names,
                            const Macs3FragWorkspaceParams& params,
                            Macs3FragPeakWorkspace* workspace) {
  if (workspace == nullptr || params.frag_pileup_max_count < 0 ||
      params.effective_genome_size < 1 || params.llocal_bp < 1) {
    return false;
  }
  workspace->params = params;
  workspace->chrom_names = chrom_names;
  workspace->events_by_chrom.clear();
  workspace->treat_tracks.clear();
  workspace->lambda_tracks.clear();
  workspace->ppois_tracks.clear();
  workspace->pos_array_positions.clear();
  workspace->pos_array_starts_capture.clear();
  workspace->pos_array_ends_capture.clear();
  workspace->total_fragment_bp = 0;
  workspace->total_fragment_count = 0;
  workspace->events_by_chrom.resize(chrom_names.size());
  if (params.byte_parity_summits &&
      params.capture_raw_pos_array_endpoints) {
    workspace->pos_array_starts_capture.resize(chrom_names.size());
    workspace->pos_array_ends_capture.resize(chrom_names.size());
  }
  return true;
}

bool EnsureMacs3FragWorkspaceChromCapacity(Macs3FragPeakWorkspace* workspace,
                                           uint32_t chrom_id) {
  if (workspace == nullptr) return false;
  if (chrom_id < workspace->events_by_chrom.size()) return true;
  workspace->events_by_chrom.resize(chrom_id + 1);
  if (workspace->params.byte_parity_summits &&
      workspace->params.capture_raw_pos_array_endpoints) {
    workspace->pos_array_starts_capture.resize(chrom_id + 1);
    workspace->pos_array_ends_capture.resize(chrom_id + 1);
  }
  return true;
}

bool AddMacs3FragWorkspaceFragment(Macs3FragPeakWorkspace* workspace,
                                   uint32_t chrom_id, int32_t start,
                                   int32_t end, int32_t count) {
  if (workspace == nullptr || chrom_id >= workspace->events_by_chrom.size() ||
      end <= start || count <= 0) {
    return false;
  }
  const Macs3FragWorkspaceParams& params = workspace->params;
  const int32_t w = EffectiveFragWeight(count, params.frag_pileup_max_count,
                                        params.macs3_uint8_counts);
  if (w == 0) {
    return true;
  }
  Macs3FragEventBuffer& ev = workspace->events_by_chrom[chrom_id];
  ev.treat_span_events.push_back(Macs3FragEvent{start, w});
  ev.treat_span_events.push_back(Macs3FragEvent{end, -w});

  if (params.build_treatment_local_lambda_events) {
    const int32_t half = params.llocal_bp / 2;
    const int32_t ll = start - half;
    const int32_t rl = ll + params.llocal_bp;
    const int32_t lr = end - half;
    const int32_t rr = lr + params.llocal_bp;
    ev.lambda_end_window_events.push_back(Macs3FragEvent{ll, w});
    ev.lambda_end_window_events.push_back(Macs3FragEvent{rl, -w});
    ev.lambda_end_window_events.push_back(Macs3FragEvent{lr, w});
    ev.lambda_end_window_events.push_back(Macs3FragEvent{rr, -w});
  }

  // Byte-parity summits: capture per-chrom count-duplicated endpoints
  // (raw count, matching MACS3 FragParser instance count).
  if (params.byte_parity_summits &&
      params.capture_raw_pos_array_endpoints) {
    if (chrom_id >= workspace->pos_array_starts_capture.size()) {
      workspace->pos_array_starts_capture.resize(chrom_id + 1);
      workspace->pos_array_ends_capture.resize(chrom_id + 1);
    }
    std::vector<int32_t>& fs = workspace->pos_array_starts_capture[chrom_id];
    std::vector<int32_t>& fe = workspace->pos_array_ends_capture[chrom_id];
    for (int32_t k = 0; k < count; ++k) {
      fs.push_back(start);
      fe.push_back(end);
    }
  }

  const int32_t c = RowCountForMacsTotals(count, params.frag_pileup_max_count,
                                          params.macs3_uint8_counts);
  if (c > 0) {
    workspace->total_fragment_count += c;
    workspace->total_fragment_bp += static_cast<int64_t>(end - start) * c;
  }
  return true;
}

bool BuildMacs3FragTracksFromSortedFragments(
    std::vector<ChromFragments>* by_chrom,
    const Macs3FragWorkspaceParams& params,
    Macs3FragPeakWorkspace* workspace, int num_threads,
    bool release_fragments) {
  if (by_chrom == nullptr || workspace == nullptr || by_chrom->empty() ||
      params.frag_pileup_max_count < 0 || params.effective_genome_size < 1 ||
      params.llocal_bp < 1) {
    return false;
  }
  std::vector<std::string> chrom_names;
  chrom_names.reserve(by_chrom->size());
  for (const ChromFragments& ch : *by_chrom) chrom_names.push_back(ch.name);
  if (!InitMacs3FragWorkspace(chrom_names, params, workspace)) return false;

  const size_t n = by_chrom->size();
  workspace->treat_tracks.resize(n);
  workspace->lambda_tracks.resize(n);
  std::vector<Macs3BdgTrack> ctrl_raw_tracks;
  if (params.build_treatment_local_lambda_events) {
    ctrl_raw_tracks.resize(n);
  }
  std::vector<int64_t> fragment_counts(n, 0);
  std::vector<int64_t> fragment_bp(n, 0);
  std::vector<char> ok(n, 1);
  const int threads = NormalizeThreads(num_threads, n);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) if(threads > 1)
  for (int64_t ii = 0; ii < static_cast<int64_t>(n); ++ii) {
    const size_t i = static_cast<size_t>(ii);
    ChromFragments& ch = by_chrom->at(i);
    int32_t previous_start = -1;
    size_t raw_endpoint_count = 0;
    for (const Fragment& f : ch.frags) {
      if (f.start < previous_start || f.end <= f.start || f.count <= 0) {
        ok[i] = 0;
        break;
      }
      previous_start = f.start;
      const int32_t w = EffectiveFragWeight(
          f.count, params.frag_pileup_max_count, params.macs3_uint8_counts);
      if (w == 0) continue;
      if (params.byte_parity_summits &&
          params.capture_raw_pos_array_endpoints) {
        const size_t c = static_cast<size_t>(f.count);
        if (raw_endpoint_count > std::numeric_limits<size_t>::max() - c) {
          ok[i] = 0;
          break;
        }
        raw_endpoint_count += c;
      }
      const int32_t row_count = RowCountForMacsTotals(
          f.count, params.frag_pileup_max_count, params.macs3_uint8_counts);
      if (row_count > 0) {
        fragment_counts[i] += row_count;
        fragment_bp[i] += static_cast<int64_t>(f.end - f.start) * row_count;
      }
    }
    if (!ok[i]) continue;

    std::vector<WeightedEndpoint> starts;
    std::vector<WeightedEndpoint> ends;
    BuildSortedControlEndpoints(ch, params.frag_pileup_max_count,
                                params.macs3_uint8_counts, &starts, &ends);
    if (!starts.empty()) {
      TrackFromSortedSpanEndpoints(starts, ends, &workspace->treat_tracks[i]);
      if (params.build_treatment_local_lambda_events) {
        TrackFromSortedEndpointWindowsUnclamped(
            starts, ends, params.llocal_bp, &ctrl_raw_tracks[i]);
      }
    }

    if (params.byte_parity_summits &&
        params.capture_raw_pos_array_endpoints) {
      std::vector<int32_t>& fs = workspace->pos_array_starts_capture[i];
      std::vector<int32_t>& fe = workspace->pos_array_ends_capture[i];
      fs.reserve(raw_endpoint_count);
      fe.reserve(raw_endpoint_count);
      for (const Fragment& f : ch.frags) {
        if (EffectiveFragWeight(f.count, params.frag_pileup_max_count,
                                params.macs3_uint8_counts) == 0) {
          continue;
        }
        fs.insert(fs.end(), static_cast<size_t>(f.count), f.start);
        fe.insert(fe.end(), static_cast<size_t>(f.count), f.end);
      }
    }
    std::vector<WeightedEndpoint>().swap(starts);
    std::vector<WeightedEndpoint>().swap(ends);
    if (release_fragments) std::vector<Fragment>().swap(ch.frags);
  }
  if (std::find(ok.begin(), ok.end(), 0) != ok.end()) return false;

  for (size_t i = 0; i < n; ++i) {
    workspace->total_fragment_count += fragment_counts[i];
    workspace->total_fragment_bp += fragment_bp[i];
  }
  if (workspace->total_fragment_count <= 0 || workspace->total_fragment_bp <= 0) {
    return false;
  }
  // --nolambda callers finalize a single genome-wide background after the
  // treatment tracks are available. Avoid both endpoint-window construction
  // above and its scaling/merge pass below.
  if (!params.build_treatment_local_lambda_events) {
    return true;
  }
  // PeakDetect.__call_peaks_wo_control assigns both lambda_bg and the local
  // scale factor to C floats before se_all_in_one_pileup multiplies them by
  // the integer endpoint depth.  Casting only after a double multiplication
  // changes the last printed decimal on some runs and, in turn, their p/q
  // score lookup keys.
  const float lambda_bg = static_cast<float>(
      static_cast<double>(workspace->total_fragment_bp) /
      static_cast<double>(params.effective_genome_size));
  const float scale = static_cast<float>(
      static_cast<double>(workspace->total_fragment_bp) /
      (static_cast<double>(params.llocal_bp) *
       static_cast<double>(workspace->total_fragment_count) * 2.0));
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) if(threads > 1)
  for (int64_t ii = 0; ii < static_cast<int64_t>(n); ++ii) {
    const size_t i = static_cast<size_t>(ii);
    Macs3BdgTrack& ctrl = ctrl_raw_tracks[i];
    for (float& v : ctrl.v) {
      float z = v * scale;
      if (z < lambda_bg) z = lambda_bg;
      v = z;
    }
    BuildLambdaTrackFromTreatAndControl(
        workspace->treat_tracks[i].p, ctrl, &workspace->lambda_tracks[i],
        false, nullptr);
    std::vector<int32_t>().swap(ctrl.p);
    std::vector<float>().swap(ctrl.v);
  }
  return true;
}

bool FinalizeMacs3FragTreatTracks(Macs3FragPeakWorkspace* workspace,
                                   int num_threads) {
  if (workspace == nullptr) {
    return false;
  }
  // Skip if tracks were pre-populated by a different builder (e.g. the
  // sweep-line workspace) AND no events are buffered. Detect by: any chrom
  // has an empty event buffer AND the corresponding treat_track is
  // non-empty. Conservative: only short-circuit when treat_tracks already
  // matches events_by_chrom in size and at least one track is non-empty.
  bool any_event = false;
  for (const auto& ev : workspace->events_by_chrom) {
    if (!ev.treat_span_events.empty()) { any_event = true; break; }
  }
  bool any_track = false;
  for (const auto& tr : workspace->treat_tracks) {
    if (!tr.p.empty()) { any_track = true; break; }
  }
  if (!any_event && any_track &&
      workspace->treat_tracks.size() == workspace->events_by_chrom.size()) {
    return true;  // already pre-finalized; nothing to do
  }
  workspace->treat_tracks.clear();
  workspace->treat_tracks.resize(workspace->events_by_chrom.size());
  const size_t n = workspace->events_by_chrom.size();
  std::vector<char> ok(n, 1);
  const int threads = NormalizeThreads(num_threads, n);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) if(threads > 1)
  for (int64_t ii = 0; ii < static_cast<int64_t>(n); ++ii) {
    const size_t i = static_cast<size_t>(ii);
    std::vector<Macs3FragEvent>& ev = workspace->events_by_chrom[i].treat_span_events;
    if (!CoalesceEventPositions(&ev)) {
      ok[i] = 0;
      continue;
    }
    TrackFromEvents(ev, &workspace->treat_tracks[i]);
    ReleaseEvents(&ev);
  }
  return std::find(ok.begin(), ok.end(), 0) == ok.end();
}

bool FinalizeMacs3FragLambdaTracks(Macs3FragPeakWorkspace* workspace,
                                    int num_threads) {
  if (workspace == nullptr || workspace->total_fragment_count <= 0 ||
      workspace->total_fragment_bp <= 0 ||
      workspace->treat_tracks.size() != workspace->events_by_chrom.size()) {
    return false;
  }
  // Skip if lambda_tracks were pre-populated and no lambda-window events are
  // buffered (sweep-line workspace path).
  bool any_lambda_event = false;
  for (const auto& ev : workspace->events_by_chrom) {
    if (!ev.lambda_end_window_events.empty()) { any_lambda_event = true; break; }
  }
  bool any_lambda_track = false;
  for (const auto& tr : workspace->lambda_tracks) {
    if (!tr.p.empty()) { any_lambda_track = true; break; }
  }
  if (!any_lambda_event && any_lambda_track &&
      workspace->lambda_tracks.size() == workspace->events_by_chrom.size()) {
    return true;  // pre-finalized
  }
  const float lambda_bg = static_cast<float>(
      static_cast<double>(workspace->total_fragment_bp) /
      static_cast<double>(workspace->params.effective_genome_size));
  const float scale = static_cast<float>(
      static_cast<double>(workspace->total_fragment_bp) /
      (static_cast<double>(workspace->params.llocal_bp) *
       static_cast<double>(workspace->total_fragment_count) * 2.0));

  workspace->lambda_tracks.clear();
  workspace->lambda_tracks.resize(workspace->events_by_chrom.size());
  const size_t n = workspace->events_by_chrom.size();
  std::vector<char> ok(n, 1);
  const int threads = NormalizeThreads(num_threads, n);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) if(threads > 1)
  for (int64_t ii = 0; ii < static_cast<int64_t>(n); ++ii) {
    const size_t i = static_cast<size_t>(ii);
    std::vector<Macs3FragEvent>& ev = workspace->events_by_chrom[i].lambda_end_window_events;
    if (!CoalesceEventPositions(&ev)) {
      ok[i] = 0;
      continue;
    }
    Macs3BdgTrack ctrl;
    TrackFromEvents(ev, &ctrl);
    ReleaseEvents(&ev);
    for (float& v : ctrl.v) {
      float z = v * scale;
      if (z < lambda_bg) {
        z = lambda_bg;
      }
      v = z;
    }
    BuildLambdaTrackFromTreatAndControl(workspace->treat_tracks[i].p, ctrl,
                                        &workspace->lambda_tracks[i], false,
                                        nullptr);
  }
  return std::find(ok.begin(), ok.end(), 0) == ok.end();
}

bool FinalizeMacs3FragGlobalLambdaTracks(Macs3FragPeakWorkspace* workspace,
                                         int num_threads) {
  if (workspace == nullptr || workspace->total_fragment_count <= 0 ||
      workspace->total_fragment_bp <= 0 ||
      workspace->treat_tracks.size() != workspace->chrom_names.size()) {
    return false;
  }
  const float lambda_bg = static_cast<float>(
      static_cast<double>(workspace->total_fragment_bp) /
      static_cast<double>(workspace->params.effective_genome_size));
  const size_t n = workspace->treat_tracks.size();
  workspace->lambda_tracks.assign(n, Macs3BdgTrack());
  if (workspace->params.byte_parity_summits) {
    workspace->pos_array_positions.assign(n, std::vector<int32_t>());
  } else {
    workspace->pos_array_positions.clear();
  }
  const int threads = NormalizeThreads(num_threads, n);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) if(threads > 1)
  for (int64_t ii = 0; ii < static_cast<int64_t>(n); ++ii) {
    const size_t i = static_cast<size_t>(ii);
    const Macs3BdgTrack& treat = workspace->treat_tracks[i];
    if (treat.p.empty()) continue;
    workspace->lambda_tracks[i].p.push_back(treat.p.back());
    workspace->lambda_tracks[i].v.push_back(lambda_bg);
    // CallerFromAlignments constructs ctrl_pv at only treat_pv's final
    // endpoint, so the merged position array is exactly treat_pv's positions.
    if (workspace->params.byte_parity_summits) {
      workspace->pos_array_positions[i] = treat.p;
    }
  }
  for (Macs3FragEventBuffer& ev : workspace->events_by_chrom) {
    ReleaseEvents(&ev.lambda_end_window_events);
  }
  std::vector<std::vector<int32_t>>().swap(
      workspace->pos_array_starts_capture);
  std::vector<std::vector<int32_t>>().swap(
      workspace->pos_array_ends_capture);
  return true;
}

bool FinalizeMacs3FragLambdaTracksWithControl(
    Macs3FragPeakWorkspace* workspace,
    const std::vector<ChromFragments>& control_by_chrom,
    const Macs3FragControlParams& control_params, int num_threads) {
  if (workspace == nullptr || workspace->total_fragment_count <= 0 ||
      workspace->total_fragment_bp <= 0 ||
      workspace->treat_tracks.size() != workspace->chrom_names.size() ||
      control_by_chrom.empty() || control_params.slocal_bp < 0) {
    return false;
  }

  int64_t control_total = 0;
  for (const ChromFragments& ch : control_by_chrom) {
    for (const Fragment& f : ch.frags) {
      control_total += RowCountForMacsTotals(
          f.count, workspace->params.frag_pileup_max_count,
          workspace->params.macs3_uint8_counts);
    }
  }
  if (control_total <= 0) return false;

  const float retained_mean = static_cast<float>(
      static_cast<double>(workspace->total_fragment_bp) /
      static_cast<double>(workspace->total_fragment_count));
  const float parser_mean = static_cast<float>(
      control_params.treatment_fragment_length_bp > 0.0
          ? control_params.treatment_fragment_length_bp
          : static_cast<double>(retained_mean));
  const int32_t d_window = static_cast<int32_t>(parser_mean);
  if (d_window < 1 || d_window > workspace->params.llocal_bp ||
      (control_params.slocal_bp > 0 &&
       (d_window > control_params.slocal_bp ||
        control_params.slocal_bp > workspace->params.llocal_bp))) {
    return false;
  }

  // These casts intentionally mirror PeakDetect.py's PE control arithmetic:
  // PETrackI.finalize computes the parser mean by dividing float32 length by
  // float32 total. filter_dup(-1) returns before recomputing it, while a
  // numeric cap recomputes the retained mean. In the `all` branch Cython
  // promotes the long fragment count for multiplication and then assigns
  // the product to a C long. Keep the established numeric-cap path intact.
  const float control_sum_mean = control_params.keep_duplicates_all
      ? static_cast<float>(
            static_cast<float>(workspace->total_fragment_bp) /
            static_cast<float>(workspace->total_fragment_count))
      : retained_mean;
  const int64_t doubled_control_total = control_total * 2;
  const int64_t control_sum = control_params.keep_duplicates_all
      ? static_cast<int64_t>(static_cast<double>(doubled_control_total) *
                             static_cast<double>(control_sum_mean))
      : static_cast<int64_t>(
            static_cast<float>(doubled_control_total) * control_sum_mean);
  if (control_sum <= 0) return false;
  const double ratio = static_cast<double>(workspace->total_fragment_bp) /
                       static_cast<double>(control_sum);

  bool scale_treatment_to_control = false;
  if (control_params.scale_to_large) {
    scale_treatment_to_control =
        workspace->total_fragment_count <= doubled_control_total;
  } else {
    scale_treatment_to_control =
        workspace->total_fragment_count > doubled_control_total;
  }
  const float treat_scale = scale_treatment_to_control
                                ? static_cast<float>(1.0 / ratio)
                                : 1.0f;
  const float ctrl_depth_scale = scale_treatment_to_control
                                     ? 1.0f
                                     : static_cast<float>(ratio);
  const float lambda_bg = static_cast<float>(
      static_cast<double>(scale_treatment_to_control
                              ? control_sum
                              : workspace->total_fragment_bp) /
      static_cast<double>(workspace->params.effective_genome_size));

  for (Macs3BdgTrack& tr : workspace->treat_tracks) {
    for (float& v : tr.v) v = static_cast<float>(v * treat_scale);
  }

  std::map<std::string, const ChromFragments*> control_index;
  for (const ChromFragments& ch : control_by_chrom) {
    control_index[ch.name] = &ch;
  }
  workspace->lambda_tracks.assign(workspace->chrom_names.size(),
                                  Macs3BdgTrack());
  if (workspace->params.byte_parity_summits) {
    workspace->pos_array_positions.assign(workspace->chrom_names.size(),
                                          std::vector<int32_t>());
  } else {
    workspace->pos_array_positions.clear();
  }
  const size_t n = workspace->chrom_names.size();
  std::vector<char> ok(n, 1);
  std::vector<char> common(n, 0);
  if (control_params.no_lambda) {
    const int threads = NormalizeThreads(num_threads, n);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) if(threads > 1)
    for (int64_t ii = 0; ii < static_cast<int64_t>(n); ++ii) {
      const size_t i = static_cast<size_t>(ii);
      const auto found = control_index.find(workspace->chrom_names[i]);
      if (found == control_index.end() || found->second->frags.empty() ||
          workspace->treat_tracks[i].p.empty()) {
        // MACS3 calls only chromosomes shared by treatment and control.
        workspace->treat_tracks[i] = Macs3BdgTrack();
        continue;
      }
      workspace->lambda_tracks[i].p.push_back(
          workspace->treat_tracks[i].p.back());
      workspace->lambda_tracks[i].v.push_back(lambda_bg);
      if (workspace->params.byte_parity_summits) {
        workspace->pos_array_positions[i] = workspace->treat_tracks[i].p;
      }
    }
    return true;
  }
  const bool lazy_control_windows =
      EnvValueIsOne("RAPIDMACS_LAZY_CONTROL_WINDOWS");
  const bool profile_control_windows =
      EnvValueIsOne("RAPIDMACS_CONTROL_WINDOWS_PROFILE");
  std::vector<double> end_sort_seconds;
  std::vector<double> window_track_seconds;
  std::vector<std::array<size_t, 3>> max_heap_depths;
  if (profile_control_windows) {
    end_sort_seconds.assign(n, 0.0);
    window_track_seconds.assign(n, 0.0);
    max_heap_depths.resize(n);
    for (size_t i = 0; i < n; ++i) max_heap_depths[i].fill(0);
  }
  const int threads = NormalizeThreads(num_threads, n);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) if(threads > 1)
  for (int64_t ii = 0; ii < static_cast<int64_t>(n); ++ii) {
    const size_t i = static_cast<size_t>(ii);
    const auto found = control_index.find(workspace->chrom_names[i]);
    if (found == control_index.end() || found->second->frags.empty() ||
        workspace->treat_tracks[i].p.empty()) {
      continue;
    }
    common[i] = 1;
    const ChromFragments& control_ch = *found->second;
    std::vector<WeightedEndpoint> control_starts;
    std::vector<WeightedEndpoint> control_ends;
    BuildSortedControlEndpoints(
        control_ch, workspace->params.frag_pileup_max_count,
        workspace->params.macs3_uint8_counts, &control_starts, &control_ends,
        !lazy_control_windows,
        profile_control_windows ? &end_sort_seconds[i] : nullptr);
    LocalBdgTrack dlocal;
    LocalBdgTrack slocal;
    LocalBdgTrack llocal;
    const bool have_slocal = control_params.slocal_bp > 0;
    const bool have_llocal =
        workspace->params.llocal_bp > control_params.slocal_bp;
    const float slocal_scale = have_slocal
        ? static_cast<float>(static_cast<double>(parser_mean) /
                             static_cast<double>(control_params.slocal_bp) *
                             ctrl_depth_scale)
        : 0.f;
    const float llocal_scale = have_llocal
        ? static_cast<float>(static_cast<double>(parser_mean) /
                             static_cast<double>(workspace->params.llocal_bp) *
                             ctrl_depth_scale)
        : 0.f;
    std::chrono::steady_clock::time_point windows_t0;
    if (profile_control_windows) {
      windows_t0 = std::chrono::steady_clock::now();
    }
    bool tracks_ok = true;
    if (lazy_control_windows) {
      int32_t windows[3];
      float scales[3];
      LocalBdgTrack* outputs[3];
      size_t heap_depth[3] = {0, 0, 0};
      size_t nw = 0;
      windows[nw] = d_window;
      scales[nw] = ctrl_depth_scale;
      outputs[nw++] = &dlocal;
      if (have_slocal) {
        windows[nw] = control_params.slocal_bp;
        scales[nw] = slocal_scale;
        outputs[nw++] = &slocal;
      }
      if (have_llocal) {
        windows[nw] = workspace->params.llocal_bp;
        scales[nw] = llocal_scale;
        outputs[nw++] = &llocal;
      }
      tracks_ok = TrackFromLazyEndpointWindowsShared(
          control_starts, control_ends, control_ch.reference_length, windows,
          scales, lambda_bg, outputs, nw, heap_depth);
      if (profile_control_windows) {
        size_t hw = 0;
        max_heap_depths[i][0] = heap_depth[hw++];
        if (have_slocal) max_heap_depths[i][1] = heap_depth[hw++];
        if (have_llocal) max_heap_depths[i][2] = heap_depth[hw++];
      }
    } else {
      tracks_ok = TrackFromSortedEndpointWindows(
          control_starts, control_ends, control_ch.reference_length,
          d_window, ctrl_depth_scale, lambda_bg, &dlocal);
      if (tracks_ok && have_slocal) {
        tracks_ok = TrackFromSortedEndpointWindows(
            control_starts, control_ends, control_ch.reference_length,
            control_params.slocal_bp, slocal_scale, lambda_bg, &slocal);
      }
      if (tracks_ok && have_llocal) {
        tracks_ok = TrackFromSortedEndpointWindows(
            control_starts, control_ends, control_ch.reference_length,
            workspace->params.llocal_bp, llocal_scale, lambda_bg, &llocal);
      }
    }
    if (profile_control_windows) {
      const auto windows_t1 = std::chrono::steady_clock::now();
      window_track_seconds[i] =
          std::chrono::duration<double>(windows_t1 - windows_t0).count();
    }
    if (!tracks_ok) {
      ok[i] = 0;
      continue;
    }

    LocalBdgTrack ctrl_max;
    if (have_slocal && have_llocal) {
      MaxOverThreePvArrays(dlocal, slocal, llocal, &ctrl_max);
    } else if (have_slocal) {
      MaxOverTwoPvArrays(dlocal, slocal, &ctrl_max);
    } else if (have_llocal) {
      MaxOverTwoPvArrays(dlocal, llocal, &ctrl_max);
    } else {
      ctrl_max = std::move(dlocal);
    }
    std::vector<WeightedEndpoint>().swap(control_starts);
    std::vector<WeightedEndpoint>().swap(control_ends);
    std::vector<LocalBdgRun>().swap(dlocal.runs);
    std::vector<LocalBdgRun>().swap(slocal.runs);
    std::vector<LocalBdgRun>().swap(llocal.runs);

    // TrackFromSortedEndpointWindows and MaxOver* preserve the exact MACS3
    // control breakpoint union. Build the final lambda and optional summit
    // position union in the same coordinate walk.
    BuildLambdaTrackFromTreatAndLocalControl(
        workspace->treat_tracks[i].p, ctrl_max,
        &workspace->lambda_tracks[i],
        workspace->params.byte_parity_summits
            ? &workspace->pos_array_positions[i]
            : nullptr);
  }
  if (profile_control_windows) {
    double sort_total = 0.0;
    double tracks_total = 0.0;
    for (size_t i = 0; i < n; ++i) {
      if (!common[i]) continue;
      sort_total += end_sort_seconds[i];
      tracks_total += window_track_seconds[i];
      std::fprintf(
          stderr,
          "[RapidMACS] control-window-profile variant=%s chrom=%s "
          "end_sort_sec=%.9f window_tracks_sec=%.9f heap_max_d=%zu "
          "heap_max_slocal=%zu heap_max_llocal=%zu\n",
          lazy_control_windows ? "lazy" : "default",
          workspace->chrom_names[i].c_str(), end_sort_seconds[i],
          window_track_seconds[i], max_heap_depths[i][0],
          max_heap_depths[i][1], max_heap_depths[i][2]);
    }
    std::fprintf(stderr,
                 "[RapidMACS] control-window-profile variant=%s chrom=ALL "
                 "end_sort_sec=%.9f window_tracks_sec=%.9f\n",
                 lazy_control_windows ? "lazy" : "default", sort_total,
                 tracks_total);
  }
  if (std::find(ok.begin(), ok.end(), 0) != ok.end() ||
      std::find(common.begin(), common.end(), 1) == common.end()) {
    return false;
  }

  std::vector<std::vector<int32_t>>().swap(
      workspace->pos_array_starts_capture);
  std::vector<std::vector<int32_t>>().swap(
      workspace->pos_array_ends_capture);

  for (Macs3FragEventBuffer& ev : workspace->events_by_chrom) {
    ReleaseEvents(&ev.lambda_end_window_events);
  }
  workspace->ppois_tracks.clear();
  return true;
}

bool FinalizeMacs3FragPpoisTracks(Macs3FragPeakWorkspace* workspace,
                                   int num_threads,
                                   bool quantize_values) {
  if (workspace == nullptr ||
      workspace->treat_tracks.size() != workspace->lambda_tracks.size()) {
    return false;
  }
  const size_t n = workspace->treat_tracks.size();
  workspace->ppois_tracks.clear();
  workspace->ppois_tracks.resize(n);
  std::vector<char> ok(n, 1);
  const int threads = NormalizeThreads(num_threads, n);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) if(threads > 1)
  for (int64_t ii = 0; ii < static_cast<int64_t>(n); ++ii) {
    const size_t i = static_cast<size_t>(ii);
    if (!BuildMacs3FragPpoisTrack(
            workspace->treat_tracks[i], workspace->lambda_tracks[i],
            workspace->params.score_pseudocount,
            &workspace->ppois_tracks[i], quantize_values)) {
      ok[i] = 0;
    }
  }
  return std::find(ok.begin(), ok.end(), 0) == ok.end();
}

bool FinalizeMacs3FragPosArrayPositions(Macs3FragPeakWorkspace* workspace,
                                        int num_threads) {
  if (workspace == nullptr) {
    return false;
  }
  if (!workspace->params.byte_parity_summits) {
    return true;
  }
  // No-op when there's nothing to consume: either the capture buffers
  // are empty (sweep workspace already populated pos_array_positions
  // directly during chrom transitions) or this function has already run.
  bool any_capture = false;
  for (const auto& v : workspace->pos_array_starts_capture) {
    if (!v.empty()) {
      any_capture = true;
      break;
    }
  }
  if (!any_capture) {
    // Drop the unused capture vector skeletons (if any) to keep memory
    // tight; leave pos_array_positions exactly as the upstream populated.
    std::vector<std::vector<int32_t>>().swap(workspace->pos_array_starts_capture);
    std::vector<std::vector<int32_t>>().swap(workspace->pos_array_ends_capture);
    return true;
  }
  const size_t n = std::max(workspace->events_by_chrom.size(),
                            workspace->pos_array_starts_capture.size());
  workspace->pos_array_positions.assign(n, std::vector<int32_t>());
  workspace->pos_array_starts_capture.resize(n);
  workspace->pos_array_ends_capture.resize(n);
  std::vector<char> ok(n, 1);
  const int threads = NormalizeThreads(num_threads, n);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) if(threads > 1)
  for (int64_t ii = 0; ii < static_cast<int64_t>(n); ++ii) {
    const size_t i = static_cast<size_t>(ii);
    BuildMacs3PosArrayPositions(
        std::move(workspace->pos_array_starts_capture[i]),
        std::move(workspace->pos_array_ends_capture[i]),
        workspace->params.llocal_bp,
        &workspace->pos_array_positions[i]);
    std::vector<int32_t>().swap(workspace->pos_array_starts_capture[i]);
    std::vector<int32_t>().swap(workspace->pos_array_ends_capture[i]);
  }
  std::vector<std::vector<int32_t>>().swap(workspace->pos_array_starts_capture);
  std::vector<std::vector<int32_t>>().swap(workspace->pos_array_ends_capture);
  (void)ok;
  return true;
}

bool BuildMacs3FragPpoisTrack(const Macs3BdgTrack& treat,
                              const Macs3BdgTrack& lambda,
                              double pseudocount, Macs3BdgTrack* out,
                              bool quantize_values) {
  const size_t merged_run_upper_bound = treat.p.size() + lambda.p.size();
  // Full ATAC data averages about one distinct (treat, lambda) pair per eight
  // merged runs. Reserve that many flat khash buckets up front. Inserts never
  // allocate individual nodes; chromosomes with unusually diverse values may
  // geometrically grow the backing arrays. The table is chromosome-local and
  // is released as soon as this call returns.
  const size_t expected_buckets =
      std::max<size_t>(16, merged_run_upper_bound / 8);
  KhashPpoisScoreLut score_cache(expected_buckets);
  PpoisCacheStats stats;
  const bool collect_stats =
      std::getenv("RAPIDMACS_PPOIS_CACHE_STATS") != nullptr ||
      std::getenv("LIBMACS3_PPOIS_CACHE_STATS") != nullptr;
  const bool ok = collect_stats
      ? BuildMacs3FragPpoisTrackCached<true>(
            treat, lambda, pseudocount, out, quantize_values, &score_cache,
            &stats)
      : BuildMacs3FragPpoisTrackCached<false>(
            treat, lambda, pseudocount, out, quantize_values, &score_cache,
            nullptr);
  if (collect_stats) {
    std::fprintf(stderr,
                 "[RapidMACS] ppois-cache segments=%zu hits=%zu misses=%zu "
                 "output_runs=%zu\n",
                 stats.segments, stats.hits, stats.misses,
                 out == nullptr ? 0 : out->p.size());
  }
  return ok;
}

bool WriteMacs3BdgTracksBedGraph(const std::string& path,
                                 const std::vector<std::string>& chrom_names,
                                 const std::vector<Macs3BdgTrack>& tracks) {
  if (path.empty() || chrom_names.size() != tracks.size()) {
    return false;
  }
  FILE* fp = std::fopen(path.c_str(), "w");
  if (fp == nullptr) {
    return false;
  }
  for (size_t i = 0; i < tracks.size(); ++i) {
    const Macs3BdgTrack& tr = tracks[i];
    int32_t start = 0;
    size_t j = 0;
    while (j < tr.p.size() && j < tr.v.size()) {
      char value[64];
      std::snprintf(value, sizeof(value), "%.5f",
                    static_cast<double>(tr.v[j]));
      int32_t end = tr.p[j];
      size_t next = j + 1;
      while (next < tr.p.size() && next < tr.v.size()) {
        char next_value[64];
        std::snprintf(next_value, sizeof(next_value), "%.5f",
                      static_cast<double>(tr.v[next]));
        if (std::strcmp(value, next_value) != 0) break;
        end = tr.p[next];
        ++next;
      }
      if (std::fprintf(fp, "%s\t%d\t%d\t%s\n", chrom_names[i].c_str(), start,
                       end, value) < 0) {
        std::fclose(fp);
        return false;
      }
      start = end;
      j = next;
    }
  }
  if (std::fclose(fp) != 0) {
    return false;
  }
  return true;
}

}  // namespace peaks
}  // namespace chromap
