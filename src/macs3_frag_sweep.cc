#include "rapidmacs/macs3_frag_sweep.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <utility>

#include "rapidmacs/macs3_pos_array.h"

namespace chromap {
namespace peaks {

namespace {

float Quantize5Float(double x) {
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%.5f", x);
  return std::strtof(buf, nullptr);
}

int32_t EffectiveFragWeight(int32_t raw_count, int32_t max_count,
                            bool macs3_uint8_counts) {
  int64_t w = raw_count;
  if (max_count > 0) w = std::min<int64_t>(w, max_count);
  if (w <= 0) return 0;
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
  if (max_count > 0) t = std::min<int64_t>(t, max_count);
  if (t <= 0) return 0;
  if (!macs3_uint8_counts) {
    if (t > std::numeric_limits<int32_t>::max()) {
      return std::numeric_limits<int32_t>::max();
    }
    return static_cast<int32_t>(t);
  }
  const auto u = static_cast<uint8_t>(static_cast<unsigned int>(t));
  return static_cast<int32_t>(u);
}

}  // namespace

Macs3FragSweepWorkspace::Macs3FragSweepWorkspace()
    : lambda_half_(0),
      cur_chrom_id_(-1),
      prev_start_(-1),
      total_fragment_count_(0),
      total_fragment_bp_(0),
      finalized_(false) {}

bool Macs3FragSweepWorkspace::Init(const Macs3FragWorkspaceParams& params) {
  if (params.frag_pileup_max_count < 0 || params.effective_genome_size < 1 ||
      params.llocal_bp < 1) {
    error_ = "invalid Macs3FragWorkspaceParams";
    return false;
  }
  params_ = params;
  lambda_half_ = params.llocal_bp / 2;
  cur_chrom_id_ = -1;
  prev_start_ = -1;
  total_fragment_count_ = 0;
  total_fragment_bp_ = 0;
  finalized_ = false;
  ResetStream(&treat_stream_);
  ResetStream(&ctrl_stream_);
  chrom_names_.clear();
  treat_tracks_.clear();
  ctrl_raw_tracks_.clear();
  lambda_tracks_.clear();
  frag_starts_.clear();
  frag_ends_.clear();
  pos_array_positions_.clear();
  error_.clear();
  return true;
}

void Macs3FragSweepWorkspace::ResetStream(EventStream* stream) {
  std::priority_queue<Event, std::vector<Event>, std::greater<Event>> empty;
  stream->heap.swap(empty);
  stream->s = 0;
  stream->z = 0.0;
  stream->pre_z = -1.0e4;
  stream->track.p.clear();
  stream->track.v.clear();
}

void Macs3FragSweepWorkspace::EmitClose(EventStream* stream, int32_t pos) {
  // Close the segment (stream->s, pos) with current value stream->z. Coalesce
  // with the previous segment if its emitted value matches stream->pre_z (so
  // adjacent same-value runs collapse to a single bedGraph row).
  if (stream->z == stream->pre_z && !stream->track.p.empty()) {
    stream->track.p.back() = pos;
  } else {
    stream->track.p.push_back(pos);
    stream->track.v.push_back(Quantize5Float(stream->z));
    stream->pre_z = stream->z;
  }
}

void Macs3FragSweepWorkspace::DrainBefore(EventStream* stream, int32_t target) {
  while (!stream->heap.empty() && stream->heap.top().pos < target) {
    const Event e = stream->heap.top();
    stream->heap.pop();
    if (e.pos != stream->s) {
      EmitClose(stream, e.pos);
    }
    stream->z += static_cast<double>(e.delta);
    stream->s = e.pos;
  }
}

void Macs3FragSweepWorkspace::DrainFully(EventStream* stream) {
  while (!stream->heap.empty()) {
    const Event e = stream->heap.top();
    stream->heap.pop();
    if (e.pos != stream->s) {
      EmitClose(stream, e.pos);
    }
    stream->z += static_cast<double>(e.delta);
    stream->s = e.pos;
  }
}

bool Macs3FragSweepWorkspace::TransitionToChrom(uint32_t new_chrom_id) {
  if (cur_chrom_id_ >= 0) {
    DrainFully(&treat_stream_);
    if (params_.build_treatment_local_lambda_events) {
      DrainFully(&ctrl_stream_);
    }
    if (cur_chrom_id_ >= static_cast<int32_t>(treat_tracks_.size())) {
      treat_tracks_.resize(cur_chrom_id_ + 1);
      ctrl_raw_tracks_.resize(cur_chrom_id_ + 1);
    }
    treat_tracks_[cur_chrom_id_] = std::move(treat_stream_.track);
    if (params_.build_treatment_local_lambda_events) {
      ctrl_raw_tracks_[cur_chrom_id_] = std::move(ctrl_stream_.track);
    }
    ResetStream(&treat_stream_);
    ResetStream(&ctrl_stream_);
    // Byte-parity summits: capture buffers are kept until MoveInto, then
    // handed to dest as pos_array_starts_capture / pos_array_ends_capture
    // so FinalizeMacs3FragPosArrayPositions can build them in parallel
    // across chroms (the per-chrom build is the second-largest serial
    // cost after the sweep itself; running it post-sweep with #pragma omp
    // saves ~15-20 s wall at PBMC scale).
  }
  cur_chrom_id_ = static_cast<int32_t>(new_chrom_id);
  prev_start_ = -1;
  return true;
}

bool Macs3FragSweepWorkspace::AddFragment(uint32_t chrom_id, int32_t start,
                                          int32_t end, int32_t count) {
  if (finalized_) {
    error_ = "AddFragment called after Finalize";
    return false;
  }
  if (end <= start || count <= 0) {
    error_ = "invalid fragment (end <= start or count <= 0)";
    return false;
  }
  const int32_t cid = static_cast<int32_t>(chrom_id);
  if (cur_chrom_id_ >= 0 && cid < cur_chrom_id_) {
    error_ = "input not chrom-grouped (chrom_id decreased after transition)";
    return false;
  }
  if (cid != cur_chrom_id_) {
    if (!TransitionToChrom(chrom_id)) return false;
  }
  if (start < prev_start_) {
    error_ = "input not start-sorted within chromosome";
    return false;
  }

  const int32_t w = EffectiveFragWeight(count, params_.frag_pileup_max_count,
                                        params_.macs3_uint8_counts);
  if (w == 0) {
    prev_start_ = start;
    return true;
  }

  // Drain past events strictly before the new fragment's earliest event for
  // each stream. The treat stream's earliest new event is at `start`; the
  // ctrl (lambda-window) stream's earliest is at `start - lambda_half_`.
  DrainBefore(&treat_stream_, start);
  if (params_.build_treatment_local_lambda_events) {
    DrainBefore(&ctrl_stream_, start - lambda_half_);
  }

  // Treat span events: +w at [start, end).
  treat_stream_.heap.push({start, w});
  treat_stream_.heap.push({end, -w});

  // Lambda fragment-end window events: +w over [start - half, start + half]
  // and over [end - half, end + half], matching MACS3
  // CallerFromAlignments::pileup_a_chromosome_c with ctrl_d_s = [llocal_bp].
  if (params_.build_treatment_local_lambda_events) {
    const int32_t ll = start - lambda_half_;
    const int32_t rl = ll + params_.llocal_bp;
    const int32_t lr = end - lambda_half_;
    const int32_t rr = lr + params_.llocal_bp;
    ctrl_stream_.heap.push({ll, w});
    ctrl_stream_.heap.push({rl, -w});
    ctrl_stream_.heap.push({lr, w});
    ctrl_stream_.heap.push({rr, -w});
  }

  // Byte-parity summits: capture per-chrom count-duplicated endpoints.
  // Each fragment with raw count=N contributes N copies of (start, end)
  // (matches MACS3 FragParser, which creates N PETrackI instances).
  if (params_.byte_parity_summits &&
      params_.capture_raw_pos_array_endpoints) {
    if (cur_chrom_id_ >= static_cast<int32_t>(frag_starts_.size())) {
      frag_starts_.resize(cur_chrom_id_ + 1);
      frag_ends_.resize(cur_chrom_id_ + 1);
    }
    std::vector<int32_t>& fs = frag_starts_[cur_chrom_id_];
    std::vector<int32_t>& fe = frag_ends_[cur_chrom_id_];
    for (int32_t k = 0; k < count; ++k) {
      fs.push_back(start);
      fe.push_back(end);
    }
  }

  prev_start_ = start;

  const int32_t row_count = RowCountForMacsTotals(
      count, params_.frag_pileup_max_count, params_.macs3_uint8_counts);
  if (row_count > 0) {
    total_fragment_count_ += row_count;
    total_fragment_bp_ +=
        static_cast<int64_t>(end - start) * static_cast<int64_t>(row_count);
  }
  return true;
}

bool Macs3FragSweepWorkspace::Finalize() {
  if (finalized_) return true;
  if (cur_chrom_id_ >= 0) {
    DrainFully(&treat_stream_);
    if (params_.build_treatment_local_lambda_events) {
      DrainFully(&ctrl_stream_);
    }
    if (cur_chrom_id_ >= static_cast<int32_t>(treat_tracks_.size())) {
      treat_tracks_.resize(cur_chrom_id_ + 1);
      ctrl_raw_tracks_.resize(cur_chrom_id_ + 1);
    }
    treat_tracks_[cur_chrom_id_] = std::move(treat_stream_.track);
    if (params_.build_treatment_local_lambda_events) {
      ctrl_raw_tracks_[cur_chrom_id_] = std::move(ctrl_stream_.track);
    }
    ResetStream(&treat_stream_);
    ResetStream(&ctrl_stream_);
    // Byte-parity capture buffers stay populated; consumed in MoveInto.
  }

  // Apply lambda_bg / scale / Build... to convert raw ctrl tracks to final
  // lambda tracks (mirrors FinalizeMacs3FragLambdaTracks).
  const size_t n = treat_tracks_.size();
  if (ctrl_raw_tracks_.size() < n) ctrl_raw_tracks_.resize(n);
  lambda_tracks_.clear();
  lambda_tracks_.resize(n);
  if (params_.build_treatment_local_lambda_events &&
      total_fragment_count_ > 0 && total_fragment_bp_ > 0) {
    // MACS3 stores both quantities as C floats before multiplying the local
    // scale by endpoint depth.  Preserve that operation order for score
    // parity at float32 rounding boundaries.
    const float lambda_bg = static_cast<float>(
        static_cast<double>(total_fragment_bp_) /
        static_cast<double>(params_.effective_genome_size));
    const float scale = static_cast<float>(
        static_cast<double>(total_fragment_bp_) /
        (static_cast<double>(params_.llocal_bp) *
         static_cast<double>(total_fragment_count_) * 2.0));
    for (size_t i = 0; i < n; ++i) {
      Macs3BdgTrack ctrl = ctrl_raw_tracks_[i];
      for (float& v : ctrl.v) {
        float z = v * scale;
        if (z < lambda_bg) z = lambda_bg;
        v = z;
      }
      // Keep MACS3's original float32 lambda for Poisson scoring.  The
      // diagnostic bedGraph writer performs its own five-decimal formatting;
      // quantizing the in-memory track here changes p/q scores at boundaries.
      BuildLambdaTrackFromTreatAndControl(treat_tracks_[i].p, ctrl,
                                          &lambda_tracks_[i], false, nullptr);
    }
  }
  ctrl_raw_tracks_.clear();
  finalized_ = true;
  return true;
}

void Macs3FragSweepWorkspace::SetChromNames(std::vector<std::string> names) {
  chrom_names_ = std::move(names);
}

bool Macs3FragSweepWorkspace::MoveInto(Macs3FragPeakWorkspace* dest) {
  if (dest == nullptr || !finalized_) {
    error_ = "MoveInto called before Finalize or with null dest";
    return false;
  }
  dest->params = params_;
  dest->chrom_names = std::move(chrom_names_);
  dest->treat_tracks = std::move(treat_tracks_);
  dest->lambda_tracks = std::move(lambda_tracks_);
  dest->ppois_tracks.clear();
  // events_by_chrom is not used by the sweep path but the downstream
  // pipeline checks size matches treat/lambda; pad with empty entries.
  dest->events_by_chrom.assign(dest->treat_tracks.size(), Macs3FragEventBuffer());
  // Byte-parity summits: transfer the per-chrom (start, end) capture
  // buffers to the events-workspace finalize slots so the downstream
  // FinalizeMacs3FragPosArrayPositions step can build pos_array
  // positions in parallel across chroms (#pragma omp parallel for).
  if (params_.byte_parity_summits &&
      params_.capture_raw_pos_array_endpoints) {
    frag_starts_.resize(dest->treat_tracks.size());
    frag_ends_.resize(dest->treat_tracks.size());
    dest->pos_array_starts_capture = std::move(frag_starts_);
    dest->pos_array_ends_capture = std::move(frag_ends_);
    dest->pos_array_positions.assign(dest->treat_tracks.size(),
                                     std::vector<int32_t>());
  } else {
    std::vector<std::vector<int32_t>>().swap(frag_starts_);
    std::vector<std::vector<int32_t>>().swap(frag_ends_);
  }
  // Drop the now-unused inline pos_array store (kept around for back-
  // compat with previous shape, but never populated under the
  // post-sweep parallel-build flow).
  std::vector<std::vector<int32_t>>().swap(pos_array_positions_);
  // Reconcile sizes between chrom_names and the per-chrom data structures.
  // chrom_names may be larger (e.g. when the caller pre-populates a full
  // genome chrom table via SetChromNames before any fragment is added; the
  // sweep only resizes treat_tracks_ up to max_chrom_id_seen+1) or smaller
  // (when the caller doesn't supply chrom_names; treat_tracks_ governs).
  // The downstream FromWorkspace lex-sort + parallel-for require equal
  // sizes across all per-chrom vectors, so pad whichever is shorter to
  // the max of the two.
  const size_t target_n = std::max(dest->chrom_names.size(),
                                   dest->treat_tracks.size());
  if (dest->chrom_names.size() < target_n) {
    dest->chrom_names.resize(target_n);
  }
  if (dest->treat_tracks.size() < target_n) {
    dest->treat_tracks.resize(target_n);
  }
  if (dest->lambda_tracks.size() < target_n) {
    dest->lambda_tracks.resize(target_n);
  }
  if (dest->events_by_chrom.size() < target_n) {
    dest->events_by_chrom.resize(target_n);
  }
  if (params_.byte_parity_summits &&
      params_.capture_raw_pos_array_endpoints) {
    if (dest->pos_array_starts_capture.size() < target_n) {
      dest->pos_array_starts_capture.resize(target_n);
    }
    if (dest->pos_array_ends_capture.size() < target_n) {
      dest->pos_array_ends_capture.resize(target_n);
    }
    if (dest->pos_array_positions.size() < target_n) {
      dest->pos_array_positions.resize(target_n);
    }
  }
  dest->total_fragment_bp = total_fragment_bp_;
  dest->total_fragment_count = total_fragment_count_;

  // Reset *this so it can be re-Init()ed.
  cur_chrom_id_ = -1;
  prev_start_ = -1;
  total_fragment_count_ = 0;
  total_fragment_bp_ = 0;
  finalized_ = false;
  return true;
}

}  // namespace peaks
}  // namespace chromap
