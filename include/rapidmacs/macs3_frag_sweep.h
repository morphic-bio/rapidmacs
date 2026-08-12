#ifndef CHROMAP_PEAK_CALLER_MACS3_FRAG_SWEEP_H_
#define CHROMAP_PEAK_CALLER_MACS3_FRAG_SWEEP_H_

#include <cstdint>
#include <queue>
#include <string>
#include <vector>

#include "rapidmacs/macs3_bdg_track.h"
#include "rapidmacs/macs3_frag_workspace.h"

namespace chromap {
namespace peaks {

// Streaming MACS3 FRAG workspace using sweep-line BdgTrack emission. Memory
// is bounded by the maximum number of "in-flight" events (open fragments
// for treat; fragment-end window events for lambda) — typically tens of KB
// total even at PBMC-3k scale, vs the ~3.8 GB events buffer used by the
// classic Macs3FragPeakWorkspace.
//
// **Input order constraint**: fragments must arrive in chrom-grouped,
// start-sorted order. The sweep maintains:
//   - chrom_id monotonically non-decreasing across calls,
//   - start monotonically non-decreasing within the same chrom_id.
// Violations return false and set Error() to a path:line:reason message.
//
// On Finalize(), produces treat_tracks and lambda_tracks bytewise identical
// to FinalizeMacs3FragTreatTracks + FinalizeMacs3FragLambdaTracks on the
// same fragments. Move the result into a Macs3FragPeakWorkspace via
// MoveInto() and feed it to RunMacs3FragPeakPipelineFromWorkspace.
class Macs3FragSweepWorkspace {
 public:
  Macs3FragSweepWorkspace();

  bool Init(const Macs3FragWorkspaceParams& params);

  // Add one fragment. See class doc for ordering constraints.
  bool AddFragment(uint32_t chrom_id, int32_t start, int32_t end, int32_t count);

  // Drain remaining heap events, scale + lambda_bg-floor the raw lambda
  // tracks, then run the BuildLambdaTrackFromTreatAndControl re-alignment.
  // After this returns true, treat_tracks() / lambda_tracks() are ready.
  bool Finalize();

  // Optional: set the chrom name table. If unset, names default to "" and
  // can be set by the caller post-Finalize.
  void SetChromNames(std::vector<std::string> names);

  // Move workspace contents into `dest` (a classic Macs3FragPeakWorkspace).
  // dest->events_by_chrom is left empty (sweep produces tracks directly).
  // After this returns, *this is reset to its post-Init() state.
  bool MoveInto(Macs3FragPeakWorkspace* dest);

  int64_t total_fragment_count() const { return total_fragment_count_; }
  int64_t total_fragment_bp() const { return total_fragment_bp_; }
  const std::string& Error() const { return error_; }

 private:
  struct Event {
    int32_t pos;
    int32_t delta;
    bool operator>(const Event& other) const { return pos > other.pos; }
  };

  struct EventStream {
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>> heap;
    int32_t s = 0;          // last event position seen
    double z = 0.0;         // current accumulated value
    double pre_z = -1.0e4;  // last emitted value (sentinel for first emit)
    Macs3BdgTrack track;
  };

  bool TransitionToChrom(uint32_t chrom_id);
  void DrainBefore(EventStream* stream, int32_t target);
  void DrainFully(EventStream* stream);
  void EmitClose(EventStream* stream, int32_t pos);
  void ResetStream(EventStream* stream);

  Macs3FragWorkspaceParams params_;
  int32_t lambda_half_;

  int32_t cur_chrom_id_;
  int32_t prev_start_;

  EventStream treat_stream_;
  EventStream ctrl_stream_;  // raw lambda-endpoint counts (pre scale + lambda_bg)

  std::vector<std::string> chrom_names_;
  std::vector<Macs3BdgTrack> treat_tracks_;
  std::vector<Macs3BdgTrack> ctrl_raw_tracks_;
  std::vector<Macs3BdgTrack> lambda_tracks_;
  // Byte-parity summit mode (params.byte_parity_summits). Per-chrom
  // count-duplicated fragment-endpoint capture. frag_starts_/frag_ends_
  // grow during AddFragment; on chrom transition / finalize we build
  // pos_array_positions_[chrom_id] from them and release the capture
  // buffers for that chrom.
  std::vector<std::vector<int32_t>> frag_starts_;
  std::vector<std::vector<int32_t>> frag_ends_;
  std::vector<std::vector<int32_t>> pos_array_positions_;

  int64_t total_fragment_count_;
  int64_t total_fragment_bp_;
  bool finalized_;
  std::string error_;
};

}  // namespace peaks
}  // namespace chromap

#endif  // CHROMAP_PEAK_CALLER_MACS3_FRAG_SWEEP_H_
