#ifndef CHROMAP_PEAK_CALLER_MACS3_FRAG_WORKSPACE_H_
#define CHROMAP_PEAK_CALLER_MACS3_FRAG_WORKSPACE_H_

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "rapidmacs/fragment_input.h"
#include "rapidmacs/macs3_bdg_track.h"

namespace chromap {
namespace peaks {

struct Macs3FragWorkspaceParams {
  int32_t frag_pileup_max_count = 0;
  bool macs3_uint8_counts = true;
  int64_t effective_genome_size = 2913022398LL;
  int32_t llocal_bp = 10000;
  double score_pseudocount = 0.0;
  // Treatment-derived llocal events are needed only when no matched control
  // is supplied. Controlled BAMPE/BEDPE calling builds lambda entirely from
  // the control endpoints, so retaining four unused events per treatment
  // fragment only increases memory.
  bool build_treatment_local_lambda_events = true;
  // Raw treatment endpoints are required by the no-control pos-array
  // reconstruction. Controlled calling reuses the already exact union of
  // finalized treatment/control breakpoints and can disable this capture.
  bool capture_raw_pos_array_endpoints = true;
  // Default: build MACS3-byte-identical pos_array positions and feed them to
  // the summit-finding tied-plateau split. The treatment-only path captures
  // count-duplicated fragment endpoints (treat_pv ∪ ctrl_pv); controlled
  // calling obtains the same exact union from finalized track breakpoints and
  // does not retain a second endpoint copy. Validated 0 / 50,274 summit
  // differences vs compiled MACS3 v3.0.3 on the full PBMC 3k benchmark; the
  // treatment-only capture costs ~1 GB RSS and ~6 s vs the collapsed shortcut.
  // Set to false to take the simple/fast shortcut (median over
  // collapsed treat-track pieces); summit positions then diverge from
  // the MACS3 binary on ~2% of peaks at PBMC scale (max ~222 bp;
  // BED3 / treat / lambda / ppois bedGraphs and pValue/qValue/signal
  // columns within 1e-5 are unaffected by this flag).
  bool byte_parity_summits = true;
};

// Parameters specific to MACS3's paired-end treatment-vs-control path.
struct Macs3FragControlParams {
  int32_t slocal_bp = 1000;
  // MACS3 uses the treatment parser's pre-dedup mean fragment length for
  // its d-sized control window. If <= 0, use the retained treatment mean.
  double treatment_fragment_length_bp = 0.0;
  bool keep_duplicates_all = false;
  // False is MACS3's default `--scale-to small`; true matches
  // `--scale-to large` / legacy `--to-large`.
  bool scale_to_large = false;
  // Use only genome-wide lambda after applying the normal treatment/control
  // depth-scaling decision; do not construct d/slocal/llocal windows.
  bool no_lambda = false;
};

struct Macs3FragEvent {
  Macs3FragEvent() : pos(0), delta(0) {}
  Macs3FragEvent(int32_t p, int32_t d) : pos(p), delta(d) {}

  int32_t pos;
  int32_t delta;
};

struct Macs3FragEventBuffer {
  std::vector<Macs3FragEvent> treat_span_events;
  std::vector<Macs3FragEvent> lambda_end_window_events;
};

struct Macs3FragPeakWorkspace {
  Macs3FragWorkspaceParams params;
  std::vector<std::string> chrom_names;
  std::vector<Macs3FragEventBuffer> events_by_chrom;
  std::vector<Macs3BdgTrack> treat_tracks;
  std::vector<Macs3BdgTrack> lambda_tracks;
  std::vector<Macs3BdgTrack> ppois_tracks;
  // Populated only when params.byte_parity_summits == true. Each
  // per-chrom entry is the sorted unique pos_array positions emitted by
  // MACS3 callpeak's pileup_treat_ctrl_a_chromosome -> __chrom_pair_treat_ctrl
  // for that chrom. Consumed by CollectMacs3FragNarrowPeakRowsFromTracks
  // to split tied summit plateaus at the same per-event boundaries
  // MACS3 sees, recovering byte-identical summit positions.
  std::vector<std::vector<int32_t>> pos_array_positions;
  // Transient capture buffers for the events-workspace path (sweep
  // workspace populates pos_array_positions directly). Filled by
  // AddMacs3FragWorkspaceFragment when byte_parity_summits is on, and
  // consumed + cleared by FinalizeMacs3FragPosArrayPositions.
  std::vector<std::vector<int32_t>> pos_array_starts_capture;
  std::vector<std::vector<int32_t>> pos_array_ends_capture;
  int64_t total_fragment_bp = 0;
  int64_t total_fragment_count = 0;
};

bool InitMacs3FragWorkspace(const std::vector<std::string>& chrom_names,
                            const Macs3FragWorkspaceParams& params,
                            Macs3FragPeakWorkspace* workspace);

bool AddMacs3FragWorkspaceFragment(Macs3FragPeakWorkspace* workspace,
                                   uint32_t chrom_id, int32_t start,
                                   int32_t end, int32_t count);

// Resize events_by_chrom up to chrom_id + 1 if needed. Used by the
// streaming entry point (RunMacs3FragPeakPipelineFromIterator) when the
// iterator interns chrom names in first-seen order rather than supplying
// the full chrom list up front. Idempotent; safe to call before each
// AddMacs3FragWorkspaceFragment.
bool EnsureMacs3FragWorkspaceChromCapacity(Macs3FragPeakWorkspace* workspace,
                                           uint32_t chrom_id);

bool BuildMacs3FragWorkspaceFromFragments(
    const std::vector<ChromFragments>& by_chrom,
    const Macs3FragWorkspaceParams& params,
    Macs3FragPeakWorkspace* workspace);

// Build the finalized treatment and treatment-derived local-lambda tracks
// directly from chromosome-grouped, start-sorted fragments.  Unlike the
// general events workspace, this path merges the already sorted start and end
// event streams without materializing and comparison-sorting six events per
// fragment.  Chromosomes are independent and are built in parallel.
//
// When release_fragments is true, each per-chrom fragment vector is released
// after its tracks and optional summit captures have been constructed.
bool BuildMacs3FragTracksFromSortedFragments(
    std::vector<ChromFragments>* by_chrom,
    const Macs3FragWorkspaceParams& params,
    Macs3FragPeakWorkspace* workspace,
    int num_threads = 1,
    bool release_fragments = false);

// Finalization is one-shot for the event-backed workspace: treat events are
// released by FinalizeMacs3FragTreatTracks, and lambda events are released by
// FinalizeMacs3FragLambdaTracks.
bool FinalizeMacs3FragTreatTracks(Macs3FragPeakWorkspace* workspace,
                                   int num_threads = 1);
bool FinalizeMacs3FragLambdaTracks(Macs3FragPeakWorkspace* workspace,
                                    int num_threads = 1);

// Build a constant genome-wide lambda track ending at each treatment track's
// last position. This is MACS3 --nolambda for treatment-only paired input and
// deliberately avoids all local-window event construction.
bool FinalizeMacs3FragGlobalLambdaTracks(
    Macs3FragPeakWorkspace* workspace, int num_threads = 1);

// Replace the treatment-derived no-control lambda with MACS3's paired-end
// control lambda: the maximum of d-, slocal-, and llocal-sized control
// endpoint pileups plus genome-wide lambda, with treatment/control depth
// normalization. Treatment tracks are scaled in-place when MACS3 selects the
// control depth as the target. Only chromosomes shared by treatment and
// control receive scoreable tracks, matching CallerFromAlignments.
bool FinalizeMacs3FragLambdaTracksWithControl(
    Macs3FragPeakWorkspace* workspace,
    const std::vector<ChromFragments>& control_by_chrom,
    const Macs3FragControlParams& control_params,
    int num_threads = 1);
bool FinalizeMacs3FragPpoisTracks(Macs3FragPeakWorkspace* workspace,
                                   int num_threads = 1,
                                   bool quantize_values = true);
// Build pos_array_positions from the count-duplicated capture buffers
// (events-workspace path) and release the buffers. No-op if
// params.byte_parity_summits is false. Idempotent.
bool FinalizeMacs3FragPosArrayPositions(Macs3FragPeakWorkspace* workspace,
                                        int num_threads = 1);

bool BuildMacs3FragPpoisTrack(const Macs3BdgTrack& treat,
                              const Macs3BdgTrack& lambda,
                              double pseudocount, Macs3BdgTrack* out,
                              bool quantize_values = true);

// Realigns a (scaled) ctrl track to the treat track's breakpoints and
// emits the final per-chrom lambda track. Used by both the events-based
// finalize and the sweep-line workspace post-processing.
void BuildLambdaTrackFromTreatAndControl(const std::vector<int32_t>& t_p,
                                         const Macs3BdgTrack& ctrl,
                                         Macs3BdgTrack* out,
                                         bool quantize_values = true,
                                         std::vector<int32_t>* merged_positions = nullptr);

bool WriteMacs3BdgTracksBedGraph(const std::string& path,
                                 const std::vector<std::string>& chrom_names,
                                 const std::vector<Macs3BdgTrack>& tracks);

}  // namespace peaks
}  // namespace chromap

#endif  // CHROMAP_PEAK_CALLER_MACS3_FRAG_WORKSPACE_H_
