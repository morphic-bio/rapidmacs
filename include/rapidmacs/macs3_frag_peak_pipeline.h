#ifndef CHROMAP_PEAK_CALLER_MACS3_FRAG_PEAK_PIPELINE_H_
#define CHROMAP_PEAK_CALLER_MACS3_FRAG_PEAK_PIPELINE_H_

#include <cstdint>
#include <string>
#include <vector>

#include "rapidmacs/bdgpeakcall.h"
#include "rapidmacs/fragment_input.h"
#include "rapidmacs/fragments.h"  // macs3::FragmentIterator

namespace chromap {
namespace peaks {

class StageProfileCollector;
struct Macs3FragPeakWorkspace;

// Parameters for the validated MACS3-compatible FRAG pileup → lambda → ppois
// → bdgpeakcall → narrowPeak/summits path (see parity harnesses).
struct Macs3FragPeakPipelineParams {
  // New public API: threshold_mode/qvalue_cutoff. The older
  // cutoff_is_qvalue field below remains supported for existing consumers.
  Macs3FragThresholdMode threshold_mode = Macs3FragThresholdMode::kPValue;
  double qvalue_cutoff = 0.0;
  // -log10(p) by default, or -log10(q) when cutoff_is_qvalue is true.
  float bdgpeakcall_cutoff = 5.f;
  bool cutoff_is_qvalue = false;
  // MACS3 CallPeakUnit stores p/q histogram keys and arithmetic results in
  // float32 for FRAG, BAMPE, and BEDPE.  False retains the historical
  // RapidMACS round-5/double mode for callers that explicitly request it.
  bool macs3_float_score_math = true;
  int32_t min_length = 200;
  int32_t max_gap = 30;
  // MACS3 --max-count: applies to all count uses, including summit
  // endpoint multiplicities; 0 disables the cap.
  int32_t frag_pileup_max_count = 0;
  bool macs3_uint8_counts = true;
  int64_t effective_genome_size = 2913022398LL;
  int32_t slocal_bp = 1000;
  int32_t llocal_bp = 10000;
  double treatment_fragment_length_bp = 0.0;
  // Select MACS3's parser-time float division for PE `--keep-dup all`
  // control-depth scaling. Numeric caps use the retained-fragment mean.
  bool keep_duplicates_all = false;
  bool scale_to_large = false;
  // Skip treatment/control local-window construction and use only MACS3's
  // genome-wide background lambda.
  bool no_lambda = false;
  // Resolve multiple summits within a peak using MACS3's smoothed-signal
  // subpeak algorithm. False preserves the established single-summit output.
  bool call_summits = false;
  double score_pseudocount = 0.0;
  // MACS3 callpeak applies PeakIO.filter_fc after region/summit calling.
  // Its default is 1.0, which removes statistically significant depletion
  // regions (treatment < local lambda) from the final peak files.
  double fold_enrichment_cutoff = 1.0;
  int peak_caller_threads = 1;
  // If non-empty, peaks are named "<name_prefix>_peak_<N>" (matches MACS3
  // callpeak's `-n NAME` convention). Otherwise "peak_<N>".
  std::string name_prefix;
  // Default: produce summit positions byte-identical to the compiled
  // MACS3 binary by capturing per-chrom count-duplicated fragment
  // endpoints and computing pos_array positions (treat_pv ∪ ctrl_pv)
  // at finalize. Set to false for a smaller-memory shortcut that
  // diverges on ~2% of summits at PBMC scale.
  bool byte_parity_summits = true;
};

// Optional explicit intermediate bedGraph paths (empty = files under work_dir).
struct Macs3FragPeakPipelinePaths {
  std::string treat_bdg;
  std::string lambda_bdg;
  std::string ppois_bdg;
};

// -log10(p) cutoff for bdgpeakcall on the ppois score track (MACS3 -p 1e-5 → 5).
float BdgPeakCallCutoffFromPValue(double p_value);

// Fragment-span pileup, control_lambda, ppois score, then narrowPeak + summits.
// keep_intermediates_dir: if non-empty, writes treat/lambda/ppois bedGraphs there
// (directory is created) and does not delete it. If empty, uses a unique
// directory under work_dir_parent (or /tmp) and removes it after success.
bool RunMacs3FragPeakPipelineFromFragments(
    std::vector<ChromFragments>* by_chrom,
    const Macs3FragPeakPipelineParams& params,
    const Macs3FragPeakPipelinePaths& user_paths,
    const std::string& narrowpeak_out, const std::string& summits_out,
    const std::string& keep_intermediates_dir,
    const std::string& work_dir_parent, std::string* work_dir_used,
    std::string* error_message,
    StageProfileCollector* stage_profile = nullptr,
    bool release_fragments_after_workspace = true,
    const std::vector<ChromFragments>* control_by_chrom = nullptr);

// Streaming entry point: drains a RapidMACS FragmentIterator directly into
// the workspace without materializing the per-chrom fragment list first.
// Same outputs as RunMacs3FragPeakPipelineFromFragments (and same
// downstream pipeline code), but callers (libchromap, file-source) don't
// need to allocate std::vector<ChromFragments> up front.
//
// The iterator is consumed; its ChromNames() at end-of-stream is used as
// the workspace chrom name table.
//
// Memory: this entry point uses the events-buffer workspace
// (~48 B/fragment). For sorted input (typical Chromap/ARC fragments TSV)
// prefer RunMacs3FragPeakPipelineFromSortedIterator below — same outputs,
// dramatically less memory (~tens of KB total instead of GB).
bool RunMacs3FragPeakPipelineFromIterator(
    macs3::FragmentIterator& iter,
    const Macs3FragPeakPipelineParams& params,
    const Macs3FragPeakPipelinePaths& user_paths,
    const std::string& narrowpeak_out, const std::string& summits_out,
    const std::string& keep_intermediates_dir,
    const std::string& work_dir_parent, std::string* work_dir_used,
    std::string* error_message,
    StageProfileCollector* stage_profile = nullptr);

// Sweep-line streaming entry point: same outputs as
// RunMacs3FragPeakPipelineFromIterator but uses the Macs3FragSweepWorkspace
// instead of the events-buffer workspace, producing treat/lambda BdgTracks
// incrementally. **Requires chrom-grouped, start-sorted input** —
// monotonically non-decreasing chrom_id across calls and monotonically
// non-decreasing start within each chrom_id. On violation returns false
// with a descriptive error. Memory peak is bounded by max in-flight
// fragment-end window events (typically tens of KB total).
//
// Use this for FileFragmentIterator (Chromap/ARC fragments TSV are
// already chrom + start sorted). The bytewise output of treat / lambda /
// ppois / narrowPeak is identical to the non-streaming path.
bool RunMacs3FragPeakPipelineFromSortedIterator(
    macs3::FragmentIterator& iter,
    const Macs3FragPeakPipelineParams& params,
    const Macs3FragPeakPipelinePaths& user_paths,
    const std::string& narrowpeak_out, const std::string& summits_out,
    const std::string& keep_intermediates_dir,
    const std::string& work_dir_parent, std::string* work_dir_used,
    std::string* error_message,
    StageProfileCollector* stage_profile = nullptr);

bool RunMacs3FragPeakPipelineFromWorkspace(
    Macs3FragPeakWorkspace* workspace,
    const Macs3FragPeakPipelineParams& params,
    const Macs3FragPeakPipelinePaths& user_paths,
    const std::string& narrowpeak_out, const std::string& summits_out,
    const std::string& keep_intermediates_dir,
    const std::string& work_dir_parent, std::string* work_dir_used,
    std::string* error_message,
    StageProfileCollector* stage_profile = nullptr);

}  // namespace peaks
}  // namespace chromap

#endif  // CHROMAP_PEAK_CALLER_MACS3_FRAG_PEAK_PIPELINE_H_
