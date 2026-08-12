#ifndef CHROMAP_PEAK_CALLER_PEAK_IO_H_
#define CHROMAP_PEAK_CALLER_PEAK_IO_H_

#include <cstdio>
#include <string>
#include <vector>

#include "rapidmacs/call_peaks.h"
#include "rapidmacs/fragment_input.h"

namespace chromap {
namespace peaks {

// Writes 10-col narrowPeak (tab); sorts peaks by (chrom, start, end) first.
// Returns false on I/O or invalid intervals.
// `name_prefix`: if non-empty, peaks are named "<name_prefix>_peak_<N>"
// (matches MACS3 callpeak's `-n NAME` naming convention). If empty, the
// shorter "peak_<N>" form is used (RapidMACS's historical default).
bool WriteNarrowPeak(const std::string& path,
                     const std::vector<Peak>& peaks,
                     const std::string& name_prefix = "");

// Simple tab-separated metrics: peak_count, total_bp, median_width_bp,
// argv, params, and (when saturated_bins >= 0) total_bins / saturated_bins
// / saturated_bin_fraction. Emit -1 for saturated_bins to suppress those
// rows when the caller didn't track them.
bool WriteRunSummaryTsv(
    const std::string& path, const std::vector<Peak>& all_peaks,
    const std::string& argv_line, const std::string& program_params,
    int64_t total_bins, int64_t saturated_bins);

// bedGraph (4 columns): chrom, start, end, signal. Merges adjacent bins with
// identical nonzero counts. Omits zero bins (no MACS3-style zero runs).
bool WritePileupBedGraph(const std::string& path,
                         const std::vector<ChromFragments>& by_chrom,
                         int32_t bin_width_bp, int32_t ext_size_bp);

// Fragment-span pileup (MACS3 `pileup -f FRAG` semantics): each row adds its
// effective weight along [start, end). Writes bedGraph like MACS3 pileup_PV:
// segments [prev, pos) with float values, including zero-coverage runs from 0
// to the first breakpoint. No UCSC track line (matches macs3 pileup_cmd).
//
// max_count: if > 0, per-row weight is min(col5, max_count) (MACS3 --max-count).
// macs3_uint8_counts: if true, coerce per-row weight through uint8_t (numpy u1)
// after max_count, matching MACS3 PETrackII storage; use for exact parity vs
// macs3 when col5 can exceed 255.
bool WriteFragSpanPileupBedGraph(const std::string& path,
                                 const std::vector<ChromFragments>& by_chrom,
                                 int32_t max_count, bool macs3_uint8_counts);

// MACS3 `callpeak -f FRAG` without control: local lambda on fragment-end windows
// (llocal) after PeakDetect / PETrackII.pileup_a_chromosome_c, plus lambda_bg
// = total_fragment_bp / effective_genome. Diagnostic only. Use
// -g hs (2913022398) to match `macs3 callpeak -g hs`.
// MACS3 `control_lambda`: `__chrom_pair_treat_ctrl` (min-merge of treat/ctrl
// pileup `p` arrays with `c_v[ic]`) + `__write_bedGraph` 1e-5 control row merge.
// No tail merge past the shorter p-array (matches published MACS3 Cython).
// Only llocal (no slocal) for no-control FRAG, per `__call_peaks_wo_control`.
bool WriteFragMacs3NoControlLambdaBedGraph(
    const std::string& path, const std::vector<ChromFragments>& by_chrom,
    int32_t max_count, bool macs3_uint8_counts, int64_t effective_genome_size,
    int32_t llocal_bp);

// MACS3 `bdgcmp` parity on `*_treat_pileup.bdg` vs `*_control_lambda.bdg`:
// merge breakpoints like BedGraph.make_ScoreTrackII_for_macs (stops when either
// track ends), apply the same pseudocount as `macs3 bdgcmp -p` (default 0),
// then emit score bedGraphs with 5-digit merge (|Δ|>1e-5) like ScoreTrackII.
// ppois: -log10 Poisson upper tail with MACS3's discrete convention
// (ScoreTrack.get_pscore / Prob.poisson_cdf log10 mode).
// FE: (treat+pseudocount)/(control+pseudocount). Empty path skips that output.
bool WriteFragMacs3BdgcmpScoreBedGraphs(
    const std::string& path_ppois, const std::string& path_fe,
    const std::vector<ChromFragments>& by_chrom, int32_t max_count,
    bool macs3_uint8_counts, int64_t effective_genome_size, int32_t llocal_bp,
    double pseudocount);

// Same scoring as WriteFragMacs3BdgcmpScoreBedGraphs but reads treat/lambda from
// existing bedGraphs (no scratch rebuild from fragments). Use when those files
// are already the MACS3-style pileup and control_lambda tracks.
bool WriteMacs3BdgcmpScoreBedGraphsFromTreatLambdaBdgs(
    const std::string& path_ppois, const std::string& path_fe,
    const std::string& path_treat_bdg, const std::string& path_lambda_bdg,
    double pseudocount);

// MACS3 bdgcmp -m ppois value at (treat, control): -log10(Poisson p). Treat/ctrl
// are the scan-quantized float pileup values as in the MACS3 bedGraph inputs.
double Macs3BdgcmpPpoisNegLog10P(float treat_f, float ctrl_f, double pseudocount);

}  // namespace peaks
}  // namespace chromap

#endif  // CHROMAP_PEAK_CALLER_PEAK_IO_H_
