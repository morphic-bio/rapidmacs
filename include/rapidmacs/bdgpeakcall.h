#ifndef CHROMAP_PEAK_CALLER_BDGPEAKCALL_H_
#define CHROMAP_PEAK_CALLER_BDGPEAKCALL_H_

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "rapidmacs/macs3_bdg_track.h"

namespace chromap {
namespace peaks {

class StageProfileCollector;

enum class Macs3FragThresholdMode {
  kPValue,
  kQValue,
};

struct Macs3FragNarrowPeakRow {
  std::string chrom;
  int32_t start = 0;
  int32_t end = 0;
  int32_t summit = 0;
  int32_t peak_off = 0;
  double fc = 0.0;
  double p_mlog = 0.0;
  double q_mlog = 0.0;
};

// MACS3 `bdgpeakcall` region calling on a score bedGraph (diagnostic).
// Reads the file like MACS3 BedGraphIO.read_bedGraph / bedGraphTrackI.add_loc
// (sequential merge; no implicit gap fill). Runs BedGraph.call_peaks with
// call_summits=False (interval endpoints from merged above-cutoff segments).
//
// Writes BED3 (chrom, start, end), one interval per line, sorted lexicographically
// by chrom then start, end.
bool RunMacs3StyleBdgPeakCallFromBedGraph(const std::string& path_bdg,
                                          const std::string& path_bed_out,
                                          float cutoff, int32_t min_length,
                                          int32_t max_gap);

// MACS3-style narrowPeak + summits: regions from ppois bdg (same as
// `bdgpeakcall` on that track). By default, summit is the midpoint of the
// treat-pileup-maximal above-cutoff subsegment. The track overload's optional
// call_summits flag instead applies MACS3's smoothed multi-subpeak algorithm.
// p/q and signal are evaluated from treat+lambda at each selected summit.
bool RunMacs3FragPpoisNarrowPeaks(
    const std::string& path_ppois_bdg, const std::string& path_treat_bdg,
    const std::string& path_lambda_bdg, float cutoff, int32_t min_length,
    int32_t max_gap, double pseudocount, const std::string& path_narrowpeak,
    const std::string& path_summits_bed, StageProfileCollector* stage_profile = nullptr,
    int64_t profile_input_frag_rows = -1,
    const std::string& name_prefix = "",
    Macs3FragThresholdMode threshold_mode = Macs3FragThresholdMode::kPValue,
    double qvalue_cutoff = 0.0);

bool RunMacs3FragPpoisNarrowPeaksFromTracks(
    const std::vector<std::string>& chrom_names,
    const std::vector<Macs3BdgTrack>& ppois_tracks,
    const std::vector<Macs3BdgTrack>& treat_tracks,
    const std::vector<Macs3BdgTrack>& lambda_tracks,
    float cutoff, int32_t min_length, int32_t max_gap, double pseudocount,
    const std::string& path_narrowpeak, const std::string& path_summits_bed,
    StageProfileCollector* stage_profile = nullptr,
    int64_t profile_input_frag_rows = -1,
    const std::string& name_prefix = "",
    bool cutoff_is_qvalue = false,
    bool macs3_float_qvalue_math = false,
    const std::vector<std::vector<int32_t>>* pos_array_positions = nullptr,
    bool call_summits = false,
    double fold_enrichment_cutoff = 0.0);

// Compatibility overload used by the MACS3 BED/Signac profile and by callers
// of the previously published q-value threshold API.
bool RunMacs3FragPpoisNarrowPeaksFromTracks(
    const std::vector<std::string>& chrom_names,
    const std::vector<Macs3BdgTrack>& ppois_tracks,
    const std::vector<Macs3BdgTrack>& treat_tracks,
    const std::vector<Macs3BdgTrack>& lambda_tracks,
    float cutoff, int32_t min_length, int32_t max_gap, double pseudocount,
    const std::string& path_narrowpeak, const std::string& path_summits_bed,
    StageProfileCollector* stage_profile,
    int64_t profile_input_frag_rows,
    const std::string& name_prefix,
    Macs3FragThresholdMode threshold_mode,
    double qvalue_cutoff);

// `pos_array_positions`: optional sorted unique pos_array positions for
// this chrom (treat_pv ∪ ctrl_pv positions, mirroring MACS3
// __chrom_pair_treat_ctrl). When non-null and non-empty, tied summit
// plateaus are split at every interior pos_array position before the
// median tie-pick, so summit positions match the compiled MACS3 binary
// byte-for-byte. May be nullptr (or empty) for the legacy collapsed-
// pieces behavior, which follows MACS3's source-documented algorithm
// against RapidMACS's collapsed treat track.
bool CollectMacs3FragNarrowPeakRowsFromTracks(
    const std::string& chrom_name, const Macs3BdgTrack& ppois,
    const Macs3BdgTrack& treat, const Macs3BdgTrack& lambda, float cutoff,
    int32_t min_length, int32_t max_gap, double pseudocount,
    std::vector<Macs3FragNarrowPeakRow>* rows,
    const std::vector<int32_t>* pos_array_positions = nullptr,
    bool macs3_float_peak_scores = false,
    const std::map<double, double>* region_p_to_q = nullptr,
    bool call_summits = false);

// Region selection may use a qscore track while summit p/q values continue to
// be evaluated from the original ppois track.
bool CollectMacs3FragNarrowPeakRowsFromThresholdTrack(
    const std::string& chrom_name, const Macs3BdgTrack& threshold_score,
    const Macs3BdgTrack& ppois, const Macs3BdgTrack& treat,
    const Macs3BdgTrack& lambda, float cutoff, int32_t min_length,
    int32_t max_gap, double pseudocount,
    std::vector<Macs3FragNarrowPeakRow>* rows,
    const std::vector<int32_t>* pos_array_positions = nullptr,
    bool macs3_float_peak_scores = false,
    bool call_summits = false);

// Optional p_to_q: genome-wide pscore -> qscore table (round-5 keys).
// If non-null and non-empty, qscore is looked up from this table (matches
// MACS3 callpeak's pqtable). Otherwise BH FDR is computed over the emitted
// peak rows only (legacy behavior).
// name_prefix: if non-empty, peaks are named "<name_prefix>_peak_<N>"
// matching MACS3 callpeak's `-n NAME` convention; otherwise "peak_<N>".
bool WriteMacs3FragNarrowPeakRows(
    std::vector<Macs3FragNarrowPeakRow>* rows,
    const std::string& path_narrowpeak, const std::string& path_summits_bed,
    const std::map<double, double>* p_to_q = nullptr,
    const std::string& name_prefix = "",
    bool narrowpeak_score_is_qvalue = false);

// MACS3 callpeak filters completed peaks by summit fold enrichment after peak
// calling. A non-positive cutoff disables filtering for lower-level
// bdgpeakcall-compatible callers.
void FilterMacs3FragNarrowPeakRowsByFoldEnrichment(
    std::vector<Macs3FragNarrowPeakRow>* rows, double cutoff);

bool WriteMacs3FragNarrowPeakRowsWithScoreCap(
    std::vector<Macs3FragNarrowPeakRow>* rows,
    const std::string& path_narrowpeak, const std::string& path_summits_bed,
    const std::map<double, double>* p_to_q,
    const std::string& name_prefix,
    bool narrowpeak_score_uses_qvalue,
    bool cap_narrowpeak_score,
    bool macs3_float_output = false);

// Build MACS3-compatible pqtable: genome-wide BH-adjusted -log10(qvalue)
// for each unique -log10(pvalue) seen across all ppois track segments.
// Paired CallPeakUnit parity uses exact float32 keys; the legacy FRAG path
// retains its established five-decimal buckets. Empty output map if input is
// empty / total bp is non-positive.
void BuildMacs3PqTable(const std::vector<Macs3BdgTrack>& ppois_tracks,
                       std::map<double, double>* p_to_q,
                       bool macs3_float_math = false);

// Convert pscore tracks to MACS3 qscore tracks using a genome-wide pqtable.
// Positions are unchanged and values use float32, matching CallPeakUnit's
// qscore arrays. Used when callpeak's default `-q` cutoff is requested.
void BuildMacs3QscoreTracks(
    const std::vector<Macs3BdgTrack>& ppois_tracks,
    const std::map<double, double>& p_to_q,
    std::vector<Macs3BdgTrack>* qscore_tracks);

bool BuildMacs3QScoreTrackFromPpois(const Macs3BdgTrack& ppois,
                                    const std::map<double, double>& p_to_q,
                                    Macs3BdgTrack* qscore,
                                    bool macs3_float_qtable = false);

}  // namespace peaks
}  // namespace chromap

#endif  // CHROMAP_PEAK_CALLER_BDGPEAKCALL_H_
