#include "rapidmacs/bdgpeakcall.h"

#include "rapidmacs/call_peaks.h"
#include "rapidmacs/peak_io.h"
#include "rapidmacs/stage_profile.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace chromap {
namespace peaks {

// (BuildMacs3PqTable is declared in rapidmacs/bdgpeakcall.h)

namespace {

bool StartsWith(const std::string& s, const char* prefix) {
  const size_t n = std::strlen(prefix);
  return s.size() >= n && s.compare(0, n, prefix) == 0;
}

void AddLocMacs3(std::map<std::string, Macs3BdgTrack>* by_chrom,
                 const std::string& chromosome, int32_t startpos, int32_t endpos,
                 float value, float baseline_value) {
  if (endpos <= 0) {
    return;
  }
  if (startpos < 0) {
    startpos = 0;
  }
  Macs3BdgTrack& ch = (*by_chrom)[chromosome];
  if (ch.p.empty()) {
    if (startpos > 0) {
      ch.p.push_back(startpos);
      ch.v.push_back(baseline_value);
    }
    ch.p.push_back(endpos);
    ch.v.push_back(value);
    return;
  }
  const float pre_v = ch.v.back();
  if (pre_v == value) {
    ch.p.back() = endpos;
  } else {
    ch.p.push_back(endpos);
    ch.v.push_back(value);
  }
}

bool ReadBedGraphMacs3Style(const std::string& path, float baseline_value,
                            std::map<std::string, Macs3BdgTrack>* out) {
  std::ifstream in(path.c_str());
  if (!in) {
    return false;
  }
  out->clear();
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    if (StartsWith(line, "track") || StartsWith(line, "#") ||
        StartsWith(line, "browse")) {
      continue;
    }
    std::istringstream iss(line);
    std::string chrom;
    std::string s1, s2, s3;
    if (!(iss >> chrom >> s1 >> s2 >> s3)) {
      continue;
    }
    const int32_t start = static_cast<int32_t>(std::strtol(s1.c_str(), nullptr, 10));
    const int32_t end = static_cast<int32_t>(std::strtol(s2.c_str(), nullptr, 10));
    const float val = static_cast<float>(std::strtod(s3.c_str(), nullptr));
    if (end <= start) {
      continue;
    }
    AddLocMacs3(out, chrom, start, end, val, baseline_value);
  }
  return in.eof() || !in.bad();
}

struct PeakPiece {
  int32_t start = 0;
  int32_t end = 0;
  float value = 0.f;
};

void ClosePeakMacs3(const std::vector<PeakPiece>& peak_content, int32_t min_length,
                    std::vector<std::pair<int32_t, int32_t>>* out_ivs) {
  if (peak_content.empty()) {
    return;
  }
  const int32_t peak_length = peak_content.back().end - peak_content.front().start;
  if (peak_length < min_length) {
    return;
  }
  out_ivs->emplace_back(peak_content.front().start, peak_content.back().end);
}

// Mirrors MACS3 Signal/BedGraph.py bedGraphTrackI.call_peaks (call_summits ignored
// for interval bounds; summits are not written in diagnostic BED3 mode).
void CallPeaksOnTrack(const Macs3BdgTrack& tr, float cutoff, int32_t min_length,
                      int32_t max_gap,
                      std::vector<std::pair<int32_t, int32_t>>* out_ivs) {
  out_ivs->clear();
  const size_t n = tr.v.size();
  if (n == 0 || tr.p.size() != n) {
    return;
  }
  std::vector<char> above(n);
  size_t n_above = 0;
  for (size_t i = 0; i < n; ++i) {
    // MACS3 BedGraph/ScoreTrack bdgpeakcall includes values equal to the
    // cutoff. This diagnostic path intentionally differs from callpeak's
    // CallPeakUnit thresholding below, which uses a strict comparison.
    above[i] = (tr.v[i] >= cutoff) ? 1 : 0;
    if (above[i]) {
      ++n_above;
    }
  }
  if (n_above == 0) {
    return;
  }
  std::vector<int32_t> endpos;
  std::vector<int32_t> startpos;
  std::vector<float> vals;
  endpos.reserve(n_above);
  startpos.reserve(n_above);
  vals.reserve(n_above);
  for (size_t i = 0; i < n; ++i) {
    if (!above[i]) {
      continue;
    }
    endpos.push_back(tr.p[i]);
    vals.push_back(tr.v[i]);
  }
  for (size_t i = 0; i < n; ++i) {
    if (!above[i]) continue;
    startpos.push_back(i == 0 ? 0 : tr.p[i - 1]);
  }
  if (startpos.size() != endpos.size()) {
    return;
  }
  std::vector<PeakPiece> peak_content;
  peak_content.reserve(endpos.size());
  {
    PeakPiece pc;
    pc.start = startpos[0];
    pc.end = endpos[0];
    pc.value = vals[0];
    peak_content.push_back(pc);
  }
  for (size_t i = 1; i < startpos.size(); ++i) {
    if (startpos[i] - peak_content.back().end <= max_gap) {
      PeakPiece pc;
      pc.start = startpos[i];
      pc.end = endpos[i];
      pc.value = vals[i];
      peak_content.push_back(pc);
    } else {
      ClosePeakMacs3(peak_content, min_length, out_ivs);
      peak_content.clear();
      PeakPiece pc;
      pc.start = startpos[i];
      pc.end = endpos[i];
      pc.value = vals[i];
      peak_content.push_back(pc);
    }
  }
  if (!peak_content.empty()) {
    ClosePeakMacs3(peak_content, min_length, out_ivs);
  }
}

}  // namespace

bool RunMacs3StyleBdgPeakCallFromBedGraph(const std::string& path_bdg,
                                          const std::string& path_bed_out,
                                          float cutoff, int32_t min_length,
                                          int32_t max_gap) {
  if (path_bdg.empty() || path_bed_out.empty() || min_length < 1 || max_gap < 0) {
    return false;
  }
  std::map<std::string, Macs3BdgTrack> tracks;
  if (!ReadBedGraphMacs3Style(path_bdg, 0.f, &tracks)) {
    return false;
  }
  std::vector<std::string> chroms;
  chroms.reserve(tracks.size());
  for (const auto& kv : tracks) {
    chroms.push_back(kv.first);
  }
  std::sort(chroms.begin(), chroms.end());

  std::vector<std::pair<std::string, std::pair<int32_t, int32_t>>> all;
  for (const std::string& c : chroms) {
    const Macs3BdgTrack& tr = tracks[c];
    std::vector<std::pair<int32_t, int32_t>> ivs;
    CallPeaksOnTrack(tr, cutoff, min_length, max_gap, &ivs);
    for (const auto& iv : ivs) {
      all.emplace_back(c, iv);
    }
  }
  std::sort(all.begin(), all.end(),
            [](const std::pair<std::string, std::pair<int32_t, int32_t>>& a,
               const std::pair<std::string, std::pair<int32_t, int32_t>>& b) {
              if (a.first != b.first) {
                return a.first < b.first;
              }
              if (a.second.first != b.second.first) {
                return a.second.first < b.second.first;
              }
              return a.second.second < b.second.second;
            });

  FILE* fp = std::fopen(path_bed_out.c_str(), "w");
  if (fp == nullptr) {
    return false;
  }
  for (const auto& e : all) {
    if (e.second.second <= e.second.first) {
      std::fclose(fp);
      return false;
    }
    if (std::fprintf(fp, "%s\t%d\t%d\n", e.first.c_str(), static_cast<int>(e.second.first),
                     static_cast<int>(e.second.second)) < 0) {
      std::fclose(fp);
      return false;
    }
  }
  if (std::fclose(fp) != 0) {
    return false;
  }
  return true;
}

// --- MACS3-style narrowPeak: bdgpeakcall regions on ppois, summit = max(treat) chunk ---

struct PpoisMergePiece {
  int32_t start = 0;
  int32_t end = 0;
};

struct PpoisMergedRegion {
  int32_t start = 0;
  int32_t end = 0;
  std::vector<PpoisMergePiece> pieces;
};

static double LookupQFromPqTable(const std::map<double, double>& p_to_q,
                                 double p_mlog);

float BdgValueAtOrZero(const Macs3BdgTrack& tr, int32_t pos) {
  if (tr.p.empty() || tr.v.size() != tr.p.size()) {
    return 0.f;
  }
  const auto it = std::upper_bound(tr.p.begin(), tr.p.end(), pos);
  if (it == tr.p.end()) {
    return 0.f;
  }
  const size_t j = static_cast<size_t>(it - tr.p.begin());
  const int32_t left = (j == 0) ? 0 : tr.p[j - 1];
  return (pos >= left) ? tr.v[j] : 0.f;
}

// Same as CallPeaksOnTrack, but also records each above-cutoff segment as a PpoisMergePiece
// (genomic [start, end) from bdg, value = ppois for that run).
void CallPpoisTrackMergedRegions(
    const Macs3BdgTrack& tr, float cutoff, int32_t min_length, int32_t max_gap,
    std::vector<PpoisMergedRegion>* out,
    const std::map<double, double>* p_to_q = nullptr,
    bool pscore_keys_are_float = false) {
  out->clear();
  const size_t n = tr.v.size();
  if (n == 0 || tr.p.size() != n) {
    return;
  }
  double pscore_cutoff = std::numeric_limits<double>::infinity();
  if (p_to_q != nullptr) {
    // qscore is monotone in pscore. Convert the q cutoff once to the lowest
    // pscore key that passes instead of materializing a duplicate q track or
    // performing a map lookup for every genomic run.
    for (const auto& kv : *p_to_q) {
      if (static_cast<float>(kv.second) > cutoff) {
        pscore_cutoff = kv.first;
        break;
      }
    }
  }
  std::vector<char> above(n);
  size_t n_above = 0;
  for (size_t i = 0; i < n; ++i) {
    bool pass = false;
    if (p_to_q == nullptr) {
      // MACS3 CallPeakUnit.apply_multiple_cutoffs uses score > cutoff for
      // both p- and q-value modes. Equality must not enter a peak region.
      pass = tr.v[i] > cutoff;
    } else if (std::isfinite(pscore_cutoff)) {
      const double key =
          pscore_keys_are_float
              ? static_cast<double>(static_cast<float>(tr.v[i]))
              : std::round(static_cast<double>(tr.v[i]) * 1e5) / 1e5;
      pass = key >= pscore_cutoff;
    }
    above[i] = pass ? 1 : 0;
    if (above[i]) {
      ++n_above;
    }
  }
  if (n_above == 0) {
    return;
  }
  std::vector<int32_t> endpos;
  std::vector<int32_t> startpos;
  endpos.reserve(n_above);
  startpos.reserve(n_above);
  for (size_t i = 0; i < n; ++i) {
    if (!above[i]) {
      continue;
    }
    endpos.push_back(tr.p[i]);
  }
  for (size_t i = 0; i < n; ++i) {
    if (!above[i]) continue;
    startpos.push_back(i == 0 ? 0 : tr.p[i - 1]);
  }
  if (startpos.size() != endpos.size()) {
    return;
  }
  struct PeakPiece {
    int32_t start = 0;
    int32_t end = 0;
  };
  auto close_peak = [&](const std::vector<PeakPiece>& pc,
                          std::vector<PpoisMergedRegion>* dst) {
    if (pc.empty()) {
      return;
    }
    const int32_t peak_length = pc.back().end - pc.front().start;
    if (peak_length < min_length) {
      return;
    }
    PpoisMergedRegion r;
    r.start = pc.front().start;
    r.end = pc.back().end;
    r.pieces.reserve(pc.size());
    for (const auto& a : pc) {
      PpoisMergePiece e;
      e.start = a.start;
      e.end = a.end;
      r.pieces.push_back(e);
    }
    dst->push_back(std::move(r));
  };
  std::vector<PeakPiece> peak_content;
  peak_content.reserve(endpos.size());
  {
    PeakPiece pc;
    pc.start = startpos[0];
    pc.end = endpos[0];
    peak_content.push_back(pc);
  }
  for (size_t i = 1; i < startpos.size(); ++i) {
    if (startpos[i] - peak_content.back().end <= max_gap) {
      PeakPiece pc;
      pc.start = startpos[i];
      pc.end = endpos[i];
      peak_content.push_back(pc);
    } else {
      close_peak(peak_content, out);
      peak_content.clear();
      PeakPiece pc;
      pc.start = startpos[i];
      pc.end = endpos[i];
      peak_content.push_back(pc);
    }
  }
  if (!peak_content.empty()) {
    close_peak(peak_content, out);
  }
}

// MACS3 narrowPeak col 5 (score): int(10 * peak[score_column]). With `-p`
// cutoff, score_column is 'pscore' (-log10 pvalue); with `-q` cutoff, it is
// 'qscore'. Python's int() truncates toward zero — NOT round-half-up.
// MACS3 does not clamp this value to the nominal narrowPeak [0,1000] range;
// highly significant peaks can therefore have scores above 1000.
static int NarrowPeakIntScoreQ(double neg_log10_q) {
  if (neg_log10_q < 0.0) {
    return 0;
  }
  int s = static_cast<int>(10.0 * neg_log10_q);
  if (s < 0) {
    s = 0;
  }
  return s;
}

bool RunMacs3FragPpoisNarrowPeaks(
    const std::string& path_ppois_bdg, const std::string& path_treat_bdg,
    const std::string& path_lambda_bdg, float cutoff, int32_t min_length,
    int32_t max_gap, double pseudocount, const std::string& path_narrowpeak,
    const std::string& path_summits_bed, StageProfileCollector* stage_profile,
    int64_t profile_input_frag_rows, const std::string& name_prefix,
    Macs3FragThresholdMode threshold_mode, double qvalue_cutoff) {
  if (path_ppois_bdg.empty() || path_treat_bdg.empty() || path_lambda_bdg.empty() ||
      (path_narrowpeak.empty() && path_summits_bed.empty()) || min_length < 1 ||
      max_gap < 0) {
    return false;
  }
  std::map<std::string, Macs3BdgTrack> m_p;
  std::map<std::string, Macs3BdgTrack> m_t;
  std::map<std::string, Macs3BdgTrack> m_l;
  if (!ReadBedGraphMacs3Style(path_ppois_bdg, 0.f, &m_p) ||
      !ReadBedGraphMacs3Style(path_treat_bdg, 0.f, &m_t) ||
      !ReadBedGraphMacs3Style(path_lambda_bdg, 0.f, &m_l)) {
    return false;
  }
  std::vector<std::string> chrom_names;
  for (const auto& kv : m_p) {
    chrom_names.push_back(kv.first);
  }
  std::sort(chrom_names.begin(), chrom_names.end());
  std::vector<Macs3BdgTrack> ppois_tracks;
  std::vector<Macs3BdgTrack> treat_tracks;
  std::vector<Macs3BdgTrack> lambda_tracks;
  ppois_tracks.reserve(chrom_names.size());
  treat_tracks.reserve(chrom_names.size());
  lambda_tracks.reserve(chrom_names.size());
  for (const std::string& chrom : chrom_names) {
    ppois_tracks.push_back(m_p[chrom]);
    const auto itt = m_t.find(chrom);
    const auto itl = m_l.find(chrom);
    treat_tracks.push_back(itt == m_t.end() ? Macs3BdgTrack() : itt->second);
    lambda_tracks.push_back(itl == m_l.end() ? Macs3BdgTrack() : itl->second);
  }
  return RunMacs3FragPpoisNarrowPeaksFromTracks(
      chrom_names, ppois_tracks, treat_tracks, lambda_tracks, cutoff, min_length,
      max_gap, pseudocount, path_narrowpeak, path_summits_bed, stage_profile,
      profile_input_frag_rows, name_prefix, threshold_mode, qvalue_cutoff);
}

bool RunMacs3FragPpoisNarrowPeaksFromTracks(
    const std::vector<std::string>& chrom_names,
    const std::vector<Macs3BdgTrack>& ppois_tracks,
    const std::vector<Macs3BdgTrack>& treat_tracks,
    const std::vector<Macs3BdgTrack>& lambda_tracks,
    float cutoff, int32_t min_length, int32_t max_gap, double pseudocount,
    const std::string& path_narrowpeak, const std::string& path_summits_bed,
    StageProfileCollector* stage_profile, int64_t profile_input_frag_rows,
    const std::string& name_prefix, Macs3FragThresholdMode threshold_mode,
    double qvalue_cutoff) {
  const bool q_mode = threshold_mode == Macs3FragThresholdMode::kQValue;
  if (threshold_mode != Macs3FragThresholdMode::kPValue && !q_mode) {
    return false;
  }
  if ((!q_mode && (std::isnan(cutoff) || cutoff < 0.f)) ||
      (q_mode && (std::isnan(qvalue_cutoff) || qvalue_cutoff <= 0.0 ||
                  qvalue_cutoff > 1.0))) {
    return false;
  }
  const float effective_cutoff =
      q_mode ? static_cast<float>(-std::log10(qvalue_cutoff)) : cutoff;
  return RunMacs3FragPpoisNarrowPeaksFromTracks(
      chrom_names, ppois_tracks, treat_tracks, lambda_tracks, effective_cutoff,
      min_length, max_gap, pseudocount, path_narrowpeak, path_summits_bed,
      stage_profile, profile_input_frag_rows, name_prefix, q_mode,
      /*macs3_float_qvalue_math=*/true,
      /*pos_array_positions=*/nullptr);
}

bool RunMacs3FragPpoisNarrowPeaksFromTracks(
    const std::vector<std::string>& chrom_names,
    const std::vector<Macs3BdgTrack>& ppois_tracks,
    const std::vector<Macs3BdgTrack>& treat_tracks,
    const std::vector<Macs3BdgTrack>& lambda_tracks,
    float cutoff, int32_t min_length, int32_t max_gap, double pseudocount,
    const std::string& path_narrowpeak, const std::string& path_summits_bed,
    StageProfileCollector* stage_profile, int64_t profile_input_frag_rows,
    const std::string& name_prefix, bool cutoff_is_qvalue,
    bool macs3_float_qvalue_math,
    const std::vector<std::vector<int32_t>>* pos_array_positions,
    bool call_summits, double fold_enrichment_cutoff) {
  if ((path_narrowpeak.empty() && path_summits_bed.empty()) || min_length < 1 ||
      max_gap < 0 || chrom_names.size() != ppois_tracks.size() ||
      chrom_names.size() != treat_tracks.size() ||
      chrom_names.size() != lambda_tracks.size() ||
      std::isnan(fold_enrichment_cutoff) || fold_enrichment_cutoff < 0.0 ||
      (pos_array_positions != nullptr &&
       chrom_names.size() != pos_array_positions->size())) {
    return false;
  }
  StageProfileCollector::Tick t_regions_start{};
  if (stage_profile) {
    t_regions_start = StageProfileCollector::Now();
  }
  // Build the genome-wide pqtable before region collection: q-value mode
  // thresholds the qscore track, while peak p/q columns are still evaluated
  // from the original pscore track at the selected summit.
  std::map<double, double> pqtable;
  BuildMacs3PqTable(ppois_tracks, &pqtable, macs3_float_qvalue_math);
  std::vector<Macs3FragNarrowPeakRow> w;
  w.reserve(256);
  for (size_t ci = 0; ci < chrom_names.size(); ++ci) {
    const std::vector<int32_t>* pap =
        pos_array_positions == nullptr ? nullptr : &pos_array_positions->at(ci);
    if (!CollectMacs3FragNarrowPeakRowsFromTracks(
            chrom_names[ci], ppois_tracks[ci], treat_tracks[ci],
            lambda_tracks[ci],
            cutoff, min_length, max_gap, pseudocount, &w, pap,
            macs3_float_qvalue_math,
            cutoff_is_qvalue ? &pqtable : nullptr,
            call_summits)) {
      return false;
    }
  }
  FilterMacs3FragNarrowPeakRowsByFoldEnrichment(
      &w, fold_enrichment_cutoff);
  StageProfileCollector::Tick t_np_start{};
  if (stage_profile) {
    const int64_t in_fr =
        profile_input_frag_rows >= 0 ? profile_input_frag_rows : 0;
    stage_profile->Record("bdgpeakcall_regions", t_regions_start,
                          StageProfileCollector::Now(), in_fr,
                          static_cast<int64_t>(w.size()), 0,
                          "in_memory_tracks");
    t_np_start = StageProfileCollector::Now();
  }
  if (!WriteMacs3FragNarrowPeakRows(&w, path_narrowpeak, path_summits_bed,
                                    &pqtable, name_prefix, cutoff_is_qvalue)) {
    return false;
  }
  if (stage_profile) {
    const StageProfileCollector::Tick t_summits_end = StageProfileCollector::Now();
    const StageProfileCollector::Tick tm0 = StageProfileCollector::Now();
    int64_t bbytes = 0, sbytes = 0;
    if (!path_narrowpeak.empty()) {
      (void)StageProfileCollector::FileMetrics(path_narrowpeak, &bbytes, nullptr);
    }
    if (!path_summits_bed.empty()) {
      (void)StageProfileCollector::FileMetrics(path_summits_bed, &sbytes, nullptr);
    }
    const StageProfileCollector::Tick tm1 = StageProfileCollector::Now();
    const int64_t in_fr = profile_input_frag_rows >= 0 ? profile_input_frag_rows : 0;
    stage_profile->Record("narrowpeak_summits", t_np_start, t_summits_end, in_fr,
                          static_cast<int64_t>(w.size()), bbytes + sbytes,
                          "in_memory_tracks");
    stage_profile->Record("profile_file_metrics_narrowpeak_summits", tm0, tm1, 0, 0,
                          bbytes + sbytes, "stat_size_only");
  }
  return true;
}

namespace {

// A MACS3 peak_content chunk after splitting score runs at every position in
// the merged treatment/control position array.
struct Macs3SummitChunk {
  int32_t start = 0;
  int32_t end = 0;
  float treat = 0.f;
  float control = 0.f;
};

int SmoothedDerivativeSign(double x) {
  // SignalProcessing.maxima rounds the Savitzky-Golay derivative to 16
  // decimals before np.sign. nearbyint follows the default ties-to-even mode.
  const double rounded = std::nearbyint(x * 1.0e16) / 1.0e16;
  return (rounded > 0.0) - (rounded < 0.0);
}

std::vector<int32_t> Macs3SignalMaxima(const std::vector<float>& signal,
                                       int32_t window_size) {
  std::vector<int32_t> maxima;
  if (signal.empty() || window_size < 1) return maxima;
  window_size = (window_size / 2) * 2 + 1;
  const int32_t half = (window_size - 1) / 2;
  const size_t n = signal.size();
  std::vector<float> padded;
  padded.reserve(n + static_cast<size_t>(half) * 2);
  // NumPy slicing truncates naturally when half exceeds signal length.
  const size_t available = n > 0 ? std::min<size_t>(half, n - 1) : 0;
  for (size_t j = available; j > 0; --j) {
    const float delta = std::fabs(signal[j] - signal[0]);
    padded.push_back(static_cast<float>(signal[0] - delta));
  }
  padded.insert(padded.end(), signal.begin(), signal.end());
  for (size_t j = 1; j <= available; ++j) {
    const size_t source = n - 1 - j;
    const float delta = std::fabs(signal[source] - signal.back());
    padded.push_back(static_cast<float>(signal.back() + delta));
  }
  const size_t width = static_cast<size_t>(window_size);
  if (padded.size() < width) return maxima;
  const size_t output_n = padded.size() - width + 1;
  // For a symmetric order-2 Savitzky-Golay first derivative, the coefficient
  // at k is k/sum(k^2). Evaluate each dot product independently: the rounded
  // derivative sign is deliberately sensitive near flat plateaus, and a
  // sliding accumulator introduces drift that NumPy's convolution does not.
  double denom = 0.0;
  for (int32_t k = -half; k <= half; ++k) {
    denom += static_cast<double>(k) * static_cast<double>(k);
  }
  if (denom == 0.0) return maxima;
  std::vector<int8_t> signs(output_n, 0);
  for (size_t i = 0; i < output_n; ++i) {
    double derivative = 0.0;
    for (int32_t k = -half; k <= half; ++k) {
      derivative +=
          (static_cast<double>(k) / denom) *
          static_cast<double>(padded[i + static_cast<size_t>(k + half)]);
    }
    signs[i] = static_cast<int8_t>(SmoothedDerivativeSign(derivative));
  }
  for (size_t i = 0; i + 1 < signs.size(); ++i) {
    if (static_cast<int>(signs[i + 1]) - static_cast<int>(signs[i]) <= -1) {
      maxima.push_back(static_cast<int32_t>(i));
    }
  }
  return maxima;
}

bool Macs3SignalPeakIsValid(const std::vector<float>& signal,
                            int32_t maximum) {
  // NumPy accepts maximum == len(signal) here; hard_clip's right scan is
  // simply empty. This occurs for adjacent maxima in the first subregion.
  if (maximum < 0 || static_cast<size_t>(maximum) > signal.size()) {
    return false;
  }
  // Match SignalProcessing.hard_clip. Its left-hand loop is empty for the
  // positive range arguments used there, so only the first negative value to
  // the right clips the signal.
  size_t right = signal.size();
  for (size_t i = static_cast<size_t>(maximum); i < signal.size(); ++i) {
    if (signal[i] < 0.f) {
      right = i;
      break;
    }
  }
  if (right < 50) return false;
  std::vector<float> unique(signal.begin(), signal.begin() + right);
  std::sort(unique.begin(), unique.end());
  unique.erase(std::unique(unique.begin(), unique.end()), unique.end());
  return unique.size() >= 6;
}

std::vector<int32_t> Macs3EnforcePeakyness(
    const std::vector<float>& signal, const std::vector<int32_t>& maxima) {
  if (maxima.size() <= 1) return maxima;
  std::vector<int32_t> minima;
  minima.reserve(maxima.size() - 1);
  for (size_t i = 0; i + 1 < maxima.size(); ++i) {
    const int32_t begin = maxima[i];
    const int32_t end = maxima[i + 1];
    int32_t best = begin;
    for (int32_t p = begin + 1; p < end; ++p) {
      if (signal[static_cast<size_t>(p)] <
          signal[static_cast<size_t>(best)]) {
        best = p;
      }
    }
    minima.push_back(best);
  }
  if (minima.empty()) return maxima;
  std::vector<int32_t> kept;
  kept.reserve(maxima.size());
  auto subtract_and_validate = [&](size_t begin, size_t end, float threshold,
                                   int32_t maximum) {
    std::vector<float> adjusted;
    adjusted.reserve(end - begin);
    for (size_t i = begin; i < end; ++i) {
      adjusted.push_back(static_cast<float>(signal[i] - threshold));
    }
    return Macs3SignalPeakIsValid(adjusted, maximum);
  };

  float threshold = signal[static_cast<size_t>(minima.front())];
  threshold = static_cast<float>(threshold + std::sqrt(threshold));
  const float first_threshold =
      static_cast<float>(threshold + std::sqrt(threshold));
  if (subtract_and_validate(0, static_cast<size_t>(minima.front()),
                            first_threshold, maxima.front())) {
    kept.push_back(maxima.front());
  }
  for (size_t i = 0; i + 1 < minima.size(); ++i) {
    threshold = std::max(signal[static_cast<size_t>(minima[i])],
                         signal[static_cast<size_t>(minima[i + 1])]);
    threshold = static_cast<float>(threshold + std::sqrt(threshold));
    const int32_t new_maximum = maxima[i + 1] - minima[i];
    if (subtract_and_validate(static_cast<size_t>(minima[i]),
                              static_cast<size_t>(minima[i + 1]), threshold,
                              new_maximum)) {
      kept.push_back(maxima[i + 1]);
    }
  }
  threshold = signal[static_cast<size_t>(minima.back())];
  threshold = static_cast<float>(threshold + std::sqrt(threshold));
  const int32_t last_maximum = maxima.back() - minima.back();
  if (subtract_and_validate(static_cast<size_t>(minima.back()), signal.size(),
                            threshold, last_maximum)) {
    kept.push_back(maxima.back());
  }
  return kept;
}

std::string Macs3SubpeakLetters(size_t i) {
  const char suffix = static_cast<char>('a' + (i % 26));
  if (i < 26) return std::string(1, suffix);
  return Macs3SubpeakLetters(i / 26) + suffix;
}

std::vector<std::string> Macs3PeakNames(
    const std::vector<Macs3FragNarrowPeakRow>& rows,
    const std::string& name_prefix) {
  std::vector<std::string> names(rows.size());
  size_t peak_number = 0;
  size_t begin = 0;
  while (begin < rows.size()) {
    size_t end = begin + 1;
    while (end < rows.size() && rows[end].chrom == rows[begin].chrom &&
           rows[end].end == rows[begin].end) {
      ++end;
    }
    ++peak_number;
    const std::string base =
        (name_prefix.empty() ? "peak_" : name_prefix + "_peak_") +
        std::to_string(peak_number);
    if (end - begin == 1) {
      names[begin] = base;
    } else {
      for (size_t i = begin; i < end; ++i) {
        names[i] = base + Macs3SubpeakLetters(i - begin);
      }
    }
    begin = end;
  }
  return names;
}

}  // namespace

static bool CollectMacs3FragNarrowPeakRowsImpl(
    const std::string& chrom_name, const Macs3BdgTrack& threshold_score,
    const Macs3BdgTrack& ppois,
    const Macs3BdgTrack& treat, const Macs3BdgTrack& lambda, float cutoff,
    int32_t min_length, int32_t max_gap, double pseudocount,
    std::vector<Macs3FragNarrowPeakRow>* rows,
    const std::vector<int32_t>* pos_array_positions,
    bool macs3_float_peak_scores,
    const std::map<double, double>* region_p_to_q,
    bool call_summits) {
  if (rows == nullptr || min_length < 1 || max_gap < 0 ||
      threshold_score.p.size() != threshold_score.v.size() ||
      ppois.p.size() != ppois.v.size() || treat.p.size() != treat.v.size() ||
      lambda.p.size() != lambda.v.size()) {
    return false;
  }
  if (threshold_score.p.empty() || ppois.p.empty() || treat.p.empty() ||
      lambda.p.empty()) {
    return true;
  }
  std::vector<PpoisMergedRegion> regs;
  CallPpoisTrackMergedRegions(threshold_score, cutoff, min_length, max_gap, &regs,
                              region_p_to_q, macs3_float_peak_scores);
  // When pos_array_positions is provided, each tied plateau is split at
  // its interior positions before the median tie-pick, recovering the
  // same per-event chunking MACS3 callpeak's peak_content sees and
  // producing byte-identical summit positions vs the MACS3 binary.
  const bool have_pap = (pos_array_positions != nullptr &&
                         !pos_array_positions->empty());
  for (const PpoisMergedRegion& r : regs) {
    if (r.pieces.empty()) {
      continue;
    }
    std::vector<Macs3SummitChunk> chunks;
    auto add_chunk = [&](int32_t s, int32_t e) {
      if (e <= s) return;
      Macs3SummitChunk chunk;
      chunk.start = s;
      chunk.end = e;
      chunk.treat = BdgValueAtOrZero(treat, s);
      chunk.control = BdgValueAtOrZero(lambda, s);
      chunks.push_back(chunk);
    };
    for (size_t i = 0; i < r.pieces.size(); ++i) {
      const PpoisMergePiece& pc = r.pieces[i];
      // BED-derived score tracks can contain repeated endpoints when two
      // tag events share a coordinate. Such zero-width pieces carry no
      // genomic signal; skipping them also keeps upper_bound(start) from
      // advancing past lower_bound(end) in the split-point walk below.
      if (pc.end <= pc.start) {
        continue;
      }
      if (!have_pap) {
        add_chunk(pc.start, pc.end);
        continue;
      }
      // A p/q-score run can span several pos_array chunks whose treatment
      // pileups differ (different treatment/control pairs can have the same
      // score). MACS3 evaluates the summit pileup for every such chunk, so
      // re-evaluate treatment after each split instead of inheriting the
      // value at the beginning of the score run.
      auto lo = std::upper_bound(pos_array_positions->begin(),
                                 pos_array_positions->end(), pc.start);
      auto hi = std::lower_bound(pos_array_positions->begin(),
                                 pos_array_positions->end(), pc.end);
      int32_t cur = pc.start;
      for (auto it = lo; it != hi; ++it) {
        const int32_t p = *it;
        add_chunk(cur, p);
        cur = p;
      }
      add_chunk(cur, pc.end);
    }
    if (chunks.empty()) {
      continue;
    }
    auto append_row = [&](int32_t summit_pos, const Macs3SummitChunk& chunk) {
      Macs3FragNarrowPeakRow row;
      row.chrom = chrom_name;
      row.start = r.start;
      row.end = r.end;
      row.summit = summit_pos;
      row.peak_off = static_cast<int32_t>(row.summit - row.start);
      // MACS3 narrowPeak col 7 (signalValue) uses pseudocount=1 hardcoded;
      // pscore (col 8) keeps the configured score pseudocount.
      const double raw_fc = (static_cast<double>(chunk.treat) + 1.0) /
                            (static_cast<double>(chunk.control) + 1.0);
      row.fc = macs3_float_peak_scores ? static_cast<float>(raw_fc) : raw_fc;
      row.p_mlog = Macs3BdgcmpPpoisNegLog10P(
          chunk.treat, chunk.control, pseudocount);
      if (macs3_float_peak_scores) {
        row.p_mlog = static_cast<float>(row.p_mlog);
      }
      row.q_mlog = 0.0;
      rows->push_back(std::move(row));
    };
    auto append_single_summit = [&]() {
      float best = 0.f;
      bool have_best = false;
      std::vector<size_t> best_indices;
      for (size_t i = 0; i < chunks.size(); ++i) {
        if (!have_best || chunks[i].treat > best) {
          best = chunks[i].treat;
          have_best = true;
          best_indices.clear();
          best_indices.push_back(i);
        } else if (chunks[i].treat == best) {
          best_indices.push_back(i);
        }
      }
      if (best_indices.empty()) return;
      const size_t selected =
          best_indices[(best_indices.size() + 1) / 2 - 1];
      const Macs3SummitChunk& chunk = chunks[selected];
      append_row((chunk.start + chunk.end) / 2, chunk);
    };

    if (!call_summits) {
      append_single_summit();
      continue;
    }

    // CallPeakUnit.__close_peak_with_subpeaks: render the treatment pileup
    // base-by-base with 10-bp boundary padding, find smoothed derivative
    // maxima, then apply MACS3's width/flatness peak-shape filter.
    const int64_t padded_end64 = static_cast<int64_t>(r.end) + 10;
    int64_t padded_start64 = static_cast<int64_t>(r.start) - 10;
    int32_t start_boundary = 10;
    if (padded_start64 < 0) {
      start_boundary = static_cast<int32_t>(10 + padded_start64);
      padded_start64 = 0;
    }
    const int64_t signal_size64 = padded_end64 - padded_start64;
    if (signal_size64 <= 0 ||
        static_cast<uint64_t>(signal_size64) >
            static_cast<uint64_t>(std::numeric_limits<size_t>::max())) {
      append_single_summit();
      continue;
    }
    std::vector<float> peakdata(static_cast<size_t>(signal_size64), 0.f);
    std::vector<int32_t> peakindices(static_cast<size_t>(signal_size64), 0);
    for (size_t i = 0; i < chunks.size(); ++i) {
      const int64_t m64 = static_cast<int64_t>(chunks[i].start) -
                          padded_start64 + start_boundary;
      const int64_t n64 = static_cast<int64_t>(chunks[i].end) -
                          padded_start64 + start_boundary;
      const size_t m = static_cast<size_t>(std::max<int64_t>(0, m64));
      const size_t n = static_cast<size_t>(
          std::min<int64_t>(signal_size64, std::max<int64_t>(m64, n64)));
      for (size_t p = m; p < n; ++p) {
        peakdata[p] = chunks[i].treat;
        peakindices[p] = static_cast<int32_t>(i);
      }
    }
    std::vector<int32_t> summit_offsets =
        Macs3SignalMaxima(peakdata, min_length);
    if (!summit_offsets.empty()) {
      const int32_t peak_length = r.end - r.start;
      const auto lo = std::lower_bound(summit_offsets.begin(),
                                       summit_offsets.end(), start_boundary);
      const auto hi = std::upper_bound(
          summit_offsets.begin(), summit_offsets.end(),
          static_cast<int32_t>(peak_length + start_boundary));
      summit_offsets = std::vector<int32_t>(lo, hi);
      summit_offsets = Macs3EnforcePeakyness(peakdata, summit_offsets);
    }
    if (summit_offsets.empty()) {
      append_single_summit();
      continue;
    }
    for (int32_t offset : summit_offsets) {
      if (offset < 0 || static_cast<size_t>(offset) >= peakindices.size()) {
        continue;
      }
      const size_t chunk_index =
          static_cast<size_t>(peakindices[static_cast<size_t>(offset)]);
      if (chunk_index >= chunks.size()) continue;
      const int64_t summit64 = padded_start64 + offset - start_boundary;
      if (summit64 < std::numeric_limits<int32_t>::min() ||
          summit64 > std::numeric_limits<int32_t>::max()) {
        continue;
      }
      append_row(static_cast<int32_t>(summit64), chunks[chunk_index]);
    }
  }
  return true;
}

bool CollectMacs3FragNarrowPeakRowsFromTracks(
    const std::string& chrom_name, const Macs3BdgTrack& ppois,
    const Macs3BdgTrack& treat, const Macs3BdgTrack& lambda, float cutoff,
    int32_t min_length, int32_t max_gap, double pseudocount,
    std::vector<Macs3FragNarrowPeakRow>* rows,
    const std::vector<int32_t>* pos_array_positions,
    bool macs3_float_peak_scores,
    const std::map<double, double>* region_p_to_q,
    bool call_summits) {
  return CollectMacs3FragNarrowPeakRowsImpl(
      chrom_name, ppois, ppois, treat, lambda, cutoff, min_length, max_gap,
      pseudocount, rows, pos_array_positions, macs3_float_peak_scores,
      region_p_to_q, call_summits);
}

bool CollectMacs3FragNarrowPeakRowsFromThresholdTrack(
    const std::string& chrom_name, const Macs3BdgTrack& threshold_score,
    const Macs3BdgTrack& ppois, const Macs3BdgTrack& treat,
    const Macs3BdgTrack& lambda, float cutoff, int32_t min_length,
    int32_t max_gap, double pseudocount,
    std::vector<Macs3FragNarrowPeakRow>* rows,
    const std::vector<int32_t>* pos_array_positions,
    bool macs3_float_peak_scores,
    bool call_summits) {
  return CollectMacs3FragNarrowPeakRowsImpl(
      chrom_name, threshold_score, ppois, treat, lambda, cutoff, min_length,
      max_gap, pseudocount, rows, pos_array_positions,
      macs3_float_peak_scores, /*region_p_to_q=*/nullptr, call_summits);
}

// Build a MACS3-style pqtable: for each unique pscore (-log10 p) value v
// observed across the genome-wide ppois track, q_mlog = v + log10(rank) -
// log10(N) with monotonic non-increasing adjustment as v decreases. rank
// starts at 1 for the highest pscore and increments by the bp coverage of
// each unique pscore. N is the total bp coverage across all chromosomes.
// Matches MACS3.Signal.CallPeakUnit.__cal_pvalue_qvalue_table().
//
// The bp count per segment is (end - start) signed, including any pre-zero
// virtual offsets MACS3 emits for chromosomes shorter than its lambda
// half-window (matches MACS3 behavior, even though it is mildly weird).
void BuildMacs3PqTable(
    const std::vector<Macs3BdgTrack>& ppois_tracks,
    std::map<double, double>* p_to_q,
    bool macs3_float_math) {
  if (p_to_q == nullptr) return;
  p_to_q->clear();
  size_t segment_count = 0;
  for (const Macs3BdgTrack& tr : ppois_tracks) {
    segment_count += std::min(tr.p.size(), tr.v.size());
  }
  std::unordered_map<double, int64_t> pscore_stat;
  // Segment count can be tens of millions while unique score count is only a
  // few million. Reserving for every segment inflated the bucket array by
  // hundreds of MB on full ATAC data; cap the eager reservation and let the
  // table grow only when the unique-key count actually requires it.
  const size_t reserve_cap = static_cast<size_t>(1) << 22;
  pscore_stat.reserve(std::min(segment_count, reserve_cap));
  for (const Macs3BdgTrack& tr : ppois_tracks) {
    int32_t prev = 0;
    for (size_t i = 0; i < tr.p.size() && i < tr.v.size(); ++i) {
      const int64_t bp = static_cast<int64_t>(tr.p[i]) - static_cast<int64_t>(prev);
      // CallPeakUnit stores get_pscore() in a C float and uses that float
      // directly as the pscore_stat key.  Keep the historical round-5 keys
      // for explicit legacy mode, but use the exact promoted float key for
      // callpeak-compatible FRAG, BAMPE, and BEDPE paths.
      const double key = macs3_float_math
                             ? static_cast<double>(static_cast<float>(tr.v[i]))
                             : std::round(static_cast<double>(tr.v[i]) * 1e5) / 1e5;
      pscore_stat[key] += bp;
      prev = tr.p[i];
    }
  }
  if (pscore_stat.empty()) {
    return;
  }
  int64_t N = 0;
  for (const auto& kv : pscore_stat) N += kv.second;
  if (N <= 0) {
    return;
  }
  const double f_double = -std::log10(static_cast<double>(N));
  // MACS3-compatible mode follows CallPeakUnit's C-float f/q assignments.
  // Explicit legacy mode retains the older RapidMACS double arithmetic.
  const float f_float = static_cast<float>(f_double);
  std::vector<std::pair<double, int64_t>> unique_v;
  unique_v.reserve(pscore_stat.size());
  for (const auto& kv : pscore_stat) unique_v.push_back(kv);
  // The accumulator and ordered output contain the same 2M+ keys on large
  // ATAC datasets. Release the hash nodes/buckets before allocating the tree
  // table so their peak lifetimes do not overlap.
  std::unordered_map<double, int64_t>().swap(pscore_stat);
  std::sort(unique_v.begin(), unique_v.end(),
            [](const std::pair<double, int64_t>& a,
               const std::pair<double, int64_t>& b) {
              return a.first > b.first;
            });

  int64_t k = 1;
  double pre_q = std::numeric_limits<double>::infinity();
  auto insert_hint = p_to_q->begin();
  for (const auto& stat : unique_v) {
    const double v = stat.first;
    const int64_t l = stat.second;
    double q = v + std::log10(static_cast<double>(k)) + f_double;
    if (macs3_float_math) {
      q = static_cast<float>(static_cast<float>(v) +
                             (std::log10(static_cast<double>(k)) + f_float));
    }
    if (q > pre_q) q = pre_q;
    if (q <= 0.0) {
      insert_hint = p_to_q->emplace_hint(insert_hint, v, 0.0);
      pre_q = 0.0;
    } else {
      insert_hint = p_to_q->emplace_hint(insert_hint, v, q);
      pre_q = q;
    }
    k += l;
  }
}

// Look up qscore for a peak's p_mlog. Round to 5 digits to match the table
// keys; if exact key missing, fall back to nearest lower pscore (which is
// safe given the monotonic non-increasing q-table).
static double LookupQFromPqTable(const std::map<double, double>& p_to_q,
                                 double p_mlog) {
  if (p_to_q.empty()) return 0.0;
  // MACS3-compatible paths use exact float32 pscore keys, as MACS3's
  // Float32to32Map does. The round-5 retry retains compatibility with an
  // explicit legacy table built by this library.
  const double float_key = static_cast<double>(static_cast<float>(p_mlog));
  auto it = p_to_q.find(float_key);
  if (it != p_to_q.end()) return it->second;
  const double rounded_key = std::round(p_mlog * 1e5) / 1e5;
  it = p_to_q.find(rounded_key);
  if (it != p_to_q.end()) return it->second;
  const double key = float_key;
  auto upper = p_to_q.upper_bound(key);
  if (upper == p_to_q.begin()) return 0.0;
  --upper;
  return upper->second;
}

void BuildMacs3QscoreTracks(
    const std::vector<Macs3BdgTrack>& ppois_tracks,
    const std::map<double, double>& p_to_q,
    std::vector<Macs3BdgTrack>* qscore_tracks) {
  if (qscore_tracks == nullptr) return;
  qscore_tracks->clear();
  qscore_tracks->resize(ppois_tracks.size());
  for (size_t i = 0; i < ppois_tracks.size(); ++i) {
    qscore_tracks->at(i).p = ppois_tracks[i].p;
    qscore_tracks->at(i).v.resize(ppois_tracks[i].v.size());
    for (size_t j = 0; j < ppois_tracks[i].v.size(); ++j) {
      qscore_tracks->at(i).v[j] = static_cast<float>(
          LookupQFromPqTable(p_to_q, ppois_tracks[i].v[j]));
    }
  }
}

bool BuildMacs3QScoreTrackFromPpois(const Macs3BdgTrack& ppois,
                                    const std::map<double, double>& p_to_q,
                                    Macs3BdgTrack* qscore,
                                    bool macs3_float_qtable) {
  if (qscore == nullptr || ppois.p.size() != ppois.v.size()) {
    return false;
  }
  qscore->p.clear();
  qscore->v.clear();
  if (ppois.p.empty()) {
    return true;
  }
  if (p_to_q.empty()) {
    return false;
  }
  // LookupQFromPqTable accepts both exact float32 and legacy round-5 keys.
  // The flag remains part of the compatibility API and documents which
  // table representation the caller constructed.
  (void)macs3_float_qtable;
  qscore->p = ppois.p;
  qscore->v.reserve(ppois.v.size());
  for (float pscore : ppois.v) {
    qscore->v.push_back(
        static_cast<float>(LookupQFromPqTable(p_to_q, pscore)));
  }
  return true;
}

bool WriteMacs3FragNarrowPeakRowsWithScoreCap(
    std::vector<Macs3FragNarrowPeakRow>* rows,
    const std::string& path_narrowpeak, const std::string& path_summits_bed,
    const std::map<double, double>* p_to_q,
    const std::string& name_prefix,
    bool narrowpeak_score_uses_qvalue,
    bool cap_narrowpeak_score,
    bool macs3_float_output) {
  if (rows == nullptr || (path_narrowpeak.empty() && path_summits_bed.empty())) {
    return false;
  }
  if (p_to_q == nullptr || p_to_q->empty()) {
    std::vector<double> p_lin;
    p_lin.reserve(rows->size());
    for (const Macs3FragNarrowPeakRow& e : *rows) {
      p_lin.push_back(
          std::min(1.0, std::max(0.0, std::pow(10.0, -e.p_mlog))));
    }
    std::vector<double> qv;
    BenjaminiHochbergFdr(p_lin, &qv);
    for (size_t i = 0; i < rows->size(); ++i) {
      const double q = std::min(1.0, std::max(0.0, qv[i]));
      double q_mlog = -std::log10(std::max(q, 1e-300));
      if (macs3_float_output) {
        q_mlog = static_cast<float>(q_mlog);
      }
      rows->at(i).q_mlog = q_mlog;
    }
  } else {
    for (Macs3FragNarrowPeakRow& row : *rows) {
      double q_mlog = LookupQFromPqTable(*p_to_q, row.p_mlog);
      if (macs3_float_output) {
        q_mlog = static_cast<float>(q_mlog);
      }
      row.q_mlog = q_mlog;
    }
  }
  std::sort(rows->begin(), rows->end(),
            [](const Macs3FragNarrowPeakRow& a, const Macs3FragNarrowPeakRow& b) {
              if (a.chrom != b.chrom) return a.chrom < b.chrom;
              if (a.start != b.start) return a.start < b.start;
              if (a.end != b.end) return a.end < b.end;
              return a.peak_off < b.peak_off;
            });
  const std::vector<std::string> names =
      Macs3PeakNames(*rows, name_prefix);
  if (!path_narrowpeak.empty()) {
    FILE* fpn = std::fopen(path_narrowpeak.c_str(), "w");
    if (fpn == nullptr) return false;
    for (size_t i = 0; i < rows->size(); ++i) {
      const Macs3FragNarrowPeakRow& p = rows->at(i);
      const double fc = macs3_float_output ? static_cast<float>(p.fc) : p.fc;
      const double p_mlog =
          macs3_float_output ? static_cast<float>(p.p_mlog) : p.p_mlog;
      const double q_mlog =
          macs3_float_output ? static_cast<float>(p.q_mlog) : p.q_mlog;
      const std::string& name = names[i];
      int score = NarrowPeakIntScoreQ(
          narrowpeak_score_uses_qvalue ? q_mlog : p_mlog);
      if (cap_narrowpeak_score && score > 1000) score = 1000;
      if (std::fprintf(fpn, "%s\t%d\t%d\t%s\t%d\t.\t%.6g\t%.6g\t%.6g\t%d\n",
                       p.chrom.c_str(), static_cast<int>(p.start),
                       static_cast<int>(p.end), name.c_str(), score, fc,
                       p_mlog, q_mlog, static_cast<int>(p.peak_off)) < 0) {
        std::fclose(fpn);
        return false;
      }
    }
    if (std::fclose(fpn) != 0) return false;
  }
  if (!path_summits_bed.empty()) {
    FILE* fps = std::fopen(path_summits_bed.c_str(), "w");
    if (fps == nullptr) return false;
    for (size_t i = 0; i < rows->size(); ++i) {
      const Macs3FragNarrowPeakRow& p = rows->at(i);
      const double q_mlog =
          macs3_float_output ? static_cast<float>(p.q_mlog) : p.q_mlog;
      const std::string& name = names[i];
      // MACS3's summits file always reports qscore, independently of which
      // metric supplies narrowPeak column 5.
      if (std::fprintf(fps, "%s\t%d\t%d\t%s\t%.6g\n", p.chrom.c_str(),
                       static_cast<int>(p.summit),
                       static_cast<int>(p.summit + 1), name.c_str(), q_mlog) < 0) {
        std::fclose(fps);
        return false;
      }
    }
    if (std::fclose(fps) != 0) return false;
  }
  return true;
}

bool WriteMacs3FragNarrowPeakRows(
    std::vector<Macs3FragNarrowPeakRow>* rows,
    const std::string& path_narrowpeak, const std::string& path_summits_bed,
    const std::map<double, double>* p_to_q,
    const std::string& name_prefix,
    bool narrowpeak_score_is_qvalue) {
  if (rows == nullptr || (path_narrowpeak.empty() && path_summits_bed.empty())) {
    return false;
  }
  if (p_to_q == nullptr || p_to_q->empty()) {
    // Fallback when no genome-wide ppois track was supplied: BH over the
    // emitted peak rows (legacy behavior; does NOT match MACS3 callpeak).
    std::vector<double> p_lin;
    p_lin.reserve(rows->size());
    for (const Macs3FragNarrowPeakRow& e : *rows) {
      const double p = std::min(1.0, std::max(0.0, std::pow(10.0, -e.p_mlog)));
      p_lin.push_back(p);
    }
    std::vector<double> qv;
    BenjaminiHochbergFdr(p_lin, &qv);
    for (size_t i = 0; i < rows->size(); ++i) {
      const double q = std::min(1.0, std::max(0.0, qv[i]));
      rows->at(i).q_mlog = -std::log10(std::max(q, 1e-300));
    }
  } else {
    for (size_t i = 0; i < rows->size(); ++i) {
      rows->at(i).q_mlog = LookupQFromPqTable(*p_to_q, rows->at(i).p_mlog);
    }
  }
  std::sort(rows->begin(), rows->end(),
            [](const Macs3FragNarrowPeakRow& a, const Macs3FragNarrowPeakRow& b) {
              if (a.chrom != b.chrom) {
                return a.chrom < b.chrom;
              }
              if (a.start != b.start) {
                return a.start < b.start;
              }
              if (a.end != b.end) {
                return a.end < b.end;
              }
              return a.peak_off < b.peak_off;
            });
  const std::vector<std::string> names =
      Macs3PeakNames(*rows, name_prefix);
  if (!path_narrowpeak.empty()) {
    FILE* fpn = std::fopen(path_narrowpeak.c_str(), "w");
    if (fpn == nullptr) {
      return false;
    }
    for (size_t i = 0; i < rows->size(); ++i) {
      const Macs3FragNarrowPeakRow& p = rows->at(i);
      const std::string& name = names[i];
      const int sc = NarrowPeakIntScoreQ(
          narrowpeak_score_is_qvalue ? p.q_mlog : p.p_mlog);
      if (std::fprintf(fpn, "%s\t%d\t%d\t%s\t%d\t.\t%.6g\t%.6g\t%.6g\t%d\n",
                       p.chrom.c_str(), static_cast<int>(p.start),
                       static_cast<int>(p.end), name.c_str(), sc, p.fc,
                       p.p_mlog, p.q_mlog, static_cast<int>(p.peak_off)) < 0) {
        std::fclose(fpn);
        return false;
      }
    }
    if (std::fclose(fpn) != 0) {
      return false;
    }
  }
  if (!path_summits_bed.empty()) {
    FILE* fps = std::fopen(path_summits_bed.c_str(), "w");
    if (fps == nullptr) {
      return false;
    }
    for (size_t i = 0; i < rows->size(); ++i) {
      const Macs3FragNarrowPeakRow& p = rows->at(i);
      const std::string& name = names[i];
      const double summit_score =
          narrowpeak_score_is_qvalue ? p.q_mlog : p.p_mlog;
      if (std::fprintf(fps, "%s\t%d\t%d\t%s\t%.6g\n", p.chrom.c_str(),
                       static_cast<int>(p.summit),
                       static_cast<int>(p.summit + 1), name.c_str(),
                       summit_score) < 0) {
        std::fclose(fps);
        return false;
      }
    }
    if (std::fclose(fps) != 0) {
      return false;
    }
  }
  return true;
}

void FilterMacs3FragNarrowPeakRowsByFoldEnrichment(
    std::vector<Macs3FragNarrowPeakRow>* rows, double cutoff) {
  if (rows == nullptr || !(cutoff > 0.0)) {
    return;
  }
  rows->erase(
      std::remove_if(rows->begin(), rows->end(),
                     [cutoff](const Macs3FragNarrowPeakRow& row) {
                       return row.fc < cutoff;
                     }),
      rows->end());
}

}  // namespace peaks
}  // namespace chromap
