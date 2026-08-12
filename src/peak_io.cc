#include "rapidmacs/peak_io.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <limits>
#include <mutex>
#include <unordered_map>
#include <tuple>
#include <utility>
#include <vector>

#include "rapidmacs/binned_signal.h"

namespace chromap {
namespace peaks {

namespace {

int ScoreFromQ(double q) {
  if (q < 0.0) {
    q = 0.0;
  }
  if (q > 1.0) {
    q = 1.0;
  }
  const double s = -10.0 * std::log10(std::max(q, 1e-300));
  int out = static_cast<int>(s + 0.5);
  if (out < 0) {
    out = 0;
  }
  if (out > 1000) {
    out = 1000;
  }
  return out;
}

}  // namespace

bool WriteNarrowPeak(const std::string& path, const std::vector<Peak>& peaks,
                     const std::string& name_prefix) {
  std::vector<Peak> sorted = peaks;
  std::sort(sorted.begin(), sorted.end(), [](const Peak& a, const Peak& b) {
    if (a.chrom != b.chrom) {
      return a.chrom < b.chrom;
    }
    if (a.start != b.start) {
      return a.start < b.start;
    }
    if (a.end != b.end) {
      return a.end < b.end;
    }
    return a.peak_offset < b.peak_offset;
  });
  for (const Peak& p : sorted) {
    if (p.end <= p.start) {
      return false;
    }
  }
  for (size_t i = 1; i < sorted.size(); ++i) {
    if (sorted[i - 1].chrom == sorted[i].chrom) {
      if (sorted[i - 1].start > sorted[i].start) {
        return false;
      }
    } else if (sorted[i - 1].chrom > sorted[i].chrom) {
      return false;
    }
  }
  FILE* fp = std::fopen(path.c_str(), "w");
  if (fp == nullptr) {
    return false;
  }
  for (size_t i = 0; i < sorted.size(); ++i) {
    const Peak& p = sorted[i];
    const std::string name = name_prefix.empty()
                                 ? ("peak_" + std::to_string(i + 1))
                                 : (name_prefix + "_peak_" + std::to_string(i + 1));
    const int sc = ScoreFromQ(p.q_value);
    const double sfg = static_cast<double>(p.max_signal);
    const double pv = std::min(1.0, std::max(0.0, p.p_value));
    const double qv = std::min(1.0, std::max(0.0, p.q_value));
    const double logp = -std::log10(std::max(pv, 1e-300));
    const double logq = -std::log10(std::max(qv, 1e-300));
    if (std::fprintf(fp, "%s\t%d\t%d\t%s\t%d\t.\t%g\t%g\t%g\t%d\n", p.chrom.c_str(),
                     p.start, p.end, name.c_str(), sc, sfg, logp, logq,
                     p.peak_offset) < 0) {
      std::fclose(fp);
      return false;
    }
  }
  if (std::fclose(fp) != 0) {
    return false;
  }
  return true;
}

bool WriteRunSummaryTsv(const std::string& path,
                        const std::vector<Peak>& all_peaks,
                        const std::string& argv_line,
                        const std::string& program_params,
                        int64_t total_bins, int64_t saturated_bins) {
  std::ofstream o(path);
  if (!o) {
    return false;
  }
  o << "metric\tvalue\n";
  o << "peak_count\t" << all_peaks.size() << "\n";
  int64_t tot = 0;
  std::vector<int32_t> widths;
  for (const Peak& p : all_peaks) {
    if (p.end > p.start) {
      const int w = p.end - p.start;
      tot += w;
      widths.push_back(w);
    }
  }
  o << "total_peak_bp\t" << tot << "\n";
  if (widths.empty()) {
    o << "median_peak_width_bp\t0\n";
  } else {
    std::sort(widths.begin(), widths.end());
    const size_t m = widths.size() / 2;
    o << "median_peak_width_bp\t" << widths[m] << "\n";
  }
  if (saturated_bins >= 0) {
    o << "total_bins\t" << total_bins << "\n";
    o << "saturated_bins\t" << saturated_bins << "\n";
    const double frac =
        total_bins > 0
            ? static_cast<double>(saturated_bins) / static_cast<double>(total_bins)
            : 0.0;
    o << "saturated_bin_fraction\t" << frac << "\n";
  }
  o << "argv\t" << argv_line << "\n";
  o << "params\t" << program_params << "\n";
  o.close();
  if (!o) {
    return false;
  }
  return true;
}

bool WritePileupBedGraph(const std::string& path,
                       const std::vector<ChromFragments>& by_chrom,
                       int32_t bin_width_bp, int32_t ext_size_bp) {
  if (path.empty() || bin_width_bp < 1 || ext_size_bp < 2 ||
      (ext_size_bp % 2) != 0) {
    return false;
  }
  std::vector<const ChromFragments*> order;
  order.reserve(by_chrom.size());
  for (const ChromFragments& ch : by_chrom) {
    order.push_back(&ch);
  }
  std::sort(order.begin(), order.end(), [](const ChromFragments* a,
                                            const ChromFragments* b) {
    return a->name < b->name;
  });
  FILE* fp = std::fopen(path.c_str(), "w");
  if (fp == nullptr) {
    return false;
  }
  for (const ChromFragments* chp : order) {
    const ChromFragments& ch = *chp;
    std::vector<uint16_t> obs;
    int32_t pad = 0;
    if (!BuildBinnedCutSignal(ch, bin_width_bp, ext_size_bp, &obs, &pad)) {
      continue;
    }
    size_t i = 0;
    while (i < obs.size()) {
      if (obs[i] == 0) {
        ++i;
        continue;
      }
      const int64_t v = obs[i];
      size_t j = i + 1;
      while (j < obs.size() && obs[j] == v) {
        ++j;
      }
      const int64_t start = static_cast<int64_t>(i) * bin_width_bp;
      const int64_t end = static_cast<int64_t>(j) * bin_width_bp;
      if (end > start &&
          std::fprintf(fp, "%s\t%lld\t%lld\t%lld\n", ch.name.c_str(),
                       static_cast<long long>(start),
                       static_cast<long long>(end),
                       static_cast<long long>(v)) < 0) {
        std::fclose(fp);
        return false;
      }
      i = j;
    }
  }
  if (std::fclose(fp) != 0) {
    return false;
  }
  return true;
}

namespace macs3_pois_detail {
double LogspaceAddMacs3(double logx, double logy) {
  if (logx > logy) {
    return logx + std::log1p(std::exp(logy - logx));
  }
  return logy + std::log1p(std::exp(logx - logy));
}

double LogFactorial(uint32_t n) {
  static std::vector<double> cache(1, 0.0);  // cache[n] = log(n!)
  static std::mutex mu;
  std::lock_guard<std::mutex> lock(mu);
  while (cache.size() <= n) {
    const size_t i = cache.size();
    cache.push_back(cache.back() + std::log(static_cast<double>(i)));
  }
  return cache[n];
}

// MACS3.Signal.Prob.log10_poisson_cdf_Q_large_lambda(k, lam) -> log10(upper tail P).
double Macs3Log10PoissonUpperQ(uint32_t k, double lam) {
  const double ln_lbd = std::log(lam);
  int32_t m = static_cast<int32_t>(k) + 1;
  const double sum_ln_m = LogFactorial(static_cast<uint32_t>(m));
  double logx = m * ln_lbd - sum_ln_m;
  double residue = logx;
  while (true) {
    m += 1;
    const double logy = logx + ln_lbd - std::log(static_cast<double>(m));
    const double pre_residue = residue;
    residue = LogspaceAddMacs3(pre_residue, logy);
    if (std::fabs(pre_residue - residue) < 1e-5) {
      break;
    }
    logx = logy;
  }
  constexpr double kLn10 = 2.30258509299404568402;
  return std::round((residue - lam) / kLn10 * 1e5) / 1e5;
}
}  // namespace macs3_pois_detail

double Macs3BdgcmpPpoisNegLog10P(float treat_f, float ctrl_f, double pseudocount) {
  const float ps = static_cast<float>(pseudocount);
  const int obs = static_cast<int>(treat_f + ps);
  const float lam_f = ctrl_f + ps;
  if (!(lam_f > 0.0f)) {
    return 0.0;
  }
  const double log10p =
      macs3_pois_detail::Macs3Log10PoissonUpperQ(
          static_cast<uint32_t>(std::max(0, obs)), static_cast<double>(lam_f));
  return -log10p;
}

namespace {

double EffectiveFragWeight(int32_t raw_count, int32_t max_count,
                           bool macs3_uint8_counts) {
  int64_t w = raw_count;
  if (max_count > 0) {
    w = std::min<int64_t>(w, max_count);
  }
  if (w <= 0) {
    return 0.0;
  }
  if (!macs3_uint8_counts) {
    return static_cast<double>(w);
  }
  const auto u = static_cast<uint8_t>(static_cast<unsigned int>(w));
  return static_cast<double>(u);
}

// Integer count per row after max_count/uint8 (for MACS3 `total` / `length`)
int32_t RowCountForMcsTotals(int32_t raw_count, int32_t max_count,
                             bool macs3_uint8) {
  int64_t t = raw_count;
  if (max_count > 0) {
    t = std::min<int64_t>(t, max_count);
  }
  if (t <= 0) {
    return 0;
  }
  if (macs3_uint8) {
    return static_cast<int32_t>(static_cast<uint8_t>(static_cast<unsigned int>(t)));
  }
  if (t > std::numeric_limits<int32_t>::max()) {
    return std::numeric_limits<int32_t>::max();
  }
  return static_cast<int32_t>(t);
}

// Raw fragment-span [start,end) segments with pileup `z` (before scaling).
void MergedDeltasForFragSpan(const std::vector<Fragment>& frags,
                            int32_t max_count, bool macs3_uint8,
                            std::vector<std::pair<int32_t, double>>* events) {
  events->clear();
  for (const Fragment& f : frags) {
    const double w = EffectiveFragWeight(f.count, max_count, macs3_uint8);
    if (w == 0.0) {
      continue;
    }
    events->push_back(std::make_pair(f.start, w));
    events->push_back(std::make_pair(f.end, -w));
  }
}

// MACS3 PETrackII `pileup_a_chromosome_c` for a single d: end-centered windows
// at left (l) and right (r) cut sites.
void MergedDeltasForEndWindows(const std::vector<Fragment>& frags, int32_t d_bp,
                             int32_t max_count, bool macs3_uint8,
                             std::vector<std::pair<int32_t, double>>* events) {
  events->clear();
  if (d_bp < 1) {
    return;
  }
  const int32_t half = d_bp / 2;
  for (const Fragment& f : frags) {
    const double w = EffectiveFragWeight(f.count, max_count, macs3_uint8);
    if (w == 0.0) {
      continue;
    }
    const int32_t ll = f.start - half;
    const int32_t rl = ll + d_bp;
    const int32_t lr = f.end - half;
    const int32_t rr = lr + d_bp;
    events->push_back(std::make_pair(ll, w));
    events->push_back(std::make_pair(rl, -w));
    events->push_back(std::make_pair(lr, w));
    events->push_back(std::make_pair(rr, -w));
  }
}

void CoalesceEventPositions(std::vector<std::pair<int32_t, double>>* ev) {
  if (ev->empty()) {
    return;
  }
  std::sort(ev->begin(), ev->end(),
            [](const std::pair<int32_t, double>& a, const std::pair<int32_t, double>& b) {
              if (a.first != b.first) {
                return a.first < b.first;
              }
              return a.second < b.second;
            });
  size_t wri = 0;
  for (size_t i = 0; i < ev->size(); ++i) {
    if (wri > 0 && ev->at(i).first == ev->at(wri - 1).first) {
      ev->at(wri - 1).second += ev->at(i).second;
    } else {
      ev->at(wri++) = ev->at(i);
    }
  }
  ev->resize(wri);
}

// Half-open [start,end) segments with value `z` (pileup_PV from MACS3 PileupV2).
void PileupSegmentsFromEvents(const std::vector<std::pair<int32_t, double>>& ev,
                            std::vector<std::tuple<int32_t, int32_t, double>>* out) {
  out->clear();
  if (ev.empty()) {
    return;
  }
  double z = 0.0;
  double pre_z = -1.0e4;
  int32_t s = 0;
  std::vector<std::pair<int32_t, double>> pv;
  pv.reserve(ev.size() + 1);
  for (const auto& e : ev) {
    const int32_t pos = e.first;
    const double v = e.second;
    if (pos != s) {
      if (z == pre_z && !pv.empty()) {
        pv.back().first = pos;
      } else {
        pv.push_back(std::make_pair(pos, z));
        pre_z = z;
      }
    }
    z += v;
    s = pos;
  }
  int32_t pre = 0;
  for (const auto& seg : pv) {
    out->push_back(std::make_tuple(pre, seg.first, seg.second));
    pre = seg.first;
  }
}

// From half-open pileup segments, build MACS3 (p, v) where p[i] is the i-th run end
// and v[i] is the value on (0, p[0]) for i=0, else (p[i-1], p[i]).
void PileupToPvArray(const std::vector<std::tuple<int32_t, int32_t, double>>& segs,
                     std::vector<int32_t>* p_out, std::vector<float>* v_out) {
  p_out->clear();
  v_out->clear();
  p_out->reserve(segs.size());
  v_out->reserve(segs.size());
  for (const std::tuple<int32_t, int32_t, double>& t : segs) {
    p_out->push_back(std::get<1>(t));
    v_out->push_back(static_cast<float>(std::get<2>(t)));
  }
}

void PileupToPArrayOnly(const std::vector<std::tuple<int32_t, int32_t, double>>& segs,
                        std::vector<int32_t>* p_out) {
  p_out->clear();
  p_out->reserve(segs.size());
  for (const std::tuple<int32_t, int32_t, double>& t : segs) {
    p_out->push_back(std::get<1>(t));
  }
}

// MACS3 `CallPeakUnit.__chrom_pair_treat_ctrl` (no tail: longer array is truncated,
// matching the published Cython). Only ret_p/ret_c are needed for `control_lambda.bdg`.
void ChromPairTreatCtrl(const std::vector<int32_t>& t_p, const std::vector<int32_t>& c_p,
                        const std::vector<float>& c_v, std::vector<int32_t>* ret_p,
                        std::vector<float>* ret_c) {
  ret_p->clear();
  ret_c->clear();
  size_t it = 0;
  size_t ic = 0;
  const size_t lt = t_p.size();
  const size_t lc = c_p.size();
  ret_p->reserve(lt + lc);
  ret_c->reserve(lt + lc);
  while (it < lt && ic < lc) {
    if (t_p[it] < c_p[ic]) {
      ret_p->push_back(t_p[it]);
      ret_c->push_back(c_v[ic]);
      it++;
    } else if (t_p[it] > c_p[ic]) {
      ret_p->push_back(c_p[ic]);
      ret_c->push_back(c_v[ic]);
      ic++;
    } else {
      ret_p->push_back(t_p[it]);
      ret_c->push_back(c_v[ic]);
      it++;
      ic++;
    }
  }
}

// `control_lambda.bdg` in MACS3 `__write_bedGraph_for_a_chromosome` for control: merge
// consecutive run boundaries when |Δv|>1e-5 (5-digit precision) like CallPeakUnit.py.
bool Mcs3EmitControlBdgFromSteps(FILE* fp, const char* chrom,
                                const std::vector<int32_t>& pos_r_end,
                                const std::vector<float>& c) {
  const int l = static_cast<int>(pos_r_end.size());
  if (l == 0) {
    return true;
  }
  int32_t pre_p = 0;
  float pre_v = c[0];
  if (l == 1) {
    const int32_t p = pos_r_end[0];
    return std::fprintf(fp, "%s\t%d\t%d\t%.5f\n", chrom, pre_p, p, static_cast<double>(pre_v)) >= 0;
  }
  for (int i = 1; i < l; i++) {
    const int32_t p = pos_r_end[static_cast<size_t>(i - 1)];
    const float v = c[static_cast<size_t>(i)];
    if (std::fabs((double)pre_v - (double)v) > 1e-5) {
      if (std::fprintf(fp, "%s\t%d\t%d\t%.5f\n", chrom, pre_p, p,
                       static_cast<double>(pre_v)) < 0) {
        return false;
      }
      pre_p = p;
      pre_v = v;
    }
  }
  const int32_t p = pos_r_end[static_cast<size_t>(l - 1)];
  return std::fprintf(fp, "%s\t%d\t%d\t%.5f\n", chrom, pre_p, p, static_cast<double>(pre_v)) >= 0;
}

inline double Round5Score(double x) {
  return std::round(x * 1e5) / 1e5;
}

double Macs3BdgcmpFeScore(float treat_f, float ctrl_f, float pseudocount_f) {
  const float lam_f = ctrl_f + pseudocount_f;
  if (!(lam_f > 0.0f)) {
    return 0.0;
  }
  const float fe = (treat_f + pseudocount_f) / lam_f;
  return static_cast<double>(fe);
}

struct ScoreRunMerge {
  bool has = false;
  int32_t lo = 0;
  int32_t hi = 0;
  double val = 0.0;
};

void ScoreRunPush(FILE* fp, const char* chrom, int32_t seg_lo, int32_t seg_hi,
                  double val, ScoreRunMerge* st) {
  if (fp == nullptr || seg_hi <= seg_lo) {
    return;
  }
  if (!st->has) {
    st->has = true;
    st->lo = seg_lo;
    st->hi = seg_hi;
    st->val = val;
    return;
  }
  if (std::fabs(val - st->val) <= 1e-5 && seg_lo == st->hi) {
    st->hi = seg_hi;
    return;
  }
  std::fprintf(fp, "%s\t%d\t%d\t%.5f\n", chrom, st->lo, st->hi, st->val);
  st->lo = seg_lo;
  st->hi = seg_hi;
  st->val = val;
}

void ScoreRunFlush(FILE* fp, const char* chrom, ScoreRunMerge* st) {
  if (fp == nullptr || !st->has) {
    return;
  }
  std::fprintf(fp, "%s\t%d\t%d\t%.5f\n", chrom, st->lo, st->hi, st->val);
  st->has = false;
}

// Pairwise min-merge of treat pileup vs control_lambda breakpoints; stop when either
// p-array is exhausted (same as MACS3 BedGraphTrackI.make_ScoreTrackII_for_macs).
void EmitMacs3BdgcmpScoresForChrom(FILE* fp_p, FILE* fp_f, const char* chrom,
                                   const std::vector<int32_t>& t_p,
                                   const std::vector<float>& t_v,
                                   const std::vector<int32_t>& m_p,
                                   const std::vector<float>& m_c, double pseudocount) {
  if (t_p.empty() || m_p.empty()) {
    return;
  }
  const float pseudocount_f = static_cast<float>(pseudocount);
  size_t i1 = 0;
  size_t i2 = 0;
  ScoreRunMerge st_p;
  ScoreRunMerge st_f;
  int32_t prev = 0;
  while (i1 < t_p.size() && i2 < m_p.size()) {
    const int32_t p1 = t_p[i1];
    const int32_t p2 = m_p[i2];
    const int32_t nx = std::min(p1, p2);
    if (nx > prev) {
      const int32_t lo = std::max(prev, static_cast<int32_t>(0));
      const int32_t hi = std::max(nx, static_cast<int32_t>(0));
      if (hi > lo) {
        const float tf = t_v[i1];
        const float cf = m_c[i2];
        const double sp =
            Round5Score(Macs3BdgcmpPpoisNegLog10P(tf, cf, static_cast<double>(pseudocount_f)));
        const double sf = Round5Score(Macs3BdgcmpFeScore(tf, cf, pseudocount_f));
        ScoreRunPush(fp_p, chrom, lo, hi, sp, &st_p);
        ScoreRunPush(fp_f, chrom, lo, hi, sf, &st_f);
      }
    }
    if (p1 <= p2) {
      ++i1;
    }
    if (p2 <= p1) {
      ++i2;
    }
    prev = nx;
  }
  ScoreRunFlush(fp_p, chrom, &st_p);
  ScoreRunFlush(fp_f, chrom, &st_f);
}

}  // namespace

bool WriteFragSpanPileupBedGraph(const std::string& path,
                                 const std::vector<ChromFragments>& by_chrom,
                                 int32_t max_count, bool macs3_uint8_counts) {
  if (path.empty()) {
    return false;
  }
  if (max_count < 0) {
    return false;
  }
  std::vector<const ChromFragments*> order;
  order.reserve(by_chrom.size());
  for (const ChromFragments& ch : by_chrom) {
    if (!ch.frags.empty()) {
      order.push_back(&ch);
    }
  }
  std::sort(order.begin(), order.end(), [](const ChromFragments* a,
                                            const ChromFragments* b) {
    return a->name < b->name;
  });
  FILE* fp = std::fopen(path.c_str(), "w");
  if (fp == nullptr) {
    return false;
  }
  std::vector<std::pair<int32_t, double>> events;
  std::vector<std::tuple<int32_t, int32_t, double>> segs;
  events.reserve(1024);
  for (const ChromFragments* chp : order) {
    const ChromFragments& ch = *chp;
    MergedDeltasForFragSpan(ch.frags, max_count, macs3_uint8_counts, &events);
    CoalesceEventPositions(&events);
    if (events.empty()) {
      continue;
    }
    PileupSegmentsFromEvents(events, &segs);
    for (const std::tuple<int32_t, int32_t, double>& t : segs) {
      if (std::fprintf(fp, "%s\t%d\t%d\t%.5f\n", ch.name.c_str(),
                       std::get<0>(t), std::get<1>(t), std::get<2>(t)) < 0) {
        std::fclose(fp);
        return false;
      }
    }
  }
  if (std::fclose(fp) != 0) {
    return false;
  }
  return true;
}

bool WriteFragMacs3NoControlLambdaBedGraph(
    const std::string& path, const std::vector<ChromFragments>& by_chrom,
    int32_t max_count, bool macs3_uint8_counts, int64_t effective_genome_size,
    int32_t llocal_bp) {
  if (path.empty() || max_count < 0 || llocal_bp < 1 || effective_genome_size < 1) {
    return false;
  }
  int64_t total_len = 0;
  int64_t tot_c = 0;
  for (const ChromFragments& ch : by_chrom) {
    for (const Fragment& f : ch.frags) {
      const int32_t c = RowCountForMcsTotals(f.count, max_count, macs3_uint8_counts);
      if (c <= 0) {
        continue;
      }
      tot_c += c;
      total_len += static_cast<int64_t>(f.end - f.start) * c;
    }
  }
  if (tot_c == 0 || total_len == 0) {
    return false;
  }
  const double lambda_bg = static_cast<double>(total_len) / static_cast<double>(effective_genome_size);
  const double scale = static_cast<double>(total_len) /
                       (static_cast<double>(llocal_bp) * static_cast<double>(tot_c) * 2.0);

  std::vector<const ChromFragments*> order;
  order.reserve(by_chrom.size());
  for (const ChromFragments& ch : by_chrom) {
    if (!ch.frags.empty()) {
      order.push_back(&ch);
    }
  }
  std::sort(order.begin(), order.end(), [](const ChromFragments* a, const ChromFragments* b) {
    return a->name < b->name;
  });
  FILE* fp = std::fopen(path.c_str(), "w");
  if (fp == nullptr) {
    return false;
  }
  std::vector<std::pair<int32_t, double>> tevents;
  std::vector<std::pair<int32_t, double>> cevents;
  std::vector<std::tuple<int32_t, int32_t, double>> treat_segs;
  std::vector<std::tuple<int32_t, int32_t, double>> ctrl_segs;
  for (const ChromFragments* chp : order) {
    const ChromFragments& ch = *chp;
    MergedDeltasForFragSpan(ch.frags, max_count, macs3_uint8_counts, &tevents);
    CoalesceEventPositions(&tevents);
    treat_segs.clear();
    if (!tevents.empty()) {
      PileupSegmentsFromEvents(tevents, &treat_segs);
    }
    MergedDeltasForEndWindows(ch.frags, llocal_bp, max_count, macs3_uint8_counts, &cevents);
    CoalesceEventPositions(&cevents);
    ctrl_segs.clear();
    if (!cevents.empty()) {
      PileupSegmentsFromEvents(cevents, &ctrl_segs);
      for (std::tuple<int32_t, int32_t, double>& t : ctrl_segs) {
        const double z = std::get<2>(t);
        double v = z * scale;
        if (v < lambda_bg) {
          v = lambda_bg;
        }
        std::get<2>(t) = v;
      }
    }
    std::vector<int32_t> t_p;
    std::vector<int32_t> c_p;
    std::vector<float> c_v_f;
    PileupToPArrayOnly(treat_segs, &t_p);
    PileupToPvArray(ctrl_segs, &c_p, &c_v_f);
    if (t_p.empty() || c_p.empty()) {
      continue;
    }
    std::vector<int32_t> m_p;
    std::vector<float> m_c;
    ChromPairTreatCtrl(t_p, c_p, c_v_f, &m_p, &m_c);
    if (m_p.empty()) {
      continue;
    }
    if (!Mcs3EmitControlBdgFromSteps(fp, ch.name.c_str(), m_p, m_c)) {
      std::fclose(fp);
      return false;
    }
  }
  if (std::fclose(fp) != 0) {
    return false;
  }
  return true;
}

bool ParseBdgTrackEndsVals(
    const std::string& path,
    std::map<std::string, std::pair<std::vector<int32_t>, std::vector<float>>>* out) {
  std::ifstream in(path.c_str());
  if (!in) {
    return false;
  }
  out->clear();
  std::unordered_map<std::string, int32_t> next_start;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    if (line.size() >= 6 && line.compare(0, 6, "track ") == 0) {
      continue;
    }
    std::istringstream ls(line);
    std::string chrom;
    int32_t s = 0;
    int32_t e = 0;
    std::string vs;
    if (!(ls >> chrom >> s >> e >> vs)) {
      return false;
    }
    const float v = std::strtof(vs.c_str(), nullptr);
    int32_t exp = 0;
    const auto it = next_start.find(chrom);
    if (it != next_start.end()) {
      exp = it->second;
    }
    if (s != exp) {
      return false;
    }
    auto& pr = (*out)[chrom];
    pr.first.push_back(e);
    pr.second.push_back(v);
    next_start[chrom] = e;
  }
  return true;
}

bool WriteFragMacs3BdgcmpScoreBedGraphs(
    const std::string& path_ppois, const std::string& path_fe,
    const std::vector<ChromFragments>& by_chrom, int32_t max_count,
    bool macs3_uint8_counts, int64_t effective_genome_size, int32_t llocal_bp,
    double pseudocount) {
  if (path_ppois.empty() && path_fe.empty()) {
    return false;
  }
  if (max_count < 0 || llocal_bp < 1 || effective_genome_size < 1) {
    return false;
  }
  const std::string base = path_ppois.empty() ? path_fe : path_ppois;
  const std::string scr_treat = base + ".scr_treat.bdg";
  const std::string scr_lam = base + ".scr_lambda.bdg";
  if (!WriteFragSpanPileupBedGraph(scr_treat, by_chrom, max_count, macs3_uint8_counts)) {
    return false;
  }
  if (!WriteFragMacs3NoControlLambdaBedGraph(scr_lam, by_chrom, max_count, macs3_uint8_counts,
                                             effective_genome_size, llocal_bp)) {
    std::remove(scr_treat.c_str());
    return false;
  }
  std::map<std::string, std::pair<std::vector<int32_t>, std::vector<float>>> treat_tr;
  std::map<std::string, std::pair<std::vector<int32_t>, std::vector<float>>> ctrl_tr;
  if (!ParseBdgTrackEndsVals(scr_treat, &treat_tr) ||
      !ParseBdgTrackEndsVals(scr_lam, &ctrl_tr)) {
    std::remove(scr_treat.c_str());
    std::remove(scr_lam.c_str());
    return false;
  }
  std::remove(scr_treat.c_str());
  std::remove(scr_lam.c_str());

  FILE* fp_p = nullptr;
  FILE* fp_f = nullptr;
  if (!path_ppois.empty()) {
    fp_p = std::fopen(path_ppois.c_str(), "w");
    if (fp_p == nullptr) {
      return false;
    }
  }
  if (!path_fe.empty()) {
    fp_f = std::fopen(path_fe.c_str(), "w");
    if (fp_f == nullptr) {
      if (fp_p != nullptr) {
        std::fclose(fp_p);
      }
      return false;
    }
  }
  for (const auto& chkv : treat_tr) {
    const std::string& ch = chkv.first;
    const auto itc = ctrl_tr.find(ch);
    if (itc == ctrl_tr.end()) {
      continue;
    }
    const std::vector<int32_t>& t_p = chkv.second.first;
    const std::vector<float>& t_v = chkv.second.second;
    const std::vector<int32_t>& c_p = itc->second.first;
    const std::vector<float>& c_v = itc->second.second;
    if (t_p.empty() || c_p.empty()) {
      continue;
    }
    EmitMacs3BdgcmpScoresForChrom(fp_p, fp_f, ch.c_str(), t_p, t_v, c_p, c_v, pseudocount);
  }
  if (fp_p != nullptr && std::fclose(fp_p) != 0) {
    if (fp_f != nullptr) {
      std::fclose(fp_f);
    }
    return false;
  }
  if (fp_f != nullptr && std::fclose(fp_f) != 0) {
    return false;
  }
  return true;
}

bool WriteMacs3BdgcmpScoreBedGraphsFromTreatLambdaBdgs(
    const std::string& path_ppois, const std::string& path_fe,
    const std::string& path_treat_bdg, const std::string& path_lambda_bdg,
    double pseudocount) {
  if (path_ppois.empty() && path_fe.empty()) {
    return false;
  }
  if (path_treat_bdg.empty() || path_lambda_bdg.empty()) {
    return false;
  }
  std::map<std::string, std::pair<std::vector<int32_t>, std::vector<float>>> treat_tr;
  std::map<std::string, std::pair<std::vector<int32_t>, std::vector<float>>> ctrl_tr;
  if (!ParseBdgTrackEndsVals(path_treat_bdg, &treat_tr) ||
      !ParseBdgTrackEndsVals(path_lambda_bdg, &ctrl_tr)) {
    return false;
  }
  FILE* fp_p = nullptr;
  FILE* fp_f = nullptr;
  if (!path_ppois.empty()) {
    fp_p = std::fopen(path_ppois.c_str(), "w");
    if (fp_p == nullptr) {
      return false;
    }
  }
  if (!path_fe.empty()) {
    fp_f = std::fopen(path_fe.c_str(), "w");
    if (fp_f == nullptr) {
      if (fp_p != nullptr) {
        std::fclose(fp_p);
      }
      return false;
    }
  }
  for (const auto& chkv : treat_tr) {
    const std::string& ch = chkv.first;
    const auto itc = ctrl_tr.find(ch);
    if (itc == ctrl_tr.end()) {
      continue;
    }
    const std::vector<int32_t>& t_p = chkv.second.first;
    const std::vector<float>& t_v = chkv.second.second;
    const std::vector<int32_t>& c_p = itc->second.first;
    const std::vector<float>& c_v = itc->second.second;
    if (t_p.empty() || c_p.empty()) {
      continue;
    }
    EmitMacs3BdgcmpScoresForChrom(fp_p, fp_f, ch.c_str(), t_p, t_v, c_p, c_v,
                                  pseudocount);
  }
  if (fp_p != nullptr && std::fclose(fp_p) != 0) {
    if (fp_f != nullptr) {
      std::fclose(fp_f);
    }
    return false;
  }
  if (fp_f != nullptr && std::fclose(fp_f) != 0) {
    return false;
  }
  return true;
}

}  // namespace peaks
}  // namespace chromap
