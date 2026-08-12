#include "rapidmacs/call_peaks.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <vector>

#include "htslib/kfunc.h"

namespace chromap {
namespace peaks {

namespace {

const double kEps = 1e-12;

// Streaming local-lambda over a per-chrom obs[] array. Computes At(i) in
// amortized O(1) by maintaining a sliding window sum across monotonically
// non-decreasing i. Replaces the materialized lambda[n] array used by the
// previous implementation, saving ~8 bytes per bin (~250 MB on hg38 chr1
// transient under OMP parallelism).
class SlidingLocalLambda {
 public:
  SlidingLocalLambda(const std::vector<uint16_t>& obs, int32_t window)
      : obs_(obs), window_(window), cur_lo_(0), cur_hi_(-1), sum_(0) {
    int64_t s = 0;
    for (int64_t v : obs) s += v;
    global_ = obs.empty() ? 0.0 : static_cast<double>(s) / static_cast<double>(obs.size());
  }

  double Global() const { return global_; }

  // Caller must invoke with i monotonically non-decreasing.
  double At(int64_t i) {
    const int64_t n = static_cast<int64_t>(obs_.size());
    int64_t want_lo = i - window_;
    if (want_lo < 0) want_lo = 0;
    int64_t want_hi = i + window_;
    if (want_hi > n - 1) want_hi = n - 1;
    while (cur_hi_ < want_hi) {
      ++cur_hi_;
      sum_ += obs_[static_cast<size_t>(cur_hi_)];
    }
    while (cur_lo_ < want_lo) {
      sum_ -= obs_[static_cast<size_t>(cur_lo_)];
      ++cur_lo_;
    }
    const int64_t cnt = want_hi - want_lo;
    const int64_t sum_wo = sum_ - obs_[static_cast<size_t>(i)];
    double local = global_;
    if (cnt > 0) {
      local = std::max(local, static_cast<double>(sum_wo) / static_cast<double>(cnt));
    }
    if (local < kEps) local = kEps;
    return std::max(global_, local);
  }

 private:
  const std::vector<uint16_t>& obs_;
  int32_t window_;
  int64_t cur_lo_;
  int64_t cur_hi_;
  int64_t sum_;
  double global_;
};

}  // namespace

double PoissonSurvivalUpper(int k, double lambda) {
  if (k <= 0) {
    return 1.0;
  }
  if (lambda < kEps) {
    return k == 0 ? 1.0 : 0.0;
  }
  // P(Poisson(λ) >= k) = P(k, λ) regularized lower incomplete = kf_gammap(k, λ)
  // Verified against R ppois(k-1, lambda, lower.tail=FALSE) for k>=1.
  double p = kf_gammap(static_cast<double>(k), lambda);
  if (p < 0.0) {
    p = 0.0;
  }
  if (p > 1.0) {
    p = 1.0;
  }
  return p;
}

void BenjaminiHochbergFdr(const std::vector<double>& p_values,
                          std::vector<double>* q_values) {
  const size_t m = p_values.size();
  q_values->assign(m, 0.0);
  if (m == 0) {
    return;
  }
  std::vector<size_t> order(m);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    if (p_values[a] == p_values[b]) {
      return a < b;
    }
    return p_values[a] < p_values[b];
  });
  std::vector<double> q(m);
  const double dmin = m + 0.0;
  for (int64_t k = static_cast<int64_t>(m) - 1; k >= 0; --k) {
    const size_t i = order[static_cast<size_t>(k)];
    const double raw = p_values[i] * dmin / static_cast<double>(k + 1);
    if (k + 1 < static_cast<int64_t>(m)) {
      const size_t nxt = order[static_cast<size_t>(k + 1)];
      q[static_cast<size_t>(i)] = std::min(std::min(raw, 1.0), q[nxt]);
    } else {
      q[static_cast<size_t>(i)] = std::min(raw, 1.0);
    }
  }
  for (size_t i = 0; i < m; ++i) {
    (*q_values)[i] = std::min(std::max(q[i], 0.0), 1.0);
  }
}

void ComputeLocalBackground(const std::vector<uint16_t>& obs,
                            int32_t local_window_bins, double* global_per_bin,
                            std::vector<double>* lambda) {
  const int64_t n = static_cast<int64_t>(obs.size());
  lambda->assign(static_cast<size_t>(n), 0.0);
  *global_per_bin = 0.0;
  if (n < 1) {
    return;
  }
  int64_t s_all = 0;
  for (int64_t v : obs) {
    s_all += v;
  }
  *global_per_bin = static_cast<double>(s_all) / static_cast<double>(n);
  if (local_window_bins < 1) {
    for (int64_t i = 0; i < n; ++i) {
      (*lambda)[static_cast<size_t>(i)] =
          std::max(*global_per_bin, kEps);
    }
    return;
  }
  std::vector<int64_t> pref(static_cast<size_t>(n) + 1, 0);
  for (int64_t i = 0; i < n; ++i) {
    pref[static_cast<size_t>(i) + 1] =
        pref[static_cast<size_t>(i)] + obs[static_cast<size_t>(i)];
  }
  const int64_t w = static_cast<int64_t>(local_window_bins);
  for (int64_t i = 0; i < n; ++i) {
    int64_t lo = i - w;
    int64_t hi = i + w;
    if (lo < 0) {
      lo = 0;
    }
    if (hi >= n) {
      hi = n - 1;
    }
    int64_t sumw = pref[static_cast<size_t>(hi) + 1] - pref[static_cast<size_t>(lo)];
    int64_t sum_wo = sumw - obs[static_cast<size_t>(i)];
    int64_t cnt = (hi - lo);  // number of j != i in [lo,hi]
    double local = *global_per_bin;
    if (cnt > 0) {
      local = std::max(local, static_cast<double>(sum_wo) / static_cast<double>(cnt));
    }
    if (local < kEps) {
      local = kEps;
    }
    (*lambda)[static_cast<size_t>(i)] =
        std::max(*global_per_bin, local);
  }
}

int CallPeaksOnBins(const std::string& chrom, int32_t bin_width,
                    const std::vector<uint16_t>& obs, int32_t local_window_bins,
                    double fdr, double p_cutoff, int32_t merge_bin_gap,
                    std::vector<Peak>* out_peaks) {
  out_peaks->clear();
  const int64_t n = static_cast<int64_t>(obs.size());
  if (n < 1) {
    return 0;
  }

  // Pass 1: streaming lambda + sparse p-value computation. Only bins with
  // obs >= 1 contribute to BH (matches MACS3 behavior of avoiding m-inflation
  // by zero-count bins). p_sub / sub_ix are parallel; sub_ix is monotonically
  // increasing so we can binary-search it later for summit p/q lookup.
  SlidingLocalLambda lam(obs, std::max<int32_t>(0, local_window_bins));
  std::vector<double> p_sub;
  std::vector<int64_t> sub_ix;
  for (int64_t i = 0; i < n; ++i) {
    if (obs[static_cast<size_t>(i)] < 1) continue;
    const int o = static_cast<int>(std::min<int64_t>(
        obs[static_cast<size_t>(i)],
        static_cast<int64_t>(std::numeric_limits<int>::max() / 2)));
    sub_ix.push_back(i);
    p_sub.push_back(PoissonSurvivalUpper(o, lam.At(i)));
  }
  if (p_sub.empty()) {
    return 0;
  }

  // Pass 2: BH FDR over the sparse p-values.
  std::vector<double> q_sub;
  BenjaminiHochbergFdr(p_sub, &q_sub);

  // Pass 3: filter to significant bins; carry their p/q forward parallel-array
  // style so we never need to re-look-up by bin index during emission.
  std::vector<int32_t> sig;
  std::vector<double> sig_p;
  std::vector<double> sig_q;
  for (size_t j = 0; j < sub_ix.size(); ++j) {
    if (q_sub[j] > fdr) continue;
    if (p_sub[j] > p_cutoff) continue;
    sig.push_back(static_cast<int32_t>(sub_ix[j]));
    sig_p.push_back(p_sub[j]);
    sig_q.push_back(q_sub[j]);
  }
  if (sig.empty()) {
    return 0;
  }

  // Pass 4: merge runs (sorted sig; cluster when next - last <= cap), pick
  // summit (max obs in [b0, b1]; ties broken by smallest bin), emit Peak.
  // Summit's p/q comes from sub_ix via lower_bound — best is guaranteed
  // present because best_v starts at -1 and b0 in sig has obs >= 1, so the
  // chosen bin always has obs >= 1 and is therefore in sub_ix.
  const int cap = 1 + std::max(0, merge_bin_gap);
  size_t a = 0;
  while (a < sig.size()) {
    size_t b = a;
    while (b + 1 < sig.size() &&
           sig[b + 1] - sig[b] <= static_cast<int32_t>(cap)) {
      ++b;
    }
    const int32_t b0 = sig[a];
    const int32_t b1 = sig[b];
    int32_t best = b0;
    int64_t best_v = -1;
    for (int32_t t = b0; t <= b1; ++t) {
      const int64_t v = obs[static_cast<size_t>(t)];
      if (v > best_v || (v == best_v && t < best)) {
        best_v = v;
        best = t;
      }
    }
    auto it = std::lower_bound(sub_ix.begin(), sub_ix.end(),
                               static_cast<int64_t>(best));
    double bp = 1.0;
    double bq = 1.0;
    if (it != sub_ix.end() && *it == best) {
      const size_t idx = static_cast<size_t>(it - sub_ix.begin());
      bp = p_sub[idx];
      bq = q_sub[idx];
    }
    const int32_t start = b0 * bin_width;
    const int32_t end = (b1 + 1) * bin_width;
    const int32_t summit_0 = best * bin_width + std::max(0, bin_width / 2);
    Peak p;
    p.chrom = chrom;
    p.start = start;
    p.end = end;
    p.thick_start = start;
    p.thick_end = end;
    p.peak_offset = summit_0 - start;
    p.max_signal_bin = best;
    p.max_signal = best_v;
    p.p_value = bp;
    p.q_value = bq;
    out_peaks->push_back(p);
    a = b + 1;
  }
  return static_cast<int>(out_peaks->size());
}

}  // namespace peaks
}  // namespace chromap
