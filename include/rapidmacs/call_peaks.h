#ifndef CHROMAP_PEAK_CALLER_CALL_PEAKS_H_
#define CHROMAP_PEAK_CALLER_CALL_PEAKS_H_

#include <cstdint>
#include <string>
#include <vector>

namespace chromap {
namespace peaks {

struct Peak {
  std::string chrom;
  int32_t start = 0;  // 0-based BED inclusive
  int32_t end = 0;    // 0-based BED exclusive
  int32_t thick_start = 0;
  int32_t thick_end = 0;  // narrowPeak: thickStart, thickEnd (unused = start/end)
  int32_t peak_offset = 0;  // narrowPeak: relative to start (0-based)
  int32_t max_signal_bin = 0;
  int64_t max_signal = 0;
  double p_value = 1.0;
  double q_value = 1.0;
};

// Poisson P(X >= k) for k >= 0 (k integer count); uses hts kfunc (gamma).
double PoissonSurvivalUpper(int k, double lambda);

// Benjamini--Hochberg; p-values in [0,1], m = size.
void BenjaminiHochbergFdr(const std::vector<double>& p_values,
                          std::vector<double>* q_values);

// Local background: for each i, take max of global mean and
// (sum_{j in window, j!=i} obs[j]) / (window_count-1) when the latter is
// defined. Window is [i - half, i + half] in bins; half = local_window_bins.
void ComputeLocalBackground(const std::vector<uint16_t>& obs,
                            int32_t local_window_bins, double* global_per_bin,
                            std::vector<double>* lambda);

// Significance, merge (gap in bins), summit per merged block.
int CallPeaksOnBins(
    const std::string& chrom, int32_t bin_width,
    const std::vector<uint16_t>& obs, int32_t local_window_bins, double fdr,
    double p_cutoff, int32_t merge_bin_gap, std::vector<Peak>* out_peaks);

}  // namespace peaks
}  // namespace chromap

#endif  // CHROMAP_PEAK_CALLER_CALL_PEAKS_H_
