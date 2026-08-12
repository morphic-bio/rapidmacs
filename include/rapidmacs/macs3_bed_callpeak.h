#ifndef RAPIDMACS_MACS3_BED_CALLPEAK_H_
#define RAPIDMACS_MACS3_BED_CALLPEAK_H_

#include <cstdint>
#include <string>

#include "rapidmacs/bdgpeakcall.h"

namespace chromap {
namespace peaks {

struct Macs3BedCallPeakParams {
  Macs3FragThresholdMode threshold_mode = Macs3FragThresholdMode::kQValue;
  double qvalue_cutoff = 0.05;
  double pvalue_cutoff = 0.0;
  int64_t effective_genome_size = 2913022398LL;
  bool nomodel = true;
  int32_t extsize = 200;
  int32_t shift = 0;
  int32_t llocal_bp = 10000;
  int32_t min_length = 0;  // 0 = MACS default: fragment length d.
  int32_t max_gap = 0;     // 0 = MACS default: detected tag size.
  std::string keep_dup = "1";
  bool broad = false;
  double broad_cutoff = 0.1;
  // Use only the genome-wide background lambda, matching MACS3 --nolambda.
  bool no_lambda = false;
  bool call_summits = false;
  std::string name_prefix;
  std::string treat_bdg;
  std::string lambda_bdg;
  std::string ppois_bdg;
};

struct Macs3BedCallPeakStats {
  int64_t total_tags = 0;
  int64_t tags_after_filtering = 0;
  int32_t detected_tag_size = 0;
  int32_t resolved_min_length = 0;
  int32_t resolved_max_gap = 0;
  int32_t keep_dup_max = -1;  // -1 means all.
  int64_t peak_count = 0;
};

struct Macs3BedCallPeakProfile {
  std::string name;
  std::string input_format = "BED";
  Macs3BedCallPeakParams params;
};

// Named BED callpeak profiles are opt-in compatibility presets for reproducing
// external workflows. They do not change the default MACS3 FRAG pipeline
// settings used by Signac.
const char* SupportedMacs3BedCallPeakProfiles();

bool ResolveMacs3BedCallPeakProfile(const std::string& profile_name,
                                    Macs3BedCallPeakProfile* profile,
                                    std::string* error_message);

bool RunMacs3BedCallPeak(const std::string& bed_input,
                         const Macs3BedCallPeakParams& params,
                         const std::string& narrowpeak_out,
                         const std::string& summits_out,
                         Macs3BedCallPeakStats* stats,
                         std::string* error_message);

}  // namespace peaks
}  // namespace chromap

#endif  // RAPIDMACS_MACS3_BED_CALLPEAK_H_
