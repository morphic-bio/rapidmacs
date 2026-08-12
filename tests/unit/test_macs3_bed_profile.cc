#include "rapidmacs/macs3_bed_callpeak.h"

#include <cassert>
#include <cstdio>
#include <string>

int main() {
  chromap::peaks::Macs3BedCallPeakProfile profile;
  std::string err;
  assert(chromap::peaks::ResolveMacs3BedCallPeakProfile("signac-atac",
                                                        &profile, &err));
  assert(err.empty());
  assert(profile.name == "signac-atac");
  assert(profile.input_format == "BED");
  assert(profile.params.threshold_mode ==
         chromap::peaks::Macs3FragThresholdMode::kQValue);
  assert(profile.params.qvalue_cutoff == 0.05);
  assert(profile.params.pvalue_cutoff == 0.0);
  assert(profile.params.effective_genome_size == 2700000000LL);
  assert(profile.params.nomodel);
  assert(profile.params.extsize == 200);
  assert(profile.params.shift == -100);
  assert(profile.params.llocal_bp == 10000);
  assert(profile.params.keep_dup == "1");
  assert(!profile.params.broad);
  assert(profile.params.broad_cutoff == 0.1);
  assert(!profile.params.call_summits);
  assert(!profile.params.no_lambda);

  err.clear();
  assert(!chromap::peaks::ResolveMacs3BedCallPeakProfile("missing-profile",
                                                         &profile, &err));
  assert(err.find("signac-atac") != std::string::npos);
  assert(std::string(chromap::peaks::SupportedMacs3BedCallPeakProfiles()) ==
         "signac-atac");

  std::printf("PASS test_macs3_bed_profile\n");
  return 0;
}
