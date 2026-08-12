#ifndef CHROMAP_PEAK_CALLER_BINNED_SIGNAL_H_
#define CHROMAP_PEAK_CALLER_BINNED_SIGNAL_H_

#include <cstdint>
#include <vector>

#include "rapidmacs/fragment_input.h"

namespace chromap {
namespace peaks {

// Per-bin weighted overlap counts: each fragment's two Tn5 cuts are extended by
// half of `ext_size_bp` to each side (ext_size_bp default 150 => +/-75) and
// the duplicate count is added to every bin the extended interval covers.
// Counts saturate at UINT16_MAX (65535) per bin — adequate for 50 bp bins on
// real ATAC unless extreme hotspots (mtDNA / rDNA / centromeric satellites).
// If saturation matters, callers should switch to a uint32 substrate.
// Returns false if the chromosome has no non-empty bins.
bool BuildBinnedCutSignal(const ChromFragments& ch, int32_t bin_width_bp,
                          int32_t ext_size_bp, std::vector<uint16_t>* obs,
                          int32_t* padded_max_coord);

}  // namespace peaks
}  // namespace chromap

#endif  // CHROMAP_PEAK_CALLER_BINNED_SIGNAL_H_
