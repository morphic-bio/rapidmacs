#ifndef CHROMAP_PEAK_CALLER_MACS3_BDG_TRACK_H_
#define CHROMAP_PEAK_CALLER_MACS3_BDG_TRACK_H_

#include <cstdint>
#include <vector>

namespace chromap {
namespace peaks {

// MACS3 bedGraph track representation: p[i] is the end coordinate of run i,
// v[i] is the value on [i == 0 ? 0 : p[i - 1], p[i]).
struct Macs3BdgTrack {
  std::vector<int32_t> p;
  std::vector<float> v;
};

}  // namespace peaks
}  // namespace chromap

#endif  // CHROMAP_PEAK_CALLER_MACS3_BDG_TRACK_H_
