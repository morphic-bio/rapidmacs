#ifndef CHROMAP_PEAK_CALLER_FRAG_COMPACT_STORE_H_
#define CHROMAP_PEAK_CALLER_FRAG_COMPACT_STORE_H_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "rapidmacs/fragment_input.h"

namespace chromap {
namespace peaks {

// Single-row fragment: start (0-based BED) + packed length and duplicate count
// in the low 32 bits (layout chosen by FragPackingPolicy).
struct CompactFrag64 {
  uint32_t start;
  uint32_t packed_len_count;
};

// Explicit length + count when 32-bit packing does not meet policy thresholds.
struct WideFrag {
  uint32_t start;
  uint32_t length;
  uint32_t count;
};

enum class FragPeakStorageClass {
  kCompact64,
  kWide
};

// Runtime packing: length_bits = ceil_log2(max_insert_size + 1);
// count_bits = 32 - length_bits. Compact mode stores length in the low
// `length_bits` and count in the upper `count_bits` of packed_len_count.
struct FragPackingPolicy {
  FragPeakStorageClass storage = FragPeakStorageClass::kWide;
  int length_bits = 0;
  int count_bits = 0;
  uint32_t max_length = 0;
  uint32_t max_count = 0;

  // Computes policy; returns false if parameters are invalid (max_insert_size < 0).
  static bool FromMaxInsertSize(int max_insert_size, int min_count_bits,
                                FragPackingPolicy* out, std::string* error);
};

uint32_t CeilLog2U32(uint32_t v);

// In-memory per-(reference id) buffer used for MACS3 FRAG peak calling. Rows are
// append-only in emission order; materialization matches LoadFragmentsFromTsv
// sort/grouping. One storage class is active; the other container stays empty.
class FragPeakMemoryAccumulator {
 public:
  FragPeakMemoryAccumulator(int max_insert_size, int min_count_bits,
                            std::string* error);

  bool UseCompact64() const {
    return policy_.storage == FragPeakStorageClass::kCompact64;
  }
  const FragPackingPolicy& policy() const { return policy_; }
  const std::string& StorageModeLabel() const { return mode_label_; }
  // False if construction failed (invalid max_insert / min_count_bits).
  bool IsValid() const { return mode_label_ != "invalid"; }

  void ResizeForReferenceNames(uint32_t num_ref_seqs,
                               const std::vector<std::string>& chrom_names);

  // Append a fragment row (same numeric semantics as a fragments TSV line:
  // [start, end) and duplicate count, no aggregation of identical rows).
  bool Add(uint32_t rid, uint32_t start, uint32_t end, uint32_t count,
           std::string* error);

  // Lexicographic chrom order + (start, end) sort within chrom, as in
  // LoadFragmentsFromTsv.
  bool MaterializeToChromFragments(std::vector<ChromFragments>* by_chrom,
                                    std::string* error) const;

 private:
  bool PackOrAssignWide(uint32_t start, uint32_t length, uint32_t count,
                        CompactFrag64* compact, WideFrag* wide,
                        std::string* error) const;

  FragPackingPolicy policy_{};
  std::string mode_label_;
  int min_count_bits_ = 16;

  std::vector<std::string> chrom_names_;
  std::vector<std::vector<CompactFrag64>> compact_by_chrom_;
  std::vector<std::vector<WideFrag>> wide_by_chrom_;
  bool storage_ready_ = false;
};

}  // namespace peaks
}  // namespace chromap

#endif  // CHROMAP_PEAK_CALLER_FRAG_COMPACT_STORE_H_
