#include "rapidmacs/frag_compact_store.h"

#include <algorithm>
#include <limits>

namespace chromap {
namespace peaks {
namespace {

uint32_t MaskLowerBitsU32(int bits) {
  if (bits >= 32) {
    return 0xFFFFFFFFu;
  }
  if (bits <= 0) {
    return 0u;
  }
  return (1u << bits) - 1u;
}

}  // namespace

uint32_t CeilLog2U32(uint32_t v) {
  if (v <= 1u) {
    return 0u;
  }
  return 32u - static_cast<uint32_t>(__builtin_clz(v - 1u));
}

bool FragPackingPolicy::FromMaxInsertSize(int max_insert_size, int min_count_bits,
                                         FragPackingPolicy* out,
                                         std::string* error) {
  if (out == nullptr) {
    return false;
  }
  *out = FragPackingPolicy{};
  if (max_insert_size < 0) {
    if (error) {
      *error = "max_insert_size must be non-negative for fragment packing";
    }
    return false;
  }
  if (min_count_bits < 1 || min_count_bits > 31) {
    if (error) {
      *error = "min_count_bits must be in [1, 31]";
    }
    return false;
  }
  const uint32_t m = static_cast<uint32_t>(max_insert_size);
  // At least one bit is required to represent non-zero fragment lengths; for
  // m==0, CeilLog2(1) is 0 in the raw helper.
  const uint32_t len_bits = std::max(1u, CeilLog2U32(m + 1u));
  const int count_bits = 32 - static_cast<int>(len_bits);
  if (count_bits < min_count_bits) {
    out->storage = FragPeakStorageClass::kWide;
    out->length_bits = 0;
    out->count_bits = 0;
    out->max_length = std::numeric_limits<uint32_t>::max();
    out->max_count = std::numeric_limits<uint32_t>::max();
    return true;
  }
  if (static_cast<int>(len_bits) < 1) {
    if (error) {
      *error = "internal: invalid length_bit split";
    }
    return false;
  }
  out->storage = FragPeakStorageClass::kCompact64;
  out->length_bits = static_cast<int>(len_bits);
  out->count_bits = count_bits;
  out->max_length = MaskLowerBitsU32(out->length_bits);
  out->max_count = MaskLowerBitsU32(out->count_bits);
  return true;
}

FragPeakMemoryAccumulator::FragPeakMemoryAccumulator(int max_insert_size,
                                                     int min_count_bits,
                                                     std::string* error) {
  min_count_bits_ = min_count_bits;
  mode_label_ = "invalid";
  if (!FragPackingPolicy::FromMaxInsertSize(max_insert_size, min_count_bits, &policy_,
                                             error)) {
    if (error && error->empty()) {
      *error = "failed to build fragment packing policy";
    }
    return;
  }
  if (policy_.storage == FragPeakStorageClass::kCompact64) {
    mode_label_ = "Compact64";
  } else {
    mode_label_ = "Wide12";
  }
}

void FragPeakMemoryAccumulator::ResizeForReferenceNames(
    uint32_t num_ref_seqs, const std::vector<std::string>& chrom_names) {
  if (!IsValid()) {
    return;
  }
  if (chrom_names.size() != static_cast<size_t>(num_ref_seqs)) {
    return;
  }
  chrom_names_ = chrom_names;
  compact_by_chrom_.assign(num_ref_seqs, std::vector<CompactFrag64>{});
  wide_by_chrom_.assign(num_ref_seqs, std::vector<WideFrag>{});
  storage_ready_ = true;
}

bool FragPeakMemoryAccumulator::PackOrAssignWide(
    uint32_t start, uint32_t length, uint32_t count, CompactFrag64* compact,
    WideFrag* wide, std::string* error) const {
  if (policy_.storage == FragPeakStorageClass::kCompact64) {
    if (length > policy_.max_length) {
      if (error) {
        *error = "fragment length " + std::to_string(length) +
                 " exceeds packed max_length " + std::to_string(policy_.max_length);
      }
      return false;
    }
    if (count > policy_.max_count) {
      if (error) {
        *error = "fragment duplicate count " + std::to_string(count) +
                 " exceeds packed max_count " + std::to_string(policy_.max_count);
      }
      return false;
    }
    const uint32_t low = length & MaskLowerBitsU32(policy_.length_bits);
    const uint32_t packed = (count << policy_.length_bits) | low;
    compact->start = start;
    compact->packed_len_count = packed;
    (void)wide;
    return true;
  }
  wide->start = start;
  wide->length = length;
  wide->count = count;
  (void)compact;
  return true;
}

bool FragPeakMemoryAccumulator::Add(uint32_t rid, uint32_t start, uint32_t end,
                                    uint32_t count, std::string* error) {
  if (!IsValid()) {
    if (error) {
      *error = "invalid fragment memory accumulator (construction failed)";
    }
    return false;
  }
  if (!storage_ready_) {
    if (error) {
      *error = "fragment memory accumulator: ResizeForReference not called";
    }
    return false;
  }
  if (end <= start) {
    if (error) {
      *error = "fragment end must be > start (got end=" + std::to_string(end) +
               " start=" + std::to_string(start) + ")";
    }
    return false;
  }
  if (count == 0) {
    if (error) {
      *error = "fragment duplicate count must be > 0";
    }
    return false;
  }
  if (rid >= wide_by_chrom_.size() || rid >= compact_by_chrom_.size()) {
    if (error) {
      *error = "reference id " + std::to_string(rid) + " out of range for accumulator";
    }
    return false;
  }
  const uint32_t len = end - start;
  if (policy_.storage == FragPeakStorageClass::kCompact64) {
    CompactFrag64 c{};
    if (!PackOrAssignWide(start, len, count, &c, nullptr, error)) {
      return false;
    }
    compact_by_chrom_[rid].push_back(c);
  } else {
    WideFrag w{};
    if (!PackOrAssignWide(start, len, count, nullptr, &w, error)) {
      return false;
    }
    wide_by_chrom_[rid].push_back(w);
  }
  return true;
}

bool FragPeakMemoryAccumulator::MaterializeToChromFragments(
    std::vector<ChromFragments>* by_chrom, std::string* error) const {
  if (by_chrom == nullptr) {
    return false;
  }
  by_chrom->clear();
  if (!storage_ready_ || !IsValid()) {
    if (error) {
      *error = "cannot materialize: accumulator not ready";
    }
    return false;
  }
  if (UseCompact64()) {
    for (uint32_t rid = 0; rid < static_cast<uint32_t>(compact_by_chrom_.size());
         ++rid) {
      if (rid >= chrom_names_.size()) {
        if (error) {
          *error = "internal: missing chrom name for reference id " +
                   std::to_string(rid);
        }
        return false;
      }
      if (compact_by_chrom_[rid].empty()) {
        continue;
      }
      ChromFragments ch;
      ch.name = chrom_names_[rid];
      ch.max_end = 0;
      for (const CompactFrag64& c : compact_by_chrom_[rid]) {
        const uint32_t low = c.packed_len_count & MaskLowerBitsU32(policy_.length_bits);
        const uint32_t packed_count = c.packed_len_count >> policy_.length_bits;
        const uint32_t len = low;
        const uint32_t co = packed_count;
        if (c.start > static_cast<uint32_t>(std::numeric_limits<int32_t>::max()) ||
            len > static_cast<uint32_t>(std::numeric_limits<int32_t>::max()) ||
            static_cast<int64_t>(c.start) + static_cast<int64_t>(len) >
                static_cast<int64_t>(std::numeric_limits<int32_t>::max())) {
          if (error) {
            *error = "fragment coordinates overflow int32 during materialize";
          }
          return false;
        }
        Fragment f;
        f.start = static_cast<int32_t>(c.start);
        f.end = static_cast<int32_t>(c.start + len);
        f.count = static_cast<int32_t>(co);
        if (f.count < 0) {
          if (error) {
            *error = "fragment count overflow int32";
          }
          return false;
        }
        ch.frags.push_back(f);
        if (f.end > ch.max_end) {
          ch.max_end = f.end;
        }
      }
      by_chrom->push_back(std::move(ch));
    }
  } else {
    for (uint32_t rid = 0; rid < static_cast<uint32_t>(wide_by_chrom_.size());
         ++rid) {
      if (rid >= chrom_names_.size()) {
        if (error) {
          *error = "internal: missing chrom name for reference id " +
                   std::to_string(rid);
        }
        return false;
      }
      if (wide_by_chrom_[rid].empty()) {
        continue;
      }
      ChromFragments ch;
      ch.name = chrom_names_[rid];
      ch.max_end = 0;
      for (const WideFrag& w : wide_by_chrom_[rid]) {
        if (w.start > static_cast<uint32_t>(std::numeric_limits<int32_t>::max()) ||
            w.length > static_cast<uint32_t>(std::numeric_limits<int32_t>::max()) ||
            static_cast<int64_t>(w.start) + static_cast<int64_t>(w.length) >
                static_cast<int64_t>(std::numeric_limits<int32_t>::max())) {
          if (error) {
            *error = "fragment coordinates overflow int32 during materialize (wide)";
          }
          return false;
        }
        Fragment f;
        f.start = static_cast<int32_t>(w.start);
        f.end = static_cast<int32_t>(w.start + w.length);
        f.count = static_cast<int32_t>(w.count);
        if (f.count < 0) {
          if (error) {
            *error = "fragment count overflow int32 (wide)";
          }
          return false;
        }
        ch.frags.push_back(f);
        if (f.end > ch.max_end) {
          ch.max_end = f.end;
        }
      }
      by_chrom->push_back(std::move(ch));
    }
  }
  for (auto& ch : *by_chrom) {
    std::sort(
        ch.frags.begin(), ch.frags.end(),
        [](const Fragment& a, const Fragment& b) {
          if (a.start != b.start) {
            return a.start < b.start;
          }
          return a.end < b.end;
        });
  }
  std::sort(by_chrom->begin(), by_chrom->end(),
            [](const ChromFragments& a, const ChromFragments& b) {
              return a.name < b.name;
            });
  if (by_chrom->empty()) {
    if (error) {
      *error = "no fragments in memory accumulator";
    }
    return false;
  }
  return true;
}

}  // namespace peaks
}  // namespace chromap

static_assert(sizeof(chromap::peaks::CompactFrag64) == 8, "CompactFrag64 size");
static_assert(sizeof(chromap::peaks::WideFrag) == 12, "WideFrag size");
