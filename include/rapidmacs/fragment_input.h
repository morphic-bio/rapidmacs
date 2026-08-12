#ifndef CHROMAP_PEAK_CALLER_FRAGMENT_INPUT_H_
#define CHROMAP_PEAK_CALLER_FRAGMENT_INPUT_H_

#include <cstdint>
#include <string>
#include <vector>

namespace chromap {
namespace peaks {

// Cell Ranger/Chromap barcoded style: chrom, start, end, barcode,
// duplicate_count; Chromap bulk style: chrom, start, end, N, mapq, strand,
// duplicate_count (0-based, [start, end) BED style).
struct Fragment {
  int32_t start = 0;
  int32_t end = 0;
  int32_t count = 1;
};

struct ChromFragments {
  std::string name;
  std::vector<Fragment> frags;
  int32_t max_end = 0;
  // Reference length from an alignment header when known (BAMPE). Zero for
  // text inputs that do not carry a sequence dictionary.
  int32_t reference_length = 0;
};

// Reader accounting used by the CLI to reproduce MACS3's paired-end default
// geometry. MACS3 estimates the mean fragment length before applying its
// coordinate duplicate filter, so both pre- and post-filter totals matter.
struct FragmentInputStats {
  int64_t records_read = 0;
  int64_t fragments_parsed = 0;
  int64_t fragment_bp_parsed = 0;
  int64_t fragments_retained = 0;
  int64_t duplicates_removed = 0;
  int64_t records_skipped = 0;
};

// Reads a 5-column barcoded or 7-column Chromap bulk tab-separated fragments
// file (plain or .gz). Lines starting with '#', empty, or malformed are
// skipped. Returns true on full success (file opened, no read errors; empty
// input is an error for callers to decide).
// Fragments are grouped by `chrom` in strict lexicographic order and sorted
// by (start, end) within each group.
bool LoadFragmentsFromTsv(
    const std::string& path, std::vector<ChromFragments>* by_chrom);

// MACS3 `-f BAMPE` reader. Keeps one record for each properly paired primary
// template (read 1), including alignments marked duplicate, then applies a
// MACS3-style coordinate duplicate cap. `max_duplicates` is -1 for all or a
// positive maximum per identical (start, end) coordinate (default 1). MACS3
// does not impose a MAPQ threshold in this parser. QC-fail, unmapped,
// secondary, and supplementary records are skipped. Output ordering matches
// the TSV reader.
bool LoadFragmentsFromBampe(
    const std::string& path, std::vector<ChromFragments>* by_chrom,
    FragmentInputStats* stats = nullptr, int num_threads = 1,
    int max_duplicates = 1);

// MACS3 BEDPE is a three-column fragment file (chrom, start, end). For
// interoperability this reader also accepts standard six-column BEDPE mate
// intervals when both mates are on the same chromosome. Applies the same
// MACS3-style coordinate duplicate cap as BAMPE.
bool LoadFragmentsFromBedpe(
    const std::string& path, std::vector<ChromFragments>* by_chrom,
    FragmentInputStats* stats = nullptr, int max_duplicates = 1);

}  // namespace peaks
}  // namespace chromap

#endif  // CHROMAP_PEAK_CALLER_FRAGMENT_INPUT_H_
