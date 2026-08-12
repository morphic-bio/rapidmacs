#ifndef RAPIDMACS_FRAGMENTS_H_
#define RAPIDMACS_FRAGMENTS_H_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace macs3 {

struct FragmentRecord {
  int32_t chrom_id = -1;
  int32_t start = 0;
  int32_t end = 0;
  uint32_t count = 1;
};

class FragmentIterator {
 public:
  virtual ~FragmentIterator() = default;

  // Fill *out with the next fragment. Returns true on success, false at EOF
  // or on unrecoverable error (use Error() to distinguish).
  virtual bool Next(FragmentRecord* out) = 0;

  // Stable id -> name map. Indices match FragmentRecord::chrom_id.
  // Names are interned in first-seen order; the vector grows monotonically.
  virtual const std::vector<std::string>& ChromNames() const = 0;

  // Empty when no error. Populated only after Next() returns false due to
  // a parse or I/O failure (vs end-of-stream).
  virtual const std::string& Error() const = 0;
};

// Reads a Cell Ranger/Chromap fragments TSV (plain or gzipped):
//   chrom \t start \t end [\t barcode [\t count]]
// and Chromap's 7-column bulk fragment BED:
//   chrom \t start \t end \t N \t mapq \t strand \t count
// Barcode, name, MAPQ, and strand are ignored after schema validation. Count
// is optional for 3-/4-column inputs (absent => 1); it is column 5 for the
// barcoded schema and column 7 for the bulk schema.
//
// Strict validation. Next() returns false and sets Error() with a
// path:line:reason message on any of:
//   - missing required field (chrom / start / end)
//   - non-numeric or non-fully-consumed start / end / count
//   - negative start
//   - end <= start
//   - count < 1
//   - line longer than 65535 bytes (no fragment line should ever be this)
// Lines starting with '#' and empty lines are skipped (header comments).
std::unique_ptr<FragmentIterator> OpenFileFragmentIterator(
    const std::string& path);

// Iterator over an in-memory vector of FragmentRecords. Useful as an adapter
// for callers (e.g. libchromap concurrent path) that buffer fragments
// during mapping rather than streaming them. Construction sorts the input
// vector by (chrom_id, start, end) so the output stream satisfies the
// chrom-grouped + start-sorted constraint required by
// RunMacs3FragPeakPipelineFromSortedIterator.
//
// `chrom_names` is the id->name table used by ChromNames(); it must cover
// every chrom_id present in `fragments` (chrom_names.size() > max chrom_id).
// The iterator takes ownership of the vector via std::move.
std::unique_ptr<FragmentIterator> WrapVectorFragmentIterator(
    std::vector<FragmentRecord> fragments,
    std::vector<std::string> chrom_names);

}  // namespace macs3

#endif  // RAPIDMACS_FRAGMENTS_H_
