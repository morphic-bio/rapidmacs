// Unit: OpenFileFragmentIterator rejects malformed inputs with a clear
// path:line:reason error rather than silently defaulting / corrupting.
#include "rapidmacs/fragments.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace {

void Write(const std::string& path, const std::string& s) {
  FILE* f = std::fopen(path.c_str(), "w");
  assert(f != nullptr);
  std::fwrite(s.data(), 1, s.size(), f);
  std::fclose(f);
}

bool DrainExpectError(const std::string& path, const std::string& want_substr) {
  auto it = macs3::OpenFileFragmentIterator(path);
  if (!it) return false;
  macs3::FragmentRecord rec;
  while (it->Next(&rec)) {}
  const std::string& err = it->Error();
  if (err.empty()) {
    std::fprintf(stderr, "  expected error containing '%s' but got success\n",
                 want_substr.c_str());
    return false;
  }
  if (err.find(want_substr) == std::string::npos) {
    std::fprintf(stderr, "  expected '%s' in error, got: %s\n",
                 want_substr.c_str(), err.c_str());
    return false;
  }
  return true;
}

bool DrainExpectOk(const std::string& path, size_t want_records) {
  auto it = macs3::OpenFileFragmentIterator(path);
  if (!it) return false;
  size_t n = 0;
  macs3::FragmentRecord rec;
  while (it->Next(&rec)) ++n;
  if (!it->Error().empty()) {
    std::fprintf(stderr, "  unexpected error: %s\n", it->Error().c_str());
    return false;
  }
  if (n != want_records) {
    std::fprintf(stderr, "  want %zu records, got %zu\n", want_records, n);
    return false;
  }
  return true;
}

}  // namespace

int main() {
  const std::string path = "/tmp/rapidmacs_strict_test.tsv";
  int failed = 0;

  // 1. Bad count: non-numeric.
  Write(path, "chr1\t100\t200\tBC\tNOPE\n");
  if (!DrainExpectError(path, "bad count")) { ++failed; std::fputs("FAIL bad-count\n", stderr); }

  // 2. Zero count.
  Write(path, "chr1\t100\t200\tBC\t0\n");
  if (!DrainExpectError(path, "non-positive count")) { ++failed; std::fputs("FAIL zero-count\n", stderr); }

  // 3. Negative count.
  Write(path, "chr1\t100\t200\tBC\t-1\n");
  if (!DrainExpectError(path, "non-positive count")) { ++failed; std::fputs("FAIL negative-count\n", stderr); }

  // 4. Missing count column (4-col input): valid, default to 1.
  Write(path, "chr1\t100\t200\tBC\n");
  if (!DrainExpectOk(path, 1)) { ++failed; std::fputs("FAIL missing-count\n", stderr); }

  // 5. Three-column input (chrom/start/end only): valid, default count=1.
  Write(path, "chr1\t100\t200\n");
  if (!DrainExpectOk(path, 1)) { ++failed; std::fputs("FAIL 3col\n", stderr); }

  // 6. Bad start.
  Write(path, "chr1\tNOPE\t200\tBC\t1\n");
  if (!DrainExpectError(path, "bad start")) { ++failed; std::fputs("FAIL bad-start\n", stderr); }

  // 7. Bad end.
  Write(path, "chr1\t100\tNOPE\tBC\t1\n");
  if (!DrainExpectError(path, "bad end")) { ++failed; std::fputs("FAIL bad-end\n", stderr); }

  // 8. Negative start.
  Write(path, "chr1\t-1\t200\tBC\t1\n");
  if (!DrainExpectError(path, "negative start")) { ++failed; std::fputs("FAIL negative-start\n", stderr); }

  // 9. end <= start.
  Write(path, "chr1\t200\t100\tBC\t1\n");
  if (!DrainExpectError(path, "end <= start")) { ++failed; std::fputs("FAIL end-le-start\n", stderr); }

  // 10. end == start.
  Write(path, "chr1\t100\t100\tBC\t1\n");
  if (!DrainExpectError(path, "end <= start")) { ++failed; std::fputs("FAIL end-eq-start\n", stderr); }

  // 11. Empty chrom.
  Write(path, "\t100\t200\tBC\t1\n");
  if (!DrainExpectError(path, "empty chrom")) { ++failed; std::fputs("FAIL empty-chrom\n", stderr); }

  // 12. Comment + blank lines + valid records mixed.
  Write(path,
        "# header\n"
        "\n"
        "chr1\t100\t200\tBC\t1\n"
        "\n"
        "chr2\t10\t90\n"
        "# trailing comment\n");
  if (!DrainExpectOk(path, 2)) { ++failed; std::fputs("FAIL comment-mix\n", stderr); }

  // 13. Trailing junk after count.
  Write(path, "chr1\t100\t200\tBC\t1junk\n");
  if (!DrainExpectError(path, "bad count")) { ++failed; std::fputs("FAIL trailing-junk\n", stderr); }

  // 14. Seven-column Chromap bulk count is column 7, not MAPQ in column 5.
  Write(path, "chr1\t100\t200\tN\t60\t+\t4\n");
  if (!DrainExpectOk(path, 1)) { ++failed; std::fputs("FAIL bulk-7col\n", stderr); }

  // 15. Malformed bulk schema is rejected instead of using column 5.
  Write(path, "chr1\t100\t200\tN\t60\t?\t4\n");
  if (!DrainExpectError(path, "bad Chromap bulk strand")) { ++failed; std::fputs("FAIL bulk-strand\n", stderr); }

  // 16. Bulk duplicate count receives the same strict validation.
  Write(path, "chr1\t100\t200\tN\t60\t-\t0\n");
  if (!DrainExpectError(path, "non-positive count")) { ++failed; std::fputs("FAIL bulk-count\n", stderr); }

  std::remove(path.c_str());
  if (failed > 0) {
    std::fprintf(stderr, "FAIL test_fragment_iterator_strict: %d cases\n", failed);
    return 1;
  }
  std::printf("PASS test_fragment_iterator_strict (16 cases)\n");
  return 0;
}
