// Unit: OpenFileFragmentIterator parses synthetic 5-col and Chromap bulk
// 7-col fragment TSV rows correctly.
#include "rapidmacs/fragments.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

void WriteTempTsv(const std::string& path, const std::string& contents) {
  FILE* f = std::fopen(path.c_str(), "w");
  assert(f != nullptr);
  std::fwrite(contents.data(), 1, contents.size(), f);
  std::fclose(f);
}

}  // namespace

int main() {
  const std::string path = "/tmp/rapidmacs_test_frag.tsv";
  WriteTempTsv(path,
               "# header comment\n"
               "\n"
               "chr1\t100\t200\tBC1\t1\n"
               "chr1\t150\t250\tBC2\t3\n"
               "chr2\t10\t90\tBC1\t1\n"
               "chr1\t300\t400\tBC3\t1\n"
               "chr3\t20\t120\tN\t60\t+\t7\n"
               "chr3\t150\t260\tN\t42\t-\t2\n");

  auto it = macs3::OpenFileFragmentIterator(path);
  assert(it != nullptr);

  std::vector<macs3::FragmentRecord> records;
  macs3::FragmentRecord rec;
  while (it->Next(&rec)) records.push_back(rec);

  assert(it->Error().empty());
  assert(records.size() == 6);

  // chr1 interned first (id 0), chr2 second (id 1).
  const auto& names = it->ChromNames();
  assert(names.size() == 3);
  assert(names[0] == "chr1");
  assert(names[1] == "chr2");
  assert(names[2] == "chr3");

  assert(records[0].chrom_id == 0 && records[0].start == 100 &&
         records[0].end == 200 && records[0].count == 1);
  assert(records[1].chrom_id == 0 && records[1].start == 150 &&
         records[1].end == 250 && records[1].count == 3);
  assert(records[2].chrom_id == 1 && records[2].start == 10 &&
         records[2].end == 90 && records[2].count == 1);
  assert(records[3].chrom_id == 0 && records[3].start == 300 &&
         records[3].end == 400 && records[3].count == 1);
  // The bulk count comes from column 7; column 5 is MAPQ and must not be
  // mistaken for duplicate count.
  assert(records[4].chrom_id == 2 && records[4].start == 20 &&
         records[4].end == 120 && records[4].count == 7);
  assert(records[5].chrom_id == 2 && records[5].start == 150 &&
         records[5].end == 260 && records[5].count == 2);

  std::remove(path.c_str());
  std::printf("PASS test_fragment_iterator (parsed %zu records, %zu chroms)\n",
              records.size(), names.size());
  return 0;
}
