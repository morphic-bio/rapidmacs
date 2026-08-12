// Unit: WrapVectorFragmentIterator sorts on construction and yields
// fragments in (chrom_id, start, end) order. Useful as the adapter that
// turns an in-memory accumulated fragment buffer (e.g. from libchromap's
// concurrent mapping output) into a FragmentIterator suitable for
// RunMacs3FragPeakPipelineFromSortedIterator.
#include "rapidmacs/fragments.h"

#include <cassert>
#include <cstdio>
#include <vector>

int main() {
  std::vector<macs3::FragmentRecord> in;
  // Deliberately out of order across both chrom_id and start.
  auto add = [&](int32_t cid, int32_t s, int32_t e, uint32_t c) {
    macs3::FragmentRecord r;
    r.chrom_id = cid; r.start = s; r.end = e; r.count = c;
    in.push_back(r);
  };
  add(1, 500, 600, 1);
  add(0, 1000, 1100, 2);
  add(1, 100, 200, 3);
  add(0, 50, 150, 1);
  add(0, 50, 200, 2);  // same start, longer end
  add(1, 500, 700, 1); // same start as another, longer end

  std::vector<std::string> names = {"chr1", "chr2"};
  auto it = macs3::WrapVectorFragmentIterator(std::move(in), std::move(names));
  assert(it != nullptr);

  std::vector<macs3::FragmentRecord> out;
  macs3::FragmentRecord rec;
  while (it->Next(&rec)) out.push_back(rec);
  assert(it->Error().empty());
  assert(out.size() == 6);

  // Expected order: chrom 0 by (start, end), then chrom 1 by (start, end).
  assert(out[0].chrom_id == 0 && out[0].start == 50 && out[0].end == 150);
  assert(out[1].chrom_id == 0 && out[1].start == 50 && out[1].end == 200);
  assert(out[2].chrom_id == 0 && out[2].start == 1000 && out[2].end == 1100);
  assert(out[3].chrom_id == 1 && out[3].start == 100 && out[3].end == 200);
  assert(out[4].chrom_id == 1 && out[4].start == 500 && out[4].end == 600);
  assert(out[5].chrom_id == 1 && out[5].start == 500 && out[5].end == 700);

  assert(it->ChromNames().size() == 2);
  assert(it->ChromNames()[0] == "chr1");
  assert(it->ChromNames()[1] == "chr2");

  std::printf("PASS test_vector_fragment_iterator (6 records, sorted by chrom+start+end)\n");
  return 0;
}
