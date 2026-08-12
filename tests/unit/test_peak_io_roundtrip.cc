// Unit: WriteNarrowPeak emits 10 tab-separated columns in the documented
// MACS3 narrowPeak order, with peaks sorted by (chrom, start, end).
#include "rapidmacs/peak_io.h"

#include <cassert>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

int main() {
  using chromap::peaks::Peak;
  using chromap::peaks::WriteNarrowPeak;

  std::vector<Peak> peaks;

  // Deliberately out-of-order; WriteNarrowPeak should sort by chrom/start/end.
  Peak b;
  b.chrom = "chr1";
  b.start = 2000;
  b.end = 2150;
  b.peak_offset = 50;
  b.max_signal = 12;
  b.p_value = 1e-5;
  b.q_value = 1e-3;
  peaks.push_back(b);

  Peak a;
  a.chrom = "chr1";
  a.start = 1000;
  a.end = 1100;
  a.peak_offset = 25;
  a.max_signal = 7;
  a.p_value = 1e-2;
  a.q_value = 1e-1;
  peaks.push_back(a);

  const std::string path = "/tmp/rapidmacs_test_peaks.narrowPeak";
  bool ok = WriteNarrowPeak(path, peaks);
  assert(ok);

  std::ifstream in(path);
  assert(in);
  std::vector<std::string> lines;
  std::string line;
  while (std::getline(in, line)) lines.push_back(line);
  assert(lines.size() == 2);

  // Sorted: chr1:1000 first, chr1:2000 second.
  // Each line: 10 tab-separated columns.
  for (size_t i = 0; i < lines.size(); ++i) {
    int tabs = 0;
    for (char c : lines[i]) if (c == '\t') ++tabs;
    assert(tabs == 9);  // 9 tabs => 10 columns
  }

  // Column 1 (chrom), 2 (start), 4 (name), 6 (strand) checks.
  std::istringstream s0(lines[0]);
  std::string chrom, start, end, name, score, strand, signal, logp, logq, off;
  std::getline(s0, chrom, '\t');
  std::getline(s0, start, '\t');
  std::getline(s0, end, '\t');
  std::getline(s0, name, '\t');
  std::getline(s0, score, '\t');
  std::getline(s0, strand, '\t');
  std::getline(s0, signal, '\t');
  std::getline(s0, logp, '\t');
  std::getline(s0, logq, '\t');
  std::getline(s0, off, '\t');

  assert(chrom == "chr1");
  assert(start == "1000");
  assert(end == "1100");
  assert(name == "peak_1");
  assert(strand == ".");
  assert(off == "25");
  assert(signal == "7");

  // Second peak: peak_2 (renumbered after sort).
  std::istringstream s1(lines[1]);
  std::string c2, st2, en2, nm2;
  std::getline(s1, c2, '\t');
  std::getline(s1, st2, '\t');
  std::getline(s1, en2, '\t');
  std::getline(s1, nm2, '\t');
  assert(c2 == "chr1" && st2 == "2000" && en2 == "2150" && nm2 == "peak_2");

  std::remove(path.c_str());
  std::printf("PASS test_peak_io_roundtrip (2 peaks, sorted, 10 cols each)\n");
  return 0;
}
