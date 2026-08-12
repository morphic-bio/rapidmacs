#include "rapidmacs/signac_bed_projection.h"

#include <cassert>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

void WriteFile(const std::string& path, const std::string& text) {
  std::ofstream out(path.c_str());
  assert(out.good());
  out << text;
}

std::string ReadFile(const std::string& path) {
  std::ifstream in(path.c_str());
  assert(in.good());
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

int CountLines(const std::string& s) {
  int n = 0;
  for (char c : s) {
    if (c == '\n') ++n;
  }
  return n;
}

}  // namespace

int main() {
  const std::string frag = "/tmp/rapidmacs_projection_fragments.tsv";
  const std::string cells = "/tmp/rapidmacs_projection_cells.tsv";
  const std::string groups = "/tmp/rapidmacs_projection_groups.tsv";
  const std::string bed = "/tmp/rapidmacs_projection_single.bed";
  const std::string dir = "/tmp/rapidmacs_projection_groups";

  WriteFile(frag,
            "# comment\n"
            "chr1\t10\t20\tBC1\t2\n"
            "chr1\t30\t40\tBC2\t1\n"
            "chr2\t5\t15\tBC3\t3\n"
            "chr2\t20\t30\tBC1\t1\n");
  WriteFile(cells, "barcode\nBC1\nBC3\n");
  WriteFile(groups, "barcode\tmacs_group\nBC1\tA\nBC3\tB\n");

  chromap::peaks::SignacBedProjectionParams p;
  p.fragments_path = frag;
  p.cells_path = cells;
  p.bed_out = bed;
  chromap::peaks::SignacBedProjectionResult r;
  std::string err;
  assert(chromap::peaks::ProjectFragmentsToSignacBed(p, &r, &err));
  assert(err.empty());
  assert(r.input_fragment_rows == 4);
  assert(r.written_bed_rows == 3);
  assert(r.selected_cell_count == 2);
  assert(ReadFile(bed) ==
         "chr1\t10\t20\tBC1\t2\n"
         "chr2\t5\t15\tBC3\t3\n"
         "chr2\t20\t30\tBC1\t1\n");

  chromap::peaks::SignacBedProjectionParams gp;
  gp.fragments_path = frag;
  gp.cells_path = cells;
  gp.barcode_group_tsv_path = groups;
  gp.group_column = "macs_group";
  gp.bed_out_dir = dir;
  gp.expand_counts = true;
  chromap::peaks::SignacBedProjectionResult gr;
  err.clear();
  assert(chromap::peaks::ProjectFragmentsToSignacBed(gp, &gr, &err));
  assert(err.empty());
  assert(gr.group_count == 2);
  assert(gr.written_bed_rows == 6);
  assert(gr.group_rows["A"] == 3);
  assert(gr.group_rows["B"] == 3);
  assert(CountLines(ReadFile(gr.group_bed_paths["A"])) == 3);
  assert(CountLines(ReadFile(gr.group_bed_paths["B"])) == 3);

  std::remove(frag.c_str());
  std::remove(cells.c_str());
  std::remove(groups.c_str());
  std::remove(bed.c_str());
  std::remove(gr.group_bed_paths["A"].c_str());
  std::remove(gr.group_bed_paths["B"].c_str());
  std::printf("PASS test_signac_bed_projection\n");
  return 0;
}
