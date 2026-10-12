#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include <unistd.h>

#include "rapidmacs/fragments.h"
#include "rapidmacs/macs3_frag_peak_pipeline.h"
#include "rapidmacs/macs3_frag_sweep.h"

namespace {

void Require(bool ok, const std::string& message) {
  if (!ok) {
    std::fprintf(stderr, "FAIL test_frag_max_count: %s\n", message.c_str());
    std::exit(1);
  }
}

std::string Read(const std::string& path) {
  std::ifstream in(path.c_str(), std::ios::binary);
  Require(in.good(), "cannot read " + path);
  return std::string(std::istreambuf_iterator<char>(in),
                     std::istreambuf_iterator<char>());
}

}  // namespace

int main() {
  using namespace chromap::peaks;
  const char* test_root = std::getenv("RAPIDMACS_TEST_ROOT");
  const std::string root = test_root ? test_root : ".";
  const std::string input = root + "/tests/fixtures/capped_counts.fragments.tsv";
  char temp[] = "/tmp/rapidmacs_count_paths_XXXXXX";
  Require(mkdtemp(temp) != nullptr, "mkdtemp failed");
  const std::string work(temp);
  const std::string np = work + "/peaks.narrowPeak";
  const std::string sm = work + "/summits.bed";
  std::vector<ChromFragments> fragments;
  Require(LoadFragmentsFromTsv(input, &fragments), "load fixture");

  for (int cap : {0, 1, 2, 5}) {
    Macs3FragPeakPipelineParams params;
    params.frag_pileup_max_count = cap;
    params.name_prefix = "capped";
    params.peak_caller_threads = 2;
    const std::string golden = root + "/tests/golden/capped_counts/" +
                               std::to_string(cap) + "/capped_";
    const std::string expected_np = Read(golden + "peaks.narrowPeak");
    const std::string expected_sm = Read(golden + "summits.bed");
    Require(!expected_np.empty() && !expected_sm.empty(), "empty golden");

    for (int path = 0; path < 4; ++path) {
      Macs3FragPeakPipelinePaths paths;
      std::string error;
      bool ok = false;
      if (path == 0) {
        ok = RunMacs3FragPeakPipelineFromFragments(
            &fragments, params, paths, np, sm, "", "", nullptr, &error,
            nullptr, false);
      } else if (path == 1 || path == 2) {
        auto iter = macs3::OpenFileFragmentIterator(input);
        Require(iter != nullptr, "open iterator");
        if (path == 1) {
          ok = RunMacs3FragPeakPipelineFromIterator(
              *iter, params, paths, np, sm, "", "", nullptr, &error);
        } else {
          ok = RunMacs3FragPeakPipelineFromSortedIterator(
              *iter, params, paths, np, sm, "", "", nullptr, &error);
        }
      } else {
        // Exercise the public sweep workspace separately: the sorted
        // iterator currently uses a batch endpoint merge instead.
        Macs3FragWorkspaceParams wp;
        wp.frag_pileup_max_count = cap;
        Macs3FragSweepWorkspace sweep;
        Require(sweep.Init(wp), "initialize sweep");
        std::vector<std::string> names;
        for (size_t i = 0; i < fragments.size(); ++i) {
          names.push_back(fragments[i].name);
          for (const Fragment& f : fragments[i].frags) {
            Require(sweep.AddFragment(i, f.start, f.end, f.count),
                    "add sweep fragment");
          }
        }
        sweep.SetChromNames(names);
        Require(sweep.Finalize(), "finalize sweep");
        Macs3FragPeakWorkspace workspace;
        Require(sweep.MoveInto(&workspace), "move sweep");
        ok = RunMacs3FragPeakPipelineFromWorkspace(
            &workspace, params, paths, np, sm, "", "", nullptr, &error);
      }
      const std::string label = "cap=" + std::to_string(cap) +
                                " path=" + std::to_string(path);
      Require(ok, label + ": " + error);
      Require(Read(np) == expected_np, label + " narrowPeak differs from MACS3");
      Require(Read(sm) == expected_sm, label + " summits differ from MACS3");
    }
  }
  std::remove(np.c_str());
  std::remove(sm.c_str());
  rmdir(work.c_str());
  std::printf("PASS test_frag_max_count (4 entry points, caps 0/1/2/5, MACS3 goldens)\n");
  return 0;
}
