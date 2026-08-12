#include "rapidmacs/bdgpeakcall.h"
#include "rapidmacs/macs3_frag_workspace.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

void Fail(const std::string& message) {
  std::fprintf(stderr, "FAIL test_nolambda_call_summits: %s\n",
               message.c_str());
  std::exit(1);
}

std::vector<std::string> ReadNames(const std::string& path) {
  std::ifstream in(path.c_str());
  std::vector<std::string> names;
  std::string line;
  while (std::getline(in, line)) {
    std::stringstream ss(line);
    std::string field;
    for (int i = 0; i < 4 && std::getline(ss, field, '\t'); ++i) {
      if (i == 3) names.push_back(field);
    }
  }
  return names;
}

}  // namespace

int main() {
  using chromap::peaks::ChromFragments;
  using chromap::peaks::Fragment;
  using chromap::peaks::Macs3BdgTrack;
  using chromap::peaks::Macs3FragNarrowPeakRow;
  using chromap::peaks::Macs3FragPeakWorkspace;
  using chromap::peaks::Macs3FragWorkspaceParams;

  // --nolambda must not allocate treatment local-window events and must emit
  // a single constant background segment ending with the treatment track.
  ChromFragments chrom;
  chrom.name = "chr1";
  Fragment first;
  first.start = 100;
  first.end = 200;
  first.count = 2;
  chrom.frags.push_back(first);
  Fragment second;
  second.start = 150;
  second.end = 250;
  second.count = 1;
  chrom.frags.push_back(second);
  std::vector<ChromFragments> fragments(1, chrom);
  Macs3FragWorkspaceParams params;
  params.effective_genome_size = 1000;
  params.build_treatment_local_lambda_events = false;
  params.capture_raw_pos_array_endpoints = false;
  params.byte_parity_summits = true;
  Macs3FragPeakWorkspace workspace;
  if (!chromap::peaks::BuildMacs3FragWorkspaceFromFragments(
          fragments, params, &workspace)) {
    Fail("workspace build failed");
  }
  if (!workspace.events_by_chrom[0].lambda_end_window_events.empty()) {
    Fail("--nolambda retained local-window events");
  }
  if (!chromap::peaks::FinalizeMacs3FragTreatTracks(&workspace) ||
      !chromap::peaks::FinalizeMacs3FragGlobalLambdaTracks(&workspace)) {
    Fail("global-lambda finalization failed");
  }
  if (workspace.lambda_tracks.size() != 1 ||
      workspace.lambda_tracks[0].p.size() != 1 ||
      workspace.lambda_tracks[0].p[0] != 250 ||
      std::fabs(workspace.lambda_tracks[0].v[0] - 0.3f) > 1e-6f) {
    Fail("global lambda is not total_fragment_bp/effective_genome_size");
  }
  if (workspace.pos_array_positions[0] != workspace.treat_tracks[0].p) {
    Fail("--nolambda position array is not the treatment breakpoint array");
  }

  // Frozen MACS3 subpeak case. The old single-summit method chooses the
  // midpoint of the 7-read plateau (offset 113); --call-summits selects the
  // smoothed derivative maximum at offset 105.
  Macs3BdgTrack ppois;
  ppois.p = {1000, 1237, 2000};
  ppois.v = {0.f, 10.f, 0.f};
  Macs3BdgTrack treat;
  treat.p = {1000, 1003, 1022, 1040, 1069, 1087, 1139,
             1150, 1153, 1163, 1172, 1190, 1219, 1237, 2000};
  treat.v = {0.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f,
             6.f, 5.f, 4.f, 5.f, 4.f, 3.f, 2.f, 0.f};
  Macs3BdgTrack lambda;
  lambda.p = {2000};
  lambda.v = {0.02f};
  const std::vector<int32_t> positions = treat.p;
  std::vector<Macs3FragNarrowPeakRow> single;
  if (!chromap::peaks::CollectMacs3FragNarrowPeakRowsFromTracks(
          "chr1", ppois, treat, lambda, 5.f, 150, 0, 0.0, &single,
          &positions, true, nullptr, false) ||
      single.size() != 1 || single[0].summit != 1113) {
    Fail("single-summit baseline changed");
  }
  std::vector<Macs3FragNarrowPeakRow> subpeaks;
  if (!chromap::peaks::CollectMacs3FragNarrowPeakRowsFromTracks(
          "chr1", ppois, treat, lambda, 5.f, 150, 0, 0.0, &subpeaks,
          &positions, true, nullptr, true) ||
      subpeaks.size() != 1 || subpeaks[0].summit != 1105) {
    Fail("MACS3 call-summits derivative maximum was not reproduced");
  }

  // Multiple rows for one region share a peak number and receive MACS3's
  // recursive subpeak suffixes; the following region increments once.
  Macs3FragNarrowPeakRow a = subpeaks[0];
  Macs3FragNarrowPeakRow b = a;
  b.summit += 20;
  b.peak_off += 20;
  Macs3FragNarrowPeakRow c = a;
  c.start = 1400;
  c.end = 1600;
  c.summit = 1500;
  c.peak_off = 100;
  std::vector<Macs3FragNarrowPeakRow> named = {a, b, c};
  const std::string np = "/tmp/rapidmacs_nolambda_subpeaks.narrowPeak";
  const std::string sm = "/tmp/rapidmacs_nolambda_subpeaks.summits.bed";
  std::remove(np.c_str());
  std::remove(sm.c_str());
  if (!chromap::peaks::WriteMacs3FragNarrowPeakRows(
          &named, np, sm, nullptr, "sample", false)) {
    Fail("subpeak writer failed");
  }
  const std::vector<std::string> np_names = ReadNames(np);
  const std::vector<std::string> sm_names = ReadNames(sm);
  const std::vector<std::string> expected = {
      "sample_peak_1a", "sample_peak_1b", "sample_peak_2"};
  if (np_names != expected || sm_names != expected) {
    Fail("subpeak names do not match MACS3 PeakIO naming");
  }
  std::remove(np.c_str());
  std::remove(sm.c_str());
  std::printf("PASS test_nolambda_call_summits\n");
  return 0;
}
