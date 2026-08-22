#include "rapidmacs/bdgpeakcall.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

void Fail(const std::string& msg) {
  std::fprintf(stderr, "FAIL test_qvalue_threshold: %s\n", msg.c_str());
  std::exit(1);
}

std::vector<std::string> SplitTabs(const std::string& line) {
  std::vector<std::string> fields;
  std::stringstream ss(line);
  std::string field;
  while (std::getline(ss, field, '\t')) {
    fields.push_back(field);
  }
  return fields;
}

std::vector<std::string> ReadFirstNarrowPeakRow(const std::string& path) {
  std::ifstream in(path.c_str());
  std::string line;
  if (!std::getline(in, line)) {
    Fail("missing narrowPeak row at " + path);
  }
  return SplitTabs(line);
}

}  // namespace

int main() {
  chromap::peaks::Macs3BdgTrack ppois;
  ppois.p = {100, 200, 300, 1000};
  ppois.v = {0.0f, 2.0f, 5.0f, 0.0f};

  chromap::peaks::Macs3BdgTrack treat;
  treat.p = {100, 200, 300, 1000};
  treat.v = {0.0f, 10.0f, 20.0f, 0.0f};

  chromap::peaks::Macs3BdgTrack lambda;
  lambda.p = {1000};
  lambda.v = {1.0f};

  std::vector<chromap::peaks::Macs3BdgTrack> all_ppois;
  all_ppois.push_back(ppois);
  std::map<double, double> pqtable;
  chromap::peaks::BuildMacs3PqTable(all_ppois, &pqtable);

  chromap::peaks::Macs3BdgTrack qscore;
  if (!chromap::peaks::BuildMacs3QScoreTrackFromPpois(ppois, pqtable, &qscore)) {
    Fail("BuildMacs3QScoreTrackFromPpois failed");
  }
  const float q_cutoff = static_cast<float>(-std::log10(0.05));
  if (qscore.v.size() != 4 || !(qscore.v[1] < q_cutoff) ||
      !(qscore.v[2] >= q_cutoff)) {
    Fail("synthetic qscore track does not straddle q=0.05 cutoff");
  }

  std::vector<chromap::peaks::Macs3FragNarrowPeakRow> p_rows;
  if (!chromap::peaks::CollectMacs3FragNarrowPeakRowsFromThresholdTrack(
          "chr1", ppois, ppois, treat, lambda, 2.0f, 1, 0, 0.0, &p_rows)) {
    Fail("p-mode row collection failed");
  }
  if (p_rows.size() != 1 || p_rows[0].start != 200 || p_rows[0].end != 300) {
    std::stringstream msg;
    msg << "p-mode did not exclude the pscore exactly equal to its cutoff"
        << " (rows=" << p_rows.size();
    if (!p_rows.empty()) {
      msg << ", first=" << p_rows[0].start << "-" << p_rows[0].end;
    }
    msg << ")";
    Fail(msg.str());
  }

  // The optimized q-value path converts the q cutoff to the first pscore key
  // whose qscore passes. MACS3 uses qscore > cutoff, not >=, so an exact
  // qscore equality must also be excluded.
  std::map<double, double> exact_q_table;
  exact_q_table[0.0] = 0.0;
  exact_q_table[2.0] = 2.0;
  exact_q_table[5.0] = 3.0;
  std::vector<chromap::peaks::Macs3FragNarrowPeakRow> exact_q_rows;
  if (!chromap::peaks::CollectMacs3FragNarrowPeakRowsFromTracks(
          "chr1", ppois, treat, lambda, 2.0f, 1, 0, 0.0, &exact_q_rows,
          nullptr, true, &exact_q_table)) {
    Fail("optimized q-mode equality row collection failed");
  }
  if (exact_q_rows.size() != 1 || exact_q_rows[0].start != 200 ||
      exact_q_rows[0].end != 300) {
    std::stringstream msg;
    msg << "q-mode did not exclude the qscore exactly equal to its cutoff"
        << " (rows=" << exact_q_rows.size();
    if (!exact_q_rows.empty()) {
      msg << ", first=" << exact_q_rows[0].start << "-"
          << exact_q_rows[0].end;
    }
    msg << ")";
    Fail(msg.str());
  }

  std::vector<chromap::peaks::Macs3FragNarrowPeakRow> q_rows;
  if (!chromap::peaks::CollectMacs3FragNarrowPeakRowsFromThresholdTrack(
          "chr1", qscore, ppois, treat, lambda, q_cutoff, 1, 0, 0.0,
          &q_rows)) {
    Fail("q-mode row collection failed");
  }
  if (q_rows.size() != 1 || q_rows[0].start != 200 || q_rows[0].end != 300) {
    std::stringstream msg;
    msg << "q-mode did not threshold on qscore-derived regions"
        << " (rows=" << q_rows.size();
    if (!q_rows.empty()) {
      msg << ", first=" << q_rows[0].start << "-" << q_rows[0].end;
    }
    msg << ")";
    Fail(msg.str());
  }

  // `macs3 callpeak` applies its default --fe-cutoff 1 after peak calling.
  // This matters at permissive p-value thresholds: a zero-treatment region
  // against a tiny positive lambda has a large upper-tail pscore but is a
  // depletion, not an enrichment peak.
  chromap::peaks::Macs3FragNarrowPeakRow depletion;
  depletion.fc = 0.9995;
  chromap::peaks::Macs3FragNarrowPeakRow enrichment;
  enrichment.fc = 1.0;
  std::vector<chromap::peaks::Macs3FragNarrowPeakRow> fc_rows;
  fc_rows.push_back(depletion);
  fc_rows.push_back(enrichment);
  chromap::peaks::FilterMacs3FragNarrowPeakRowsByFoldEnrichment(
      &fc_rows, 1.0);
  if (fc_rows.size() != 1 || fc_rows[0].fc != 1.0) {
    Fail("MACS3 default fold-enrichment filter was not reproduced");
  }

  // BED tag events can yield repeated score-track endpoints. Region
  // collection must ignore the resulting zero-width piece when summit
  // plateaus are split at paired treatment/lambda positions.
  chromap::peaks::Macs3BdgTrack repeated_ppois;
  repeated_ppois.p = {100, 105, 105, 200};
  repeated_ppois.v = {1.0f, 2.0f, 2.0f, 1.0f};
  chromap::peaks::Macs3BdgTrack repeated_treat;
  repeated_treat.p = {100, 105, 200};
  repeated_treat.v = {1.0f, 2.0f, 1.0f};
  chromap::peaks::Macs3BdgTrack repeated_lambda;
  repeated_lambda.p = {200};
  repeated_lambda.v = {0.5f};
  const std::vector<int32_t> repeated_positions = {100, 105, 200};
  std::vector<chromap::peaks::Macs3FragNarrowPeakRow> repeated_rows;
  if (!chromap::peaks::CollectMacs3FragNarrowPeakRowsFromTracks(
          "chr1", repeated_ppois, repeated_treat, repeated_lambda, 0.0f,
          1, 0, 0.0, &repeated_rows, &repeated_positions)) {
    Fail("repeated-endpoint BED region collection failed");
  }

  chromap::peaks::Macs3FragNarrowPeakRow row;
  row.chrom = "chr1";
  row.start = 100;
  row.end = 200;
  row.summit = 150;
  row.peak_off = 50;
  row.fc = 3.0;
  row.p_mlog = 5.0;

  std::map<double, double> score_table;
  score_table[5.0] = 2.5;

  const std::string p_path = "/tmp/rapidmacs_qvalue_threshold_p.narrowPeak";
  const std::string q_path = "/tmp/rapidmacs_qvalue_threshold_q.narrowPeak";
  std::remove(p_path.c_str());
  std::remove(q_path.c_str());

  std::vector<chromap::peaks::Macs3FragNarrowPeakRow> write_p;
  write_p.push_back(row);
  if (!chromap::peaks::WriteMacs3FragNarrowPeakRows(
          &write_p, p_path, "", &score_table, "", false)) {
    Fail("p-mode narrowPeak write failed");
  }
  std::vector<chromap::peaks::Macs3FragNarrowPeakRow> write_q;
  write_q.push_back(row);
  if (!chromap::peaks::WriteMacs3FragNarrowPeakRows(
          &write_q, q_path, "", &score_table, "", true)) {
    Fail("q-mode narrowPeak write failed");
  }

  const std::vector<std::string> p_fields = ReadFirstNarrowPeakRow(p_path);
  const std::vector<std::string> q_fields = ReadFirstNarrowPeakRow(q_path);
  if (p_fields.size() != 10 || q_fields.size() != 10) {
    Fail("unexpected narrowPeak field count");
  }
  if (p_fields[4] != "50") {
    Fail("p-mode narrowPeak score column did not use pscore");
  }
  if (q_fields[4] != "25") {
    Fail("q-mode narrowPeak score column did not use qscore");
  }

  std::remove(p_path.c_str());
  std::remove(q_path.c_str());
  std::printf("PASS test_qvalue_threshold\n");
  return 0;
}
