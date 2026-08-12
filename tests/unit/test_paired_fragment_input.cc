#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <unistd.h>

#include "rapidmacs/fragment_input.h"

namespace {

bool Check(bool condition, const char* message) {
  if (!condition) std::cerr << "FAIL: " << message << "\n";
  return condition;
}

}  // namespace

int main() {
  const std::string base =
      "/tmp/rapidmacs_paired_input_" + std::to_string(getpid());
  const std::string sam = base + ".sam";
  const std::string bedpe = base + ".bedpe";
  const std::string fragments = base + ".fragments.tsv";
  {
    std::ofstream o(sam.c_str());
    o << "@HD\tVN:1.6\tSO:coordinate\n"
      << "@SQ\tSN:chr2\tLN:1000000\n"
      << "@SQ\tSN:chr1\tLN:2000000\n"
      << "q1\t99\tchr2\t101\t60\t50M\t=\t151\t100\t*\t*\n"
      << "q1\t147\tchr2\t151\t60\t50M\t=\t101\t-100\t*\t*\n"
      // Duplicate-flagged records are parsed, then coordinate-deduplicated.
      << "q2\t1123\tchr2\t101\t60\t50M\t=\t151\t100\t*\t*\n"
      << "q2\t1171\tchr2\t151\t60\t50M\t=\t101\t-100\t*\t*\n"
      << "q9\t99\tchr2\t101\t60\t50M\t=\t151\t100\t*\t*\n"
      << "q9\t147\tchr2\t151\t60\t50M\t=\t101\t-100\t*\t*\n"
      // MAPQ zero is retained, matching MACS3 BAMPEParser.
      << "q3\t99\tchr2\t301\t0\t50M\t=\t351\t100\t*\t*\n"
      << "q3\t147\tchr2\t351\t0\t50M\t=\t301\t-100\t*\t*\n"
      << "q4\t611\tchr2\t501\t60\t50M\t=\t551\t100\t*\t*\n"
      << "q5\t355\tchr2\t601\t60\t50M\t=\t651\t100\t*\t*\n"
      << "q6\t2147\tchr2\t701\t60\t50M\t=\t751\t100\t*\t*\n"
      << "q7\t97\tchr2\t801\t60\t50M\t=\t851\t100\t*\t*\n"
      << "q8\t99\tchr1\t501\t60\t40M\t=\t541\t80\t*\t*\n"
      << "q8\t147\tchr1\t541\t60\t40M\t=\t501\t-80\t*\t*\n";
  }

  std::vector<chromap::peaks::ChromFragments> bam_chroms;
  chromap::peaks::FragmentInputStats bam_stats;
  bool pass = Check(chromap::peaks::LoadFragmentsFromBampe(
                        sam, &bam_chroms, &bam_stats),
                    "BAMPE/SAM load") &&
              Check(bam_chroms.size() == 2, "BAMPE chromosome count") &&
              Check(bam_chroms[0].name == "chr1" &&
                        bam_chroms[1].name == "chr2",
                    "BAMPE lexicographic chromosome order") &&
              Check(bam_chroms[0].reference_length == 2000000 &&
                        bam_chroms[1].reference_length == 1000000,
                    "BAMPE reference lengths") &&
              Check(bam_chroms[0].frags.size() == 1 &&
                        bam_chroms[0].frags[0].start == 500 &&
                        bam_chroms[0].frags[0].end == 580,
                    "BAMPE chr1 template interval") &&
              Check(bam_chroms[1].frags.size() == 2 &&
                        bam_chroms[1].frags[0].start == 100 &&
                        bam_chroms[1].frags[0].end == 200 &&
                        bam_chroms[1].frags[1].start == 300 &&
                        bam_chroms[1].frags[1].end == 400,
                    "BAMPE filtering, MAPQ, sorting, and dedup") &&
              Check(bam_stats.fragments_parsed == 5 &&
                        bam_stats.fragments_retained == 3 &&
                        bam_stats.duplicates_removed == 2,
                    "BAMPE accounting");

  std::vector<chromap::peaks::ChromFragments> bam_keep_two;
  chromap::peaks::FragmentInputStats bam_keep_two_stats;
  pass = Check(chromap::peaks::LoadFragmentsFromBampe(
                   sam, &bam_keep_two, &bam_keep_two_stats, 1, 2),
               "BAMPE --keep-dup 2 load") &&
         Check(bam_keep_two_stats.fragments_retained == 4 &&
                   bam_keep_two_stats.duplicates_removed == 1,
               "BAMPE --keep-dup 2 accounting") &&
         pass;

  std::vector<chromap::peaks::ChromFragments> bam_keep_all;
  chromap::peaks::FragmentInputStats bam_keep_all_stats;
  pass = Check(chromap::peaks::LoadFragmentsFromBampe(
                   sam, &bam_keep_all, &bam_keep_all_stats, 1, -1),
               "BAMPE --keep-dup all load") &&
         Check(bam_keep_all_stats.fragments_retained == 5 &&
                   bam_keep_all_stats.duplicates_removed == 0,
               "BAMPE --keep-dup all accounting") &&
         pass;

  {
    std::ofstream o(bedpe.c_str());
    o << "track name=test\n"
      << "chr2\t100\t200\n"
      << "chr2\t100\t200\n"
      << "chr2\t100\t200\n"
      << "chr1\t500\t550\tchr1\t580\t600\tq1\t60\t+\t-\n"
      << "chr1\t1\t20\tchr2\t30\t50\tcross\n"
      << "malformed\n";
  }
  std::vector<chromap::peaks::ChromFragments> bedpe_chroms;
  chromap::peaks::FragmentInputStats bedpe_stats;
  pass = Check(chromap::peaks::LoadFragmentsFromBedpe(
                   bedpe, &bedpe_chroms, &bedpe_stats),
               "BEDPE load") &&
         Check(bedpe_chroms.size() == 2 &&
                   bedpe_chroms[0].name == "chr1" &&
                   bedpe_chroms[1].name == "chr2",
               "BEDPE chromosome order") &&
         Check(bedpe_chroms[0].frags[0].start == 500 &&
                   bedpe_chroms[0].frags[0].end == 600,
               "standard six-column BEDPE conversion") &&
         Check(bedpe_stats.fragments_parsed == 4 &&
                   bedpe_stats.fragments_retained == 2 &&
                   bedpe_stats.duplicates_removed == 2 &&
                   bedpe_stats.records_skipped == 2,
               "BEDPE filtering and accounting") &&
         pass;

  std::vector<chromap::peaks::ChromFragments> bedpe_keep_two;
  chromap::peaks::FragmentInputStats bedpe_keep_two_stats;
  pass = Check(chromap::peaks::LoadFragmentsFromBedpe(
                   bedpe, &bedpe_keep_two, &bedpe_keep_two_stats, 2),
               "BEDPE --keep-dup 2 load") &&
         Check(bedpe_keep_two_stats.fragments_retained == 3 &&
                   bedpe_keep_two_stats.duplicates_removed == 1,
               "BEDPE --keep-dup 2 accounting") &&
         pass;

  std::vector<chromap::peaks::ChromFragments> bedpe_keep_all;
  chromap::peaks::FragmentInputStats bedpe_keep_all_stats;
  pass = Check(chromap::peaks::LoadFragmentsFromBedpe(
                   bedpe, &bedpe_keep_all, &bedpe_keep_all_stats, -1),
               "BEDPE --keep-dup all load") &&
         Check(bedpe_keep_all_stats.fragments_retained == 4 &&
                   bedpe_keep_all_stats.duplicates_removed == 0 &&
                   bedpe_keep_all[1].frags.size() == 3,
               "BEDPE --keep-dup all retains coordinate copies") &&
         pass;

  {
    std::ofstream o(fragments.c_str());
    o << "chr2\t100\t200\tBC1\t3\n"
      << "chr1\t500\t650\tN\t60\t+\t7\n"
      << "chr1\t700\t850\tN\t42\t-\t2\n";
  }
  std::vector<chromap::peaks::ChromFragments> fragment_chroms;
  pass = Check(chromap::peaks::LoadFragmentsFromTsv(fragments,
                                                     &fragment_chroms),
               "FRAG TSV load") &&
         Check(fragment_chroms.size() == 2 &&
                   fragment_chroms[0].name == "chr1" &&
                   fragment_chroms[1].name == "chr2",
               "FRAG TSV chromosome order") &&
         Check(fragment_chroms[0].frags.size() == 2 &&
                   fragment_chroms[0].frags[0].count == 7 &&
                   fragment_chroms[0].frags[1].count == 2 &&
                   fragment_chroms[1].frags[0].count == 3,
               "five-column and seven-column duplicate counts") &&
         pass;

  std::remove(sam.c_str());
  std::remove(bedpe.c_str());
  std::remove(fragments.c_str());
  if (!pass) return 1;
  std::cout << "PASS test_paired_fragment_input (BAMPE + BEDPE + FRAG schemas)\n";
  return 0;
}
