// RapidMACS: standalone MACS3-compatible paired-fragment narrow peak caller.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cerrno>
#include <cctype>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <sys/stat.h>
#include <sys/types.h>

#include "rapidmacs/binned_signal.h"
#include "rapidmacs/bdgpeakcall.h"
#include "rapidmacs/call_peaks.h"
#include "rapidmacs/fragment_input.h"
#include "rapidmacs/fragments.h"
#include "rapidmacs/macs3_bed_callpeak.h"
#include "rapidmacs/macs3_frag_peak_pipeline.h"
#include "rapidmacs/peak_io.h"
#include "rapidmacs/signac_bed_projection.h"
#include "rapidmacs/stage_profile.h"
#include "rapidmacs/streaming.h"

namespace {

void Usage() {
  std::cerr
      << "Usage: rapidmacs -t <treatment> [-c <control>] "
         "[--out-prefix <path/prefix>] [options]\n"
         "  --treatment|-t, --input|-i  Treatment fragments/BAMPE/BEDPE input\n"
         "  --control|-c PATH            Optional matched input/IgG control\n"
         "  --format|-f {FRAG,BAMPE,BEDPE}  Paired-end input format (default: FRAG;\n"
         "                accepts 5-col barcoded or 7-col Chromap bulk fragments)\n"
         "  --out-prefix   Writes <prefix>.narrowPeak and companion outputs\n"
         "  --name|-n NAME Peaks are named \"<NAME>_peak_<i>\". "
         "Callpeak output defaults to \"NA\", matching MACS3.\n"
         "  --input-fragments PATH  Project fragments to Signac/MACS BED input\n"
         "  --cells PATH   Optional barcode allow-list for --input-fragments\n"
         "  --barcode-group-tsv PATH  Optional barcode/group table for per-group BED projection\n"
         "  --group-column NAME  Group column when --barcode-group-tsv has a header\n"
         "  --macs-profile signac-atac  Resolve BED callpeak defaults: BED, --nomodel,\n"
         "                --extsize 200 --shift -100 -g 2700000000 -q 0.05 --llocal 10000\n"
         "  --macs-bed-input PATH  Existing BED/tag input for BED callpeak compatibility\n"
         "  --macs-input-format BED  Required format for the BED compatibility path\n"
         "  --macs-bed-out PATH  Projected BED output for single-sample projection\n"
         "  --macs-bed-out-dir DIR  Projected per-group BED output directory\n"
         "  --macs-bed-peaks-out PATH  BED callpeak narrowPeak output\n"
         "  --macs-bed-summits-out PATH  BED callpeak summits BED output\n"
         "  --effective-genome-size N  BED callpeak -g effective genome size\n"
         "  --qvalue Q / --pvalue P  BED callpeak q/p threshold (mutually exclusive)\n"
         "  --nomodel --extsize N --shift N --llocal N  BED callpeak options\n"
         "  --nolambda  Use genome-wide background only; skip local-window work\n"
         "  --call-summits  Call multiple subpeak summits (paired-end and BED)\n"
         "  --keep-dup {all,N}  Keep all or at most positive N duplicate coordinates\n"
         "                (BAMPE/BEDPE and BED; default: 1). FRAG always keeps all.\n"
         "  --macs-bed-summary PATH  Resolved BED callpeak parameters and projection stats\n"
         "  --pileup-bdg   Binned Tn5-extended *cut* pileup bedGraph (diagnostic; "
         "differs from MACS3 `pileup -f FRAG`)\n"
         "  --frag-span-pileup-bdg  Fragment-span pileup bedGraph ([start,end) weighted by "
         "col5; matches MACS3 `pileup -f FRAG` when --frag-pileup-macs3-uint8-counts is set)\n"
         "  --frag-pileup-max-count N  Cap col5 like MACS3 --max-count (0 = no cap; applies "
         "to --frag-span-pileup-bdg only)\n"
         "  --frag-pileup-macs3-uint8-counts  Coerce fragment weights through uint8 (MACS3 "
         "PETrackII); recommended for byte-for-byte parity vs macs3 pileup\n"
         "  --no-macs3-frag-summit-parity  Disable byte-identical-summit mode (default ON);\n"
         "                takes a smaller-memory shortcut that diverges from the MACS3 binary on\n"
         "                ~2% of summit positions at PBMC scale (max ~222 bp; BED3/bdg unaffected)\n"
         "  --frag-lambda-bdg  MACS3 `control_lambda` diagnostic (matched control, or "
         "treatment-derived llocal when no control is supplied)\n"
         "  --frag-lambda-effective-genome N  -g size (default 2913022398 = hs GRCh38)\n"
         "  --gsize|-g {hs,mm,ce,dm,N}  Effective genome size\n"
         "  --frag-lambda-llocal  Large local window in bp (default: 10000, MACS3 llocal)\n"
         "  --frag-lambda-slocal N  Small local control window (default: 1000)\n"
         "  --frag-score-ppois-bdg PATH  MACS3 `bdgcmp -m ppois` score bedGraph (diagnostic)\n"
         "  --frag-score-fe-bdg PATH     MACS3 `bdgcmp -m FE` score bedGraph (diagnostic)\n"
         "  --frag-score-pseudocount X   Same as `macs3 bdgcmp -p` (default: 0)\n"
         "  --bdgpeakcall-input PATH     Score bedGraph; run MACS3-style bdgpeakcall (diagnostic)\n"
         "  --bdgpeakcall-output PATH    BED3 regions (requires --bdgpeakcall-input)\n"
         "  --preset NAME                Apply a named parameter preset, then continue\n"
         "                parsing flags. Place this BEFORE per-flag overrides if you\n"
         "                want explicit flags to win. Supported: shiftTAG, chip, cutrun.\n"
         "                chip/cutrun select BAMPE, -q 0.05, and MACS3's dynamic\n"
         "                paired-end peak geometry; chip requires -c. shiftTAG sets\n"
         "                --bdgpeakcall-cutoff 2 (= -log10(0.01)) so the FRAG\n"
         "                narrowPeak path matches MACS3 callpeak with -p 0.01\n"
         "                --extsize 150 --shift -75 --keep-dup all --nomodel\n"
         "                --call-summits -g hs (other defaults already align).\n"
         "  --bdgpeakcall-cutoff X       -c cutoff (default: 5)\n"
         "  --macs3-frag-pvalue P        Raw p-value threshold for --macs3-frag-* "
         "(sets --bdgpeakcall-cutoff=-log10(P); default: 1e-5)\n"
         "  --macs3-frag-qvalue Q        q-value/FDR threshold for --macs3-frag-* "
         "(mutually exclusive with explicit --macs3-frag-pvalue)\n"
         "  --bdgpeakcall-min-len N      -l min length (default: 200)\n"
         "  --bdgpeakcall-max-gap N      -g max gap (default: 30)\n"
         "  --macs3-frag-narrowpeak PATH  Diagnostic narrowPeak (FRAG: ppois+lambda+span;\n"
         "                pileup/lambda/ppois bedGraphs are written automatically if omitted)\n"
         "  --macs3-frag-summits PATH    Summits BED5 (one row per called summit)\n"
         "  --macs3-frag-ppois-bdg PATH  Override ppois bedGraph for --macs3-frag-* (else --frag-score-ppois-bdg)\n"
         "  --bin-width    Bin width in bp (default: 50)\n"
         "  --ext-size     Tn5 extension in bp, symmetric; default 150 (2x75; shiftTAG-like)\n"
         "  --local-window  Half-window in bins for local lambda (default: 200; ~20kb for 50bp)\n"
         "  --fdr          FDR (BH) threshold (default: 0.05)\n"
         "  --p-value      Raw p upper tail threshold (default: 0.01)\n"
         "  --merge-gap     Merge significant bins separated by up to this many bins (0)\n"
         "  --min-peak-bp   Drop peaks narrower than this width in bp (0 = off)\n"
         "  --min-summit-cuts  Drop peaks whose summit bin count is below this (0 = off)\n"
         "  --profile-stages PATH  Optional stage timing TSV (wall/user/sys/RSS/rows/bytes)\n"
         "  --peak-caller-threads N  Threads for per-chromosome MACS3 FRAG stages (default: 1)\n"
         "  -q|--qvalue Q  MACS3 q-value cutoff (default in chip/cutrun: 0.05)\n"
         "  -p|--pvalue P  MACS3 p-value cutoff\n"
         "  --scale-to {small,large}  Treatment/control depth target (default: small)\n"
         "  --to-large     Alias for --scale-to large\n"
         "  -h|--help      This message\n";
}

std::string JoinArgv(int argc, char** argv) {
  std::ostringstream o;
  for (int i = 0; i < argc; ++i) {
    if (i) {
      o << ' ';
    }
    o << argv[i];
  }
  return o.str();
}

std::string UpperAscii(std::string s) {
  for (char& c : s) {
    if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
  }
  return s;
}

bool LoadInput(const std::string& format, const std::string& path,
               std::vector<chromap::peaks::ChromFragments>* out,
               chromap::peaks::FragmentInputStats* stats,
               int num_threads, int max_duplicates) {
  if (format == "BAMPE") {
    return chromap::peaks::LoadFragmentsFromBampe(path, out, stats,
                                                   num_threads, max_duplicates);
  }
  if (format == "BEDPE") {
    return chromap::peaks::LoadFragmentsFromBedpe(path, out, stats,
                                                   max_duplicates);
  }
  if (format == "FRAG") {
    const bool ok = chromap::peaks::LoadFragmentsFromTsv(path, out);
    if (stats) {
      *stats = chromap::peaks::FragmentInputStats();
      for (const chromap::peaks::ChromFragments& ch : *out) {
        for (const chromap::peaks::Fragment& f : ch.frags) {
          stats->records_read += 1;
          stats->fragments_parsed += f.count;
          stats->fragments_retained += f.count;
          stats->fragment_bp_parsed +=
              static_cast<int64_t>(f.end - f.start) * f.count;
        }
      }
    }
    return ok;
  }
  return false;
}

bool ParseKeepDupMax(const std::string& value, int* max_duplicates) {
  if (max_duplicates == nullptr) return false;
  if (value == "all") {
    *max_duplicates = -1;
    return true;
  }
  if (value.empty()) return false;
  for (char c : value) {
    if (c < '0' || c > '9') return false;
  }
  char* end = nullptr;
  errno = 0;
  const long parsed = std::strtol(value.c_str(), &end, 10);
  if (errno != 0 || end == value.c_str() || *end != '\0' || parsed < 1 ||
      parsed > std::numeric_limits<int>::max()) {
    return false;
  }
  *max_duplicates = static_cast<int>(parsed);
  return true;
}

bool EnsureDir(const std::string& path) {
  if (path.empty() || path == ".") {
    return true;
  }
  struct stat st;
  if (stat(path.c_str(), &st) == 0) {
    return S_ISDIR(st.st_mode);
  }
  if (mkdir(path.c_str(), 0775) == 0) {
    return true;
  }
  if (errno == EEXIST && stat(path.c_str(), &st) == 0) {
    return S_ISDIR(st.st_mode);
  }
  return false;
}

std::string PathJoin(const std::string& dir, const std::string& name) {
  if (dir.empty() || dir == ".") {
    return name;
  }
  if (!dir.empty() && dir.back() == '/') {
    return dir + name;
  }
  return dir + "/" + name;
}

std::string SanitizeForPath(const std::string& s) {
  std::string out;
  out.reserve(s.size());
  for (char ch : s) {
    const unsigned char c = static_cast<unsigned char>(ch);
    if (std::isalnum(c) || ch == '_' || ch == '-' || ch == '.') {
      out.push_back(ch);
    } else {
      out.push_back('_');
    }
  }
  return out.empty() ? "group" : out;
}

int64_t ParseGenomeSize(const std::string& s) {
  if (s == "hs") return 2913022398LL;
  if (s == "mm") return 2652783500LL;
  if (s == "ce") return 100286401LL;
  if (s == "dm") return 142573017LL;
  const double v = std::atof(s.c_str());
  return static_cast<int64_t>(v + 0.5);
}

const char* ThresholdModeName(chromap::peaks::Macs3FragThresholdMode mode) {
  return mode == chromap::peaks::Macs3FragThresholdMode::kQValue ? "qvalue"
                                                                 : "pvalue";
}

bool WriteMacsBedSummary(
    const std::string& path, const std::string& argv_line,
    const std::string& profile, const std::string& bed_input,
    const std::string& narrowpeak_out, const std::string& summits_out,
    const chromap::peaks::Macs3BedCallPeakParams& params,
    const chromap::peaks::Macs3BedCallPeakStats& stats,
    const chromap::peaks::SignacBedProjectionResult* projection) {
  if (path.empty()) {
    return true;
  }
  std::ofstream out(path.c_str());
  if (!out) {
    return false;
  }
  out << "argv\t" << argv_line << "\n";
  out << "macs_profile\t" << (profile.empty() ? "none" : profile) << "\n";
  out << "macs_input_format\tBED\n";
  out << "bed_input\t" << bed_input << "\n";
  out << "narrowpeak_out\t" << narrowpeak_out << "\n";
  out << "summits_out\t" << summits_out << "\n";
  out << "nomodel\t" << (params.nomodel ? "true" : "false") << "\n";
  out << "extsize\t" << params.extsize << "\n";
  out << "shift\t" << params.shift << "\n";
  out << "effective_genome_size\t" << params.effective_genome_size << "\n";
  out << "threshold_mode\t" << ThresholdModeName(params.threshold_mode) << "\n";
  out << "qvalue\t" << params.qvalue_cutoff << "\n";
  out << "pvalue\t" << params.pvalue_cutoff << "\n";
  out << "keep_dup\t" << params.keep_dup << "\n";
  out << "keep_dup_max\t" << stats.keep_dup_max << "\n";
  out << "llocal\t" << params.llocal_bp << "\n";
  out << "broad\t" << (params.broad ? "true" : "false") << "\n";
  out << "broad_cutoff\t" << params.broad_cutoff << "\n";
  out << "nolambda\t" << (params.no_lambda ? "true" : "false") << "\n";
  out << "call_summits\t" << (params.call_summits ? "true" : "false") << "\n";
  out << "detected_tag_size\t" << stats.detected_tag_size << "\n";
  out << "resolved_min_length\t" << stats.resolved_min_length << "\n";
  out << "resolved_max_gap\t" << stats.resolved_max_gap << "\n";
  out << "total_tags\t" << stats.total_tags << "\n";
  out << "tags_after_filtering\t" << stats.tags_after_filtering << "\n";
  out << "peak_count\t" << stats.peak_count << "\n";
  if (projection != nullptr) {
    out << "projection_input_fragment_rows\t" << projection->input_fragment_rows << "\n";
    out << "projection_written_bed_rows\t" << projection->written_bed_rows << "\n";
    out << "projection_selected_cell_count\t" << projection->selected_cell_count << "\n";
    out << "projection_group_count\t" << projection->group_count << "\n";
    out << "projection_count_expansion_used\t"
        << (projection->count_expansion_used ? "true" : "false") << "\n";
    for (const auto& kv : projection->group_rows) {
      out << "projection_group_rows." << kv.first << "\t" << kv.second << "\n";
    }
  }
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  std::string input;
  std::string control;
  std::string input_format = "FRAG";
  std::string input_fragments;
  std::string cells_path;
  std::string barcode_group_tsv;
  std::string group_column;
  std::string macs_profile;
  std::string macs_input_format;
  std::string macs_bed_input;
  std::string macs_bed_out;
  std::string macs_bed_out_dir;
  std::string macs_bed_peaks_out;
  std::string macs_bed_summits_out;
  std::string macs_bed_summary;
  std::string outdir = ".";
  bool macs_bed_expand_counts = false;
  bool macs_bed_nomodel = false;
  bool macs_bed_nomodel_explicit = false;
  bool macs_bed_extsize_explicit = false;
  bool macs_bed_shift_explicit = false;
  bool macs_bed_egs_explicit = false;
  bool macs_bed_qvalue_explicit = false;
  bool macs_bed_pvalue_explicit = false;
  bool macs_bed_keep_dup_explicit = false;
  bool macs_bed_llocal_explicit = false;
  bool macs_bed_broad_cutoff_explicit = false;
  chromap::peaks::Macs3FragThresholdMode macs_bed_threshold_mode =
      chromap::peaks::Macs3FragThresholdMode::kQValue;
  int32_t macs_bed_extsize = 200;
  int32_t macs_bed_shift = 0;
  int64_t macs_bed_effective_genome_size = 2913022398LL;
  double macs_bed_qvalue = 0.05;
  double macs_bed_pvalue = 0.0;
  std::string macs_bed_keep_dup = "1";
  int32_t macs_bed_llocal = 10000;
  bool macs_bed_broad = false;
  double macs_bed_broad_cutoff = 0.1;
  bool macs_no_lambda = false;
  bool macs_bed_call_summits = false;
  std::string macs_bed_treat_bdg;
  std::string macs_bed_lambda_bdg;
  std::string macs_bed_ppois_bdg;
  std::string out_prefix;
  std::string name_prefix;
  std::string pileup_bdg;
  std::string frag_span_pileup_bdg;
  std::string frag_lambda_bdg;
  std::string frag_score_ppois_bdg;
  std::string frag_score_fe_bdg;
  int32_t frag_pileup_max_count = 0;
  bool frag_pileup_macs3_uint8 = false;
  bool macs3_frag_summit_parity = true;  // default-on as of 2026-04-27
  int64_t frag_lambda_effective_genome = 2913022398LL;  // MACS3 Constants.EFFECTIVEGS["hs"]
  int32_t frag_lambda_slocal = 1000;
  int32_t frag_lambda_llocal = 10000;
  bool scale_to_large = false;
  bool cutoff_is_qvalue = false;
  bool cutoff_configured = false;
  bool preset_requires_control = false;
  bool auto_min_length = false;
  bool auto_max_gap = false;
  bool min_length_explicit = false;
  bool max_gap_explicit = false;
  double frag_score_pseudocount = 0.0;
  std::string bdgpeakcall_input;
  std::string bdgpeakcall_output;
  float bdgpeakcall_cutoff = 5.f;
  double macs3_frag_qvalue = 0.0;
  bool macs3_frag_pvalue_explicit = false;
  bool macs3_frag_qvalue_explicit = false;
  bool bdgpeakcall_cutoff_explicit = false;
  int32_t bdgpeakcall_min_len = 200;
  int32_t bdgpeakcall_max_gap = 30;
  std::string macs3_frag_narrowpeak_out;
  std::string macs3_frag_summits_out;
  std::string macs3_frag_ppois_override;
  int32_t bin_width = 50;
  int32_t ext_size = 150;
  int32_t local_window = 200;
  double fdr = 0.05;
  double pval = 0.01;
  int32_t merge_gap = 0;
  int32_t min_peak_bp = 0;
  int64_t min_summit_cuts = 0;
  std::string profile_stages_path;
  int peak_caller_threads = 1;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "-h" || a == "--help") {
      Usage();
      return 0;
    }
    if (a == "--input" || a == "-i" || a == "--treatment" || a == "-t") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      input = argv[++i];
    } else if (a == "--control" || a == "-c") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      control = argv[++i];
    } else if (a == "--format" || a == "-f") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      input_format = UpperAscii(argv[++i]);
      if (input_format != "FRAG" && input_format != "BAMPE" &&
          input_format != "BEDPE") {
        std::cerr << "Unsupported format: " << input_format
                  << " (expected FRAG, BAMPE, or BEDPE)\n";
        return 2;
      }
    } else if (a == "--input-fragments") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      input_fragments = argv[++i];
    } else if (a == "--cells") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      cells_path = argv[++i];
    } else if (a == "--barcode-group-tsv") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      barcode_group_tsv = argv[++i];
    } else if (a == "--group-column") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      group_column = argv[++i];
    } else if (a == "--macs-profile") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs_profile = argv[++i];
    } else if (a == "--macs-input-format") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs_input_format = argv[++i];
    } else if (a == "--macs-bed-input") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs_bed_input = argv[++i];
    } else if (a == "--macs-bed-out") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs_bed_out = argv[++i];
    } else if (a == "--macs-bed-out-dir") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs_bed_out_dir = argv[++i];
    } else if (a == "--macs-bed-expand-counts") {
      macs_bed_expand_counts = true;
    } else if (a == "--macs-bed-peaks-out") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs_bed_peaks_out = argv[++i];
    } else if (a == "--macs-bed-summits-out") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs_bed_summits_out = argv[++i];
    } else if (a == "--macs-bed-summary") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs_bed_summary = argv[++i];
    } else if (a == "--outdir") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      outdir = argv[++i];
    } else if (a == "--effective-genome-size") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs_bed_effective_genome_size = ParseGenomeSize(argv[++i]);
      macs_bed_egs_explicit = true;
    } else if (a == "--qvalue") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs_bed_qvalue = std::atof(argv[++i]);
      if (!(macs_bed_qvalue > 0.0 && macs_bed_qvalue <= 1.0)) {
        std::cerr << "--qvalue must be in (0,1]\n";
        return 2;
      }
      macs_bed_qvalue_explicit = true;
      macs_bed_threshold_mode = chromap::peaks::Macs3FragThresholdMode::kQValue;
      bdgpeakcall_cutoff =
          chromap::peaks::BdgPeakCallCutoffFromPValue(macs_bed_qvalue);
      cutoff_is_qvalue = true;
      cutoff_configured = true;
    } else if (a == "--pvalue") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs_bed_pvalue = std::atof(argv[++i]);
      if (!(macs_bed_pvalue > 0.0 && macs_bed_pvalue <= 1.0)) {
        std::cerr << "--pvalue must be in (0,1]\n";
        return 2;
      }
      macs_bed_pvalue_explicit = true;
      macs_bed_threshold_mode = chromap::peaks::Macs3FragThresholdMode::kPValue;
      bdgpeakcall_cutoff =
          chromap::peaks::BdgPeakCallCutoffFromPValue(macs_bed_pvalue);
      cutoff_is_qvalue = false;
      cutoff_configured = true;
    } else if (a == "--nomodel") {
      macs_bed_nomodel = true;
      macs_bed_nomodel_explicit = true;
    } else if (a == "--extsize") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs_bed_extsize = std::atoi(argv[++i]);
      macs_bed_extsize_explicit = true;
    } else if (a == "--shift") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs_bed_shift = std::atoi(argv[++i]);
      macs_bed_shift_explicit = true;
    } else if (a == "--keep-dup") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs_bed_keep_dup = argv[++i];
      macs_bed_keep_dup_explicit = true;
    } else if (a == "--llocal") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs_bed_llocal = std::atoi(argv[++i]);
      macs_bed_llocal_explicit = true;
    } else if (a == "--broad") {
      macs_bed_broad = true;
    } else if (a == "--broad-cutoff") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs_bed_broad_cutoff = std::atof(argv[++i]);
      macs_bed_broad_cutoff_explicit = true;
    } else if (a == "--nolambda") {
      macs_no_lambda = true;
    } else if (a == "--call-summits") {
      macs_bed_call_summits = true;
    } else if (a == "--macs-bed-treat-bdg") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs_bed_treat_bdg = argv[++i];
    } else if (a == "--macs-bed-lambda-bdg") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs_bed_lambda_bdg = argv[++i];
    } else if (a == "--macs-bed-ppois-bdg") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs_bed_ppois_bdg = argv[++i];
    } else if (a == "--out-prefix") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      out_prefix = argv[++i];
    } else if (a == "--name" || a == "-n") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      name_prefix = argv[++i];
    } else if (a == "--pileup-bdg") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      pileup_bdg = argv[++i];
    } else if (a == "--frag-span-pileup-bdg") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      frag_span_pileup_bdg = argv[++i];
    } else if (a == "--frag-pileup-max-count") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      frag_pileup_max_count = std::atoi(argv[++i]);
    } else if (a == "--frag-pileup-macs3-uint8-counts") {
      frag_pileup_macs3_uint8 = true;
    } else if (a == "--macs3-frag-summit-parity") {
      macs3_frag_summit_parity = true;  // accepted for back-compat; default-on
    } else if (a == "--no-macs3-frag-summit-parity") {
      macs3_frag_summit_parity = false;
    } else if (a == "--frag-lambda-bdg") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      frag_lambda_bdg = argv[++i];
    } else if (a == "--frag-lambda-effective-genome") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      frag_lambda_effective_genome = std::atoll(argv[++i]);
    } else if (a == "--gsize" || a == "-g") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      const std::string g = argv[++i];
      if (g == "hs") {
        frag_lambda_effective_genome = 2913022398LL;
      } else if (g == "mm") {
        frag_lambda_effective_genome = 2652783500LL;
      } else if (g == "ce") {
        frag_lambda_effective_genome = 100286401LL;
      } else if (g == "dm") {
        frag_lambda_effective_genome = 142573017LL;
      } else {
        frag_lambda_effective_genome = std::atoll(g.c_str());
      }
    } else if (a == "--frag-lambda-llocal") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      frag_lambda_llocal = std::atoi(argv[++i]);
    } else if (a == "--frag-lambda-slocal") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      frag_lambda_slocal = std::atoi(argv[++i]);
    } else if (a == "--qvalue" || a == "-q") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      const double q = std::atof(argv[++i]);
      if (!(q > 0.0 && q <= 1.0)) {
        std::cerr << "-q/--qvalue must be in (0,1]\n";
        return 2;
      }
      bdgpeakcall_cutoff = chromap::peaks::BdgPeakCallCutoffFromPValue(q);
      cutoff_is_qvalue = true;
      cutoff_configured = true;
      macs_bed_qvalue = q;
      macs_bed_qvalue_explicit = true;
      macs_bed_threshold_mode = chromap::peaks::Macs3FragThresholdMode::kQValue;
    } else if (a == "--pvalue" || a == "-p") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      const double p = std::atof(argv[++i]);
      if (!(p > 0.0 && p <= 1.0)) {
        std::cerr << "-p/--pvalue must be in (0,1]\n";
        return 2;
      }
      bdgpeakcall_cutoff = chromap::peaks::BdgPeakCallCutoffFromPValue(p);
      cutoff_is_qvalue = false;
      cutoff_configured = true;
      macs_bed_pvalue = p;
      macs_bed_pvalue_explicit = true;
      macs_bed_threshold_mode = chromap::peaks::Macs3FragThresholdMode::kPValue;
    } else if (a == "--scale-to") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      const std::string target = argv[++i];
      if (target == "small") {
        scale_to_large = false;
      } else if (target == "large") {
        scale_to_large = true;
      } else {
        std::cerr << "--scale-to must be small or large\n";
        return 2;
      }
    } else if (a == "--to-large") {
      scale_to_large = true;
    } else if (a == "--frag-score-ppois-bdg") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      frag_score_ppois_bdg = argv[++i];
    } else if (a == "--frag-score-fe-bdg") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      frag_score_fe_bdg = argv[++i];
    } else if (a == "--frag-score-pseudocount") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      frag_score_pseudocount = std::atof(argv[++i]);
    } else if (a == "--bdgpeakcall-input") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      bdgpeakcall_input = argv[++i];
    } else if (a == "--bdgpeakcall-output") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      bdgpeakcall_output = argv[++i];
    } else if (a == "--preset") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      const std::string preset = argv[++i];
      if (preset == "shiftTAG" || preset == "shifttag" ||
          preset == "shift-tag" || preset == "shift_tag") {
        // MACS3 callpeak shiftTAG settings: -p 0.01 → bdgpeakcall cutoff
        // 2.0 in -log10 space. Other shiftTAG defaults (--extsize 150,
        // -g hs, llocal=10000, -q 0.05, --nomodel, --keep-dup all) are
        // already the RapidMACS defaults.
        bdgpeakcall_cutoff = 2.f;
        cutoff_is_qvalue = false;
        cutoff_configured = true;
        input_format = "FRAG";
        preset_requires_control = false;
      } else if (preset == "chip") {
        input_format = "BAMPE";
        bdgpeakcall_cutoff =
            chromap::peaks::BdgPeakCallCutoffFromPValue(0.05);
        cutoff_is_qvalue = true;
        cutoff_configured = true;
        preset_requires_control = true;
        auto_min_length = true;
        auto_max_gap = true;
      } else if (preset == "cutrun" || preset == "cut-and-run" ||
                 preset == "cut_and_run") {
        input_format = "BAMPE";
        bdgpeakcall_cutoff =
            chromap::peaks::BdgPeakCallCutoffFromPValue(0.05);
        cutoff_is_qvalue = true;
        cutoff_configured = true;
        preset_requires_control = false;
        auto_min_length = true;
        auto_max_gap = true;
      } else {
        std::cerr << "Unknown --preset \"" << preset
                  << "\" (supported: shiftTAG, chip, cutrun)\n";
        return 2;
      }
    } else if (a == "--bdgpeakcall-cutoff") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      bdgpeakcall_cutoff = static_cast<float>(std::atof(argv[++i]));
      cutoff_is_qvalue = false;
      cutoff_configured = true;
      bdgpeakcall_cutoff_explicit = true;
    } else if (a == "--macs3-frag-pvalue") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs3_frag_pvalue_explicit = true;
      const double p = std::atof(argv[++i]);
      bdgpeakcall_cutoff = chromap::peaks::BdgPeakCallCutoffFromPValue(p);
    } else if (a == "--macs3-frag-qvalue") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs3_frag_qvalue_explicit = true;
      macs3_frag_qvalue = std::atof(argv[++i]);
    } else if (a == "--bdgpeakcall-min-len") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      bdgpeakcall_min_len = std::atoi(argv[++i]);
      min_length_explicit = true;
      auto_min_length = false;
    } else if (a == "--bdgpeakcall-max-gap") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      bdgpeakcall_max_gap = std::atoi(argv[++i]);
      max_gap_explicit = true;
      auto_max_gap = false;
    } else if (a == "--macs3-frag-narrowpeak") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs3_frag_narrowpeak_out = argv[++i];
    } else if (a == "--macs3-frag-summits") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs3_frag_summits_out = argv[++i];
    } else if (a == "--macs3-frag-ppois-bdg") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs3_frag_ppois_override = argv[++i];
    } else if (a == "--bin-width") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      bin_width = std::atoi(argv[++i]);
    } else if (a == "--ext-size") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      ext_size = std::atoi(argv[++i]);
    } else if (a == "--local-window") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      local_window = std::atoi(argv[++i]);
    } else if (a == "--fdr") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      fdr = std::atof(argv[++i]);
    } else if (a == "--p-value") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      pval = std::atof(argv[++i]);
    } else if (a == "--merge-gap") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      merge_gap = std::atoi(argv[++i]);
    } else if (a == "--min-peak-bp") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      min_peak_bp = std::atoi(argv[++i]);
    } else if (a == "--min-summit-cuts") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      min_summit_cuts = std::atoll(argv[++i]);
    } else if (a == "--profile-stages") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      profile_stages_path = argv[++i];
    } else if (a == "--peak-caller-threads") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      peak_caller_threads = std::atoi(argv[++i]);
    } else {
      std::cerr << "Unknown argument: " << a << "\n";
      Usage();
      return 2;
    }
  }
  const bool paired_alignment_input =
      input_format == "BAMPE" || input_format == "BEDPE";
  if (paired_alignment_input) {
    if (!min_length_explicit) auto_min_length = true;
    if (!max_gap_explicit) auto_max_gap = true;
    if (!cutoff_configured) {
      bdgpeakcall_cutoff =
          chromap::peaks::BdgPeakCallCutoffFromPValue(0.05);
      cutoff_is_qvalue = true;
    }
  }
  if (preset_requires_control && control.empty()) {
    std::cerr << "--preset chip requires -c/--control\n";
    return 2;
  }
  if (!control.empty() && input.empty()) {
    std::cerr << "A treatment input is required with -c/--control\n";
    return 2;
  }
  if (!control.empty() && !frag_score_fe_bdg.empty()) {
    std::cerr << "--frag-score-fe-bdg is currently a no-control diagnostic; "
                 "it cannot be combined with -c/--control\n";
    return 2;
  }
  // For paired-end formats, --out-prefix denotes the validated MACS3 path,
  // not the historical binned ATAC caller.
  if ((paired_alignment_input || !control.empty()) && !out_prefix.empty()) {
    if (macs3_frag_narrowpeak_out.empty()) {
      macs3_frag_narrowpeak_out = out_prefix + ".narrowPeak";
    }
    if (macs3_frag_summits_out.empty()) {
      macs3_frag_summits_out = out_prefix + ".summits.bed";
    }
  }

  if (!macs_profile.empty()) {
    chromap::peaks::Macs3BedCallPeakProfile profile;
    std::string err;
    if (!chromap::peaks::ResolveMacs3BedCallPeakProfile(macs_profile, &profile,
                                                        &err)) {
      std::cerr << err << "\n";
      return 2;
    }
    if (macs_input_format.empty()) {
      macs_input_format = profile.input_format;
    }
    const chromap::peaks::Macs3BedCallPeakParams& pp = profile.params;
    if (!macs_bed_nomodel_explicit) {
      macs_bed_nomodel = pp.nomodel;
    }
    if (!macs_bed_extsize_explicit) {
      macs_bed_extsize = pp.extsize;
    }
    if (!macs_bed_shift_explicit) {
      macs_bed_shift = pp.shift;
    }
    if (!macs_bed_egs_explicit) {
      macs_bed_effective_genome_size = pp.effective_genome_size;
    }
    if (!macs_bed_qvalue_explicit && !macs_bed_pvalue_explicit) {
      macs_bed_threshold_mode = pp.threshold_mode;
      macs_bed_qvalue = pp.qvalue_cutoff;
      macs_bed_pvalue = pp.pvalue_cutoff;
    }
    if (!macs_bed_llocal_explicit) {
      macs_bed_llocal = pp.llocal_bp;
    }
    if (!macs_bed_keep_dup_explicit) {
      macs_bed_keep_dup = pp.keep_dup;
    }
    if (!macs_bed_broad_cutoff_explicit) {
      macs_bed_broad_cutoff = pp.broad_cutoff;
    }
  }

  // Shared MACS options such as --qvalue also belong to BAMPE/BEDPE. Enter
  // the BED compatibility surface only when a BED/projection-specific option
  // anchors that mode.
  const bool macs_bed_surface =
      !macs_profile.empty() || !input_fragments.empty() || !cells_path.empty() ||
      !barcode_group_tsv.empty() || !group_column.empty() ||
      !macs_input_format.empty() || !macs_bed_input.empty() ||
      !macs_bed_out.empty() || !macs_bed_out_dir.empty() ||
      !macs_bed_peaks_out.empty() || !macs_bed_summits_out.empty() ||
      !macs_bed_summary.empty() ||
      !macs_bed_treat_bdg.empty() || !macs_bed_lambda_bdg.empty() ||
      !macs_bed_ppois_bdg.empty();

  int max_duplicates = 1;
  if (!macs_bed_surface && input_format == "FRAG") {
    // MACS3's FRAG parser treats the count column as the already-collapsed
    // duplicate multiplicity. It accepts every --keep-dup spelling, warns,
    // and forces `all` instead of trying to filter rows again.
    if (macs_bed_keep_dup_explicit && macs_bed_keep_dup != "all") {
      std::cerr << "Since the format is 'FRAG', `--keep-dup` will be set as "
                   "'all'.\n";
    }
    max_duplicates = -1;
  } else if (!ParseKeepDupMax(macs_bed_keep_dup, &max_duplicates)) {
    std::cerr << "--keep-dup must be 'all' or a positive integer\n";
    return 2;
  }

  if (macs_bed_surface) {
    if (!macs_input_format.empty() && macs_input_format != "BED") {
      std::cerr << "--macs-input-format currently supports only BED for this "
                   "compatibility path\n";
      return 2;
    }
    if (macs_bed_qvalue_explicit && macs_bed_pvalue_explicit) {
      std::cerr << "--qvalue and --pvalue are mutually exclusive\n";
      return 2;
    }
    if (macs_bed_effective_genome_size < 1 || macs_bed_extsize < 1 ||
        macs_bed_llocal < 1) {
      std::cerr << "Invalid BED callpeak genome size, extsize, or llocal\n";
      return 2;
    }
    if (!EnsureDir(outdir)) {
      std::cerr << "Failed to create --outdir: " << outdir << "\n";
      return 1;
    }
    const std::string default_name = name_prefix.empty() ? "macs_bed" : name_prefix;
    const bool wants_call =
        !macs_profile.empty() || !macs_bed_input.empty() ||
        !macs_bed_peaks_out.empty() || !macs_bed_summits_out.empty() ||
        macs_bed_qvalue_explicit || macs_bed_pvalue_explicit ||
        macs_bed_nomodel_explicit || macs_bed_extsize_explicit ||
        macs_bed_shift_explicit || macs_bed_egs_explicit ||
        macs_bed_keep_dup_explicit || macs_bed_llocal_explicit ||
        macs_bed_broad || macs_bed_broad_cutoff_explicit ||
        macs_no_lambda || macs_bed_call_summits ||
        !macs_bed_treat_bdg.empty() ||
        !macs_bed_lambda_bdg.empty() || !macs_bed_ppois_bdg.empty();

    chromap::peaks::SignacBedProjectionResult projection;
    bool have_projection = false;
    std::string bed_input = macs_bed_input;
    const bool grouped = !barcode_group_tsv.empty();
    if (!input_fragments.empty()) {
      chromap::peaks::SignacBedProjectionParams pp;
      pp.fragments_path = input_fragments;
      pp.cells_path = cells_path;
      pp.barcode_group_tsv_path = barcode_group_tsv;
      pp.group_column = group_column;
      pp.expand_counts = macs_bed_expand_counts;
      if (grouped) {
        pp.bed_out_dir =
            macs_bed_out_dir.empty()
                ? PathJoin(outdir, SanitizeForPath(default_name) + "_signac_beds")
                : macs_bed_out_dir;
      } else {
        pp.bed_out =
            macs_bed_out.empty()
                ? (wants_call ? PathJoin(outdir, SanitizeForPath(default_name) +
                                                    "_signac_input.bed")
                              : "")
                : macs_bed_out;
      }
      std::string err;
      if (!chromap::peaks::ProjectFragmentsToSignacBed(pp, &projection, &err)) {
        std::cerr << err << "\n";
        return 1;
      }
      have_projection = true;
      if (!grouped) {
        bed_input = projection.bed_path;
      }
    }
    if (!wants_call) {
      return 0;
    }
    if (grouped) {
      if (!macs_bed_peaks_out.empty() || !macs_bed_summits_out.empty()) {
        std::cerr << "--macs-bed-peaks-out/--macs-bed-summits-out are single-output "
                     "paths; omit them for grouped calls\n";
        return 2;
      }
      for (const auto& kv : projection.group_bed_paths) {
        const std::string group = kv.first;
        const std::string group_base = SanitizeForPath(group);
        chromap::peaks::Macs3BedCallPeakParams bp;
        bp.threshold_mode = macs_bed_threshold_mode;
        bp.qvalue_cutoff = macs_bed_qvalue;
        bp.pvalue_cutoff = macs_bed_pvalue;
        bp.effective_genome_size = macs_bed_effective_genome_size;
        bp.nomodel = macs_bed_nomodel;
        bp.extsize = macs_bed_extsize;
        bp.shift = macs_bed_shift;
        bp.llocal_bp = macs_bed_llocal;
        bp.keep_dup = macs_bed_keep_dup;
        bp.broad = macs_bed_broad;
        bp.broad_cutoff = macs_bed_broad_cutoff;
        bp.no_lambda = macs_no_lambda;
        bp.call_summits = macs_bed_call_summits;
        bp.name_prefix = group;
        const std::string np = PathJoin(outdir, group_base + "_peaks.narrowPeak");
        const std::string sm = PathJoin(outdir, group_base + "_summits.bed");
        const std::string summary =
            macs_bed_summary.empty()
                ? PathJoin(outdir, group_base + "_chromap_macs_bed_summary.tsv")
                : macs_bed_summary + "." + group_base + ".tsv";
        chromap::peaks::Macs3BedCallPeakStats st;
        std::string err;
        if (!chromap::peaks::RunMacs3BedCallPeak(kv.second, bp, np, sm, &st, &err)) {
          std::cerr << err << "\n";
          return 1;
        }
        if (!WriteMacsBedSummary(summary, JoinArgv(argc, argv), macs_profile,
                                 kv.second, np, sm, bp, st,
                                 have_projection ? &projection : nullptr)) {
          std::cerr << "Failed to write BED callpeak summary: " << summary << "\n";
          return 1;
        }
      }
      return 0;
    }
    if (bed_input.empty()) {
      std::cerr << "BED callpeak requires --macs-bed-input or --input-fragments\n";
      return 2;
    }
    chromap::peaks::Macs3BedCallPeakParams bp;
    bp.threshold_mode = macs_bed_threshold_mode;
    bp.qvalue_cutoff = macs_bed_qvalue;
    bp.pvalue_cutoff = macs_bed_pvalue;
    bp.effective_genome_size = macs_bed_effective_genome_size;
    bp.nomodel = macs_bed_nomodel;
    bp.extsize = macs_bed_extsize;
    bp.shift = macs_bed_shift;
    bp.llocal_bp = macs_bed_llocal;
    bp.keep_dup = macs_bed_keep_dup;
    bp.broad = macs_bed_broad;
    bp.broad_cutoff = macs_bed_broad_cutoff;
    bp.no_lambda = macs_no_lambda;
    bp.call_summits = macs_bed_call_summits;
    bp.name_prefix = name_prefix;
    bp.treat_bdg = macs_bed_treat_bdg;
    bp.lambda_bdg = macs_bed_lambda_bdg;
    bp.ppois_bdg = macs_bed_ppois_bdg;
    const std::string np =
        macs_bed_peaks_out.empty()
            ? PathJoin(outdir, SanitizeForPath(default_name) + "_peaks.narrowPeak")
            : macs_bed_peaks_out;
    const std::string sm =
        macs_bed_summits_out.empty()
            ? PathJoin(outdir, SanitizeForPath(default_name) + "_summits.bed")
            : macs_bed_summits_out;
    const std::string summary =
        macs_bed_summary.empty()
            ? PathJoin(outdir, SanitizeForPath(default_name) +
                                  "_chromap_macs_bed_summary.tsv")
            : macs_bed_summary;
    chromap::peaks::Macs3BedCallPeakStats st;
    std::string err;
    if (!chromap::peaks::RunMacs3BedCallPeak(bed_input, bp, np, sm, &st, &err)) {
      std::cerr << err << "\n";
      return 1;
    }
    if (!WriteMacsBedSummary(summary, JoinArgv(argc, argv), macs_profile,
                             bed_input, np, sm, bp, st,
                             have_projection ? &projection : nullptr)) {
      std::cerr << "Failed to write BED callpeak summary: " << summary << "\n";
      return 1;
    }
    return 0;
  }
  const bool bdgpeakcall_pair =
      !bdgpeakcall_input.empty() || !bdgpeakcall_output.empty();
  if (bdgpeakcall_pair &&
      (bdgpeakcall_input.empty() || bdgpeakcall_output.empty())) {
    std::cerr << "Both --bdgpeakcall-input and --bdgpeakcall-output are required "
                 "together.\n";
    return 2;
  }
  const bool need_macs3_frag_np =
      !macs3_frag_narrowpeak_out.empty() || !macs3_frag_summits_out.empty();
  const bool need_fragment_input =
      !out_prefix.empty() || !pileup_bdg.empty() || !frag_span_pileup_bdg.empty() ||
      !frag_lambda_bdg.empty() || !frag_score_ppois_bdg.empty() ||
      !frag_score_fe_bdg.empty() || need_macs3_frag_np;
  if (need_fragment_input && input.empty()) {
    Usage();
    return 2;
  }
  if (input.empty() && !bdgpeakcall_pair) {
    Usage();
    return 2;
  }
  if (!bdgpeakcall_pair && out_prefix.empty() && pileup_bdg.empty() &&
      frag_span_pileup_bdg.empty() && frag_lambda_bdg.empty() &&
      frag_score_ppois_bdg.empty() && frag_score_fe_bdg.empty() &&
      !need_macs3_frag_np) {
    Usage();
    return 2;
  }
  const bool need_bins = !out_prefix.empty() || !pileup_bdg.empty();
  if (need_bins && (bin_width < 1 || ext_size < 2 || (ext_size % 2) != 0)) {
    std::cerr
        << "Invalid --bin-width or --ext-size (ext-size must be even, >=2)\n";
    return 2;
  }
  if (frag_pileup_max_count < 0) {
    std::cerr << "Invalid --frag-pileup-max-count\n";
    return 2;
  }
  if (min_peak_bp < 0 || min_summit_cuts < 0) {
    std::cerr << "Invalid --min-peak-bp or --min-summit-cuts\n";
    return 2;
  }
  if (frag_lambda_slocal < 0 || frag_lambda_llocal < 1 ||
      (!control.empty() && frag_lambda_slocal > frag_lambda_llocal) ||
      frag_lambda_effective_genome < 1) {
    std::cerr << "Invalid --frag-lambda-slocal/--frag-lambda-llocal or "
                 "--frag-lambda-effective-genome\n";
    return 2;
  }
  if (bdgpeakcall_min_len < 1 || bdgpeakcall_max_gap < 0) {
    std::cerr << "Invalid --bdgpeakcall-min-len or --bdgpeakcall-max-gap\n";
    return 2;
  }
  if (macs3_frag_pvalue_explicit && macs3_frag_qvalue_explicit) {
    std::cerr << "--macs3-frag-pvalue and --macs3-frag-qvalue are mutually exclusive\n";
    return 2;
  }
  if (bdgpeakcall_cutoff_explicit && macs3_frag_qvalue_explicit) {
    std::cerr << "--bdgpeakcall-cutoff and --macs3-frag-qvalue are mutually exclusive\n";
    return 2;
  }
  if (std::isnan(bdgpeakcall_cutoff) || bdgpeakcall_cutoff < 0.f) {
    std::cerr << "Invalid --macs3-frag-pvalue or --bdgpeakcall-cutoff\n";
    return 2;
  }
  if (macs3_frag_qvalue_explicit &&
      (std::isnan(macs3_frag_qvalue) || macs3_frag_qvalue <= 0.0 ||
       macs3_frag_qvalue > 1.0)) {
    std::cerr << "--macs3-frag-qvalue must be in (0, 1]\n";
    return 2;
  }
  if (peak_caller_threads < 1) {
    std::cerr << "Invalid --peak-caller-threads\n";
    return 2;
  }
  if (fdr < 0.0 || fdr > 1.0 || pval < 0.0 || pval > 1.0) {
    return 2;
  }
  std::unique_ptr<chromap::peaks::StageProfileCollector> profile_coll;
  std::unique_ptr<chromap::peaks::StageProfileTotalGuard> profile_total_guard;
  int64_t profile_frags = 0;
  int64_t profile_out_peaks = 0;
  if (!profile_stages_path.empty()) {
    profile_coll.reset(new chromap::peaks::StageProfileCollector(profile_stages_path));
    if (!profile_coll->IsOpen()) {
      std::cerr << "Failed to open --profile-stages output: " << profile_stages_path << "\n";
      return 1;
    }
    profile_total_guard.reset(new chromap::peaks::StageProfileTotalGuard(
        profile_coll.get(), chromap::peaks::StageProfileCollector::Now(), &profile_frags,
        &profile_out_peaks));
  }
  if (!bdgpeakcall_input.empty()) {
    if (!chromap::peaks::RunMacs3StyleBdgPeakCallFromBedGraph(
            bdgpeakcall_input, bdgpeakcall_output, bdgpeakcall_cutoff,
            bdgpeakcall_min_len, bdgpeakcall_max_gap)) {
      std::cerr << "Failed MACS3-style bdgpeakcall diagnostic\n";
      return 1;
    }
    if (!need_fragment_input) {
      return 0;
    }
  }
  // Streaming fast path: when only narrowPeak is requested (no diagnostic
  // bedGraph emission, no MACS3-style FRAG narrowPeak, no profile-stages
  // accounting), drain the fragments file through StreamingPeakCaller without
  // ever materializing the per-chrom fragment list. This is the default for
  // standalone rapidmacs invocations.
  const bool can_stream =
      input_format == "FRAG" && control.empty() && !out_prefix.empty() &&
      pileup_bdg.empty() && frag_span_pileup_bdg.empty() &&
      frag_lambda_bdg.empty() && frag_score_ppois_bdg.empty() &&
      frag_score_fe_bdg.empty() && !need_macs3_frag_np && !profile_coll;
  if (can_stream) {
    auto it = macs3::OpenFileFragmentIterator(input);
    if (!it || !it->Error().empty()) {
      std::cerr << "Failed to open fragments: " << input;
      if (it && !it->Error().empty()) std::cerr << ": " << it->Error();
      std::cerr << "\n";
      return 1;
    }
    macs3::StreamingParams sp;
    sp.bin_width = bin_width;
    sp.ext_size = ext_size;
    sp.local_window_bins = local_window;
    sp.fdr = fdr;
    sp.p_value = pval;
    sp.merge_bin_gap = merge_gap;
    macs3::StreamingPeakCaller scc(sp);
    std::string serr;
    std::vector<chromap::peaks::Peak> all;
    macs3::StreamingStats stream_stats;
    if (!scc.CallPeaks(*it, &all, &serr, &stream_stats)) {
      std::cerr << "fragments stream error: " << serr << "\n";
      return 1;
    }
    if (min_peak_bp > 0 || min_summit_cuts > 0) {
      std::vector<chromap::peaks::Peak> keep;
      keep.reserve(all.size());
      for (const auto& p : all) {
        if (p.end <= p.start) continue;
        if (min_peak_bp > 0 && (p.end - p.start) < min_peak_bp) continue;
        if (min_summit_cuts > 0 && p.max_signal < min_summit_cuts) continue;
        keep.push_back(p);
      }
      all.swap(keep);
    }
    const std::string np = out_prefix + ".narrowPeak";
    const std::string sm = out_prefix + ".summary.tsv";
    if (!chromap::peaks::WriteNarrowPeak(np, all, name_prefix)) {
      std::cerr << "Failed to write " << np << "\n";
      return 1;
    }
    std::ostringstream pstr;
    pstr.setf(std::ios::fixed);
    pstr << "bin_width=" << bin_width << " ext_size=" << ext_size
         << " local_window_bins=" << local_window
         << " fdr=" << std::setprecision(4) << fdr << " p_value=" << pval
         << " merge_bin_gap=" << merge_gap << " min_peak_bp=" << min_peak_bp
         << " min_summit_cuts=" << min_summit_cuts << " input=" << input;
    if (!chromap::peaks::WriteRunSummaryTsv(sm, all, JoinArgv(argc, argv),
                                            pstr.str(), stream_stats.total_bins,
                                            stream_stats.saturated_bins)) {
      std::cerr << "Failed to write " << sm << "\n";
      return 1;
    }
    return 0;
  }

  // MACS3 FRAG streaming fast path: when only --macs3-frag-narrowpeak is
  // requested (no flags that need a materialized fragment list — the legacy
  // custom narrowPeak via --out-prefix, --pileup-bdg, or --frag-score-fe-bdg
  // each need direct fragment access), drain the FragmentIterator straight
  // into the workspace via RunMacs3FragPeakPipelineFromIterator. The
  // intermediate treat/lambda/ppois bedGraphs (--frag-span-pileup-bdg,
  // --frag-lambda-bdg, --frag-score-ppois-bdg) still emit because the
  // pipeline writes them from the workspace tracks.
  const bool can_iter_macs3_frag =
      input_format == "FRAG" && control.empty() && need_macs3_frag_np &&
      out_prefix.empty() && pileup_bdg.empty() &&
      frag_score_fe_bdg.empty();
  if (can_iter_macs3_frag) {
    auto it = macs3::OpenFileFragmentIterator(input);
    if (!it || !it->Error().empty()) {
      std::cerr << "Failed to open fragments: " << input;
      if (it && !it->Error().empty()) std::cerr << ": " << it->Error();
      std::cerr << "\n";
      return 1;
    }
    chromap::peaks::Macs3FragPeakPipelineParams pr;
    pr.threshold_mode = macs3_frag_qvalue_explicit
                            ? chromap::peaks::Macs3FragThresholdMode::kQValue
                            : chromap::peaks::Macs3FragThresholdMode::kPValue;
    pr.bdgpeakcall_cutoff = bdgpeakcall_cutoff;
    pr.cutoff_is_qvalue = cutoff_is_qvalue;
    // CallPeakUnit stores pscore histogram keys and qscore values as float32
    // for FRAG as well as BAMPE/BEDPE.  Use that representation for every
    // callpeak-compatible path; the old round-5/double mode changed ATAC
    // q-values despite identical peak coordinates.
    pr.macs3_float_score_math = true;
    pr.qvalue_cutoff = macs3_frag_qvalue;
    pr.min_length = bdgpeakcall_min_len;
    pr.max_gap = bdgpeakcall_max_gap;
    pr.frag_pileup_max_count = frag_pileup_max_count;
    pr.macs3_uint8_counts = frag_pileup_macs3_uint8;
    pr.byte_parity_summits = macs3_frag_summit_parity;
    pr.effective_genome_size = frag_lambda_effective_genome;
    pr.slocal_bp = frag_lambda_slocal;
    pr.llocal_bp = frag_lambda_llocal;
    pr.scale_to_large = scale_to_large;
    pr.no_lambda = macs_no_lambda;
    pr.call_summits = macs_bed_call_summits;
    pr.score_pseudocount = frag_score_pseudocount;
    pr.peak_caller_threads = peak_caller_threads;
    pr.name_prefix = name_prefix.empty() ? "NA" : name_prefix;
    chromap::peaks::Macs3FragPeakPipelinePaths io;
    io.treat_bdg = frag_span_pileup_bdg;
    io.lambda_bdg = frag_lambda_bdg;
    io.ppois_bdg = !macs3_frag_ppois_override.empty() ? macs3_frag_ppois_override
                                                      : frag_score_ppois_bdg;
    std::string err;
    std::string work_dir;
    // FileFragmentIterator emits chrom-grouped + start-sorted output for
    // typical Chromap/ARC fragments TSVs, so we try the sweep-line workspace
    // first (~tens of KB peak memory vs the events-buffer's ~48 B/fragment).
    // If the input violates sort order, fall back to the events-buffer path
    // by re-opening the file. Other errors (open / parse / I/O) propagate.
    if (!chromap::peaks::RunMacs3FragPeakPipelineFromSortedIterator(
            *it, pr, io, macs3_frag_narrowpeak_out, macs3_frag_summits_out, "",
            "", &work_dir, &err, profile_coll.get())) {
      const bool sort_violation =
          err.find("not chrom-grouped") != std::string::npos ||
          err.find("not start-sorted") != std::string::npos;
      if (!sort_violation) {
        std::cerr << err << "\n";
        return 1;
      }
      const std::string sweep_err = err;
      err.clear();
      auto fallback = macs3::OpenFileFragmentIterator(input);
      if (!fallback || !fallback->Error().empty()) {
        std::cerr << "[rapidmacs] sweep failed (" << sweep_err
                  << "); fallback re-open failed";
        if (fallback && !fallback->Error().empty())
          std::cerr << ": " << fallback->Error();
        std::cerr << "\n";
        return 1;
      }
      std::cerr << "[rapidmacs] sweep-line workspace rejected input ("
                << sweep_err
                << "); falling back to events-buffer workspace\n";
      if (!chromap::peaks::RunMacs3FragPeakPipelineFromIterator(
              *fallback, pr, io, macs3_frag_narrowpeak_out,
              macs3_frag_summits_out, "", "", &work_dir, &err,
              profile_coll.get())) {
        std::cerr << err << "\n";
        return 1;
      }
    }
    return 0;
  }

  std::vector<chromap::peaks::ChromFragments> chs;
  std::vector<chromap::peaks::ChromFragments> control_chs;
  chromap::peaks::FragmentInputStats treatment_stats;
  chromap::peaks::FragmentInputStats control_stats;
  if (profile_coll) {
    const chromap::peaks::StageProfileCollector::Tick t0 =
        chromap::peaks::StageProfileCollector::Now();
    if (!LoadInput(input_format, input, &chs, &treatment_stats,
                   peak_caller_threads, max_duplicates)) {
      std::cerr << "Failed to read input: " << input << "\n";
      return 1;
    }
    if (!control.empty() &&
        !LoadInput(input_format, control, &control_chs, &control_stats,
                   peak_caller_threads, max_duplicates)) {
      std::cerr << "Failed to read control: " << control << "\n";
      return 1;
    }
    const chromap::peaks::StageProfileCollector::Tick t_load_end =
        chromap::peaks::StageProfileCollector::Now();
    for (const chromap::peaks::ChromFragments& ch : chs) {
      profile_frags += static_cast<int64_t>(ch.frags.size());
    }
    profile_coll->Record("load_fragments", t0, t_load_end, profile_frags,
                         static_cast<int64_t>(chs.size()), 0, "");
    const chromap::peaks::StageProfileCollector::Tick ts =
        chromap::peaks::StageProfileCollector::Now();
    profile_coll->Record("sort_fragments_by_chrom", ts, ts, profile_frags,
                         static_cast<int64_t>(chs.size()), 0,
                         "included_in_load_fragments;not_timed_separately");
  } else {
    if (!LoadInput(input_format, input, &chs, &treatment_stats,
                   peak_caller_threads, max_duplicates)) {
      std::cerr << "Failed to read input: " << input << "\n";
      return 1;
    }
    if (!control.empty() &&
        !LoadInput(input_format, control, &control_chs, &control_stats,
                   peak_caller_threads, max_duplicates)) {
      std::cerr << "Failed to read control: " << control << "\n";
      return 1;
    }
  }
  double treatment_fragment_length = 0.0;
  if (treatment_stats.fragments_parsed > 0) {
    if (input_format == "BEDPE") {
      // BEDPEParser.build_petrack explicitly casts the accumulated length to
      // float32 before division; BAMPEParser divides the integer totals and
      // only rounds on assignment. The one-ULP distinction propagates into
      // paired-control score formatting on deep inputs.
      treatment_fragment_length = static_cast<float>(
          static_cast<float>(treatment_stats.fragment_bp_parsed) /
          static_cast<float>(treatment_stats.fragments_parsed));
    } else {
      treatment_fragment_length =
          static_cast<double>(treatment_stats.fragment_bp_parsed) /
          static_cast<double>(treatment_stats.fragments_parsed);
    }
  }
  if (treatment_fragment_length <= 0.0) {
    std::cerr << "Treatment contains no usable paired-end fragments\n";
    return 1;
  }
  const int32_t macs_default_peak_geometry =
      std::max<int32_t>(1, static_cast<int32_t>(treatment_fragment_length));
  if (auto_min_length) bdgpeakcall_min_len = macs_default_peak_geometry;
  if (auto_max_gap) bdgpeakcall_max_gap = macs_default_peak_geometry;
  if (!control.empty() &&
      treatment_fragment_length > static_cast<double>(frag_lambda_llocal)) {
    std::cerr << "Treatment mean fragment length exceeds --frag-lambda-llocal\n";
    return 2;
  }
  if (!pileup_bdg.empty()) {
    if (!chromap::peaks::WritePileupBedGraph(pileup_bdg, chs, bin_width, ext_size)) {
      std::cerr << "Failed to write pileup bedGraph: " << pileup_bdg << "\n";
      return 1;
    }
  }
  if (need_macs3_frag_np) {
    chromap::peaks::Macs3FragPeakPipelineParams pr;
    pr.threshold_mode = macs3_frag_qvalue_explicit
                            ? chromap::peaks::Macs3FragThresholdMode::kQValue
                            : chromap::peaks::Macs3FragThresholdMode::kPValue;
    pr.bdgpeakcall_cutoff = bdgpeakcall_cutoff;
    pr.cutoff_is_qvalue = cutoff_is_qvalue;
    pr.macs3_float_score_math = true;
    pr.qvalue_cutoff = macs3_frag_qvalue;
    pr.min_length = bdgpeakcall_min_len;
    pr.max_gap = bdgpeakcall_max_gap;
    pr.frag_pileup_max_count = frag_pileup_max_count;
    pr.macs3_uint8_counts = frag_pileup_macs3_uint8;
    pr.byte_parity_summits = macs3_frag_summit_parity;
    pr.effective_genome_size = frag_lambda_effective_genome;
    pr.slocal_bp = frag_lambda_slocal;
    pr.llocal_bp = frag_lambda_llocal;
    pr.treatment_fragment_length_bp = treatment_fragment_length;
    pr.keep_duplicates_all = max_duplicates < 0;
    pr.scale_to_large = scale_to_large;
    pr.no_lambda = macs_no_lambda;
    pr.call_summits = macs_bed_call_summits;
    pr.score_pseudocount = frag_score_pseudocount;
    pr.peak_caller_threads = peak_caller_threads;
    pr.name_prefix = name_prefix.empty() ? "NA" : name_prefix;
    chromap::peaks::Macs3FragPeakPipelinePaths io;
    io.treat_bdg = frag_span_pileup_bdg;
    io.lambda_bdg = frag_lambda_bdg;
    io.ppois_bdg = !macs3_frag_ppois_override.empty() ? macs3_frag_ppois_override
                                                         : frag_score_ppois_bdg;
    std::string err;
    std::string work_dir;
    const bool release_fragments_after_workspace = frag_score_fe_bdg.empty();
    const std::vector<chromap::peaks::ChromFragments>* control_ptr =
        control.empty() ? nullptr : &control_chs;
    if (!chromap::peaks::RunMacs3FragPeakPipelineFromFragments(
            &chs, pr, io, macs3_frag_narrowpeak_out, macs3_frag_summits_out, "", "",
            &work_dir, &err, profile_coll.get(), release_fragments_after_workspace,
            control_ptr)) {
      std::cerr << err << "\n";
      return 1;
    }
    if (profile_coll && !macs3_frag_narrowpeak_out.empty()) {
      int64_t nb = 0, nl = 0;
      (void)chromap::peaks::StageProfileCollector::FileMetrics(macs3_frag_narrowpeak_out, &nb,
                                                               &nl);
      profile_out_peaks = nl;
    }
    if (!frag_score_fe_bdg.empty()) {
      if (!chromap::peaks::WriteFragMacs3BdgcmpScoreBedGraphs(
              "", frag_score_fe_bdg, chs, frag_pileup_max_count, frag_pileup_macs3_uint8,
              frag_lambda_effective_genome, frag_lambda_llocal, frag_score_pseudocount)) {
        std::cerr << "Failed to write MACS3-style FE score bedGraph\n";
        return 1;
      }
    }
  } else {
    if (!frag_span_pileup_bdg.empty()) {
      if (!chromap::peaks::WriteFragSpanPileupBedGraph(
              frag_span_pileup_bdg, chs, frag_pileup_max_count, frag_pileup_macs3_uint8)) {
        std::cerr << "Failed to write fragment-span pileup bedGraph: "
                  << frag_span_pileup_bdg << "\n";
        return 1;
      }
    }
    if (!frag_lambda_bdg.empty()) {
      if (!chromap::peaks::WriteFragMacs3NoControlLambdaBedGraph(
              frag_lambda_bdg, chs, frag_pileup_max_count, frag_pileup_macs3_uint8,
              frag_lambda_effective_genome, frag_lambda_llocal)) {
        std::cerr << "Failed to write frag lambda bedGraph: " << frag_lambda_bdg << "\n";
        return 1;
      }
    }
    if (!frag_score_ppois_bdg.empty() || !frag_score_fe_bdg.empty()) {
      if (!chromap::peaks::WriteFragMacs3BdgcmpScoreBedGraphs(
              frag_score_ppois_bdg, frag_score_fe_bdg, chs, frag_pileup_max_count,
              frag_pileup_macs3_uint8, frag_lambda_effective_genome, frag_lambda_llocal,
              frag_score_pseudocount)) {
        std::cerr << "Failed to write MACS3-style score bedGraph(s)\n";
        return 1;
      }
    }
  }
  if (paired_alignment_input || !control.empty()) {
    return 0;
  }
  if (out_prefix.empty()) {
    return 0;
  }
  std::vector<chromap::peaks::Peak> all;
  int64_t in_mem_total_bins = 0;
  int64_t in_mem_saturated_bins = 0;
  for (const chromap::peaks::ChromFragments& ch : chs) {
    std::vector<uint16_t> obs;
    int32_t max_c = 0;
    if (!chromap::peaks::BuildBinnedCutSignal(ch, bin_width, ext_size, &obs, &max_c)) {
      continue;
    }
    for (uint16_t v : obs) {
      ++in_mem_total_bins;
      if (v == 65535u) ++in_mem_saturated_bins;
    }
    std::vector<chromap::peaks::Peak> out;
    chromap::peaks::CallPeaksOnBins(ch.name, bin_width, obs, local_window, fdr, pval,
                                merge_gap, &out);
    for (const auto& p : out) {
      all.push_back(p);
    }
  }
  if (min_peak_bp > 0 || min_summit_cuts > 0) {
    std::vector<chromap::peaks::Peak> keep;
    keep.reserve(all.size());
    for (const auto& p : all) {
      if (p.end <= p.start) {
        continue;
      }
      if (min_peak_bp > 0 && (p.end - p.start) < min_peak_bp) {
        continue;
      }
      if (min_summit_cuts > 0 && p.max_signal < min_summit_cuts) {
        continue;
      }
      keep.push_back(p);
    }
    all.swap(keep);
  }
  const std::string np = out_prefix + ".narrowPeak";
  const std::string sm = out_prefix + ".summary.tsv";
  if (!chromap::peaks::WriteNarrowPeak(np, all)) {
    std::cerr << "Failed to write " << np << "\n";
    return 1;
  }
  std::ostringstream pstr;
  pstr.setf(std::ios::fixed);
  pstr << "bin_width=" << bin_width << " ext_size=" << ext_size
       << " local_window_bins=" << local_window << " fdr=" << std::setprecision(4) << fdr
       << " p_value=" << pval << " merge_bin_gap=" << merge_gap
       << " min_peak_bp=" << min_peak_bp << " min_summit_cuts=" << min_summit_cuts
       << " input=" << input;
  if (!chromap::peaks::WriteRunSummaryTsv(sm, all, JoinArgv(argc, argv), pstr.str(),
                                          in_mem_total_bins, in_mem_saturated_bins)) {
    std::cerr << "Failed to write " << sm << "\n";
    return 1;
  }
  return 0;
}
