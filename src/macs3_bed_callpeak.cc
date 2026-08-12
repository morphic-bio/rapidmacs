#include "rapidmacs/macs3_bed_callpeak.h"

#include <zlib.h>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "rapidmacs/macs3_bdg_track.h"
#include "rapidmacs/peak_io.h"

namespace chromap {
namespace peaks {
namespace {

struct BedChromTags {
  std::string name;
  std::vector<int32_t> plus;
  std::vector<int32_t> minus;
};

struct BedTagTrack {
  std::map<std::string, BedChromTags> chroms;
  int64_t total_tags = 0;
  int32_t detected_tag_size = 0;
};

std::vector<std::string> SplitTab(const std::string& line) {
  std::vector<std::string> out;
  size_t b = 0;
  while (true) {
    const size_t e = line.find('\t', b);
    if (e == std::string::npos) {
      out.push_back(line.substr(b));
      break;
    }
    out.push_back(line.substr(b, e - b));
    b = e + 1;
  }
  return out;
}

std::string TrimCr(std::string s) {
  if (!s.empty() && s.back() == '\r') {
    s.pop_back();
  }
  return s;
}

bool StartsWith(const std::string& s, const char* prefix) {
  const size_t n = std::strlen(prefix);
  return s.size() >= n && s.compare(0, n, prefix) == 0;
}

bool ParseInt32(const std::string& s, int32_t* out) {
  if (out == nullptr || s.empty()) {
    return false;
  }
  char* end = nullptr;
  errno = 0;
  const long v = std::strtol(s.c_str(), &end, 10);
  if (errno != 0 || end == s.c_str() || *end != '\0' ||
      v < std::numeric_limits<int32_t>::min() ||
      v > std::numeric_limits<int32_t>::max()) {
    return false;
  }
  *out = static_cast<int32_t>(v);
  return true;
}

float QuantizeBedGraphFloat5(float x) {
  char buf[64];
  const int n = std::snprintf(buf, sizeof(buf), "%.5f", static_cast<double>(x));
  if (n <= 0 || n >= static_cast<int>(sizeof(buf))) {
    return x;
  }
  char* end = nullptr;
  const float y = std::strtof(buf, &end);
  return (end != buf && *end == '\0') ? y : x;
}

bool IsIgnorableBedLine(const std::string& line) {
  return line.empty() || line[0] == '#' || StartsWith(line, "track") ||
         StartsWith(line, "browser");
}

bool LoadBedTags(const std::string& path, BedTagTrack* out, std::string* error) {
  if (out == nullptr) {
    if (error) *error = "BED tag output pointer is null";
    return false;
  }
  *out = BedTagTrack();
  gzFile gz = gzopen(path.c_str(), "rb");
  if (gz == nullptr) {
    if (error) *error = "failed to open BED input: " + path;
    return false;
  }
  static const int kBufLen = 1 << 20;
  std::vector<char> buf(kBufLen);
  int64_t line_no = 0;
  int64_t tag_len_sum = 0;
  int64_t tag_len_n = 0;
  while (gzgets(gz, buf.data(), kBufLen) != nullptr) {
    ++line_no;
    std::string line(buf.data());
    if (!line.empty() && line.back() == '\n') line.pop_back();
    line = TrimCr(line);
    if (IsIgnorableBedLine(line)) {
      continue;
    }
    std::vector<std::string> f = SplitTab(line);
    if (f.size() < 3) {
      if (error) {
        *error = path + ":" + std::to_string(line_no) +
                 ": BED input requires at least 3 columns";
      }
      gzclose(gz);
      return false;
    }
    int32_t start = 0;
    int32_t end = 0;
    if (!ParseInt32(f[1], &start) || !ParseInt32(f[2], &end)) {
      if (error) {
        *error = path + ":" + std::to_string(line_no) + ": invalid BED start/end";
      }
      gzclose(gz);
      return false;
    }
    if (end <= start) {
      if (error) {
        *error = path + ":" + std::to_string(line_no) + ": BED end <= start";
      }
      gzclose(gz);
      return false;
    }
    int strand = 0;
    int32_t pos = start;
    if (f.size() >= 6) {
      if (f[5] == "+") {
        strand = 0;
        pos = start;
      } else if (f[5] == "-") {
        strand = 1;
        pos = end;
      } else {
        if (error) {
          *error = path + ":" + std::to_string(line_no) +
                   ": BED strand column must be + or - when present";
        }
        gzclose(gz);
        return false;
      }
    }
    if (pos < 0) {
      continue;
    }
    if (tag_len_n < 10) {
      tag_len_sum += static_cast<int64_t>(end) - static_cast<int64_t>(start);
      ++tag_len_n;
    }
    BedChromTags& ch = out->chroms[f[0]];
    if (ch.name.empty()) {
      ch.name = f[0];
    }
    if (strand == 0) {
      ch.plus.push_back(pos);
    } else {
      ch.minus.push_back(pos);
    }
    ++out->total_tags;
  }
  int zerr = Z_OK;
  if (!gzeof(gz)) {
    (void)gzerror(gz, &zerr);
    if (error) *error = "error reading BED input: " + path;
    gzclose(gz);
    return false;
  }
  if (gzclose(gz) != Z_OK) {
    if (error) *error = "error closing BED input: " + path;
    return false;
  }
  if (tag_len_n > 0) {
    out->detected_tag_size = static_cast<int32_t>(tag_len_sum / tag_len_n);
  }
  for (auto& kv : out->chroms) {
    std::sort(kv.second.plus.begin(), kv.second.plus.end());
    std::sort(kv.second.minus.begin(), kv.second.minus.end());
  }
  return true;
}

bool ResolveKeepDup(const std::string& keep_dup, int* keep, std::string* error) {
  if (keep_dup == "all") {
    *keep = -1;
    return true;
  }
  if (keep_dup == "auto") {
    if (error) {
      *error = "--keep-dup auto is not implemented in the BED compatibility path; "
               "use an integer or all";
    }
    return false;
  }
  char* end = nullptr;
  errno = 0;
  const long v = std::strtol(keep_dup.c_str(), &end, 10);
  if (errno != 0 || end == keep_dup.c_str() || *end != '\0' || v < 1 ||
      v > std::numeric_limits<int>::max()) {
    if (error) *error = "invalid --keep-dup value: " + keep_dup;
    return false;
  }
  *keep = static_cast<int>(v);
  return true;
}

int64_t FilterDupSorted(std::vector<int32_t>* v, int keep) {
  if (v == nullptr) {
    return 0;
  }
  if (keep < 0 || v->size() <= 1) {
    return static_cast<int64_t>(v->size());
  }
  size_t w = 0;
  size_t i = 0;
  while (i < v->size()) {
    size_t j = i + 1;
    while (j < v->size() && (*v)[j] == (*v)[i]) {
      ++j;
    }
    const size_t n = std::min(static_cast<size_t>(keep), j - i);
    for (size_t k = 0; k < n; ++k) {
      (*v)[w++] = (*v)[i];
    }
    i = j;
  }
  v->resize(w);
  return static_cast<int64_t>(w);
}

int64_t ApplyKeepDup(BedTagTrack* track, int keep) {
  int64_t total = 0;
  for (auto& kv : track->chroms) {
    total += FilterDupSorted(&kv.second.plus, keep);
    total += FilterDupSorted(&kv.second.minus, keep);
  }
  return total;
}

int32_t ClampCoord(int64_t x, int32_t rlength) {
  if (x < 0) {
    return 0;
  }
  if (x > static_cast<int64_t>(rlength)) {
    return rlength;
  }
  return static_cast<int32_t>(x);
}

Macs3BdgTrack BuildSingleEndPileupTrack(const std::vector<int32_t>& plus,
                                        const std::vector<int32_t>& minus,
                                        int64_t five_shift,
                                        int64_t three_shift,
                                        int32_t rlength,
                                        float scale_factor,
                                        float baseline_value) {
  Macs3BdgTrack out;
  const size_t n = plus.size() + minus.size();
  if (n == 0) {
    return out;
  }
  std::vector<int32_t> starts;
  std::vector<int32_t> ends;
  starts.reserve(n);
  ends.reserve(n);
  for (int32_t p : plus) {
    starts.push_back(ClampCoord(static_cast<int64_t>(p) - five_shift, rlength));
    ends.push_back(ClampCoord(static_cast<int64_t>(p) + three_shift, rlength));
  }
  for (int32_t p : minus) {
    starts.push_back(ClampCoord(static_cast<int64_t>(p) - three_shift, rlength));
    ends.push_back(ClampCoord(static_cast<int64_t>(p) + five_shift, rlength));
  }
  std::sort(starts.begin(), starts.end());
  std::sort(ends.begin(), ends.end());
  size_t is = 0;
  size_t ie = 0;
  int32_t pileup = 0;
  int32_t pre_p = std::min(starts[0], ends[0]);
  if (pre_p != 0) {
    out.p.push_back(pre_p);
    out.v.push_back(std::max(0.0f, baseline_value));
  }
  auto emit = [&](int32_t p) {
    if (p != pre_p) {
      out.p.push_back(p);
      out.v.push_back(std::max(static_cast<float>(pileup) * scale_factor,
                               baseline_value));
      pre_p = p;
    }
  };
  while (is < starts.size() && ie < ends.size()) {
    if (starts[is] < ends[ie]) {
      emit(starts[is]);
      ++pileup;
      ++is;
    } else if (starts[is] > ends[ie]) {
      emit(ends[ie]);
      --pileup;
      ++ie;
    } else {
      ++is;
      ++ie;
    }
  }
  while (ie < ends.size()) {
    emit(ends[ie]);
    --pileup;
    ++ie;
  }
  return out;
}

bool BuildPpoisTrackMinMerge(const Macs3BdgTrack& treat,
                             const Macs3BdgTrack& lambda,
                             double pseudocount, Macs3BdgTrack* ppois,
                             bool quantize_inputs_to_bedgraph = false) {
  if (ppois == nullptr || treat.p.size() != treat.v.size() ||
      lambda.p.size() != lambda.v.size()) {
    return false;
  }
  ppois->p.clear();
  ppois->v.clear();
  if (treat.p.empty() || lambda.p.empty()) {
    return true;
  }
  size_t i1 = 0;
  size_t i2 = 0;
  int32_t prev = 0;
  bool have = false;
  float prev_v = 0.0f;
  while (i1 < treat.p.size() && i2 < lambda.p.size()) {
    const int32_t p1 = treat.p[i1];
    const int32_t p2 = lambda.p[i2];
    const int32_t nx = std::min(p1, p2);
    if (nx > prev) {
      float treat_v = treat.v[i1];
      float lambda_v = lambda.v[i2];
      if (quantize_inputs_to_bedgraph) {
        treat_v = QuantizeBedGraphFloat5(treat_v);
        lambda_v = QuantizeBedGraphFloat5(lambda_v);
      }
      const float v = static_cast<float>(Macs3BdgcmpPpoisNegLog10P(
          treat_v, lambda_v, pseudocount));
      if (have && std::fabs(static_cast<double>(v) - prev_v) <= 1e-6) {
        ppois->p.back() = nx;
      } else {
        ppois->p.push_back(nx);
        ppois->v.push_back(v);
        prev_v = v;
        have = true;
      }
    }
    if (p1 <= p2) {
      ++i1;
    }
    if (p2 <= p1) {
      ++i2;
    }
    prev = nx;
  }
  return true;
}

void BuildMacs3PairedPositions(const Macs3BdgTrack& treat,
                               const Macs3BdgTrack& lambda,
                               std::vector<int32_t>* positions) {
  positions->clear();
  if (treat.p.empty() || lambda.p.empty()) {
    return;
  }
  positions->reserve(treat.p.size() + lambda.p.size());
  size_t it = 0;
  size_t ic = 0;
  while (it < treat.p.size() && ic < lambda.p.size()) {
    if (treat.p[it] < lambda.p[ic]) {
      positions->push_back(treat.p[it]);
      ++it;
    } else if (treat.p[it] > lambda.p[ic]) {
      positions->push_back(lambda.p[ic]);
      ++ic;
    } else {
      positions->push_back(treat.p[it]);
      ++it;
      ++ic;
    }
  }
}

bool WriteTrackBedGraph(const std::string& path,
                        const std::vector<std::string>& chrom_names,
                        const std::vector<Macs3BdgTrack>& tracks) {
  if (path.empty()) {
    return true;
  }
  FILE* fp = std::fopen(path.c_str(), "w");
  if (fp == nullptr) {
    return false;
  }
  for (size_t ci = 0; ci < chrom_names.size() && ci < tracks.size(); ++ci) {
    const Macs3BdgTrack& tr = tracks[ci];
    if (tr.p.empty() || tr.v.empty()) {
      continue;
    }
    int32_t pre = 0;
    float pre_v = tr.v[0];
    for (size_t i = 1; i < tr.p.size() && i < tr.v.size(); ++i) {
      const int32_t p = tr.p[i - 1];
      const float v = tr.v[i];
      if (std::fabs(pre_v - v) > 1e-5f) {
        if (std::fprintf(fp, "%s\t%d\t%d\t%.5f\n", chrom_names[ci].c_str(),
                         pre, p, static_cast<double>(pre_v)) < 0) {
          std::fclose(fp);
          return false;
        }
        pre_v = v;
        pre = p;
      }
    }
    if (std::fprintf(fp, "%s\t%d\t%d\t%.5f\n", chrom_names[ci].c_str(),
                     pre, tr.p.back(), static_cast<double>(pre_v)) < 0) {
      std::fclose(fp);
      return false;
    }
  }
  return std::fclose(fp) == 0;
}

bool WritePairedTreatLambdaBedGraphs(
    const std::string& treat_path, const std::string& lambda_path,
    const std::vector<std::string>& chrom_names,
    const std::vector<Macs3BdgTrack>& treat_tracks,
    const std::vector<Macs3BdgTrack>& lambda_tracks) {
  if (treat_path.empty() && lambda_path.empty()) {
    return true;
  }
  FILE* fpt = nullptr;
  FILE* fpc = nullptr;
  if (!treat_path.empty()) {
    fpt = std::fopen(treat_path.c_str(), "w");
    if (fpt == nullptr) {
      return false;
    }
  }
  if (!lambda_path.empty()) {
    fpc = std::fopen(lambda_path.c_str(), "w");
    if (fpc == nullptr) {
      if (fpt != nullptr) std::fclose(fpt);
      return false;
    }
  }
  for (size_t ci = 0; ci < chrom_names.size() && ci < treat_tracks.size() &&
                      ci < lambda_tracks.size(); ++ci) {
    const Macs3BdgTrack& treat = treat_tracks[ci];
    const Macs3BdgTrack& lambda = lambda_tracks[ci];
    if (treat.p.empty() || lambda.p.empty()) {
      continue;
    }
    size_t it = 0;
    size_t ic = 0;
    bool first = true;
    int32_t pre_p_t = 0;
    int32_t pre_p_c = 0;
    int32_t pending_p = 0;
    float pre_v_t = 0.0f;
    float pre_v_c = 0.0f;
    while (it < treat.p.size() && ic < lambda.p.size()) {
      const int32_t p1 = treat.p[it];
      const int32_t p2 = lambda.p[ic];
      const int32_t p = std::min(p1, p2);
      const float v_t = treat.v[it];
      const float v_c = lambda.v[ic];
      if (first) {
        pre_v_t = v_t;
        pre_v_c = v_c;
        pending_p = p;
        first = false;
      } else {
        if (fpt != nullptr && std::fabs(pre_v_t - v_t) > 1e-5f) {
          if (std::fprintf(fpt, "%s\t%d\t%d\t%.5f\n", chrom_names[ci].c_str(),
                           pre_p_t, pending_p, static_cast<double>(pre_v_t)) < 0) {
            if (fpt != nullptr) std::fclose(fpt);
            if (fpc != nullptr) std::fclose(fpc);
            return false;
          }
          pre_v_t = v_t;
          pre_p_t = pending_p;
        }
        if (fpc != nullptr && std::fabs(pre_v_c - v_c) > 1e-5f) {
          if (std::fprintf(fpc, "%s\t%d\t%d\t%.5f\n", chrom_names[ci].c_str(),
                           pre_p_c, pending_p, static_cast<double>(pre_v_c)) < 0) {
            if (fpt != nullptr) std::fclose(fpt);
            if (fpc != nullptr) std::fclose(fpc);
            return false;
          }
          pre_v_c = v_c;
          pre_p_c = pending_p;
        }
        pending_p = p;
      }
      if (p1 <= p2) {
        ++it;
      }
      if (p2 <= p1) {
        ++ic;
      }
    }
    if (!first) {
      if (fpt != nullptr &&
          std::fprintf(fpt, "%s\t%d\t%d\t%.5f\n", chrom_names[ci].c_str(),
                       pre_p_t, pending_p, static_cast<double>(pre_v_t)) < 0) {
        if (fpt != nullptr) std::fclose(fpt);
        if (fpc != nullptr) std::fclose(fpc);
        return false;
      }
      if (fpc != nullptr &&
          std::fprintf(fpc, "%s\t%d\t%d\t%.5f\n", chrom_names[ci].c_str(),
                       pre_p_c, pending_p, static_cast<double>(pre_v_c)) < 0) {
        if (fpt != nullptr) std::fclose(fpt);
        if (fpc != nullptr) std::fclose(fpc);
        return false;
      }
    }
  }
  const bool ok_t = (fpt == nullptr || std::fclose(fpt) == 0);
  const bool ok_c = (fpc == nullptr || std::fclose(fpc) == 0);
  return ok_t && ok_c;
}

int64_t CountLines(const std::string& path) {
  FILE* fp = std::fopen(path.c_str(), "r");
  if (fp == nullptr) {
    return 0;
  }
  int64_t n = 0;
  char buf[65536];
  while (std::fgets(buf, sizeof(buf), fp) != nullptr) {
    if (buf[0] != '\0' && buf[0] != '#') {
      ++n;
    }
  }
  std::fclose(fp);
  return n;
}

}  // namespace

const char* SupportedMacs3BedCallPeakProfiles() {
  return "signac-atac";
}

bool ResolveMacs3BedCallPeakProfile(const std::string& profile_name,
                                    Macs3BedCallPeakProfile* profile,
                                    std::string* error_message) {
  if (profile == nullptr) {
    if (error_message) *error_message = "BED callpeak profile output pointer is null";
    return false;
  }
  if (profile_name != "signac-atac") {
    if (error_message) {
      *error_message = "Unknown MACS BED callpeak profile \"" + profile_name +
                       "\" (supported: " + SupportedMacs3BedCallPeakProfiles() + ")";
    }
    return false;
  }

  *profile = Macs3BedCallPeakProfile();
  profile->name = "signac-atac";
  profile->input_format = "BED";
  profile->params.threshold_mode = Macs3FragThresholdMode::kQValue;
  profile->params.qvalue_cutoff = 0.05;
  profile->params.pvalue_cutoff = 0.0;
  profile->params.effective_genome_size = 2700000000LL;
  profile->params.nomodel = true;
  profile->params.extsize = 200;
  profile->params.shift = -100;
  profile->params.llocal_bp = 10000;
  profile->params.keep_dup = "1";
  profile->params.broad = false;
  profile->params.broad_cutoff = 0.1;
  profile->params.no_lambda = false;
  profile->params.call_summits = false;
  return true;
}

bool RunMacs3BedCallPeak(const std::string& bed_input,
                         const Macs3BedCallPeakParams& params,
                         const std::string& narrowpeak_out,
                         const std::string& summits_out,
                         Macs3BedCallPeakStats* stats,
                         std::string* error_message) {
  if (stats) {
    *stats = Macs3BedCallPeakStats();
  }
  if (bed_input.empty()) {
    if (error_message) *error_message = "--macs-bed-input is required";
    return false;
  }
  if (narrowpeak_out.empty() || summits_out.empty()) {
    if (error_message) *error_message = "narrowPeak and summits outputs are required";
    return false;
  }
  if (!params.nomodel) {
    if (error_message) *error_message = "BED callpeak compatibility currently requires --nomodel";
    return false;
  }
  if (params.broad) {
    if (error_message) *error_message = "--broad is parsed but not implemented in the BED compatibility path";
    return false;
  }
  if (params.extsize < 1 || params.effective_genome_size < 1 ||
      params.llocal_bp < 1) {
    if (error_message) {
      *error_message = "invalid --extsize, --effective-genome-size, or --llocal";
    }
    return false;
  }
  const bool q_mode = params.threshold_mode == Macs3FragThresholdMode::kQValue;
  if (q_mode) {
    if (std::isnan(params.qvalue_cutoff) || params.qvalue_cutoff <= 0.0 ||
        params.qvalue_cutoff > 1.0) {
      if (error_message) *error_message = "invalid --qvalue";
      return false;
    }
  } else {
    if (std::isnan(params.pvalue_cutoff) || params.pvalue_cutoff <= 0.0 ||
        params.pvalue_cutoff > 1.0) {
      if (error_message) *error_message = "invalid --pvalue";
      return false;
    }
  }
  int keep = 1;
  if (!ResolveKeepDup(params.keep_dup, &keep, error_message)) {
    return false;
  }

  BedTagTrack tags;
  if (!LoadBedTags(bed_input, &tags, error_message)) {
    return false;
  }
  const int64_t total_tags = tags.total_tags;
  const int64_t filtered_tags = ApplyKeepDup(&tags, keep);
  const int32_t min_length =
      params.min_length > 0 ? params.min_length : params.extsize;
  const int32_t max_gap =
      params.max_gap > 0
          ? params.max_gap
          : (tags.detected_tag_size > 0 ? tags.detected_tag_size : params.extsize);

  std::vector<std::string> chrom_names;
  std::vector<Macs3BdgTrack> treat_tracks;
  std::vector<Macs3BdgTrack> lambda_tracks;
  std::vector<Macs3BdgTrack> ppois_tracks;
  chrom_names.reserve(tags.chroms.size());
  treat_tracks.reserve(tags.chroms.size());
  lambda_tracks.reserve(tags.chroms.size());
  ppois_tracks.reserve(tags.chroms.size());

  const int32_t rlength = std::numeric_limits<int32_t>::max();
  const int64_t treat_five_shift = -static_cast<int64_t>(params.shift);
  const int64_t treat_three_shift =
      static_cast<int64_t>(params.shift) + static_cast<int64_t>(params.extsize);
  const int64_t lambda_five_shift = static_cast<int64_t>(params.llocal_bp / 2);
  const int64_t lambda_three_shift =
      static_cast<int64_t>(params.llocal_bp) - static_cast<int64_t>(params.llocal_bp / 2);
  const float lambda_bg = static_cast<float>(
      static_cast<double>(params.extsize) * static_cast<double>(filtered_tags) /
      static_cast<double>(params.effective_genome_size));
  const float lambda_scale = static_cast<float>(
      static_cast<double>(params.extsize) / static_cast<double>(params.llocal_bp));

  for (const auto& kv : tags.chroms) {
    const BedChromTags& ch = kv.second;
    chrom_names.push_back(ch.name);
    treat_tracks.push_back(BuildSingleEndPileupTrack(
        ch.plus, ch.minus, treat_five_shift, treat_three_shift, rlength, 1.0f, 0.0f));
    if (params.no_lambda) {
      Macs3BdgTrack global_lambda;
      if (!treat_tracks.back().p.empty()) {
        global_lambda.p.push_back(treat_tracks.back().p.back());
        global_lambda.v.push_back(lambda_bg);
      }
      lambda_tracks.push_back(std::move(global_lambda));
    } else {
      lambda_tracks.push_back(BuildSingleEndPileupTrack(
          ch.plus, ch.minus, lambda_five_shift, lambda_three_shift, rlength,
          lambda_scale, lambda_bg));
    }
    Macs3BdgTrack ppois;
    if (!BuildPpoisTrackMinMerge(treat_tracks.back(), lambda_tracks.back(), 0.0,
                                 &ppois)) {
      if (error_message) *error_message = "failed to build BED ppois track";
      return false;
    }
    ppois_tracks.push_back(std::move(ppois));
  }

  if (!WritePairedTreatLambdaBedGraphs(params.treat_bdg, params.lambda_bdg,
                                       chrom_names, treat_tracks, lambda_tracks)) {
    if (error_message) *error_message = "failed to write BED callpeak bedGraph output";
    return false;
  }
  if (!params.ppois_bdg.empty()) {
    std::vector<Macs3BdgTrack> ppois_bdg_tracks;
    ppois_bdg_tracks.reserve(chrom_names.size());
    for (size_t i = 0; i < chrom_names.size(); ++i) {
      Macs3BdgTrack ppois_bdg;
      if (!BuildPpoisTrackMinMerge(treat_tracks[i], lambda_tracks[i], 0.0,
                                   &ppois_bdg, true)) {
        if (error_message) *error_message = "failed to build BED ppois bedGraph track";
        return false;
      }
      ppois_bdg_tracks.push_back(std::move(ppois_bdg));
    }
    if (!WriteTrackBedGraph(params.ppois_bdg, chrom_names, ppois_bdg_tracks)) {
      if (error_message) *error_message = "failed to write BED ppois bedGraph output";
      return false;
    }
  }

  std::map<double, double> pqtable;
  BuildMacs3PqTable(ppois_tracks, &pqtable, true);
  std::vector<Macs3FragNarrowPeakRow> rows;
  const float threshold =
      q_mode ? static_cast<float>(-std::log10(params.qvalue_cutoff))
             : static_cast<float>(-std::log10(params.pvalue_cutoff));
  for (size_t i = 0; i < chrom_names.size(); ++i) {
    std::vector<int32_t> paired_positions;
    BuildMacs3PairedPositions(treat_tracks[i], lambda_tracks[i], &paired_positions);
    if (q_mode) {
      Macs3BdgTrack qscore;
      if (!BuildMacs3QScoreTrackFromPpois(ppois_tracks[i], pqtable, &qscore,
                                          true) ||
          !CollectMacs3FragNarrowPeakRowsFromThresholdTrack(
              chrom_names[i], qscore, ppois_tracks[i], treat_tracks[i],
              lambda_tracks[i], threshold, min_length, max_gap, 0.0, &rows,
              &paired_positions, false, params.call_summits)) {
        if (error_message) *error_message = "failed BED q-value peak collection";
        return false;
      }
    } else {
      if (!CollectMacs3FragNarrowPeakRowsFromTracks(
              chrom_names[i], ppois_tracks[i], treat_tracks[i], lambda_tracks[i],
              threshold, min_length, max_gap, 0.0, &rows, &paired_positions,
              false, nullptr, params.call_summits)) {
        if (error_message) *error_message = "failed BED p-value peak collection";
        return false;
      }
    }
  }
  if (!WriteMacs3FragNarrowPeakRowsWithScoreCap(
          &rows, narrowpeak_out, summits_out, &pqtable, params.name_prefix,
          q_mode, false, true)) {
    if (error_message) *error_message = "failed to write BED narrowPeak/summits";
    return false;
  }
  if (stats) {
    stats->total_tags = total_tags;
    stats->tags_after_filtering = filtered_tags;
    stats->detected_tag_size = tags.detected_tag_size;
    stats->resolved_min_length = min_length;
    stats->resolved_max_gap = max_gap;
    stats->keep_dup_max = keep;
    stats->peak_count = CountLines(narrowpeak_out);
  }
  return true;
}

}  // namespace peaks
}  // namespace chromap
