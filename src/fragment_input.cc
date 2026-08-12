#include "rapidmacs/fragment_input.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <vector>

#include <htslib/sam.h>
#include <zlib.h>

namespace chromap {
namespace peaks {

namespace {

void TrimInPlace(std::string* s) {
  while (!s->empty() && (s->back() == '\r' || s->back() == '\n')) {
    s->pop_back();
  }
  size_t a = 0;
  while (a < s->size() && std::isspace(static_cast<unsigned char>((*s)[a]))) {
    ++a;
  }
  if (a > 0) {
    s->erase(0, a);
  }
  while (!s->empty() && std::isspace(static_cast<unsigned char>(s->back()))) {
    s->pop_back();
  }
}

bool ParseInt32(const char* p, int32_t* out) {
  char* endp = nullptr;
  long v = std::strtol(p, &endp, 10);
  if (endp == p || (endp && *endp != '\0')) {
    return false;
  }
  if (v < std::numeric_limits<int32_t>::min() ||
      v > std::numeric_limits<int32_t>::max()) {
    return false;
  }
  *out = static_cast<int32_t>(v);
  return true;
}

int64_t FilterCoordinateDuplicates(std::vector<Fragment>* fragments,
                                   int max_duplicates) {
  if (fragments == nullptr || max_duplicates < 0 || fragments->size() <= 1) {
    return 0;
  }
  const size_t before = fragments->size();
  size_t read = 0;
  size_t write = 0;
  while (read < before) {
    size_t run_end = read + 1;
    while (run_end < before &&
           (*fragments)[run_end].start == (*fragments)[read].start &&
           (*fragments)[run_end].end == (*fragments)[read].end) {
      ++run_end;
    }
    const size_t keep = std::min(
        run_end - read, static_cast<size_t>(max_duplicates));
    for (size_t i = 0; i < keep; ++i) {
      if (write != read + i) {
        (*fragments)[write] = (*fragments)[read + i];
      }
      ++write;
    }
    read = run_end;
  }
  fragments->resize(write);
  return static_cast<int64_t>(before - write);
}

void FinalizeFragmentMap(
    std::map<std::string, ChromFragments>* grouped, int max_duplicates,
    std::vector<ChromFragments>* by_chrom, FragmentInputStats* stats) {
  by_chrom->clear();
  by_chrom->reserve(grouped->size());
  for (auto& kv : *grouped) {
    ChromFragments& ch = kv.second;
    std::sort(ch.frags.begin(), ch.frags.end(),
              [](const Fragment& a, const Fragment& b) {
                if (a.start != b.start) return a.start < b.start;
                return a.end < b.end;
              });
    if (stats) {
      stats->duplicates_removed +=
          FilterCoordinateDuplicates(&ch.frags, max_duplicates);
    } else {
      FilterCoordinateDuplicates(&ch.frags, max_duplicates);
    }
    ch.max_end = 0;
    for (const Fragment& f : ch.frags) {
      ch.max_end = std::max(ch.max_end, f.end);
    }
    if (!ch.frags.empty()) by_chrom->push_back(std::move(ch));
  }
  if (stats) {
    stats->fragments_retained = 0;
    for (const ChromFragments& ch : *by_chrom) {
      stats->fragments_retained += static_cast<int64_t>(ch.frags.size());
    }
  }
}

void FinalizeBamFragmentVector(
    std::vector<ChromFragments>* grouped, int max_duplicates,
    std::vector<ChromFragments>* by_chrom, FragmentInputStats* stats,
    int num_threads) {
  const size_t n = grouped->size();
  std::vector<int64_t> duplicate_counts(n, 0);
  const int threads = std::max(
      1, std::min(num_threads, static_cast<int>(std::max<size_t>(1, n))));
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) if(threads > 1)
  for (int64_t ii = 0; ii < static_cast<int64_t>(n); ++ii) {
    ChromFragments& ch = grouped->at(static_cast<size_t>(ii));
    std::sort(ch.frags.begin(), ch.frags.end(),
              [](const Fragment& a, const Fragment& b) {
                if (a.start != b.start) return a.start < b.start;
                return a.end < b.end;
              });
    duplicate_counts[static_cast<size_t>(ii)] =
        FilterCoordinateDuplicates(&ch.frags, max_duplicates);
    ch.max_end = 0;
    for (const Fragment& f : ch.frags) {
      ch.max_end = std::max(ch.max_end, f.end);
    }
  }

  by_chrom->clear();
  by_chrom->reserve(n);
  for (size_t i = 0; i < n; ++i) {
    if (!grouped->at(i).frags.empty()) {
      by_chrom->push_back(std::move(grouped->at(i)));
    }
  }
  std::sort(by_chrom->begin(), by_chrom->end(),
            [](const ChromFragments& a, const ChromFragments& b) {
              return a.name < b.name;
            });
  if (stats) {
    stats->duplicates_removed = 0;
    stats->fragments_retained = 0;
    for (int64_t d : duplicate_counts) stats->duplicates_removed += d;
    for (const ChromFragments& ch : *by_chrom) {
      stats->fragments_retained += static_cast<int64_t>(ch.frags.size());
    }
  }
}

void AddFragment(std::map<std::string, ChromFragments>* grouped,
                 const std::string& chrom, int32_t start, int32_t end,
                 int32_t reference_length = 0) {
  ChromFragments& ch = (*grouped)[chrom];
  if (ch.name.empty()) ch.name = chrom;
  if (reference_length > 0) ch.reference_length = reference_length;
  Fragment f;
  f.start = start;
  f.end = end;
  f.count = 1;
  ch.frags.push_back(f);
}

template <typename ProcessLine>
bool ReadTextLines(const std::string& path, ProcessLine process_line) {
  const bool use_gz =
      path.size() >= 3 && path.compare(path.size() - 3, 3, ".gz") == 0;
  if (use_gz) {
    gzFile gz = gzopen(path.c_str(), "rb");
    if (gz == nullptr) return false;
    char buf[1 << 16];
    std::string pending;
    for (;;) {
      char* got = gzgets(gz, buf, sizeof(buf));
      if (got == nullptr) break;
      pending.append(got);
      if (!pending.empty() && pending.back() == '\n') {
        process_line(pending);
        pending.clear();
      }
    }
    if (!pending.empty()) process_line(pending);
    const int errnum = gzclose(gz);
    return errnum == Z_OK;
  }
  std::ifstream in(path.c_str());
  if (!in) return false;
  std::string line;
  while (std::getline(in, line)) process_line(line);
  return !in.bad();
}

}  // namespace

bool LoadFragmentsFromTsv(const std::string& path,
                          std::vector<ChromFragments>* by_chrom) {
  by_chrom->clear();
  std::unordered_map<std::string, size_t> index_of;
  const bool use_gz = path.size() >= 3 && path.compare(path.size() - 3, 3, ".gz") == 0;
  gzFile gz = nullptr;
  std::ifstream fin;
  if (use_gz) {
    gz = gzopen(path.c_str(), "rb");
    if (gz == nullptr) {
      return false;
    }
  } else {
    fin.open(path);
    if (!fin) {
      return false;
    }
  }
  auto process_line = [&](const std::string& raw) {
    std::string s = raw;
    TrimInPlace(&s);
    if (s.empty() || s[0] == '#') {
      return;
    }
    std::vector<std::string> fields;
    std::istringstream is(s);
    std::string col;
    while (std::getline(is, col, '\t')) {
      fields.push_back(std::move(col));
    }
    if (fields.size() < 5) {
      return;
    }
    int32_t start = 0, end = 0, dups = 1;
    size_t count_column = 4;
    if (fields.size() == 7) {
      int32_t mapq = 0;
      if (fields[3] != "N" ||
          !ParseInt32(fields[4].c_str(), &mapq) || mapq < 0 || mapq > 60 ||
          (fields[5] != "+" && fields[5] != "-")) {
        return;
      }
      count_column = 6;
    }
    if (!ParseInt32(fields[1].c_str(), &start) ||
        !ParseInt32(fields[2].c_str(), &end) ||
        !ParseInt32(fields[count_column].c_str(), &dups)) {
      return;
    }
    if (end <= start || dups < 0) {
      return;
    }
    if (dups == 0) {
      return;
    }
    const std::string& chrom = fields[0];
    size_t c = 0;
    auto it = index_of.find(chrom);
    if (it == index_of.end()) {
      c = by_chrom->size();
      index_of.emplace(chrom, c);
      ChromFragments block;
      block.name = chrom;
      by_chrom->push_back(std::move(block));
    } else {
      c = it->second;
    }
    Fragment f;
    f.start = start;
    f.end = end;
    f.count = dups;
    ChromFragments& ch = (*by_chrom)[c];
    ch.frags.push_back(f);
    if (end > ch.max_end) {
      ch.max_end = end;
    }
  };
  if (use_gz) {
    char buf[1 << 16];
    for (;;) {
      if (gzgets(gz, buf, sizeof(buf)) == nullptr) {
        break;
      }
      process_line(std::string(buf));
    }
    if (gzclose(gz) != Z_OK) {
      return false;
    }
  } else {
    std::string line;
    while (std::getline(fin, line)) {
      process_line(line);
    }
  }
  if (by_chrom->empty()) {
    return false;
  }
  for (auto& ch : *by_chrom) {
    std::sort(
        ch.frags.begin(), ch.frags.end(),
        [](const Fragment& a, const Fragment& b) {
          if (a.start != b.start) {
            return a.start < b.start;
          }
          return a.end < b.end;
        });
  }
  std::sort(by_chrom->begin(), by_chrom->end(),
            [](const ChromFragments& a, const ChromFragments& b) {
              return a.name < b.name;
            });
  return true;
}

bool LoadFragmentsFromBampe(const std::string& path,
                            std::vector<ChromFragments>* by_chrom,
                            FragmentInputStats* stats, int num_threads,
                            int max_duplicates) {
  if (by_chrom == nullptr || max_duplicates == 0 || max_duplicates < -1) {
    return false;
  }
  by_chrom->clear();
  FragmentInputStats local_stats;
  FragmentInputStats* st = stats ? stats : &local_stats;
  *st = FragmentInputStats();

  samFile* fp = sam_open(path.c_str(), "r");
  if (fp == nullptr) return false;
  if (num_threads > 1 && hts_set_threads(fp, num_threads - 1) < 0) {
    sam_close(fp);
    return false;
  }
  bam_hdr_t* hdr = sam_hdr_read(fp);
  if (hdr == nullptr) {
    sam_close(fp);
    return false;
  }
  bam1_t* rec = bam_init1();
  if (rec == nullptr) {
    bam_hdr_destroy(hdr);
    sam_close(fp);
    return false;
  }

  std::vector<ChromFragments> grouped(static_cast<size_t>(hdr->n_targets));
  for (int32_t tid = 0; tid < hdr->n_targets; ++tid) {
    ChromFragments& ch = grouped[static_cast<size_t>(tid)];
    if (hdr->target_name[tid] != nullptr) ch.name = hdr->target_name[tid];
    if (hdr->target_len[tid] <=
        static_cast<uint32_t>(std::numeric_limits<int32_t>::max())) {
      ch.reference_length = static_cast<int32_t>(hdr->target_len[tid]);
    }
  }
  int read_status = 0;
  while ((read_status = sam_read1(fp, hdr, rec)) >= 0) {
    ++st->records_read;
    const uint16_t flag = rec->core.flag;
    // Mirrors MACS3 3.0.3 bampe_pe_binary_parse: filter unmapped,
    // secondary, QC-fail, supplementary, mate-unmapped, read 2, and
    // non-proper pairs. The duplicate flag and MAPQ are intentionally kept.
    const uint16_t rejected = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL |
                              BAM_FSUPPLEMENTARY;
    if ((flag & rejected) != 0 || (flag & BAM_FPAIRED) == 0 ||
        (flag & BAM_FPROPER_PAIR) == 0 || (flag & BAM_FMUNMAP) != 0 ||
        (flag & BAM_FREAD2) != 0) {
      ++st->records_skipped;
      continue;
    }
    const int32_t tid = rec->core.tid;
    if (tid < 0 || tid >= hdr->n_targets || rec->core.pos < 0 ||
        rec->core.mpos < 0) {
      ++st->records_skipped;
      continue;
    }
    const int64_t tlen64 = rec->core.isize < 0
                               ? -static_cast<int64_t>(rec->core.isize)
                               : static_cast<int64_t>(rec->core.isize);
    const int64_t start64 =
        std::min<int64_t>(rec->core.pos, rec->core.mpos);
    const int64_t end64 = start64 + tlen64;
    if (tlen64 <= 0 || start64 < 0 ||
        end64 > std::numeric_limits<int32_t>::max()) {
      ++st->records_skipped;
      continue;
    }
    const int32_t start = static_cast<int32_t>(start64);
    const int32_t end = static_cast<int32_t>(end64);
    const char* target_name = hdr->target_name[tid];
    if (target_name == nullptr || target_name[0] == '\0') {
      ++st->records_skipped;
      continue;
    }
    Fragment f;
    f.start = start;
    f.end = end;
    f.count = 1;
    grouped[static_cast<size_t>(tid)].frags.push_back(f);
    ++st->fragments_parsed;
    st->fragment_bp_parsed += tlen64;
  }
  bam_destroy1(rec);
  bam_hdr_destroy(hdr);
  const int close_status = sam_close(fp);
  if (read_status < -1 || close_status != 0) {
    by_chrom->clear();
    return false;
  }
  FinalizeBamFragmentVector(&grouped, max_duplicates, by_chrom, st,
                            num_threads);
  std::fprintf(stderr,
               "[RapidMACS] BAMPE: %lld fragments parsed, %lld retained, "
               "%lld coordinate duplicates removed, %lld records skipped\n",
               static_cast<long long>(st->fragments_parsed),
               static_cast<long long>(st->fragments_retained),
               static_cast<long long>(st->duplicates_removed),
               static_cast<long long>(st->records_skipped));
  return !by_chrom->empty();
}

bool LoadFragmentsFromBedpe(const std::string& path,
                            std::vector<ChromFragments>* by_chrom,
                            FragmentInputStats* stats, int max_duplicates) {
  if (by_chrom == nullptr || max_duplicates == 0 || max_duplicates < -1) {
    return false;
  }
  by_chrom->clear();
  FragmentInputStats local_stats;
  FragmentInputStats* st = stats ? stats : &local_stats;
  *st = FragmentInputStats();
  std::map<std::string, ChromFragments> grouped;
  const bool ok = ReadTextLines(path, [&](const std::string& raw) {
    ++st->records_read;
    std::string s = raw;
    TrimInPlace(&s);
    if (s.empty() || s[0] == '#' || s.compare(0, 5, "track") == 0 ||
        s.compare(0, 7, "browser") == 0) {
      return;
    }
    std::vector<std::string> fields;
    std::istringstream is(s);
    std::string col;
    while (std::getline(is, col, '\t')) fields.push_back(col);
    std::string chrom;
    int32_t start = 0;
    int32_t end = 0;
    bool parsed = false;
    if (fields.size() >= 6) {
      int32_t s1 = 0, e1 = 0, s2 = 0, e2 = 0;
      if (fields[0] == fields[3] &&
          ParseInt32(fields[1].c_str(), &s1) &&
          ParseInt32(fields[2].c_str(), &e1) &&
          ParseInt32(fields[4].c_str(), &s2) &&
          ParseInt32(fields[5].c_str(), &e2)) {
        chrom = fields[0];
        start = std::min(s1, s2);
        end = std::max(e1, e2);
        parsed = true;
      }
    } else if (fields.size() >= 3 &&
               ParseInt32(fields[1].c_str(), &start) &&
               ParseInt32(fields[2].c_str(), &end)) {
      chrom = fields[0];
      parsed = true;
    }
    if (!parsed || chrom.empty() || start < 0 || end <= start) {
      ++st->records_skipped;
      return;
    }
    AddFragment(&grouped, chrom, start, end);
    ++st->fragments_parsed;
    st->fragment_bp_parsed += static_cast<int64_t>(end) - start;
  });
  if (!ok) return false;
  FinalizeFragmentMap(&grouped, max_duplicates, by_chrom, st);
  std::fprintf(stderr,
               "[RapidMACS] BEDPE: %lld fragments parsed, %lld retained, "
               "%lld coordinate duplicates removed, %lld records skipped\n",
               static_cast<long long>(st->fragments_parsed),
               static_cast<long long>(st->fragments_retained),
               static_cast<long long>(st->duplicates_removed),
               static_cast<long long>(st->records_skipped));
  return !by_chrom->empty();
}

}  // namespace peaks
}  // namespace chromap
