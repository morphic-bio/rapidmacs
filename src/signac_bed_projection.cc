#include "rapidmacs/signac_bed_projection.h"

#include <zlib.h>

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <sys/stat.h>
#include <sys/types.h>

namespace chromap {
namespace peaks {
namespace {

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

bool ParseInt64(const std::string& s, int64_t* out) {
  if (out == nullptr || s.empty()) {
    return false;
  }
  char* end = nullptr;
  errno = 0;
  const long long v = std::strtoll(s.c_str(), &end, 10);
  if (errno != 0 || end == s.c_str() || *end != '\0') {
    return false;
  }
  *out = static_cast<int64_t>(v);
  return true;
}

bool EnsureDir(const std::string& path) {
  if (path.empty()) {
    return false;
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

bool ReadCells(const std::string& path, std::unordered_set<std::string>* cells,
               std::string* error) {
  cells->clear();
  if (path.empty()) {
    return true;
  }
  std::ifstream in(path.c_str());
  if (!in) {
    if (error) *error = "failed to open --cells: " + path;
    return false;
  }
  std::string line;
  while (std::getline(in, line)) {
    line = TrimCr(line);
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::vector<std::string> f = SplitTab(line);
    if (f.empty() || f[0].empty()) {
      continue;
    }
    if (f[0] == "barcode" || f[0] == "cell" || f[0] == "cell_barcode") {
      continue;
    }
    cells->insert(f[0]);
  }
  return true;
}

int FindColumn(const std::vector<std::string>& header, const std::string& name) {
  for (size_t i = 0; i < header.size(); ++i) {
    if (header[i] == name) {
      return static_cast<int>(i);
    }
  }
  return -1;
}

int FindBarcodeColumn(const std::vector<std::string>& header) {
  static const char* kNames[] = {"barcode", "cell", "cell_barcode", "CB"};
  for (const char* name : kNames) {
    const int idx = FindColumn(header, name);
    if (idx >= 0) {
      return idx;
    }
  }
  return 0;
}

bool ReadGroups(const std::string& path, const std::string& group_column,
                std::unordered_map<std::string, std::string>* barcode_to_group,
                std::set<std::string>* groups, std::string* error) {
  barcode_to_group->clear();
  groups->clear();
  if (path.empty()) {
    return true;
  }
  std::ifstream in(path.c_str());
  if (!in) {
    if (error) *error = "failed to open --barcode-group-tsv: " + path;
    return false;
  }
  std::string line;
  int barcode_col = 0;
  int group_col = 1;
  bool header_consumed = group_column.empty();
  while (std::getline(in, line)) {
    line = TrimCr(line);
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::vector<std::string> f = SplitTab(line);
    if (f.size() < 2) {
      if (error) *error = "barcode-group TSV requires at least two columns";
      return false;
    }
    if (!header_consumed) {
      group_col = FindColumn(f, group_column);
      if (group_col < 0) {
        if (error) *error = "group column not found in --barcode-group-tsv: " + group_column;
        return false;
      }
      barcode_col = FindBarcodeColumn(f);
      header_consumed = true;
      continue;
    }
    if (group_column.empty() &&
        (f[0] == "barcode" || f[0] == "cell" || f[0] == "cell_barcode")) {
      continue;
    }
    const int need = std::max(barcode_col, group_col);
    if (static_cast<int>(f.size()) <= need) {
      if (error) *error = "barcode-group TSV row has fewer columns than header";
      return false;
    }
    if (f[barcode_col].empty() || f[group_col].empty()) {
      continue;
    }
    (*barcode_to_group)[f[barcode_col]] = f[group_col];
    groups->insert(f[group_col]);
  }
  return true;
}

std::string JoinFirstFields(const std::vector<std::string>& f, size_t n) {
  std::ostringstream out;
  const size_t lim = std::min(n, f.size());
  for (size_t i = 0; i < lim; ++i) {
    if (i) out << '\t';
    out << f[i];
  }
  out << '\n';
  return out.str();
}

bool WriteRepeated(FILE* fp, const std::string& row, int64_t count) {
  for (int64_t i = 0; i < count; ++i) {
    if (std::fwrite(row.data(), 1, row.size(), fp) != row.size()) {
      return false;
    }
  }
  return true;
}

}  // namespace

bool ProjectFragmentsToSignacBed(const SignacBedProjectionParams& params,
                                 SignacBedProjectionResult* result,
                                 std::string* error_message) {
  if (result == nullptr) {
    if (error_message) *error_message = "projection result pointer is null";
    return false;
  }
  *result = SignacBedProjectionResult();
  result->count_expansion_used = params.expand_counts;
  if (params.fragments_path.empty()) {
    if (error_message) *error_message = "--input-fragments is required for projection";
    return false;
  }

  std::unordered_set<std::string> cells;
  if (!ReadCells(params.cells_path, &cells, error_message)) {
    return false;
  }
  result->selected_cell_count = static_cast<int64_t>(cells.size());

  std::unordered_map<std::string, std::string> barcode_to_group;
  std::set<std::string> groups;
  if (!ReadGroups(params.barcode_group_tsv_path, params.group_column,
                  &barcode_to_group, &groups, error_message)) {
    return false;
  }

  const bool grouped = !params.barcode_group_tsv_path.empty();
  std::map<std::string, FILE*> group_files;
  FILE* single = nullptr;
  if (grouped) {
    std::string dir = params.bed_out_dir;
    if (dir.empty()) {
      if (error_message) *error_message = "--macs-bed-out-dir is required with --barcode-group-tsv";
      return false;
    }
    if (!EnsureDir(dir)) {
      if (error_message) *error_message = "failed to create --macs-bed-out-dir: " + dir;
      return false;
    }
    for (const std::string& group : groups) {
      const std::string path = dir + "/" + SanitizeForPath(group) + ".bed";
      FILE* fp = std::fopen(path.c_str(), "w");
      if (fp == nullptr) {
        if (error_message) *error_message = "failed to open projected BED: " + path;
        for (auto& kv : group_files) std::fclose(kv.second);
        return false;
      }
      group_files[group] = fp;
      result->group_bed_paths[group] = path;
      result->group_rows[group] = 0;
    }
  } else {
    if (params.bed_out.empty()) {
      if (error_message) *error_message = "--macs-bed-out is required for single projected BED output";
      return false;
    }
    single = std::fopen(params.bed_out.c_str(), "w");
    if (single == nullptr) {
      if (error_message) *error_message = "failed to open --macs-bed-out: " + params.bed_out;
      return false;
    }
    result->bed_path = params.bed_out;
  }

  gzFile gz = gzopen(params.fragments_path.c_str(), "rb");
  if (gz == nullptr) {
    if (error_message) *error_message = "failed to open fragments: " + params.fragments_path;
    if (single) std::fclose(single);
    for (auto& kv : group_files) std::fclose(kv.second);
    return false;
  }

  static const int kBufLen = 1 << 20;
  std::vector<char> buf(kBufLen);
  bool ok = true;
  while (gzgets(gz, buf.data(), kBufLen) != nullptr) {
    std::string line(buf.data());
    if (!line.empty() && line.back() == '\n') line.pop_back();
    line = TrimCr(line);
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::vector<std::string> f = SplitTab(line);
    if (f.size() < 4) {
      ok = false;
      if (error_message) *error_message = "fragment row requires chrom/start/end/barcode";
      break;
    }
    ++result->input_fragment_rows;
    const std::string& barcode = f[3];
    if (!cells.empty() && cells.find(barcode) == cells.end()) {
      continue;
    }
    int64_t repeat = 1;
    if (params.expand_counts && f.size() >= 5) {
      if (!ParseInt64(f[4], &repeat) || repeat < 1) {
        ok = false;
        if (error_message) *error_message = "invalid fragment count during projection";
        break;
      }
    }
    const std::string row = JoinFirstFields(f, std::min<size_t>(5, f.size()));
    if (grouped) {
      auto git = barcode_to_group.find(barcode);
      if (git == barcode_to_group.end()) {
        continue;
      }
      FILE* fp = group_files[git->second];
      if (!WriteRepeated(fp, row, repeat)) {
        ok = false;
        if (error_message) *error_message = "failed writing projected group BED";
        break;
      }
      result->group_rows[git->second] += repeat;
    } else {
      if (!WriteRepeated(single, row, repeat)) {
        ok = false;
        if (error_message) *error_message = "failed writing projected BED";
        break;
      }
    }
    result->written_bed_rows += repeat;
  }
  int zerr = Z_OK;
  if (!gzeof(gz)) {
    (void)gzerror(gz, &zerr);
    ok = false;
    if (error_message && error_message->empty()) {
      *error_message = "error reading fragments: " + params.fragments_path;
    }
  }
  if (gzclose(gz) != Z_OK) {
    ok = false;
    if (error_message && error_message->empty()) {
      *error_message = "error closing fragments: " + params.fragments_path;
    }
  }
  if (single && std::fclose(single) != 0) {
    ok = false;
    if (error_message && error_message->empty()) *error_message = "failed closing projected BED";
  }
  for (auto& kv : group_files) {
    if (std::fclose(kv.second) != 0) {
      ok = false;
      if (error_message && error_message->empty()) {
        *error_message = "failed closing projected group BED";
      }
    }
  }
  result->group_count = static_cast<int64_t>(result->group_bed_paths.size());
  if (!grouped) {
    result->group_count = 0;
  }
  return ok;
}

}  // namespace peaks
}  // namespace chromap
