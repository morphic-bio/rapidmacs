#include "rapidmacs/fragments.h"

#include <zlib.h>

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <unordered_map>
#include <utility>

namespace macs3 {
namespace {

class FileFragmentIterator final : public FragmentIterator {
 public:
  explicit FileFragmentIterator(const std::string& path) : path_(path) {
    fp_ = gzopen(path.c_str(), "rb");
    if (fp_ == nullptr) {
      error_ = "failed to open fragments file: " + path;
    }
  }

  ~FileFragmentIterator() override {
    if (fp_ != nullptr) gzclose(fp_);
  }

  bool Next(FragmentRecord* out) override {
    if (fp_ == nullptr) return false;
    while (true) {
      char* line = nullptr;
      size_t len = 0;
      if (!ReadLine(&line, &len)) return false;
      ++line_no_;
      if (len > 0 && line[len - 1] == '\r') line[--len] = '\0';
      if (len == 0 || line[0] == '#') continue;

      // Cell Ranger/Chromap barcoded FRAG:
      //   chrom \t start \t end [\t barcode [\t count]]
      // Chromap bulk fragment BED:
      //   chrom \t start \t end \t N \t mapq \t strand \t count
      char* tab1 = std::strchr(line, '\t');
      if (tab1 == nullptr) {
        error_ = ErrPrefix() + "missing chrom tab"; return false;
      }
      *tab1 = '\0';
      if (line[0] == '\0') {
        error_ = ErrPrefix() + "empty chrom"; return false;
      }
      const char* chrom = line;

      char* tab2 = std::strchr(tab1 + 1, '\t');
      if (tab2 == nullptr) {
        error_ = ErrPrefix() + "missing end"; return false;
      }
      int32_t start = 0;
      if (!ParseInt32(tab1 + 1, tab2, &start)) {
        error_ = ErrPrefix() + "bad start"; return false;
      }
      if (start < 0) {
        error_ = ErrPrefix() + "negative start"; return false;
      }

      char* tab3 = std::strchr(tab2 + 1, '\t');
      const char* end_field_end = tab3 == nullptr ? line + len : tab3;
      int32_t end = 0;
      if (!ParseInt32(tab2 + 1, end_field_end, &end)) {
        error_ = ErrPrefix() + "bad end"; return false;
      }
      if (end <= start) {
        error_ = ErrPrefix() + "end <= start"; return false;
      }

      uint32_t count = 1;
      if (tab3 != nullptr) {
        char* tab4 = std::strchr(tab3 + 1, '\t');
        if (tab4 != nullptr) {
          const char* count_begin = tab4 + 1;
          const char* count_end = line + len;
          char* tab5 = std::strchr(tab4 + 1, '\t');
          char* tab6 = tab5 == nullptr ? nullptr : std::strchr(tab5 + 1, '\t');
          char* tab7 = tab6 == nullptr ? nullptr : std::strchr(tab6 + 1, '\t');
          if (tab6 != nullptr && tab7 == nullptr) {
            // Exactly seven columns are the Chromap bulk schema. Validate the
            // identifying BED fields so a malformed/extended five-column
            // FRAG row cannot silently reinterpret MAPQ as duplicate count.
            if (tab4 != tab3 + 2 || tab3[1] != 'N') {
              error_ = ErrPrefix() +
                       "bad Chromap bulk name (expected N)";
              return false;
            }
            int64_t mapq = 0;
            if (!ParseInt64(tab4 + 1, tab5, &mapq) || mapq < 0 || mapq > 60) {
              error_ = ErrPrefix() + "bad Chromap bulk MAPQ";
              return false;
            }
            if (tab6 != tab5 + 2 || (tab5[1] != '+' && tab5[1] != '-')) {
              error_ = ErrPrefix() + "bad Chromap bulk strand";
              return false;
            }
            count_begin = tab6 + 1;
          } else {
            // Historical 5-column FRAG behavior: column 5 is count and any
            // extra non-Chromap columns are ignored.
            count_end = tab5 == nullptr ? line + len : tab5;
          }
          int64_t c = 0;
          if (!ParseInt64(count_begin, count_end, &c)) {
            error_ = ErrPrefix() + "bad count (trailing junk)"; return false;
          }
          if (c < 1) {
            error_ = ErrPrefix() + "non-positive count"; return false;
          }
          if (static_cast<uint64_t>(c) >
              static_cast<uint64_t>(std::numeric_limits<uint32_t>::max())) {
            error_ = ErrPrefix() + "bad count"; return false;
          }
          count = static_cast<uint32_t>(c);
        }
      }

      out->chrom_id = InternChrom(chrom);
      out->start = start;
      out->end = end;
      out->count = count;
      return true;
    }
  }

  const std::vector<std::string>& ChromNames() const override { return chrom_names_; }
  const std::string& Error() const override { return error_; }

 private:
  static constexpr size_t kReadBufLen = 1U << 20;
  static constexpr size_t kMaxLineLen = 65535;

  static bool ParseInt64(const char* begin, const char* end, int64_t* out) {
    if (begin == end) return false;
    bool negative = false;
    if (*begin == '-') {
      negative = true;
      if (++begin == end) return false;
    }
    uint64_t value = 0;
    const uint64_t limit = negative
                               ? static_cast<uint64_t>(
                                     std::numeric_limits<int64_t>::max()) + 1U
                               : static_cast<uint64_t>(
                                     std::numeric_limits<int64_t>::max());
    for (const char* p = begin; p != end; ++p) {
      if (*p < '0' || *p > '9') return false;
      const uint64_t digit = static_cast<uint64_t>(*p - '0');
      if (value > (limit - digit) / 10U) return false;
      value = value * 10U + digit;
    }
    if (negative) {
      if (value == limit) {
        *out = std::numeric_limits<int64_t>::min();
      } else {
        *out = -static_cast<int64_t>(value);
      }
    } else {
      *out = static_cast<int64_t>(value);
    }
    return true;
  }

  static bool ParseInt32(const char* begin, const char* end, int32_t* out) {
    int64_t value = 0;
    if (!ParseInt64(begin, end, &value) ||
        value < std::numeric_limits<int32_t>::min() ||
        value > std::numeric_limits<int32_t>::max()) {
      return false;
    }
    *out = static_cast<int32_t>(value);
    return true;
  }

  bool ReadLine(char** line, size_t* len) {
    size_t carried = 0;
    for (;;) {
      if (read_pos_ == read_end_) {
        const int got = gzread(fp_, read_buf_, static_cast<unsigned int>(kReadBufLen));
        if (got < 0) {
          int errnum = 0;
          const char* msg = gzerror(fp_, &errnum);
          error_ = std::string("gz read error: ") +
                   (msg == nullptr ? "unknown error" : msg);
          return false;
        }
        if (got == 0) {
          if (carried == 0) return false;
          line_buf_[carried] = '\0';
          *line = line_buf_;
          *len = carried;
          return true;
        }
        read_pos_ = 0;
        read_end_ = static_cast<size_t>(got);
      }

      char* begin = read_buf_ + read_pos_;
      const size_t available = read_end_ - read_pos_;
      char* newline = static_cast<char*>(std::memchr(begin, '\n', available));
      const size_t piece = newline == nullptr
                               ? available
                               : static_cast<size_t>(newline - begin);
      if (carried + piece > kMaxLineLen) {
        error_ = path_ + ":" + std::to_string(line_no_ + 1) +
                 ": line exceeds buffer (65536 bytes)";
        return false;
      }
      if (newline != nullptr && carried == 0) {
        *newline = '\0';
        *line = begin;
        *len = piece;
        read_pos_ += piece + 1;
        return true;
      }
      if (piece > 0) {
        std::memcpy(line_buf_ + carried, begin, piece);
        carried += piece;
      }
      read_pos_ += piece;
      if (newline != nullptr) {
        ++read_pos_;
        line_buf_[carried] = '\0';
        *line = line_buf_;
        *len = carried;
        return true;
      }
    }
  }

  std::string ErrPrefix() const {
    return path_ + ":" + std::to_string(line_no_) + ": ";
  }

  int32_t InternChrom(const char* name) {
    if (last_chrom_id_ >= 0 && last_chrom_name_ == name) {
      return last_chrom_id_;
    }
    auto it = chrom_index_.find(name);
    if (it == chrom_index_.end()) {
      const int32_t id = static_cast<int32_t>(chrom_names_.size());
      chrom_names_.emplace_back(name);
      it = chrom_index_.emplace(chrom_names_.back(), id).first;
    }
    last_chrom_name_ = name;
    last_chrom_id_ = it->second;
    return last_chrom_id_;
  }

  std::string path_;
  gzFile fp_ = nullptr;
  char read_buf_[kReadBufLen];
  char line_buf_[kMaxLineLen + 1];
  size_t read_pos_ = 0;
  size_t read_end_ = 0;
  int64_t line_no_ = 0;
  std::vector<std::string> chrom_names_;
  std::unordered_map<std::string, int32_t> chrom_index_;
  std::string last_chrom_name_;
  int32_t last_chrom_id_ = -1;
  std::string error_;
};

}  // namespace

std::unique_ptr<FragmentIterator> OpenFileFragmentIterator(
    const std::string& path) {
  return std::unique_ptr<FragmentIterator>(new FileFragmentIterator(path));
}

namespace {

class VectorFragmentIterator final : public FragmentIterator {
 public:
  VectorFragmentIterator(std::vector<FragmentRecord> fragments,
                         std::vector<std::string> chrom_names)
      : fragments_(std::move(fragments)),
        chrom_names_(std::move(chrom_names)),
        cursor_(0) {
    std::sort(fragments_.begin(), fragments_.end(),
              [](const FragmentRecord& a, const FragmentRecord& b) {
                if (a.chrom_id != b.chrom_id) return a.chrom_id < b.chrom_id;
                if (a.start != b.start) return a.start < b.start;
                return a.end < b.end;
              });
  }

  bool Next(FragmentRecord* out) override {
    if (cursor_ >= fragments_.size()) return false;
    *out = fragments_[cursor_++];
    return true;
  }

  const std::vector<std::string>& ChromNames() const override {
    return chrom_names_;
  }

  const std::string& Error() const override { return error_; }

 private:
  std::vector<FragmentRecord> fragments_;
  std::vector<std::string> chrom_names_;
  size_t cursor_;
  std::string error_;
};

}  // namespace

std::unique_ptr<FragmentIterator> WrapVectorFragmentIterator(
    std::vector<FragmentRecord> fragments,
    std::vector<std::string> chrom_names) {
  return std::unique_ptr<FragmentIterator>(new VectorFragmentIterator(
      std::move(fragments), std::move(chrom_names)));
}

}  // namespace macs3
