#include "rapidmacs/stage_profile.h"

#include <cstdio>
#include <cstring>
#include <fstream>
#include <vector>

#include <sys/stat.h>

namespace chromap {
namespace peaks {

namespace {

double TimevalSec(const struct timeval& tv) {
  return static_cast<double>(tv.tv_sec) +
         static_cast<double>(tv.tv_usec) * 1e-6;
}

}  // namespace

double StageProfileCollector::RusageUserSec(const struct rusage& ru) {
  return TimevalSec(ru.ru_utime);
}

double StageProfileCollector::RusageSysSec(const struct rusage& ru) {
  return TimevalSec(ru.ru_stime);
}

StageProfileCollector::StageProfileCollector(const std::string& path)
    : path_(path), fp_(nullptr), total_written_(false) {
  if (path_.empty()) {
    return;
  }
  fp_ = std::fopen(path.c_str(), "w");
  if (fp_ == nullptr) {
    return;
  }
  WriteHeader();
}

StageProfileCollector::~StageProfileCollector() {
  if (fp_ != nullptr) {
    std::fclose(fp_);
    fp_ = nullptr;
  }
}

void StageProfileCollector::WriteHeader() {
  if (fp_ == nullptr) {
    return;
  }
  std::fprintf(fp_,
               "stage\twall_sec\tuser_sec\tsys_sec\tmax_rss_kb\tinput_rows\t"
               "output_rows\toutput_bytes\tnotes\n");
  std::fflush(fp_);
}

StageProfileCollector::Tick StageProfileCollector::Now() {
  Tick t;
  t.wall = std::chrono::steady_clock::now();
  std::memset(&t.ru, 0, sizeof(t.ru));
  (void)getrusage(RUSAGE_SELF, &t.ru);
  return t;
}

bool StageProfileCollector::FileMetrics(const std::string& path, int64_t* out_bytes,
                                        int64_t* out_lines) {
  if (out_bytes) {
    *out_bytes = 0;
  }
  if (out_lines) {
    *out_lines = 0;
  }
  if (path.empty()) {
    return false;
  }
  struct stat st;
  if (stat(path.c_str(), &st) != 0) {
    return false;
  }
  if (out_bytes) {
    *out_bytes = static_cast<int64_t>(st.st_size);
  }
  if (!out_lines) {
    return true;
  }
  std::ifstream in(path.c_str(), std::ios::binary);
  if (!in) {
    return false;
  }
  constexpr size_t kBuf = 1 << 20;
  std::vector<char> buf(kBuf);
  int64_t lines = 0;
  for (;;) {
    in.read(buf.data(), static_cast<std::streamsize>(buf.size()));
    const std::streamsize n = in.gcount();
    if (n <= 0) {
      break;
    }
    for (std::streamsize i = 0; i < n; ++i) {
      if (buf[static_cast<size_t>(i)] == '\n') {
        ++lines;
      }
    }
  }
  *out_lines = lines;
  return true;
}

void StageProfileCollector::WriteRow(const std::string& stage, double wall_sec,
                                     double user_sec, double sys_sec,
                                     int64_t max_rss_kb, int64_t input_rows,
                                     int64_t output_rows, int64_t output_bytes,
                                     const std::string& notes) {
  if (fp_ == nullptr) {
    return;
  }
  std::fprintf(fp_, "%s\t%.9f\t%.9f\t%.9f\t%lld\t%lld\t%lld\t%lld\t", stage.c_str(),
               wall_sec, user_sec, sys_sec, static_cast<long long>(max_rss_kb),
               static_cast<long long>(input_rows),
               static_cast<long long>(output_rows),
               static_cast<long long>(output_bytes));
  for (size_t i = 0; i < notes.size(); ++i) {
    const char c = notes[i];
    if (c == '\t' || c == '\n' || c == '\r') {
      std::fputc(' ', fp_);
    } else {
      std::fputc(c, fp_);
    }
  }
  std::fputc('\n', fp_);
  std::fflush(fp_);
}

void StageProfileCollector::Record(const std::string& stage, const Tick& t0,
                                   const Tick& t1, int64_t input_rows,
                                   int64_t output_rows, int64_t output_bytes,
                                   const std::string& notes) {
  if (fp_ == nullptr) {
    return;
  }
  const std::chrono::duration<double> wd = t1.wall - t0.wall;
  const double user = RusageUserSec(t1.ru) - RusageUserSec(t0.ru);
  const double sys = RusageSysSec(t1.ru) - RusageSysSec(t0.ru);
  const int64_t rss = static_cast<int64_t>(t1.ru.ru_maxrss);
  WriteRow(stage, wd.count(), user, sys, rss, input_rows, output_rows, output_bytes,
           notes);
}

void StageProfileCollector::RecordTotal(const Tick& t0, const Tick& t1,
                                        int64_t input_rows, int64_t output_rows,
                                        int64_t output_bytes,
                                        const std::string& notes) {
  if (fp_ == nullptr || total_written_) {
    return;
  }
  total_written_ = true;
  Record("total", t0, t1, input_rows, output_rows, output_bytes, notes);
}

StageProfileTotalGuard::StageProfileTotalGuard(StageProfileCollector* coll,
                                               StageProfileCollector::Tick run_start,
                                               int64_t* input_rows,
                                               int64_t* output_rows)
    : coll_(coll), run_start_(run_start), input_rows_(input_rows),
      output_rows_(output_rows) {}

StageProfileTotalGuard::~StageProfileTotalGuard() {
  if (coll_ == nullptr) {
    return;
  }
  const StageProfileCollector::Tick t1 = StageProfileCollector::Now();
  const int64_t in_r = input_rows_ ? *input_rows_ : 0;
  const int64_t out_r = output_rows_ ? *output_rows_ : 0;
  coll_->RecordTotal(run_start_, t1, in_r, out_r, 0, "");
}

}  // namespace peaks
}  // namespace chromap
