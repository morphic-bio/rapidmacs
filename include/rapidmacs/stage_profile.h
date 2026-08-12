#ifndef CHROMAP_PEAK_CALLER_STAGE_PROFILE_H_
#define CHROMAP_PEAK_CALLER_STAGE_PROFILE_H_

#include <chrono>
#include <cstdint>
#include <string>

#include <sys/resource.h>

namespace chromap {
namespace peaks {

// Opt-in stage timing for MACS3 FRAG peak pipeline diagnostics. Disabled when
// no collector is passed. Uses steady_clock + getrusage(RUSAGE_SELF); no extra
// dependencies.
class StageProfileCollector {
 public:
  struct Tick {
    std::chrono::steady_clock::time_point wall;
    struct rusage ru;
  };

  explicit StageProfileCollector(const std::string& path);
  ~StageProfileCollector();

  StageProfileCollector(const StageProfileCollector&) = delete;
  StageProfileCollector& operator=(const StageProfileCollector&) = delete;

  // False if a non-empty path was given but the TSV could not be opened.
  bool IsOpen() const { return fp_ != nullptr; }

  static Tick Now();
  // Optional line count reads the whole file; capture the stage end Tick (t1)
  // before calling this if it should not affect Record() wall_sec.
  static bool FileMetrics(const std::string& path, int64_t* out_bytes,
                          int64_t* out_lines);

  void Record(const std::string& stage, const Tick& t0, const Tick& t1,
              int64_t input_rows, int64_t output_rows, int64_t output_bytes,
              const std::string& notes);

  void RecordTotal(const Tick& t0, const Tick& t1, int64_t input_rows,
                   int64_t output_rows, int64_t output_bytes,
                   const std::string& notes);

 private:
  void WriteHeader();
  void WriteRow(const std::string& stage, double wall_sec, double user_sec,
                double sys_sec, int64_t max_rss_kb, int64_t input_rows,
                int64_t output_rows, int64_t output_bytes,
                const std::string& notes);
  static double RusageUserSec(const struct rusage& ru);
  static double RusageSysSec(const struct rusage& ru);

  std::string path_;
  FILE* fp_;
  bool total_written_;
};

// Writes a final "total" row when leaving main (success or early exit).
class StageProfileTotalGuard {
 public:
  StageProfileTotalGuard(StageProfileCollector* coll, StageProfileCollector::Tick run_start,
                         int64_t* input_rows, int64_t* output_rows);
  ~StageProfileTotalGuard();

 private:
  StageProfileCollector* coll_;
  StageProfileCollector::Tick run_start_;
  int64_t* input_rows_;
  int64_t* output_rows_;
};

}  // namespace peaks
}  // namespace chromap

#endif  // CHROMAP_PEAK_CALLER_STAGE_PROFILE_H_
