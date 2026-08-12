#ifndef RAPIDMACS_SIGNAC_BED_PROJECTION_H_
#define RAPIDMACS_SIGNAC_BED_PROJECTION_H_

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace chromap {
namespace peaks {

struct SignacBedProjectionParams {
  std::string fragments_path;
  std::string cells_path;
  std::string barcode_group_tsv_path;
  std::string group_column;
  std::string bed_out;
  std::string bed_out_dir;
  bool expand_counts = false;
};

struct SignacBedProjectionResult {
  int64_t input_fragment_rows = 0;
  int64_t written_bed_rows = 0;
  int64_t selected_cell_count = 0;
  int64_t group_count = 0;
  bool count_expansion_used = false;
  std::map<std::string, int64_t> group_rows;
  std::map<std::string, std::string> group_bed_paths;
  std::string bed_path;
};

bool ProjectFragmentsToSignacBed(const SignacBedProjectionParams& params,
                                 SignacBedProjectionResult* result,
                                 std::string* error_message);

}  // namespace peaks
}  // namespace chromap

#endif  // RAPIDMACS_SIGNAC_BED_PROJECTION_H_
