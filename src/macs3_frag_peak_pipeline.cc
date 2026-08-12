#include "rapidmacs/macs3_frag_peak_pipeline.h"

#include <algorithm>
#include <cmath>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <limits>

#include <ftw.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "rapidmacs/bdgpeakcall.h"
#include "rapidmacs/macs3_frag_sweep.h"
#include "rapidmacs/macs3_frag_workspace.h"
#include "rapidmacs/peak_io.h"
#include "rapidmacs/stage_profile.h"

namespace chromap {
namespace peaks {

namespace {

bool DirExists(const std::string& path) {
  struct stat st;
  return stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
}

bool MakeDir(const std::string& path) {
  if (path.empty()) {
    return false;
  }
  if (mkdir(path.c_str(), 0755) == 0) {
    return true;
  }
  return errno == EEXIST && DirExists(path);
}

int UnlinkPath(const char* fpath, const struct stat* sb, int typeflag,
               struct FTW* ftwbuf) {
  (void)sb;
  (void)typeflag;
  (void)ftwbuf;
  return remove(fpath);
}

bool RemoveTree(const std::string& root) {
  if (root.empty()) {
    return true;
  }
  return nftw(root.c_str(), UnlinkPath, 64, FTW_DEPTH | FTW_PHYS) == 0;
}

std::string DefaultWorkParent(const std::string& work_dir_parent) {
  if (!work_dir_parent.empty()) {
    return work_dir_parent;
  }
  const char* tmp = getenv("TMPDIR");
  if (tmp && tmp[0] != '\0') {
    return std::string(tmp);
  }
  return "/tmp";
}

int NormalizeThreads(int requested, size_t n_tasks) {
  if (requested < 1 || n_tasks <= 1) {
    return 1;
  }
  const size_t rt = static_cast<size_t>(requested);
  return static_cast<int>(std::min(rt, n_tasks));
}

bool ResolvePeakThreshold(const Macs3FragPeakPipelineParams& params,
                          float* cutoff, bool* qvalue_mode) {
  if (cutoff == nullptr || qvalue_mode == nullptr) return false;
  if (params.threshold_mode == Macs3FragThresholdMode::kQValue) {
    if (std::isnan(params.qvalue_cutoff) || params.qvalue_cutoff <= 0.0 ||
        params.qvalue_cutoff > 1.0) {
      return false;
    }
    *cutoff = static_cast<float>(-std::log10(params.qvalue_cutoff));
    *qvalue_mode = true;
    return true;
  }
  if (params.threshold_mode != Macs3FragThresholdMode::kPValue ||
      std::isnan(params.bdgpeakcall_cutoff) ||
      params.bdgpeakcall_cutoff <= 0.f) {
    return false;
  }
  *cutoff = params.bdgpeakcall_cutoff;
  *qvalue_mode = params.cutoff_is_qvalue;
  return true;
}

}  // namespace

float BdgPeakCallCutoffFromPValue(double p_value) {
  if (p_value <= 0.0 || p_value > 1.0 || std::isnan(p_value)) {
    return -1.f;
  }
  return static_cast<float>(-std::log10(p_value));
}

bool RunMacs3FragPeakPipelineFromFragments(
    std::vector<ChromFragments>* by_chrom,
    const Macs3FragPeakPipelineParams& params,
    const Macs3FragPeakPipelinePaths& user_paths,
    const std::string& narrowpeak_out, const std::string& summits_out,
    const std::string& keep_intermediates_dir,
    const std::string& work_dir_parent, std::string* work_dir_used,
    std::string* error_message, StageProfileCollector* stage_profile,
    bool release_fragments_after_workspace,
    const std::vector<ChromFragments>* control_by_chrom) {
  auto fail = [&](const char* msg) {
    if (error_message) {
      *error_message = msg;
    }
    return false;
  };
  if (by_chrom == nullptr || narrowpeak_out.empty() || summits_out.empty()) {
    return fail("narrowPeak and summits output paths are required");
  }
  if (params.min_length < 1 || params.max_gap < 0 || params.slocal_bp < 0 ||
      params.llocal_bp < 1 || params.slocal_bp > params.llocal_bp ||
      params.effective_genome_size < 1 || params.peak_caller_threads < 1) {
    return fail("invalid MACS3 FRAG peak pipeline numeric parameters");
  }
  float threshold_cutoff = 0.f;
  bool qvalue_mode = false;
  if (!ResolvePeakThreshold(params, &threshold_cutoff, &qvalue_mode)) {
    return fail("invalid p-value/q-value peak threshold");
  }
  const int peak_threads =
      NormalizeThreads(params.peak_caller_threads, by_chrom->size());

  int64_t n_frag_rows = 0;
  if (stage_profile) {
    for (const ChromFragments& ch : *by_chrom) {
      n_frag_rows += static_cast<int64_t>(ch.frags.size());
    }
  }

  const bool export_treat = !keep_intermediates_dir.empty() || !user_paths.treat_bdg.empty();
  const bool export_lambda = !keep_intermediates_dir.empty() || !user_paths.lambda_bdg.empty();
  const bool export_ppois = !keep_intermediates_dir.empty() || !user_paths.ppois_bdg.empty();
  const bool any_intermediate_export = export_treat || export_lambda || export_ppois;

  std::string work;
  bool own_workdir = false;
  if (!keep_intermediates_dir.empty()) {
    work = keep_intermediates_dir;
    if (!MakeDir(work)) {
      return fail("failed to create --macs3-frag-keep-intermediates directory");
    }
  } else if (any_intermediate_export) {
    const std::string parent = DefaultWorkParent(work_dir_parent);
    if (!DirExists(parent) && !MakeDir(parent)) {
      return fail("failed to create temp parent directory for MACS3 FRAG workdir");
    }
    work = parent + "/chromap_macs3_frag_XXXXXX";
    std::vector<char> tmpl(work.begin(), work.end() + 1);
    if (mkdtemp(tmpl.data()) == nullptr) {
      return fail("mkdtemp failed for MACS3 FRAG peak pipeline");
    }
    work.assign(tmpl.data());
    own_workdir = true;
  }
  if (work_dir_used) {
    *work_dir_used = work;
  }

  const std::string treat =
      user_paths.treat_bdg.empty() ? (work.empty() ? "" : (work + "/macs3_frag_treat_pileup.bdg"))
                                   : user_paths.treat_bdg;
  const std::string lambda =
      user_paths.lambda_bdg.empty() ? (work.empty() ? "" : (work + "/macs3_frag_control_lambda.bdg"))
                                    : user_paths.lambda_bdg;
  const std::string ppois =
      user_paths.ppois_bdg.empty() ? (work.empty() ? "" : (work + "/macs3_frag_score_ppois.bdg"))
                                   : user_paths.ppois_bdg;

  Macs3FragWorkspaceParams workspace_params;
  workspace_params.frag_pileup_max_count = params.frag_pileup_max_count;
  workspace_params.macs3_uint8_counts = params.macs3_uint8_counts;
  workspace_params.effective_genome_size = params.effective_genome_size;
  workspace_params.llocal_bp = params.llocal_bp;
  workspace_params.score_pseudocount = params.score_pseudocount;
  workspace_params.byte_parity_summits = params.byte_parity_summits;
  workspace_params.build_treatment_local_lambda_events =
      control_by_chrom == nullptr && !params.no_lambda;
  workspace_params.capture_raw_pos_array_endpoints =
      control_by_chrom == nullptr && !params.no_lambda;
  Macs3FragPeakWorkspace workspace;
  {
    auto build_workspace = [&]() {
      const bool ok =
          BuildMacs3FragWorkspaceFromFragments(*by_chrom, workspace_params, &workspace);
      if (ok && release_fragments_after_workspace) {
        std::vector<ChromFragments>().swap(*by_chrom);
      }
      return ok;
    };
    if (stage_profile) {
      const StageProfileCollector::Tick t0 = StageProfileCollector::Now();
      if (!build_workspace()) {
        if (own_workdir) {
          RemoveTree(work);
        }
        return fail("BuildMacs3FragWorkspaceFromFragments failed");
      }
      stage_profile->Record("build_workspace", t0, StageProfileCollector::Now(),
                            n_frag_rows, static_cast<int64_t>(workspace.chrom_names.size()),
                            0, "in_memory_events");
    } else {
      if (!build_workspace()) {
        if (own_workdir) {
          RemoveTree(work);
        }
        return fail("BuildMacs3FragWorkspaceFromFragments failed");
      }
    }
  }

  // The treatment-only path can finalize its reconstructed MACS3 position
  // arrays immediately. Controlled PE mode consumes the captured treatment
  // endpoints together with raw control-window positions during lambda
  // finalization below.
  if (control_by_chrom == nullptr && !params.no_lambda &&
      workspace.params.byte_parity_summits &&
      !FinalizeMacs3FragPosArrayPositions(&workspace, peak_threads)) {
    if (own_workdir) {
      RemoveTree(work);
    }
    return fail("FinalizeMacs3FragPosArrayPositions failed");
  }

  {
    auto finalize_treat = [&]() {
      return FinalizeMacs3FragTreatTracks(&workspace, peak_threads);
    };
    if (stage_profile) {
      const StageProfileCollector::Tick t0 = StageProfileCollector::Now();
      if (!finalize_treat()) {
        if (own_workdir) {
          RemoveTree(work);
        }
        return fail("FinalizeMacs3FragTreatTracks failed");
      }
      const StageProfileCollector::Tick t1 = StageProfileCollector::Now();
      stage_profile->Record("finalize_treat_tracks", t0, t1, n_frag_rows,
                            static_cast<int64_t>(workspace.treat_tracks.size()), 0,
                            "in_memory_tracks;threads=" + std::to_string(peak_threads));
    } else {
      if (!finalize_treat()) {
        if (own_workdir) {
          RemoveTree(work);
        }
        return fail("FinalizeMacs3FragTreatTracks failed");
      }
    }
  }
  {
    auto finalize_lam = [&]() {
      if (control_by_chrom != nullptr) {
        Macs3FragControlParams cp;
        cp.slocal_bp = params.slocal_bp;
        cp.treatment_fragment_length_bp =
            params.treatment_fragment_length_bp;
        cp.keep_duplicates_all = params.keep_duplicates_all;
        cp.scale_to_large = params.scale_to_large;
        cp.no_lambda = params.no_lambda;
        return FinalizeMacs3FragLambdaTracksWithControl(
            &workspace, *control_by_chrom, cp, peak_threads);
      }
      if (params.no_lambda) {
        return FinalizeMacs3FragGlobalLambdaTracks(&workspace, peak_threads);
      }
      return FinalizeMacs3FragLambdaTracks(&workspace, peak_threads);
    };
    if (stage_profile) {
      const StageProfileCollector::Tick t0 = StageProfileCollector::Now();
      if (!finalize_lam()) {
        if (own_workdir) {
          RemoveTree(work);
        }
        return fail("FinalizeMacs3FragLambdaTracks failed");
      }
      const StageProfileCollector::Tick t1 = StageProfileCollector::Now();
      stage_profile->Record("finalize_lambda_tracks", t0, t1, n_frag_rows,
                            static_cast<int64_t>(workspace.lambda_tracks.size()), 0,
                            "in_memory_tracks;threads=" + std::to_string(peak_threads));
    } else {
      if (!finalize_lam()) {
        if (own_workdir) {
          RemoveTree(work);
        }
        return fail("FinalizeMacs3FragLambdaTracks failed");
      }
    }
  }
  if (export_ppois) {
    auto finalize_score = [&]() {
      return FinalizeMacs3FragPpoisTracks(
          &workspace, peak_threads, !params.macs3_float_score_math);
    };
    if (stage_profile) {
      const StageProfileCollector::Tick t0 = StageProfileCollector::Now();
      if (!finalize_score()) {
        if (own_workdir) {
          RemoveTree(work);
        }
        return fail("FinalizeMacs3FragPpoisTracks failed");
      }
      const StageProfileCollector::Tick t1 = StageProfileCollector::Now();
      stage_profile->Record("finalize_ppois_tracks", t0, t1, n_frag_rows,
                            static_cast<int64_t>(workspace.ppois_tracks.size()), 0,
                            "in_memory_tracks;threads=" + std::to_string(peak_threads));
    } else {
      if (!finalize_score()) {
        if (own_workdir) {
          RemoveTree(work);
        }
        return fail("FinalizeMacs3FragPpoisTracks failed");
      }
    }
  }
  auto export_track = [&](const char* stage, const std::string& path,
                          const std::vector<Macs3BdgTrack>& tracks) {
    if (path.empty()) {
      return true;
    }
    if (stage_profile) {
      const StageProfileCollector::Tick t0 = StageProfileCollector::Now();
      if (!WriteMacs3BdgTracksBedGraph(path, workspace.chrom_names, tracks)) {
        return false;
      }
      const StageProfileCollector::Tick t1 = StageProfileCollector::Now();
      const StageProfileCollector::Tick tm0 = StageProfileCollector::Now();
      int64_t bytes = 0, lines = 0;
      (void)StageProfileCollector::FileMetrics(path, &bytes, &lines);
      const StageProfileCollector::Tick tm1 = StageProfileCollector::Now();
      stage_profile->Record(stage, t0, t1, n_frag_rows, lines, bytes, "");
      stage_profile->Record(std::string("profile_file_metrics_") + stage, tm0, tm1,
                            n_frag_rows, lines, bytes, "line_count_and_size_scan");
      return true;
    }
    return WriteMacs3BdgTracksBedGraph(path, workspace.chrom_names, tracks);
  };
  if (!export_track("export_treat_pileup", treat, workspace.treat_tracks)) {
    if (own_workdir) {
      RemoveTree(work);
    }
    return fail("WriteFragSpanPileupBedGraph failed");
  }
  if (!export_track("export_control_lambda", lambda, workspace.lambda_tracks)) {
    if (own_workdir) {
      RemoveTree(work);
    }
    return fail("WriteFragMacs3NoControlLambdaBedGraph failed");
  }
  if (export_ppois) {
    if (!export_track("export_score_ppois", ppois, workspace.ppois_tracks)) {
      if (own_workdir) {
        RemoveTree(work);
      }
      return fail("WriteMacs3BdgTracksBedGraph ppois failed");
    }
    if (!RunMacs3FragPpoisNarrowPeaksFromTracks(
            workspace.chrom_names, workspace.ppois_tracks, workspace.treat_tracks,
            workspace.lambda_tracks, threshold_cutoff, params.min_length,
            params.max_gap, params.score_pseudocount, narrowpeak_out, summits_out,
            stage_profile, n_frag_rows, params.name_prefix,
            qvalue_mode,
            params.macs3_float_score_math,
            workspace.params.byte_parity_summits
                ? &workspace.pos_array_positions
                : nullptr,
            params.call_summits)) {
      if (own_workdir) {
        RemoveTree(work);
      }
      return fail("RunMacs3FragPpoisNarrowPeaks failed");
    }
  } else {
    StageProfileCollector::Tick t_ppois_start{};
    StageProfileCollector::Tick t_regions_start{};
    if (stage_profile) {
      t_ppois_start = StageProfileCollector::Now();
    }
    const size_t n_chrom = workspace.chrom_names.size();
    std::vector<std::vector<Macs3FragNarrowPeakRow>> rows_by_chrom(n_chrom);
    // Retain per-chrom ppois tracks so we can build a genome-wide pqtable
    // for MACS3-compatible qscore lookup. Discarded after the table builds.
    std::vector<Macs3BdgTrack> ppois_tracks_all(n_chrom);
    std::vector<char> ok(n_chrom, 1);
    const int threads = NormalizeThreads(peak_threads, n_chrom);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) if(threads > 1)
    for (int64_t ii = 0; ii < static_cast<int64_t>(n_chrom); ++ii) {
      const size_t i = static_cast<size_t>(ii);
      if (!BuildMacs3FragPpoisTrack(workspace.treat_tracks[i],
                                    workspace.lambda_tracks[i],
                                    params.score_pseudocount,
                                    &ppois_tracks_all[i],
                                    !params.macs3_float_score_math)) {
        ok[i] = 0;
      }
    }
    if (std::find(ok.begin(), ok.end(), 0) != ok.end()) {
      if (own_workdir) RemoveTree(work);
      return fail("chromosome-at-a-time ppois construction failed");
    }
    std::map<double, double> pqtable;
    BuildMacs3PqTable(ppois_tracks_all, &pqtable,
                      params.macs3_float_score_math);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) if(threads > 1)
    for (int64_t ii = 0; ii < static_cast<int64_t>(n_chrom); ++ii) {
      const size_t i = static_cast<size_t>(ii);
      const std::vector<int32_t>* pap =
          (workspace.params.byte_parity_summits &&
           i < workspace.pos_array_positions.size())
              ? &workspace.pos_array_positions[i]
              : nullptr;
      if (!CollectMacs3FragNarrowPeakRowsFromTracks(
              workspace.chrom_names[i], ppois_tracks_all[i],
              workspace.treat_tracks[i],
              workspace.lambda_tracks[i], threshold_cutoff, params.min_length,
              params.max_gap, params.score_pseudocount, &rows_by_chrom[i], pap,
              params.macs3_float_score_math,
              qvalue_mode ? &pqtable : nullptr,
              params.call_summits)) {
        ok[i] = 0;
        continue;
      }
      if (i < workspace.pos_array_positions.size()) {
        std::vector<int32_t>().swap(workspace.pos_array_positions[i]);
      }
      if (!export_treat) {
        std::vector<int32_t>().swap(workspace.treat_tracks[i].p);
        std::vector<float>().swap(workspace.treat_tracks[i].v);
      }
      if (!export_lambda) {
        std::vector<int32_t>().swap(workspace.lambda_tracks[i].p);
        std::vector<float>().swap(workspace.lambda_tracks[i].v);
      }
    }
    if (std::find(ok.begin(), ok.end(), 0) != ok.end()) {
      if (own_workdir) {
        RemoveTree(work);
      }
      return fail("chromosome-at-a-time ppois peak collection failed");
    }
    std::vector<Macs3FragNarrowPeakRow> rows;
    size_t row_count = 0;
    for (const auto& chrom_rows : rows_by_chrom) {
      row_count += chrom_rows.size();
    }
    rows.reserve(row_count);
    for (auto& chrom_rows : rows_by_chrom) {
      rows.insert(rows.end(), chrom_rows.begin(), chrom_rows.end());
    }
    if (stage_profile) {
      t_regions_start = StageProfileCollector::Now();
      stage_profile->Record("finalize_ppois_tracks", t_ppois_start, t_regions_start,
                            n_frag_rows, static_cast<int64_t>(workspace.chrom_names.size()),
                            0, "chromosome_at_a_time;includes_region_collection;threads=" +
                                   std::to_string(threads));
      stage_profile->Record("bdgpeakcall_regions", t_regions_start, t_regions_start,
                            n_frag_rows, static_cast<int64_t>(rows.size()), 0,
                            "included_in_finalize_ppois_tracks");
    }
    StageProfileCollector::Tick t_np_start{};
    if (stage_profile) {
      t_np_start = StageProfileCollector::Now();
    }
    if (!WriteMacs3FragNarrowPeakRows(&rows, narrowpeak_out, summits_out,
                                      &pqtable, params.name_prefix,
                                      qvalue_mode)) {
      if (own_workdir) {
        RemoveTree(work);
      }
      return fail("WriteMacs3FragNarrowPeakRows failed");
    }
    if (stage_profile) {
      const StageProfileCollector::Tick t_summits_end = StageProfileCollector::Now();
      const StageProfileCollector::Tick tm0 = StageProfileCollector::Now();
      int64_t bbytes = 0, sbytes = 0;
      if (!narrowpeak_out.empty()) {
        (void)StageProfileCollector::FileMetrics(narrowpeak_out, &bbytes, nullptr);
      }
      if (!summits_out.empty()) {
        (void)StageProfileCollector::FileMetrics(summits_out, &sbytes, nullptr);
      }
      const StageProfileCollector::Tick tm1 = StageProfileCollector::Now();
      stage_profile->Record("narrowpeak_summits", t_np_start, t_summits_end, n_frag_rows,
                            static_cast<int64_t>(rows.size()), bbytes + sbytes,
                            "chromosome_at_a_time");
      stage_profile->Record("profile_file_metrics_narrowpeak_summits", tm0, tm1, 0, 0,
                            bbytes + sbytes, "stat_size_only");
    }
  }
  {
    auto do_cleanup = [&]() {
      if (own_workdir) {
        RemoveTree(work);
      }
    };
    if (stage_profile) {
      const StageProfileCollector::Tick t0 = StageProfileCollector::Now();
      do_cleanup();
      const std::string notes = own_workdir ? "removed_temp_workdir" : "no_temp_workdir";
      stage_profile->Record("cleanup", t0, StageProfileCollector::Now(), 0, 0, 0, notes);
    } else {
      do_cleanup();
    }
  }
  return true;
}

bool RunMacs3FragPeakPipelineFromIterator(
    macs3::FragmentIterator& iter,
    const Macs3FragPeakPipelineParams& params,
    const Macs3FragPeakPipelinePaths& user_paths,
    const std::string& narrowpeak_out, const std::string& summits_out,
    const std::string& keep_intermediates_dir,
    const std::string& work_dir_parent, std::string* work_dir_used,
    std::string* error_message, StageProfileCollector* stage_profile) {
  auto fail = [&](const char* msg) {
    if (error_message) *error_message = msg;
    return false;
  };
  if (narrowpeak_out.empty() || summits_out.empty()) {
    return fail("narrowPeak and summits output paths are required");
  }
  Macs3FragWorkspaceParams workspace_params;
  workspace_params.frag_pileup_max_count = params.frag_pileup_max_count;
  workspace_params.macs3_uint8_counts = params.macs3_uint8_counts;
  workspace_params.effective_genome_size = params.effective_genome_size;
  workspace_params.llocal_bp = params.llocal_bp;
  workspace_params.score_pseudocount = params.score_pseudocount;
  workspace_params.byte_parity_summits = params.byte_parity_summits;
  workspace_params.build_treatment_local_lambda_events = !params.no_lambda;
  workspace_params.capture_raw_pos_array_endpoints = !params.no_lambda;
  Macs3FragPeakWorkspace workspace;
  // Init with empty chrom_names; we grow events_by_chrom on demand as the
  // iterator interns new chrom ids in first-seen order, then sync the name
  // table from iter.ChromNames() once the stream is drained.
  if (!InitMacs3FragWorkspace({}, workspace_params, &workspace)) {
    return fail("InitMacs3FragWorkspace failed");
  }

  StageProfileCollector::Tick t_drain_start{};
  if (stage_profile) t_drain_start = StageProfileCollector::Now();
  int64_t n_frag_rows = 0;
  macs3::FragmentRecord rec;
  while (iter.Next(&rec)) {
    if (rec.chrom_id < 0 || rec.end <= rec.start) continue;
    if (!EnsureMacs3FragWorkspaceChromCapacity(
            &workspace, static_cast<uint32_t>(rec.chrom_id))) {
      return fail("EnsureMacs3FragWorkspaceChromCapacity failed");
    }
    if (!AddMacs3FragWorkspaceFragment(
            &workspace, static_cast<uint32_t>(rec.chrom_id),
            rec.start, rec.end, static_cast<int32_t>(rec.count))) {
      return fail("AddMacs3FragWorkspaceFragment failed");
    }
    ++n_frag_rows;
  }
  if (!iter.Error().empty()) {
    if (error_message) *error_message = iter.Error();
    return false;
  }
  // Sync chrom name table — workspace owns its names from here on.
  workspace.chrom_names = iter.ChromNames();
  // Pad chrom_names if any chrom_ids were observed beyond the iterator's
  // current name table (shouldn't happen with FileFragmentIterator but
  // guards against partial-stream scenarios).
  if (workspace.chrom_names.size() < workspace.events_by_chrom.size()) {
    workspace.chrom_names.resize(workspace.events_by_chrom.size());
  }
  if (stage_profile) {
    stage_profile->Record("drain_iterator_to_workspace", t_drain_start,
                          StageProfileCollector::Now(), n_frag_rows,
                          static_cast<int64_t>(workspace.chrom_names.size()),
                          0, "in_memory_events");
  }

  return RunMacs3FragPeakPipelineFromWorkspace(
      &workspace, params, user_paths, narrowpeak_out, summits_out,
      keep_intermediates_dir, work_dir_parent, work_dir_used, error_message,
      stage_profile);
}

bool RunMacs3FragPeakPipelineFromSortedIterator(
    macs3::FragmentIterator& iter,
    const Macs3FragPeakPipelineParams& params,
    const Macs3FragPeakPipelinePaths& user_paths,
    const std::string& narrowpeak_out, const std::string& summits_out,
    const std::string& keep_intermediates_dir,
    const std::string& work_dir_parent, std::string* work_dir_used,
    std::string* error_message, StageProfileCollector* stage_profile) {
  auto fail = [&](const std::string& msg) {
    if (error_message) *error_message = msg;
    return false;
  };
  if (narrowpeak_out.empty() || summits_out.empty()) {
    return fail("narrowPeak and summits output paths are required");
  }
  Macs3FragWorkspaceParams workspace_params;
  workspace_params.frag_pileup_max_count = params.frag_pileup_max_count;
  workspace_params.macs3_uint8_counts = params.macs3_uint8_counts;
  workspace_params.effective_genome_size = params.effective_genome_size;
  workspace_params.llocal_bp = params.llocal_bp;
  workspace_params.score_pseudocount = params.score_pseudocount;
  workspace_params.byte_parity_summits = params.byte_parity_summits;
  workspace_params.build_treatment_local_lambda_events = !params.no_lambda;
  workspace_params.capture_raw_pos_array_endpoints = !params.no_lambda;

  StageProfileCollector::Tick t_load_start{};
  if (stage_profile) t_load_start = StageProfileCollector::Now();
  int64_t n_frag_rows = 0;
  std::vector<ChromFragments> by_chrom;
  int32_t previous_chrom_id = -1;
  int32_t previous_start = -1;
  macs3::FragmentRecord rec;
  while (iter.Next(&rec)) {
    if (rec.chrom_id < 0 || rec.end <= rec.start) continue;
    if (rec.chrom_id < previous_chrom_id) {
      return fail("input not chrom-grouped (chrom_id decreased after transition)");
    }
    if (rec.chrom_id != previous_chrom_id) {
      previous_chrom_id = rec.chrom_id;
      previous_start = -1;
      if (static_cast<size_t>(rec.chrom_id) >= by_chrom.size()) {
        by_chrom.resize(static_cast<size_t>(rec.chrom_id) + 1);
      }
    }
    if (rec.start < previous_start) {
      return fail("input not start-sorted within chromosome");
    }
    if (rec.count >
        static_cast<uint32_t>(std::numeric_limits<int32_t>::max())) {
      return fail("fragment count exceeds int32 range");
    }
    previous_start = rec.start;
    Fragment f;
    f.start = rec.start;
    f.end = rec.end;
    f.count = static_cast<int32_t>(rec.count);
    ChromFragments& ch = by_chrom[static_cast<size_t>(rec.chrom_id)];
    ch.frags.push_back(f);
    ch.max_end = std::max(ch.max_end, f.end);
    ++n_frag_rows;
  }
  if (!iter.Error().empty()) return fail(iter.Error());
  const std::vector<std::string>& chrom_names = iter.ChromNames();
  if (chrom_names.size() < by_chrom.size()) {
    return fail("iterator chromosome table is smaller than observed chrom ids");
  }
  for (size_t i = 0; i < by_chrom.size(); ++i) by_chrom[i].name = chrom_names[i];
  if (stage_profile) {
    stage_profile->Record("load_fragments", t_load_start,
                          StageProfileCollector::Now(), n_frag_rows,
                          static_cast<int64_t>(by_chrom.size()), 0,
                          "block_parser;chrom_grouped;start_sorted");
  }

  Macs3FragPeakWorkspace workspace;
  workspace.params = workspace_params;
  StageProfileCollector::Tick t_build_start{};
  if (stage_profile) t_build_start = StageProfileCollector::Now();
  if (!BuildMacs3FragTracksFromSortedFragments(
          &by_chrom, workspace_params, &workspace,
          params.peak_caller_threads, true)) {
    return fail("BuildMacs3FragTracksFromSortedFragments failed");
  }
  if (stage_profile) {
    stage_profile->Record("build_workspace", t_build_start,
                          StageProfileCollector::Now(), n_frag_rows,
                          static_cast<int64_t>(workspace.chrom_names.size()), 0,
                          "parallel_sorted_endpoint_merge");
  }

  return RunMacs3FragPeakPipelineFromWorkspace(
      &workspace, params, user_paths, narrowpeak_out, summits_out,
      keep_intermediates_dir, work_dir_parent, work_dir_used, error_message,
      stage_profile);
}

bool RunMacs3FragPeakPipelineFromWorkspace(
    Macs3FragPeakWorkspace* workspace,
    const Macs3FragPeakPipelineParams& params,
    const Macs3FragPeakPipelinePaths& user_paths,
    const std::string& narrowpeak_out, const std::string& summits_out,
    const std::string& keep_intermediates_dir,
    const std::string& work_dir_parent, std::string* work_dir_used,
    std::string* error_message,
    StageProfileCollector* stage_profile) {
  auto fail = [&](const char* msg) {
    if (error_message) {
      *error_message = msg;
    }
    return false;
  };
  if (workspace == nullptr) {
    return fail("MACS3 FRAG workspace is null");
  }
  if (narrowpeak_out.empty() || summits_out.empty()) {
    return fail("narrowPeak and summits output paths are required");
  }
  if (params.min_length < 1 || params.max_gap < 0 || params.slocal_bp < 0 ||
      params.llocal_bp < 1 || params.slocal_bp > params.llocal_bp ||
      params.effective_genome_size < 1 || params.peak_caller_threads < 1) {
    return fail("invalid MACS3 FRAG peak pipeline numeric parameters");
  }
  float threshold_cutoff = 0.f;
  bool qvalue_mode = false;
  if (!ResolvePeakThreshold(params, &threshold_cutoff, &qvalue_mode)) {
    return fail("invalid p-value/q-value peak threshold");
  }
  const int peak_threads =
      NormalizeThreads(params.peak_caller_threads, workspace->chrom_names.size());

  const bool export_treat = !keep_intermediates_dir.empty() || !user_paths.treat_bdg.empty();
  const bool export_lambda = !keep_intermediates_dir.empty() || !user_paths.lambda_bdg.empty();
  const bool export_ppois = !keep_intermediates_dir.empty() || !user_paths.ppois_bdg.empty();
  const bool any_intermediate_export = export_treat || export_lambda || export_ppois;

  std::string work;
  bool own_workdir = false;
  if (!keep_intermediates_dir.empty()) {
    work = keep_intermediates_dir;
    if (!MakeDir(work)) {
      return fail("failed to create --macs3-frag-keep-intermediates directory");
    }
  } else if (any_intermediate_export) {
    const std::string parent = DefaultWorkParent(work_dir_parent);
    if (!DirExists(parent) && !MakeDir(parent)) {
      return fail("failed to create temp parent directory for MACS3 FRAG workdir");
    }
    work = parent + "/chromap_macs3_frag_XXXXXX";
    std::vector<char> tmpl(work.begin(), work.end() + 1);
    if (mkdtemp(tmpl.data()) == nullptr) {
      return fail("mkdtemp failed for MACS3 FRAG peak pipeline");
    }
    work.assign(tmpl.data());
    own_workdir = true;
  }
  if (work_dir_used) {
    *work_dir_used = work;
  }

  const std::string treat =
      user_paths.treat_bdg.empty() ? (work.empty() ? "" : (work + "/macs3_frag_treat_pileup.bdg"))
                                   : user_paths.treat_bdg;
  const std::string lambda =
      user_paths.lambda_bdg.empty() ? (work.empty() ? "" : (work + "/macs3_frag_control_lambda.bdg"))
                                    : user_paths.lambda_bdg;
  const std::string ppois =
      user_paths.ppois_bdg.empty() ? (work.empty() ? "" : (work + "/macs3_frag_score_ppois.bdg"))
                                   : user_paths.ppois_bdg;

  auto cleanup_on_fail = [&]() {
    if (own_workdir) {
      RemoveTree(work);
    }
  };
  // Sort chrom_names + parallel track vectors lexicographically so the
  // exported bedGraphs (treat / lambda / ppois) are in a deterministic
  // chrom order regardless of producer (events buffer => already lex via
  // std::map; sweep workspace => first-seen via FragmentIterator). This
  // makes the bdg outputs byte-identical between paths on the same input.
  // No-op if already sorted; a few-microsecond pass otherwise.
  {
    const size_t n = workspace->chrom_names.size();
    if (n > 1) {
      bool already_sorted = true;
      for (size_t i = 1; i < n; ++i) {
        if (workspace->chrom_names[i - 1] > workspace->chrom_names[i]) {
          already_sorted = false;
          break;
        }
      }
      if (!already_sorted) {
        std::vector<size_t> order(n);
        for (size_t i = 0; i < n; ++i) order[i] = i;
        std::sort(order.begin(), order.end(),
                  [&](size_t a, size_t b) {
                    return workspace->chrom_names[a] < workspace->chrom_names[b];
                  });
        std::vector<std::string> new_names(n);
        std::vector<Macs3BdgTrack> new_treat(n), new_lambda(n), new_ppois(n);
        std::vector<Macs3FragEventBuffer> new_events(n);
        std::vector<std::vector<int32_t>> new_pos(n);
        std::vector<std::vector<int32_t>> new_cap_s(n);
        std::vector<std::vector<int32_t>> new_cap_e(n);
        const bool have_treat = workspace->treat_tracks.size() == n;
        const bool have_lambda = workspace->lambda_tracks.size() == n;
        const bool have_ppois = workspace->ppois_tracks.size() == n;
        const bool have_events = workspace->events_by_chrom.size() == n;
        const bool have_pos = workspace->pos_array_positions.size() == n;
        const bool have_cap_s = workspace->pos_array_starts_capture.size() == n;
        const bool have_cap_e = workspace->pos_array_ends_capture.size() == n;
        for (size_t i = 0; i < n; ++i) {
          new_names[i] = std::move(workspace->chrom_names[order[i]]);
          if (have_treat) new_treat[i] = std::move(workspace->treat_tracks[order[i]]);
          if (have_lambda) new_lambda[i] = std::move(workspace->lambda_tracks[order[i]]);
          if (have_ppois) new_ppois[i] = std::move(workspace->ppois_tracks[order[i]]);
          if (have_events) new_events[i] = std::move(workspace->events_by_chrom[order[i]]);
          if (have_pos) new_pos[i] = std::move(workspace->pos_array_positions[order[i]]);
          if (have_cap_s) new_cap_s[i] = std::move(workspace->pos_array_starts_capture[order[i]]);
          if (have_cap_e) new_cap_e[i] = std::move(workspace->pos_array_ends_capture[order[i]]);
        }
        workspace->chrom_names = std::move(new_names);
        if (have_treat) workspace->treat_tracks = std::move(new_treat);
        if (have_lambda) workspace->lambda_tracks = std::move(new_lambda);
        if (have_ppois) workspace->ppois_tracks = std::move(new_ppois);
        if (have_events) workspace->events_by_chrom = std::move(new_events);
        if (have_pos) workspace->pos_array_positions = std::move(new_pos);
        if (have_cap_s) workspace->pos_array_starts_capture = std::move(new_cap_s);
        if (have_cap_e) workspace->pos_array_ends_capture = std::move(new_cap_e);
      }
    }
  }
  const int64_t n_frag_rows = workspace->total_fragment_count;
  StageProfileCollector::Tick tk0{};
  if (stage_profile) tk0 = StageProfileCollector::Now();
  // Build pos_array positions if requested (no-op for sweep workspace,
  // which already populated; events workspace consumes the capture
  // buffers here).
  if (workspace->params.byte_parity_summits && !params.no_lambda) {
    if (!FinalizeMacs3FragPosArrayPositions(workspace, peak_threads)) {
      cleanup_on_fail();
      return fail("FinalizeMacs3FragPosArrayPositions failed");
    }
    if (stage_profile) {
      const auto tk1 = StageProfileCollector::Now();
      stage_profile->Record("finalize_pos_array_positions", tk0, tk1, n_frag_rows,
                            static_cast<int64_t>(workspace->pos_array_positions.size()), 0,
                            "byte_parity_summits");
      tk0 = tk1;
    }
  }
  if (!FinalizeMacs3FragTreatTracks(workspace, peak_threads)) {
    cleanup_on_fail();
    return fail("FinalizeMacs3FragTreatTracks failed");
  }
  if (stage_profile) {
    const auto tk1 = StageProfileCollector::Now();
    stage_profile->Record("finalize_treat_tracks", tk0, tk1, n_frag_rows,
                          static_cast<int64_t>(workspace->treat_tracks.size()), 0,
                          "from_workspace");
    tk0 = tk1;
  }
  const bool lambda_ok = params.no_lambda
                             ? FinalizeMacs3FragGlobalLambdaTracks(
                                   workspace, peak_threads)
                             : FinalizeMacs3FragLambdaTracks(workspace,
                                                             peak_threads);
  if (!lambda_ok) {
    cleanup_on_fail();
    return fail("FinalizeMacs3FragLambdaTracks failed");
  }
  if (stage_profile) {
    const auto tk1 = StageProfileCollector::Now();
    stage_profile->Record("finalize_lambda_tracks", tk0, tk1, n_frag_rows,
                          static_cast<int64_t>(workspace->lambda_tracks.size()), 0,
                          "from_workspace");
    tk0 = tk1;
  }
  if (!treat.empty()) {
    if (!WriteMacs3BdgTracksBedGraph(treat, workspace->chrom_names,
                                     workspace->treat_tracks)) {
      cleanup_on_fail();
      return fail("WriteFragSpanPileupBedGraph failed");
    }
    if (stage_profile) {
      const auto tk1 = StageProfileCollector::Now();
      stage_profile->Record("export_treat_pileup", tk0, tk1, n_frag_rows,
                            static_cast<int64_t>(workspace->treat_tracks.size()), 0,
                            "");
      tk0 = tk1;
    }
  }
  if (!lambda.empty()) {
    if (!WriteMacs3BdgTracksBedGraph(lambda, workspace->chrom_names,
                                     workspace->lambda_tracks)) {
      cleanup_on_fail();
      return fail("WriteFragMacs3NoControlLambdaBedGraph failed");
    }
    if (stage_profile) {
      const auto tk1 = StageProfileCollector::Now();
      stage_profile->Record("export_control_lambda", tk0, tk1, n_frag_rows,
                            static_cast<int64_t>(workspace->lambda_tracks.size()), 0,
                            "");
      tk0 = tk1;
    }
  }
  if (export_ppois) {
    if (!FinalizeMacs3FragPpoisTracks(
            workspace, peak_threads, !params.macs3_float_score_math)) {
      cleanup_on_fail();
      return fail("FinalizeMacs3FragPpoisTracks failed");
    }
    if (stage_profile) {
      const auto tk1 = StageProfileCollector::Now();
      stage_profile->Record("finalize_ppois_tracks", tk0, tk1, n_frag_rows,
                            static_cast<int64_t>(workspace->ppois_tracks.size()), 0,
                            "from_workspace");
      tk0 = tk1;
    }
    if (!WriteMacs3BdgTracksBedGraph(ppois, workspace->chrom_names,
                                     workspace->ppois_tracks)) {
      cleanup_on_fail();
      return fail("WriteMacs3BdgTracksBedGraph ppois failed");
    }
    if (stage_profile) {
      const auto tk1 = StageProfileCollector::Now();
      stage_profile->Record("export_score_ppois", tk0, tk1, n_frag_rows,
                            static_cast<int64_t>(workspace->ppois_tracks.size()), 0,
                            "");
      tk0 = tk1;
    }
    if (!RunMacs3FragPpoisNarrowPeaksFromTracks(
            workspace->chrom_names, workspace->ppois_tracks, workspace->treat_tracks,
            workspace->lambda_tracks, threshold_cutoff, params.min_length,
            params.max_gap, params.score_pseudocount, narrowpeak_out, summits_out,
            stage_profile, n_frag_rows, params.name_prefix,
            qvalue_mode, params.macs3_float_score_math,
            workspace->params.byte_parity_summits
                ? &workspace->pos_array_positions
                : nullptr,
            params.call_summits)) {
      cleanup_on_fail();
      return fail("RunMacs3FragPpoisNarrowPeaks failed");
    }
  } else {
    const size_t n_chrom = workspace->chrom_names.size();
    std::vector<std::vector<Macs3FragNarrowPeakRow>> rows_by_chrom(n_chrom);
    std::vector<Macs3BdgTrack> ppois_tracks_all(n_chrom);
    std::vector<char> ok(n_chrom, 1);
    const int threads = NormalizeThreads(peak_threads, n_chrom);
    StageProfileCollector::Tick t_ppois_start{};
    if (stage_profile) t_ppois_start = StageProfileCollector::Now();
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) if(threads > 1)
    for (int64_t ii = 0; ii < static_cast<int64_t>(n_chrom); ++ii) {
      const size_t i = static_cast<size_t>(ii);
      if (!BuildMacs3FragPpoisTrack(workspace->treat_tracks[i],
                                    workspace->lambda_tracks[i],
                                    params.score_pseudocount,
                                    &ppois_tracks_all[i],
                                    !params.macs3_float_score_math)) {
        ok[i] = 0;
      }
    }
    if (std::find(ok.begin(), ok.end(), 0) != ok.end()) {
      cleanup_on_fail();
      return fail("chromosome-at-a-time ppois construction failed");
    }
    StageProfileCollector::Tick t_ppois_end{};
    if (stage_profile) {
      t_ppois_end = StageProfileCollector::Now();
      stage_profile->Record("score_build_ppois", t_ppois_start, t_ppois_end,
                            n_frag_rows, static_cast<int64_t>(n_chrom), 0,
                            "chromosome_parallel");
    }
    std::map<double, double> pqtable;
    BuildMacs3PqTable(ppois_tracks_all, &pqtable,
                      params.macs3_float_score_math);
    StageProfileCollector::Tick t_pqtable_end{};
    if (stage_profile) {
      t_pqtable_end = StageProfileCollector::Now();
      stage_profile->Record("score_build_pqtable", t_ppois_end, t_pqtable_end,
                            n_frag_rows, static_cast<int64_t>(pqtable.size()), 0,
                            "global_serial");
    }
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) if(threads > 1)
    for (int64_t ii = 0; ii < static_cast<int64_t>(n_chrom); ++ii) {
      const size_t i = static_cast<size_t>(ii);
      const std::vector<int32_t>* pap =
          (workspace->params.byte_parity_summits &&
           i < workspace->pos_array_positions.size())
              ? &workspace->pos_array_positions[i]
              : nullptr;
      if (!CollectMacs3FragNarrowPeakRowsFromTracks(
              workspace->chrom_names[i], ppois_tracks_all[i],
              workspace->treat_tracks[i],
              workspace->lambda_tracks[i], threshold_cutoff, params.min_length,
              params.max_gap, params.score_pseudocount, &rows_by_chrom[i], pap,
              params.macs3_float_score_math,
              qvalue_mode ? &pqtable : nullptr,
              params.call_summits)) {
        ok[i] = 0;
        continue;
      }
      if (i < workspace->pos_array_positions.size()) {
        std::vector<int32_t>().swap(workspace->pos_array_positions[i]);
      }
      if (!export_treat) {
        std::vector<int32_t>().swap(workspace->treat_tracks[i].p);
        std::vector<float>().swap(workspace->treat_tracks[i].v);
      }
      if (!export_lambda) {
        std::vector<int32_t>().swap(workspace->lambda_tracks[i].p);
        std::vector<float>().swap(workspace->lambda_tracks[i].v);
      }
    }
    if (std::find(ok.begin(), ok.end(), 0) != ok.end()) {
      cleanup_on_fail();
      return fail("chromosome-at-a-time ppois peak collection failed");
    }
    if (stage_profile) {
      stage_profile->Record("score_collect_regions", t_pqtable_end,
                            StageProfileCollector::Now(), n_frag_rows,
                            static_cast<int64_t>(n_chrom), 0,
                            "chromosome_parallel");
    }
    std::vector<Macs3FragNarrowPeakRow> rows;
    size_t row_count = 0;
    for (const auto& chrom_rows : rows_by_chrom) {
      row_count += chrom_rows.size();
    }
    rows.reserve(row_count);
    for (auto& chrom_rows : rows_by_chrom) {
      rows.insert(rows.end(), chrom_rows.begin(), chrom_rows.end());
    }
    if (stage_profile) {
      const auto tk1 = StageProfileCollector::Now();
      stage_profile->Record("bdgpeakcall_regions", tk0, tk1, n_frag_rows,
                            static_cast<int64_t>(rows.size()), 0,
                            "from_workspace;in_memory_tracks");
      tk0 = tk1;
    }
    if (!WriteMacs3FragNarrowPeakRows(&rows, narrowpeak_out, summits_out,
                                      &pqtable, params.name_prefix,
                                      qvalue_mode)) {
      cleanup_on_fail();
      return fail("WriteMacs3FragNarrowPeakRows failed");
    }
    if (stage_profile) {
      const auto tk1 = StageProfileCollector::Now();
      stage_profile->Record("narrowpeak_summits", tk0, tk1, n_frag_rows,
                            static_cast<int64_t>(rows.size()), 0, "from_workspace");
      tk0 = tk1;
    }
  }
  if (own_workdir) {
    RemoveTree(work);
  }
  return true;
}

}  // namespace peaks
}  // namespace chromap
