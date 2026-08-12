# Library API

RapidMACS builds to `lib/librapidmacs.a`, a static archive whose only
dependencies are htslib and zlib. There is no Python or Cython runtime, so a
host program — an aligner, a pipeline binary — can call peaks in the same
process on fragments it already holds, instead of writing a temporary file and
starting `macs3`.

```sh
g++ -std=c++11 -Iinclude my_tool.cc lib/librapidmacs.a -lhts -lm -lz -lpthread -fopenmp
```

Public headers live in `include/rapidmacs/`. Two namespaces are in play:
`macs3` for the fragment-source types, `chromap::peaks` for the pipeline.

## The fragment source

Everything enters through one interface, in
[`rapidmacs/fragments.h`](../include/rapidmacs/fragments.h):

```cpp
namespace macs3 {

struct FragmentRecord {
  int32_t chrom_id = -1;
  int32_t start = 0;
  int32_t end = 0;
  uint32_t count = 1;   // duplicate multiplicity
};

class FragmentIterator {
 public:
  virtual ~FragmentIterator() = default;

  // Fill *out with the next fragment. false at end of stream or on error.
  virtual bool Next(FragmentRecord* out) = 0;

  // Stable id -> name map; indices match FragmentRecord::chrom_id.
  virtual const std::vector<std::string>& ChromNames() const = 0;

  // Empty unless Next() returned false because of a parse or I/O failure,
  // as opposed to reaching the end of the stream.
  virtual const std::string& Error() const = 0;
};

}  // namespace macs3
```

`Next()` returning `false` means either end-of-stream or failure; check
`Error()` to tell them apart. Chromosome names are interned in first-seen
order and the table only grows, so ids stay valid for the life of the iterator.

Two implementations ship:

```cpp
// Reads a fragments TSV, plain or gzipped, 3- to 7-column schemas.
std::unique_ptr<FragmentIterator> OpenFileFragmentIterator(
    const std::string& path);

// Adapts an in-memory vector. Construction sorts by (chrom_id, start, end),
// so the result satisfies the sorted-input contract below. Takes ownership.
// chrom_names must cover every chrom_id present.
std::unique_ptr<FragmentIterator> WrapVectorFragmentIterator(
    std::vector<FragmentRecord> fragments,
    std::vector<std::string> chrom_names);
```

Implementing the interface yourself is the point of it: an aligner emits
`FragmentRecord`s as it maps, and peaks are called without the data ever
reaching disk.

## Parameters

`chromap::peaks::Macs3FragPeakPipelineParams`, from
[`rapidmacs/macs3_frag_peak_pipeline.h`](../include/rapidmacs/macs3_frag_peak_pipeline.h).
The defaults reproduce MACS3 v3.0.3 for fragment input; change only what you
mean to.

| Field | Default | Meaning |
|---|---|---|
| `threshold_mode` | `kPValue` | p-value or q-value thresholding |
| `bdgpeakcall_cutoff` | `5.f` | `-log10` cutoff; 5 is `-p 1e-5` |
| `qvalue_cutoff` | `0.0` | used when `threshold_mode` is the q-value form |
| `min_length` | `200` | minimum peak length in bp |
| `max_gap` | `30` | maximum gap to link in bp |
| `effective_genome_size` | `2913022398` | hs GRCh38 |
| `slocal_bp` / `llocal_bp` | `1000` / `10000` | local background windows |
| `keep_duplicates_all` | `false` | selects MACS3's parser-time float division for PE `--keep-dup all` control-depth scaling |
| `scale_to_large` | `false` | depth target |
| `no_lambda` | `false` | genome-wide background only |
| `call_summits` | `false` | MACS3 subpeak algorithm |
| `peak_caller_threads` | `1` | chromosome-level concurrency |
| `name_prefix` | empty | peaks named `<prefix>_peak_<N>`; otherwise `peak_<N>` |
| `byte_parity_summits` | `true` | reproduce MACS3's per-event `pos_array` summit tie-break |
| `macs3_float_score_math` | `true` | MACS3's float32 p/q-table arithmetic |
| `macs3_uint8_counts` | `true` | coerce fragment weights through uint8, as PETrackII does |
| `frag_pileup_max_count` | `0` | cap the count column; 0 disables |
| `score_pseudocount` | `0.0` | `bdgcmp -p` |
| `treatment_fragment_length_bp` | `0.0` | override the derived mean fragment length |

Leave `byte_parity_summits`, `macs3_float_score_math`, and
`macs3_uint8_counts` alone unless you have a specific reason: together they are
what makes the output byte-identical.

`Macs3FragPeakPipelinePaths` optionally names the intermediate bedGraphs
(`treat_bdg`, `lambda_bdg`, `ppois_bdg`). Leave the fields empty to have them
written into a working directory and cleaned up.

## Entry points

Four, differing in where fragments come from and how much memory that costs.

### `RunMacs3FragPeakPipelineFromSortedIterator` — the one to use

```cpp
bool RunMacs3FragPeakPipelineFromSortedIterator(
    macs3::FragmentIterator& iter,
    const Macs3FragPeakPipelineParams& params,
    const Macs3FragPeakPipelinePaths& user_paths,
    const std::string& narrowpeak_out, const std::string& summits_out,
    const std::string& keep_intermediates_dir,
    const std::string& work_dir_parent, std::string* work_dir_used,
    std::string* error_message,
    StageProfileCollector* stage_profile = nullptr);
```

Streams straight into the sweep-line workspace, building the treatment and
background tracks incrementally. Peak memory is bounded by the fragments
spanning the current position — tens of kilobytes — rather than by the size of
the input.

**It requires chrom-grouped, start-sorted input**: `chrom_id` non-decreasing
across calls, and `start` non-decreasing within each `chrom_id`. Fragments
files from Cell Ranger and Chromap already satisfy this. On violation it
returns `false` with a descriptive `error_message` rather than producing wrong
output.

### `RunMacs3FragPeakPipelineFromIterator`

Same signature, no ordering requirement. Buffers events at roughly 48 bytes per
fragment, so a 70 M-fragment dataset costs gigabytes where the sorted entry
point costs kilobytes. Use it only when you genuinely cannot present sorted
input.

### `RunMacs3FragPeakPipelineFromFragments`

Takes `std::vector<ChromFragments>*` that you have already grouped by
chromosome, plus an optional `control_by_chrom`. This is the entry point for
matched-control workflows.

```cpp
bool RunMacs3FragPeakPipelineFromFragments(
    std::vector<ChromFragments>* by_chrom,
    const Macs3FragPeakPipelineParams& params,
    const Macs3FragPeakPipelinePaths& user_paths,
    const std::string& narrowpeak_out, const std::string& summits_out,
    const std::string& keep_intermediates_dir,
    const std::string& work_dir_parent, std::string* work_dir_used,
    std::string* error_message,
    StageProfileCollector* stage_profile = nullptr,
    bool release_fragments_after_workspace = true,
    const std::vector<ChromFragments>* control_by_chrom = nullptr);
```

`release_fragments_after_workspace` frees the input once the workspace is
built, which matters at scale.

To fill those vectors from files, [`rapidmacs/fragment_input.h`](../include/rapidmacs/fragment_input.h)
provides `LoadFragmentsFromTsv`, `LoadFragmentsFromBampe`, and
`LoadFragmentsFromBedpe`, each reporting a `FragmentInputStats`.

### `RunMacs3FragPeakPipelineFromWorkspace`

Takes an already-built `Macs3FragPeakWorkspace`. For callers driving workspace
construction themselves.

## Return and error handling

All four return `bool`. On `false`, `*error_message` explains why. Nothing
throws. `work_dir_used` receives the working directory actually used, which
matters when `keep_intermediates_dir` is empty and a temporary one was chosen.

Pass a `StageProfileCollector` (see
[`rapidmacs/stage_profile.h`](../include/rapidmacs/stage_profile.h)) to get the
per-stage timings the paper reports.

## A minimal caller

```cpp
#include "rapidmacs/fragments.h"
#include "rapidmacs/macs3_frag_peak_pipeline.h"

using chromap::peaks::Macs3FragPeakPipelineParams;
using chromap::peaks::Macs3FragPeakPipelinePaths;

int main() {
  auto iter = macs3::OpenFileFragmentIterator("fragments.tsv.gz");
  if (!iter) return 1;

  Macs3FragPeakPipelineParams params;
  params.name_prefix = "atac";
  params.peak_caller_threads = 8;
  // Defaults are already MACS3 -p 1e-5 --min-length 200 --max-gap 30.

  std::string work_dir, error;
  const bool ok = chromap::peaks::RunMacs3FragPeakPipelineFromSortedIterator(
      *iter, params, Macs3FragPeakPipelinePaths{},
      "atac.narrowPeak", "atac.summits.bed",
      /*keep_intermediates_dir=*/"", /*work_dir_parent=*/"",
      &work_dir, &error);
  if (!ok) { fprintf(stderr, "rapidmacs: %s\n", error.c_str()); return 1; }
  return 0;
}
```

The narrowPeak this writes is byte-identical to what `macs3 callpeak -t
fragments.tsv.gz -f FRAG -g hs -n atac -p 1e-5 --min-length 200 --max-gap 30`
would have written.

## Other headers

| Header | Contents |
|---|---|
| `macs3_bed_callpeak.h` | Single-end BED callpeak, and `ResolveMacs3BedCallPeakProfile("signac-atac", ...)`. Use the resolver rather than copying profile constants. |
| `signac_bed_projection.h` | Projecting a fragments file to BED, optionally per barcode group |
| `peak_io.h` | narrowPeak and summit readers and writers |
| `macs3_bdg_track.h` | The run-length bedGraph track type |
| `macs3_pos_array.h` | MACS3's per-event `pos_array` summit tie-break |
| `bdgpeakcall.h` | Peak calling on an existing score track; `BdgPeakCallCutoffFromPValue` |
| `streaming.h`, `binned_signal.h`, `call_peaks.h` | The legacy binned caller, retained for the Chromap-compatible path |
