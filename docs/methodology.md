# RapidMACS technical methodology

This document describes the implemented paired-end narrow-peak caller, its
relationship to MACS3 3.0.3, and the validation and benchmarking protocol used
for the ATAC-seq, ChIP-seq, and CUT&RUN study. It is intended to be the
technical source for the paper Methods section and for reproducing the command
line results. Historical development plans are not normative; the source code
and the procedures below describe the current implementation.

## 1. Scope and reproducibility target

RapidMACS implements the MACS3 paired-end narrow-peak path for:

- Chromap/10x-style fragment files (`FRAG`);
- paired-end BAM (`BAMPE`); and
- MACS3 three-column or conventional six-column paired-end BED (`BEDPE`).

Treatment-only and treatment-plus-control analyses are supported. A control is
an input library for ChIP-seq or, typically, IgG for CUT&RUN. Version 1 is
deliberately limited to paired-end fragments and narrow peaks. Single-end model
building, cross-correlation-based fragment-length inference, and broad-peak
calling are outside the implemented scope.

The reference behavior is `macs3 callpeak` version 3.0.3 with the same
post-alignment input, effective genome size, statistical cutoff, local-window
settings, duplicate policy, and scaling direction. Validation is performed on
the final CLI products, not only on internal library objects.

RapidMACS is an independent implementation of that behavior, not a port of
MACS3's Cython sources. It was written clean-room from the published
algorithmic specification and validated by output comparison; MACS3's source
was consulted only for the low-level details required to reach byte-for-byte
equivalence, chiefly the summit tie-break of §8. This distinction matters for
interpreting the performance results: because the codebase is independent, it
selects the representations and control flow described in §§4-9, and the
single-threaded advantage follows from those choices rather than from
compiling an equivalent design in a faster language. See `NOTICE`.

## 2. End-to-end computation

The production path is:

```text
FRAG, BAMPE, or BEDPE
        |
        v
sorted [start,end) fragments grouped by chromosome
        |
        +--> treatment fragment-span pileup T(x)
        |
        +--> treatment-derived or matched-control local lambda L(x)
        |
        v
Poisson score S(x) = -log10 P[X >= T(x) | L(x)]
        |
        +--> genome-wide, base-pair-weighted p-to-q table
        |
        v
threshold + gap merge + minimum-length filter
        |
        v
MACS3-compatible summit selection and narrowPeak/summits output
```

Tracks are piecewise constant. A `Macs3BdgTrack` stores parallel endpoint and
value arrays: value `v[i]` applies from the preceding endpoint through `p[i]`.
The implementation therefore scales with the number of signal changes rather
than allocating a value for every genomic base.

Coordinates are signed 32-bit integers and input fragments use zero-based,
half-open `[start,end)` intervals. Counts and intermediate totals use wider
integer types where accumulation can exceed a chromosome coordinate.

## 3. Input normalization

### 3.1 FRAG

The FRAG reader accepts plain or gzip-compressed records of the form:

```text
chrom  start  end  barcode  count
```

The barcode is ignored and `count` is the fragment multiplicity. The optimized
reader uses 1 MiB `gzread` blocks, manual checked integer conversion, and a
last-chromosome lookup fast path. Malformed coordinates, non-positive counts,
and lines longer than 65,535 bytes are rejected by the strict iterator.

The MACS3-compatible CLI path requires chromosomes to be grouped and fragments
to be non-decreasing by start within a chromosome. It validates this contract,
groups fragments by chromosome, constructs tracks in parallel, and releases
each chromosome's fragment vector after its tracks have been produced. If the
ordering contract is violated, the CLI reopens the file and falls back to the
general event-buffer implementation.

### 3.2 BAMPE

BAMPE is read with htslib. The parser mirrors MACS3 3.0.3
`bampe_pe_binary_parse`:

- retain one record per template by accepting read 1 and rejecting read 2;
- require paired and proper-pair flags;
- reject unmapped, mate-unmapped, secondary, supplementary, and QC-fail
  records;
- retain duplicate-marked alignments; and
- do not impose a MAPQ threshold.

The fragment is the leftmost mate coordinate plus the absolute template length.
After parsing, fragments are sorted by `(chromosome,start,end)` and identical
coordinate intervals are collapsed to reproduce MACS3's default
`--keep-dup 1`. The treatment mean fragment length is calculated before this
coordinate deduplication because MACS3 uses the parser-level mean for its
paired-end `d` window and automatic peak geometry.

htslib BGZF workers are enabled when more than one peak-caller thread is
requested. Per-chromosome sorting and coordinate deduplication are also
parallelized.

### 3.3 BEDPE

Two text representations are accepted:

- MACS3 BEDPE: `chrom,start,end`; and
- conventional BEDPE: two mate intervals in six or more columns.

For six-column input, both mates must be on the same chromosome. The emitted
fragment spans the minimum mate start through the maximum mate end. BEDPE uses
the same lexicographic chromosome ordering, `(start,end)` sorting, and
`--keep-dup 1` coordinate deduplication as BAMPE.

## 4. Treatment pileup

For fragment `j`, with interval `[a_j,b_j)` and effective weight `w_j`, the
treatment pileup is

```text
T(x) = sum_j w_j I(a_j <= x < b_j).
```

Sorted start and end streams are merged directly. Equal-coordinate events are
accumulated before the running pileup changes, and adjacent runs with equal
values are coalesced. This avoids materializing and comparison-sorting two
events for every endpoint when the input is already sorted.

FRAG multiplicities follow MACS3 `PETrackII` behavior, including the optional
maximum-count cap and the compatibility conversion through `uint8_t`. BAMPE
and BEDPE fragments have weight one after default coordinate deduplication.

## 5. Local background lambda

### 5.1 Treatment-only FRAG

Treatment-only analysis follows MACS3's paired-end no-control path. The large
local window defaults to `llocal = 10,000` bp. Endpoint-window events are
constructed from treatment fragments and scaled by the treatment fragment
statistics. The genome-wide floor is

```text
lambda_bg = total treatment fragment bp / effective genome size.
```

At each position, local lambda is the scaled large-window estimate bounded
below by `lambda_bg`. The final lambda track is aligned to the union of
treatment and local-background breakpoints.

### 5.2 Matched-control BAMPE/BEDPE

Controlled analysis follows MACS3's paired-end
`__call_peaks_w_control` calculation. Three endpoint-window pileups are built
from the control sample:

- `d`: the pre-dedup treatment mean fragment length;
- `slocal`: 1,000 bp by default; and
- `llocal`: 10,000 bp by default.

The `slocal` and `llocal` estimates are normalized by `d/window_size`. Control
and treatment depths are normalized using MACS3's paired-end arithmetic. The
default is scale-to-small: the larger library is scaled toward the smaller
library. `--scale-to large` and `--to-large` select the opposite direction.

For every genomic run,

```text
L(x) = max(lambda_bg, lambda_d(x), lambda_slocal(x), lambda_llocal(x)).
```

The four already-sorted endpoint streams for each window are merged directly.
Temporary controlled-lambda tracks use an eight-byte array-of-structures run
representation (`int32_t end`, `float value`). Only the final maximum lambda is
converted into the public track representation. Treatment-only lambda event
buffers and redundant raw endpoint captures are not allocated in controlled
mode.

## 6. Poisson scores and chromosome-local memoization

Treatment and lambda tracks are traversed together over their breakpoint
union. Each merged run is assigned MACS3's Poisson upper-tail score,
`-log10(p)`, using the htslib mathematical routines already linked by the
library.

Many runs reuse the same treatment and lambda values. Recalculating the
Poisson tail for every run is unnecessarily expensive, so each chromosome
calculation creates a private khash table:

- the key is one 64-bit integer formed from the exact 32-bit floating-point
  bit patterns of `(treatment,lambda)`;
- the value is the resulting score as a `double`;
- the flat table is initially sized to one bucket per eight merged input runs;
- there are no per-entry node allocations, locks, atomics, or cross-thread
  cache sharing; and
- the table is destroyed when that chromosome finishes.

On the full ATAC benchmark, 101,027,560 of 113,684,789 scored runs were cache
hits (88.87%). Diagnostic counters are selected through
`RAPIDMACS_PPOIS_CACHE_STATS` (the former
`LIBMACS3_PPOIS_CACHE_STATS` name remains accepted for compatibility); the
counter-free template specialization is used
for normal runs so instrumentation does not remain in the hot loop.

All callpeak-compatible FRAG, BAMPE, and BEDPE paths use the float assignments
and lookup keys made by MACS3's `CallPeakUnit`, including conversions at
rounding boundaries. Local lambda remains float32 in memory for Poisson
scoring; five-decimal conversion is confined to diagnostic bedGraph output.
The lower-level API retains an explicit legacy round-five/double mode for
callers that need old RapidMACS behavior.

## 7. Genome-wide q values and peak construction

The p-to-q conversion is genome-wide, not per chromosome. The implementation
first accumulates the base-pair coverage of every unique p-score, sorts the
unique scores from highest to lowest, and applies MACS3's base-pair-weighted
Benjamini-Hochberg calculation:

```text
-log10(q) = -log10(p) + log10(rank) - log10(total covered bp),
```

followed by MACS3's monotonic adjustment. The eager hash reservation is capped
because the number of signal runs is much larger than the number of unique
scores. The accumulator is released before the ordered lookup table is built.
A duplicate full-genome q-score track is not retained; peak collection uses
the p-to-q lookup table directly when q-value thresholding is requested.

Above-threshold runs are merged when their separation does not exceed
`max-gap`, then filtered by `min-length`. The paired-end presets use the
pre-dedup mean treatment fragment length for MACS3's dynamic minimum length and
maximum gap. ATAC paper runs explicitly use 200 bp and 30 bp, respectively.

For each peak, signal value, p-score, and q-score are evaluated with MACS3's
paired-end float behavior. Output consists of narrowPeak and one BED summit
record per peak.

## 8. Summit parity

MACS3 selects the median member of a tied maximum after subdividing a peak at
positions in its internal `pos_array`. Exported bedGraph tracks can collapse
adjacent equal-valued runs and therefore do not contain every one of these
boundaries.

Default parity mode reconstructs `pos_array` per chromosome from
count-duplicated treatment endpoints and the relevant local/control endpoint
windows, takes their sorted union, and splits tied treatment plateaus at every
interior position before applying MACS3's median tie rule. This mode produced
zero summit-coordinate differences across 50,274 peaks in the full PBMC ATAC
validation and byte-identical controlled-assay summit files.

`--no-macs3-frag-summit-parity` disables this reconstruction and saves memory,
but it is not used for parity or paper benchmarks.

## 9. Parallel execution and scaling limit

OpenMP uses dynamic scheduling with one work item per chromosome for:

- BAM per-chromosome sorting and deduplication;
- treatment and local-lambda track construction;
- summit-position reconstruction;
- Poisson-score construction; and
- peak-region collection.

The khash score cache belongs to the chromosome work item executing on a
thread. It is not a persistent per-worker or process-global cache.

FRAG decompression/parsing and the genome-wide p-to-q table are serial. BAM
decompression has independent htslib workers. Consequently, chromosome-level
work is effectively bounded by the roughly 24 large human chromosomes; small
alternate contigs add tasks but little useful load balance. At high thread
counts, the final large chromosome and the serial parsing/q-table stages form
the scaling floor. Coordinate-region sharding was evaluated but is not used:
duplicated cache/output setup and repeated score work made it slower than
chromosome scheduling.

## 10. Standalone executable

Build the static library and CLI with:

```bash
git submodule update --init --recursive
make -j
```

The executable is `bin/rapidmacs`. The optimized code paths require no hidden
or experimental flag; thread count is selected with
`--peak-caller-threads`.

### ATAC-seq

For FRAG input, `--macs3-frag-narrowpeak` selects the MACS3-compatible path.
The older FRAG `--out-prefix`-only interface selects a different legacy binned
caller and must not be used for the paper comparison.

```bash
bin/rapidmacs --preset shiftTAG \
  -f FRAG -t fragments.tsv.gz -g hs -n atac_full \
  --peak-caller-threads 24 \
  --macs3-frag-narrowpeak out/atac_full.narrowPeak \
  --macs3-frag-summits out/atac_full.summits.bed \
  -p 1e-5 --bdgpeakcall-min-len 200 --bdgpeakcall-max-gap 30
```

### ChIP-seq

```bash
bin/rapidmacs --preset chip \
  -t treatment.bam -c input_control.bam -g hs -n chipseq_full \
  --peak-caller-threads 24 --out-prefix out/chipseq_full
```

### CUT&RUN

```bash
bin/rapidmacs --preset cutrun \
  -t treatment.bam -c igg_control.bam -g hs -n cutrun_full \
  --peak-caller-threads 24 --out-prefix out/cutrun_full
```

The ChIP-seq and CUT&RUN presets select BAMPE, `-q 0.05`, MACS3 float
score behavior, and dynamic paired-end peak geometry. CUT&RUN may also be run
without a control, in which case the treatment-only background path is used.

## 11. Validation datasets

### 11.1 ATAC-seq

The full ATAC dataset is the 10x Multiome PBMC 3k GRCh38 analysis used by the
existing RapidMACS work. Its post-Chromap input contains 53,969,811 fragment
rows:

```text
data
  chromap_full_dual_20260424_194803/fragments.tsv.gz
```

The committed fast regression fixture contains 100,000 input fragments. Full
validation produces 50,274 peaks and BED3 checksum
`b711b2d62ae08620d480df1b964f98e0`.

### 11.2 ChIP-seq

The full public dataset is ENCODE CTCF ChIP-seq in NCI-H929 with its declared
input-library control:

| Role | Experiment | Alignment | GRCh38 BAM MD5 |
|---|---|---|---|
| Treatment | ENCSR634OAQ | ENCFF686KKV | `bedc2707c0be9f6d61c09edf811d655e` |
| Control | ENCSR409FGC | ENCFF181CXT | `ffec8961e15783bafb161fb7ab534568` |

The smoke fixture is `chr21:20,000,000-30,000,000`. The full comparison
produces 33,985 peaks with BED3 checksum
`129a2a508b271f8f34333147c7de218e` and summit checksum
`82ae9573bb99733e30635cef96da1f69`.

### 11.3 CUT&RUN

The full public dataset is K562 CTCF CUT&RUN with its matched IgG control from
GEO series GSE172236:

| Role | GEO/SRA identifiers | Description |
|---|---|---|
| Treatment | GSM5244531 / SRX10617507 / SRR14255054 | K562 sgLUC CTCF CUT&RUN, 150-bp PE |
| Control | GSM5244530 / SRX10617506 / SRR14255053 | K562 IgG, 150-bp PE |

Both samples are aligned identically with Chromap using the commands recorded
in `tests/fixtures/README.md`. The smoke fixture is
`chr21:20,000,000-30,000,000`. The full comparison produces 9,097 peaks with
BED3 checksum `c7cf0e3b6e46202460a717473bcfec2e` and summit checksum
`1cd24d0c399a5ea0b131f8c6e5fe9c7e`.

Download URLs, input checksums, alignment options, subset commands, and golden
generation commands are maintained in `tests/fixtures/README.md`.

## 12. Correctness validation

Small regression tests are run with:

```bash
make test
```

This executes C++ unit tests, the legacy FRAG regression, the MACS3-compatible
FRAG regression, ChIP-seq treatment/control parity, and CUT&RUN
treatment/control parity. The paired fixtures are compared with goldens made
by stock MACS3 3.0.3.

For every measured execution, the full benchmark records:

1. peak count;
2. full narrowPeak MD5;
3. BED3-only MD5;
4. full summits MD5; and
5. summit-coordinate-only MD5.

For every assay, the benchmark aborts unless the complete narrowPeak and
summits files match the corresponding MACS3 reference. The final full-data
ATAC comparison is byte-identical across all 50,274 records, including names,
intervals, summit offsets and coordinates, signal values, p-scores, q-scores,
and the `-p`-selected summit BED score. Optional bedGraph diagnostics are used
during implementation, but paper performance is measured through complete
`callpeak`-style CLI runs.

## 13. Paper benchmark methodology

The benchmark unit is the standalone executable, including input parsing,
peak calling, and narrowPeak/summits output. Internal library-only timings are
diagnostic and are not used for the primary paper comparison.

### 13.1 Compared commands

MACS3 receives the same post-alignment input and assay parameters:

```bash
# ATAC
macs3 callpeak -t fragments.tsv.gz -f FRAG -g hs -n atac_full \
  -p 1e-5 --min-length 200 --max-gap 30 --outdir OUT

# ChIP-seq or CUT&RUN
macs3 callpeak -t treatment.bam -c control.bam -f BAMPE -g hs \
  -n NAME -q 0.05 --outdir OUT
```

RapidMACS uses the commands in Section 10 at 1, 4, and 24 threads. Alignment is
excluded: both tools consume the exact same FRAG or BAM files.

### 13.2 Repetition and cache policy

For every assay/tool/thread configuration:

1. run one unrecorded warm-up;
2. run three recorded repetitions serially; and
3. report the median wall time and median maximum RSS, retaining the observed
   minimum and maximum as the run-to-run range.

Warm-cache measurement is intentional. Dropping the Linux page cache requires
privileged, host-wide mutation and would disturb unrelated jobs. A full warm-up
for every configuration gives both implementations the same cache state and
avoids treating cold-cache variation as peak-caller performance.

GNU `/usr/bin/time` records elapsed wall time and maximum resident set size.
Each run writes into a separate directory. Standard output, standard error,
timing, outputs, checksums, and optional stage profiles are retained. Stage
profiling is enabled consistently for RapidMACS benchmark runs and its small
overhead is included in the reported RapidMACS wall time.

### 13.3 MACS3 timing reference

The reference version is pinned to MACS3 3.0.3. At the end of
`PETrackI.pileup_a_chromosome_c`, the release writes its internal control
pileup to `tmp_b_c.txt`, one line per breakpoint, and then returns that same
array to its caller. Nothing appears to read the file back and we found no
option that disables it, so it is produced on every run that has a control and
can write gigabytes. `paper/macs3-3.0.3-no-control-dump.patch` removes the
four lines, and `paper/scripts/setup_macs3_reference.sh` builds both the released
and the patched installation.

**The released build is what the paper reports.** Table 1 and Supplementary S3
give released-build times for both configurations with a control, and
Supplementary S4 places the two builds side by side, so the cost of the write
is disclosed rather than removed from the comparison. The three treatment-only
configurations were measured against the patched build, which is sound because
a run without a control never reaches that code; Supplementary S2 states this,
and scATAC 10x was measured both ways to bound the residual at 0.13%.

Released MACS3 3.0.3 remains the source of parity goldens. Both builds must be
verified to produce the same narrowPeak and summits files before either is
timed.

### 13.4 Reproducible harness

```bash
paper/scripts/setup_macs3_reference.sh
make -j
paper/scripts/run_paper_benchmarks.sh
```

`paper/scripts/run_paper_benchmarks.sh` writes:

- `raw.tsv`: every run, wall time, RSS, peak count, and four output hashes;
- `summary.tsv`: medians and ranges by assay/tool/thread count;
- `stages_raw.tsv`: diagnostic stage timings; and
- `stages_summary.tsv`: summarized diagnostic stage timings.

The primary paper table and figures must be generated from `summary.tsv`, not
from one-run optimization diagnostics. After a code change affecting the peak
path, the complete three-repeat matrix must be rerun. In particular, the
single-run khash measurements used during optimization are not final paper
measurements.

Before the final run, archive the source commit, compiler and flags, htslib
version, MACS3 environment, CPU model, RAM, kernel, benchmark input checksums,
and the complete raw result directory. No unrelated CPU- or I/O-intensive jobs
should share the host during measurement.

## 14. Published results

The six input configurations reported in Table 1 of the Applications Note.
Complete-program wall time in seconds, medians of three recorded warm-cache
runs after one unrecorded warm-up, on one 13th-generation Intel Core i9-13900KF
host with every input, output, and temporary file on the same NVMe device.
Parenthesised values are speedups against that configuration's own MACS3 time.
Within every row and every run, both programs produced byte-identical
narrowPeak and summit files.

| Configuration | MACS3 v3.0.3 | T=1 | T=4 | T=24 | Peaks |
|---|---:|---:|---:|---:|---:|
| scATAC-seq, 10x fragments | 747.55 | 57.87 (12.9x) | 24.13 (31.0x) | 16.37 (45.7x) | 50,003 |
| scATAC-seq, Chromap fragments | 758.62 | 58.88 (12.9x) | 25.38 (29.9x) | 16.12 (47.1x) | 50,274 |
| Bulk ATAC-seq, ENCODE BAM | 146.63 | 41.49 (3.5x) | 13.14 (11.2x) | 7.28 (20.1x) | 74,000 |
| Bulk ATAC-seq, Chromap BAM | 122.92 | 32.10 (3.8x) | 10.74 (11.4x) | 6.77 (18.2x) | 78,574 |
| ChIP-seq, ENCODE BAM + input control | 616.90 | 96.39 (6.4x) | 27.23 (22.7x) | 12.11 (50.9x) | 33,985 |
| CUT&RUN, GEO BAM + IgG control | 64.74 | 9.81 (6.6x) | 2.98 (21.7x) | 1.99 (32.5x) | 9,097 |

Peak resident set size spans 0.46-5.44 GiB against 0.35-2.74 GiB for MACS3.
Internal T=1 to T=24 scaling spans 3.5-8.0x.

All 72 measured invocations, every repetition with its per-run output
checksums, the stage profiles, host metadata, and the SHA-256 of each MACS3
executable and input are in [`paper/results/`](../paper/results/).
[`paper/results/README.md`](../paper/results/README.md) maps each recorded run
to the table rows it feeds, including which MACS3 build produced which.

## 15. Interpretation and limitations

- Scaling is primarily across chromosomes. The 1/4/24-thread comparison is a
  scaling study, not an assertion of unbounded parallelism.
- Maximum RSS includes input readers, tracks, caches, q-table construction, and
  output buffers for the complete CLI process.
- MACS3 is single-threaded for peak calling, so its reference row is reported
  at one thread. RapidMACS additionally reports one thread for a direct
  algorithm/runtime comparison and 4/24 threads for tool scaling.
- Results apply to paired-end narrow-peak calling with the pinned flags and
  datasets. They do not establish parity for single-end model building or
  broad peaks.
- Performance claims must use the final repeated CLI benchmark. Stage profiles
  may explain behavior but should not be presented as the primary comparison.

## 16. Source map

| Concern | Primary implementation |
|---|---|
| CLI and presets | `cli/rapidmacs.cc` |
| FRAG iterator/parser | `src/fragments.cc` |
| BAMPE/BEDPE normalization | `src/fragment_input.cc` |
| Track, lambda, and p-score construction | `src/macs3_frag_workspace.cc` |
| Pipeline scheduling and release policy | `src/macs3_frag_peak_pipeline.cc` |
| q-table, peak merging, and summit selection | `src/bdgpeakcall.cc` |
| MACS3 `pos_array` reconstruction | `src/macs3_pos_array.cc` |
| Full benchmark driver | `paper/scripts/run_paper_benchmarks.sh` |
| Dataset provenance and fixture creation | `tests/fixtures/README.md` |
