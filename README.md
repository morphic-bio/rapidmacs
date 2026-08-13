# RapidMACS

A C++ narrow peak caller whose output is byte-identical to MACS3 v3.0.3, and
3.5–51× faster depending on assay, input format, and thread count.

RapidMACS calls narrow peaks, with or without a matched input or IgG control,
from fragment, BAMPE, BEDPE, or single-end BED input. It ships both as a
standalone `rapidmacs` executable that replaces `macs3 callpeak`, and as a
static library with no Python or Cython runtime, so an aligner can call peaks
in-process on fragments it already holds.

Identical means literally identical: the same bytes, not the same peaks within
a tolerance. Anything computed from a RapidMACS peak or summit file returns
exactly what it would have returned from the MACS3 file, because it is handed
the same file.

## Scope

RapidMACS implements one of MACS3's fourteen subcommands. It is not a
replacement for MACS3.

- **Implemented:** narrow `callpeak` on paired-end fragment, BAMPE, and BEDPE
  input, and on single-end BED with `--nomodel`; treatment-only or matched
  control; `--nolambda`; `--call-summits`.
- **Not implemented:** broad peak calling, which is a separate algorithm with
  its own pair of cutoffs and region-linking pass; single-end fragment-size
  model building; every other subcommand.

Workflows that need any of those should keep using MACS3. Nothing prevents both
from being installed side by side.

## Building

```sh
git submodule update --init --recursive   # first time only
make                                       # lib/librapidmacs.a + bin/rapidmacs
```

The build links against system htslib by default. To build and link the pinned
submodule copy instead:

```sh
(cd third_party/htslib && autoreconf -i && ./configure && make)
make HTSLIB_DIR=third_party/htslib
```

## Usage

`-f FRAG` accepts five-column Cell Ranger / Chromap barcoded fragments
(`chrom,start,end,barcode,count`) and seven-column Chromap bulk fragment BED
(`chrom,start,end,N,mapq,strand,count`), reading the duplicate count from
column 5 or 7 according to the detected schema. As in MACS3, FRAG duplicate
counts are already collapsed, so an explicit non-`all` `--keep-dup` warns and
is overridden to `all`.

```sh
# Single-cell or bulk ATAC from a fragments file, no control.
bin/rapidmacs --preset shiftTAG -f FRAG -t fragments.tsv.gz -g hs \
  -n atac -p 1e-5 --bdgpeakcall-min-len 200 --bdgpeakcall-max-gap 30 \
  --macs3-frag-narrowpeak out/atac.narrowPeak \
  --macs3-frag-summits out/atac.summits.bed

# ChIP-seq: the chip preset requires a control.
bin/rapidmacs --preset chip -t ctcf.bam -c input.bam \
  -n ctcf --out-prefix out/ctcf

# CUT&RUN: matched IgG is optional.
bin/rapidmacs --preset cutrun -t ctcf_cutrun.bam -c igg.bam \
  -n ctcf_cutrun --out-prefix out/ctcf_cutrun

# Standard BEDPE (six mate columns).
bin/rapidmacs -f BEDPE -t fragments.bedpe -q 0.05 \
  -n bedpe --out-prefix out/bedpe
```

`--peak-caller-threads N` processes chromosomes concurrently. Presets are
applied first and then overridden by any explicit flag, so place `--preset`
before per-flag overrides.

As a library, the entry point is
`RunMacs3FragPeakPipelineFromSortedIterator`, which pulls fragments one at a
time through the `FragmentIterator` interface in
[include/rapidmacs/fragments.h](include/rapidmacs/fragments.h), so they can come
from a file or straight out of an array the host already holds.

## Output parity with MACS3

Every narrowPeak and summit field matches MACS3 v3.0.3 — intervals, peak names,
summit offsets and coordinates, signal values, p-scores and q-scores — and the
full-file MD5 sums agree. This holds across single-cell ATAC-seq, bulk ATAC-seq,
ChIP-seq, and CUT&RUN, on fragment and BAMPE input, with and without a control.

Summit BED scores follow MACS3's selected score column: p-score for `-p`,
q-score for `-q`.

Two low-level details had to be matched exactly and are worth knowing about.

**The summit tie-break.** MACS3's internal `pos_array`, the array iterated to
find the median-tied summit, is per-event rather than per-distinct-value. When
a fragment endpoint lands inside a plateau without changing the running pileup
— one fragment ending where another starts, net Δ of zero — MACS3 still emits
an entry. The bedGraph writer collapses those adjacent same-value runs, so an
exported `.bdg` agrees either way, but the median is taken over more
sub-pieces than a collapsed track exposes.
[src/macs3_pos_array.cc](src/macs3_pos_array.cc) reproduces this directly and
is on by default. `--no-macs3-frag-summit-parity` takes a lower-memory
shortcut that diverges on roughly 2% of summits (max ~222 bp at PBMC scale);
BED3, bedGraph, and score columns are unaffected.

**Score-table arithmetic.** The callpeak path retains MACS3's float32 lambda
for Poisson scoring and uses its float32 p/q-table arithmetic. Five-decimal
rounding is confined to diagnostic bedGraph output.

## Container

```sh
docker build -t rapidmacs .

docker run --rm -v "$PWD:/data" rapidmacs \
  --preset shiftTAG -f FRAG -t /data/fragments.tsv.gz -g hs -n atac \
  -p 1e-5 --bdgpeakcall-min-len 200 --bdgpeakcall-max-gap 30 \
  --macs3-frag-narrowpeak /data/atac.narrowPeak \
  --macs3-frag-summits /data/atac.summits.bed
```

The default image is about 89 MB and contains only the binary and its runtime
libraries; `rapidmacs` is the entrypoint, so arguments pass straight through.
It links Debian's htslib rather than the pinned submodule, which affects BAM
reading only and not peak output.

Tagged releases publish a multi-arch image (amd64 and arm64) to
`ghcr.io/morphic-bio/rapidmacs`, alongside per-architecture tarballs of the
binary, static library, and headers.

`--target reproduce` builds a larger image that also carries released MACS3
v3.0.3 and the `paper/` scripts, so the published comparison can be run inside
the container:

```sh
docker build -t rapidmacs:reproduce --target reproduce .
docker run --rm -it -v "$PWD:/data" rapidmacs:reproduce
# $MACS3 and $RAPIDMACS are set; paper/scripts is on disk
```

## Reproducing the published benchmark

`paper/scripts/run_10x_pbmc_example.sh` reproduces the single-cell ATAC row of
Table 1 from a public dataset. It downloads the 10x Cell Ranger ARC PBMC 3k
ATAC fragments, verifies the file against the SHA-256 the paper used, calls
peaks with released MACS3 3.0.3 and with RapidMACS at 1, 4, and 24 threads, and
exits non-zero unless every output file carries the published narrowPeak and
summits MD5 sums:

```sh
make -j
paper/scripts/setup_macs3_reference.sh VARIANT=released
MACS3=.macs3-reference/macs3-3.0.3-venv/bin/macs3 \
  paper/scripts/run_10x_pbmc_example.sh
```

Expect roughly 25 minutes, nearly all of it the one MACS3 run. On the reference
host this configuration is 747.55 s for MACS3 against 57.87, 24.13, and 16.37 s
for RapidMACS, all four producing 50,003 identical peaks. `SKIP_MACS3=1` times
RapidMACS alone against the published checksums; `WARMUPS=1 REPEATS=3`
reproduces the paper's protocol rather than the single timed run the example
defaults to.

[paper/README.md](paper/README.md) covers the complete
six-configuration matrix and both MACS3 reference builds.
[paper/results/](paper/results/) holds every recorded measurement:
medians, all three repetitions with per-run output checksums, stage profiles,
and the SHA-256 of each MACS3 executable and input.

## Workflow compatibility

The three tools that most often re-call peaks on data a pipeline has already
produced were checked against their exact command lines: Signac, ArchR (BED,
`--nolambda --call-summits`), and SCENIC+/pycisTopic (BEDPE, likewise). All
three default to `macs2` and can be pointed at `macs3`; against the MACS3
v3.0.3 commands, RapidMACS reproduces the narrowPeak and summit files byte for
byte — 12,923 and 9,181 rows respectively on the full PBMC 100K workflow input.
Commands and checksums are in
[paper/WORKFLOW-COMPATIBILITY.md](paper/WORKFLOW-COMPATIBILITY.md).

Because the files are identical, the substitution can be made part-way through
an established analysis without putting earlier comparisons in question.

## Testing

```sh
make test                    # unit + ATAC + ChIP-seq + CUT&RUN smokes
make test-unit               # C++ unit tests
make test-macs3-frag-100k    # FRAG golden
make test-chipseq            # public ENCODE BAMPE treatment/input fixture
make test-cutrun             # public GEO BAMPE CTCF/IgG fixture
```

Tests that need a fragments file take one from `$DATA_DIR` (default `./data`)
or from an explicit override:

```sh
FIXTURE_FRAGS=/path/to/fragments.tsv.gz make test-macs3-frag-100k
```

## Documentation

- [docs/walkthroughs.md](docs/walkthroughs.md) — five worked examples, from a
  first run to embedding the library.
- [docs/cli.md](docs/cli.md) — complete `rapidmacs` command-line reference.
- [docs/api.md](docs/api.md) — library API: the fragment interface, the
  parameter struct, and the four entry points.
- [docs/methodology.md](docs/methodology.md) — implemented algorithms, parity
  criteria, and the benchmark protocol.
- [paper/](paper/) — everything behind the published results: reproduction
  scripts, both MACS3 reference builds, and the recorded measurements.
- [NOTICE](NOTICE) — third-party attributions.

## Layout

```
include/rapidmacs/  public headers
src/                library implementation
cli/                standalone CLI
tests/              unit tests and parity smokes
docs/               CLI reference, API reference, walkthroughs, methodology
paper/
  scripts/          every script behind a published number
  results/          every recorded measurement
Dockerfile          slim runtime image; --target reproduce adds MACS3 3.0.3
```

## License

RapidMACS is distributed under the **MIT** license — see [LICENSE](LICENSE).

RapidMACS is an **independent implementation, not a port**. It was written
clean-room from MACS3's published algorithmic specification and validated by
output comparison against MACS3 v3.0.3; no MACS3 code is included, in original
or translated form. MACS3's source was consulted only for the handful of
low-level details needed to reach byte-for-byte equivalence, chiefly the
per-event `pos_array` summit tie-break. That independence is also why RapidMACS
is faster: a specialized paired-end narrow-peak codebase can pick
representations a general caller cannot, rather than inheriting MACS3's
structure at a faster language's speed.

MACS was originally developed in X. Shirley Liu's group and has been maintained
and extended since by Tao Liu and the MACS3 contributors. This project would not
have been possible without that work, or without its release under a permissive
open-source license: being able to read and run the reference implementation is
what made byte-level parity achievable.
