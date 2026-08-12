# Paper benchmark protocol

Everything behind the published results: the scripts that produced every
reported number, both MACS3 reference builds, and the recorded measurements.

## Paths and host guards

The scripts here are the ones that were actually run. Machine-specific paths
have been replaced by environment variables, all defaulting to directories
under the repository so nothing outside it is touched:

| Variable | Default | Holds |
|---|---|---|
| `DATA_DIR` | `./data` | input datasets |
| `WORK_ROOT` | `./work` | run outputs |
| `ARCHIVE_ROOT` | `./archive` | completed runs, copied on success |
| `MACS3_REF_ROOT` | `./.macs3-reference` | MACS3 reference installations |
| `CHROMAP_SUITE` | `./chromap-suite` | Chromap Suite, for the realignment scripts |
| `CELLRANGER_REF` | `$DATA_DIR/refdata-cellranger-arc-GRCh38-2020-A-2.0.0` | Cell Ranger ARC reference |

The paper's protocol also required the reference host and NVMe placement. Those
checks are **advisory** here, so the harnesses run anywhere: they warn rather
than abort, since timings from another machine are not comparable to the paper
anyway. Set `STRICT_HOST=1` to restore them as hard failures.

Input checksum verification is *not* relaxed. A harness still refuses to start
against the wrong MACS3 build or an input whose SHA-256 does not match, because
those would silently invalidate the comparison.

## What produced what

`results/README.md` maps each recorded run to the table rows it feeds and
records which MACS3 build produced it. In outline:

| Script | Produces |
|---|---|
| `scripts/run_10x_pbmc_example.sh` | The single-dataset reproduction. Start here. |
| `scripts/setup_macs3_reference.sh` | Both MACS3 reference builds |
| `scripts/run_paper_benchmarks.sh` | The full assay matrix |
| `scripts/run_stock_macs3_controlled_ssd.sh` | Released-build ChIP-seq and CUT&RUN baselines |
| `scripts/run_stock_macs3_atac_frag_ssd.sh` | Released-build scATAC baselines |
| `scripts/run_encode_bulk_atac_*.sh` | Bulk ATAC rows |
| `scripts/run_treatment_only_frag_benchmarks.sh` | Treatment-only fragment runs |
| `scripts/stage_*.sh`, `scripts/prepare_fastq_fragments.sh` | Input staging and Chromap realignment |
| `scripts/run_workflow_*.sh` | The Signac / ArchR / pycisTopic audit in `WORKFLOW-COMPATIBILITY.md` |
| `scripts/summarize_*.py` | Medians and ranges from raw run output |

`fastq_sources.tsv` and `fragment_sources.tsv` record the accessions and
checksums of every input.

The real-world Signac/ArchR/pycisTopic and paired-end duplicate-retention audit
is recorded in `WORKFLOW-COMPATIBILITY.md`; its three reusable command harnesses
are named `run_workflow_*_parity.sh` and `run_workflow_tool_audit.sh`.
The tool-audit harness covers the exact ArchR and pycisTopic
`--nolambda --call-summits` commands and accepts `MACS_REFERENCE_ROOT` to reuse
an already-recorded MACS3 reference without rerunning it.

`run_paper_benchmarks.sh` compares standalone RapidMACS with a clean MACS3
3.0.3 installation on the same post-alignment input.  Mapping/alignment is not
timed.  The default protocol runs one unrecorded warm-up followed by three
recorded repetitions for MACS3 and for RapidMACS at 1, 4, and 24 peak-caller
threads.  `summary.tsv` reports medians and the observed range; `raw.tsv`
retains every measurement and full narrowPeak, BED3, full-summit, and
summit-coordinate checksums. Every assay, including ATAC FRAG, requires full
narrowPeak and summit byte identity with MACS3 3.0.3.

The recorded runs are warm-cache by design.  Dropping the Linux page cache
requires root and would disrupt other work on the shared benchmark host.  An
unrecorded full run for each tool/thread setting both warms the input and avoids
mixing cold-cache noise into the tool comparison.  Wall time and maximum RSS
come from GNU `/usr/bin/time`; each run has its own output directory.

## Reproducing one row: 10x PBMC 3k

The quickest check that the headline claim holds. `run_10x_pbmc_example.sh`
fetches the public 10x Cell Ranger ARC PBMC 3k ATAC fragments, verifies the
download against the SHA-256 the paper used, runs both callers, and exits
non-zero unless every output file carries the published narrowPeak and summits
MD5 sums:

```bash
make -j
paper/scripts/setup_macs3_reference.sh VARIANT=released
MACS3=.macs3-reference/macs3-3.0.3-venv/bin/macs3 \
  paper/scripts/run_10x_pbmc_example.sh
```

Budget about 25 minutes, most of it the single MACS3 run. `SKIP_MACS3=1` times
RapidMACS alone and checks it against the published checksums instead.
`WARMUPS=1 REPEATS=3` reproduces the paper's protocol rather than the single
timed run the example defaults to.

## Building the MACS3 references

```bash
paper/scripts/setup_macs3_reference.sh
```

This builds two installations from the same SHA-256-pinned PyPI sdist:
`macs3-3.0.3-venv`, the unmodified release, and
`macs3-3.0.3-no-control-dump-venv`, the release with four lines removed.

Those four lines, at the end of `PETrackI.pileup_a_chromosome_c`, write the
internal PE control pileup to `tmp_b_c.txt`; they are present in both the PyPI
sdist and the upstream `v3.0.3` tag, and we found no option that disables them.
The removal is checked in as `macs3-3.0.3-no-control-dump.patch`, and both
builds produce identical narrowPeak and summits output.

**The released build is what Table 1 and Supplementary S3 report.** The
no-control-dump build appears in Supplementary S4, which gives both side by
side, and produced the three treatment-only configurations, where that code is
never reached. `paper/results/README.md` records which build produced which row.

Run the complete six-configuration matrix:

```bash
make -j
paper/scripts/run_paper_benchmarks.sh
```

Useful controlled reductions are available through environment variables,
for example `ASSAYS='chipseq cutrun'`, `THREADS='1 4'`, `REPEATS=3`, or
`WARMUPS=0`.  The script refuses missing datasets or a MACS3 version other than
3.0.3. After every assay it requires controlled RapidMACS runs to have the same
full narrowPeak and summit MD5 as MACS3. ATAC runs must have the same peak and
summit coordinate MD5, preventing performance numbers from being reported for
a non-parity build.

When only RapidMACS changed, reuse an existing raw MACS3 result as the parity
oracle instead of rerunning the unchanged reference caller:

```bash
TOOLS=rapidmacs \
REFERENCE_RAW_TSV=/path/to/previous/raw.tsv \
ASSAYS='chipseq cutrun' \
paper/scripts/run_paper_benchmarks.sh
```

This mode still performs the configured warm-ups and repetitions and applies
the same assay-specific parity checks. Reused ATAC references must include the
`summits_bed3_md5` column.

## What was actually measured

`results/` holds every recorded measurement behind the paper: medians, all
three repetitions, per-run output checksums, stage profiles, host metadata, and
the SHA-256 of each MACS3 executable and input. `paper/results/README.md` maps each
run directory to the table rows it feeds.

Each recorded run was driven by a harness that refuses to start against the
wrong MACS3 build, a checksum mismatch on any input, a non-NVMe output
directory, or a host other than the one the paper reports.
