# Recorded benchmark results

Every measurement reported in the paper, as the harness wrote it.

**The RapidMACS timings in the paper are v1.0.1**, from
`20260822T000000Z_rapidmacs_v1.0.1_refresh`. The MACS3 comparators are
unchanged by that release and were not relaunched, so the MACS3 medians below
remain the reported baselines. Peak counts and every output MD5 are identical
between v1.0.0 and v1.0.1, so the identity table in Supplementary S1 covers
both. The v1.0.0 RapidMACS directories are kept as the record of that release
and are no longer the reported timings.

Files are copies from the run archives, with absolute paths reduced to their
basenames; no measurement, checksum or count has been altered.

`summary.tsv` / `comparison.tsv` give medians with the observed range;
`raw.tsv` retains all three recorded repetitions with the per-run narrowPeak,
BED3, summit, and summit-coordinate MD5 sums. Together these are the 72
measured invocations Supplementary S5 refers to. `metadata.tsv` records host,
dates, tool versions, and run identifiers; `macs3_sha256s.txt` and
`input_sha256s.txt` pin the exact executable and inputs.

The per-run output files themselves (narrowPeak, summits, bedGraph, stage
profiles for every repetition) run to about 16 GB and are not in git. The
checksums here identify them exactly.

| Directory | Configurations | MACS3 build | Feeds |
|---|---|---|---|
| `20260810T163407Z_ssd_paper_matrix` | scATAC Chromap; ChIP-seq; CUT&RUN | no-control-dump | Table 1 row 2, MACS3 column; the MACS3 walls of S4's "extra file removed" rows. Its RapidMACS timings are superseded by the v1.0.1 refresh |
| `20260810T201105Z_encode_bulk_atac_bampe` | Bulk ATAC, ENCODE BAM and Chromap BAM | no-control-dump (treatment-only, so the write is unreachable) | Table 1 rows 3–4, MACS3 column. RapidMACS timings superseded |
| `20260810T212726Z_10x_pbmc_supplied_frag` | scATAC, 10x fragments | no-control-dump (treatment-only) | Historical only. Its MACS3 time of 748.49 s is superseded by the released-build run below, and its RapidMACS timings by the v1.0.1 refresh |
| `20260811T194844Z_released_macs3_controlled` | ChIP-seq; CUT&RUN | **released** | Table 1 rows 5–6, MACS3 column; the "as released" rows of S4 |
| `20260812T152652Z_released_macs3_atac_frag` | scATAC 10x; bulk ATAC Chromap fragments | **released** | Table 1 row 1, MACS3 column; the 0.13% released-vs-patched bound quoted in S2 |
| `20260822T000000Z_rapidmacs_v1.0.1_refresh` | All six, RapidMACS only | n/a (MACS3 not relaunched) | **Every RapidMACS timing in the paper: the RapidMACS columns of Table 1, S3's RapidMACS rows, S4's speedups, and S5.** `results.tsv` carries per-row parity verdicts; `stages.tsv` carries the S5 stage medians |

Two notes on what is and is not here.

The three treatment-only configurations were measured against the
no-control-dump build. That build differs from the release only in four lines
that a run without a control never reaches, so the executed path is identical.
Running scATAC 10x against both builds bounds the difference at 0.13%: 747.55 s
released against 748.49 s. Supplementary S2 states this.

The 2026-08-12 run also measured released MACS3 on the Chromap-realigned bulk
ATAC **fragments** file: 342.94 s, 77,316 peaks. There is no RapidMACS arm for
that input, so it is not a row in any table. It is kept because it was measured.
