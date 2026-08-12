# `rapidmacs` command-line reference

`rapidmacs` replaces `macs3 callpeak` for narrow peak calling. It reads
fragments, BAMPE, BEDPE, or single-end BED and writes narrowPeak and summit
files byte-identical to MACS3 v3.0.3.

`rapidmacs --help` prints the full flag list. This page explains what the
flags mean and which ones you actually need.

## Choosing an invocation

There are two output surfaces, and which flags apply depends on which you use.

| You have | Use | Outputs |
|---|---|---|
| A fragments TSV (10x / Cell Ranger ARC / Chromap) | `-f FRAG` with `--macs3-frag-narrowpeak` / `--macs3-frag-summits` | narrowPeak + summits at the paths you name |
| A coordinate-sorted paired-end BAM | `--preset chip` or `--preset cutrun` with `--out-prefix` | `<prefix>.narrowPeak` and companions |
| A BEDPE file | `-f BEDPE` with `--out-prefix` | `<prefix>.narrowPeak` and companions |
| A single-end BED tag file | `--macs-profile signac-atac` with `--macs-bed-*` | narrowPeak + summits at the paths you name |

## Presets

`--preset NAME` applies a group of MACS3-equivalent defaults and then keeps
parsing, so **any explicit flag placed after it wins**. Put `--preset` first.

| Preset | Sets | Requires |
|---|---|---|
| `shiftTAG` | `--bdgpeakcall-cutoff 2` (= `-p 0.01`), matching `macs3 callpeak --nomodel --extsize 150 --shift -75 --keep-dup all --call-summits -g hs` | — |
| `chip` | BAMPE input, `-q 0.05`, MACS3's dynamic paired-end peak geometry | `-c` |
| `cutrun` | as `chip` | control optional |

## Input

| Flag | Meaning |
|---|---|
| `-t`, `--treatment`, `-i`, `--input` | Treatment input. Required. |
| `-c`, `--control` | Matched input or IgG control. Required by `--preset chip`. |
| `-f`, `--format {FRAG,BAMPE,BEDPE}` | Input format. Default `FRAG`. |
| `-n`, `--name NAME` | Peaks are named `<NAME>_peak_<i>`. Defaults to `NA`, as MACS3 does. |
| `--keep-dup {all,N}` | Keep all duplicate coordinates, or at most `N`. Default `1`. BAMPE/BEDPE/BED only — FRAG counts are already collapsed, so FRAG always keeps all and a non-`all` value warns and is overridden. |

`-f FRAG` auto-detects two schemas and reads the duplicate count from the
right column:

- five-column barcoded: `chrom, start, end, barcode, count` (column 5)
- seven-column Chromap bulk: `chrom, start, end, N, mapq, strand, count` (column 7)

Three- and four-column inputs are accepted with an implied count of 1. Plain or
gzipped. Parsing is strict: a missing or non-numeric field, a negative start,
`end <= start`, `count < 1`, or a line over 65535 bytes is an error reported as
`path:line:reason`. Lines starting with `#` and blank lines are skipped.

## Thresholds and peak geometry

| Flag | Default | Meaning |
|---|---|---|
| `-q`, `--qvalue Q` | 0.05 in `chip`/`cutrun` | q-value cutoff |
| `-p`, `--pvalue P` | — | p-value cutoff; mutually exclusive with `-q` |
| `-g`, `--gsize {hs,mm,ce,dm,N}` | `hs` (2913022398) | effective genome size |
| `--bdgpeakcall-cutoff X` | 5 | `-log10` score cutoff on the FRAG path |
| `--macs3-frag-pvalue P` | 1e-5 | sets the cutoff from a raw p-value |
| `--macs3-frag-qvalue Q` | — | q-value form; mutually exclusive with an explicit p-value |
| `--bdgpeakcall-min-len N` | 200 | minimum peak length, MACS3 `-l` |
| `--bdgpeakcall-max-gap N` | 30 | maximum gap to link, MACS3 `-g` |
| `--frag-lambda-slocal N` | 1000 | small local background window |
| `--frag-lambda-llocal N` | 10000 | large local background window |
| `--scale-to {small,large}` | `small` | treatment/control depth target |
| `--to-large` | — | alias for `--scale-to large` |
| `--nolambda` | off | genome-wide background only; skips local-window work |
| `--call-summits` | off | MACS3's smoothed-signal subpeak algorithm; emits `1a`, `1b`, … names |

## Output

| Flag | Meaning |
|---|---|
| `--out-prefix PATH` | Writes `<PATH>.narrowPeak` and companions. The BAMPE/BEDPE surface. |
| `--macs3-frag-narrowpeak PATH` | narrowPeak output on the FRAG path. |
| `--macs3-frag-summits PATH` | Summits BED5 on the FRAG path, one row per summit. |

## Performance

| Flag | Default | Meaning |
|---|---|---|
| `--peak-caller-threads N` | 1 | Process chromosomes concurrently. The published scaling is 3.5–8.0× from 1 to 24 threads, bounded by the ~24 large human chromosomes and by serial input parsing. |
| `--profile-stages PATH` | — | Per-stage timing TSV (wall, user, sys, RSS, rows, bytes). |

Threading is over chromosomes, so more threads than large chromosomes buys
nothing. BAMPE input scales better than gzipped fragments: a BAM is a series of
independently inflatable BGZF blocks, whereas a plain gzip file is one deflate
stream that must be inflated serially.

## Parity controls

| Flag | Default | Meaning |
|---|---|---|
| `--no-macs3-frag-summit-parity` | parity **on** | Disables byte-identical summit mode in exchange for lower memory. Diverges from MACS3 on ~2% of summit positions at PBMC scale (max ~222 bp). BED3, bedGraph, and score columns are unaffected. |
| `--frag-pileup-macs3-uint8-counts` | — | Coerce fragment weights through uint8 as MACS3's PETrackII does. Recommended when comparing against `macs3 pileup`. |
| `--frag-pileup-max-count N` | 0 | Cap the count column like MACS3 `--max-count`. Applies to `--frag-span-pileup-bdg` only. |

Leave summit parity on unless you are memory-constrained and do not need
byte-identical summits.

## Single-end BED compatibility

For workflows built around Signac's BED input:

| Flag | Meaning |
|---|---|
| `--macs-profile signac-atac` | Resolves BED, `--nomodel`, `--extsize 200`, `--shift -100`, `-g 2700000000`, `-q 0.05`, `--llocal 10000` |
| `--macs-bed-input PATH` | Existing BED tag input |
| `--macs-input-format BED` | Required on this path |
| `--macs-bed-peaks-out PATH` | narrowPeak output |
| `--macs-bed-summits-out PATH` | Summits output |
| `--macs-bed-summary PATH` | Resolved parameters and projection statistics |
| `--nomodel`, `--extsize N`, `--shift N`, `--llocal N` | Individual overrides |
| `--effective-genome-size N` | `-g` on this path |

Single-end input requires `--nomodel`: MACS3's fragment-size model building is
not implemented.

To project a fragments file into BED first:

| Flag | Meaning |
|---|---|
| `--input-fragments PATH` | Fragments to project |
| `--cells PATH` | Barcode allow-list |
| `--barcode-group-tsv PATH` | Barcode-to-group table for per-group projection |
| `--group-column NAME` | Group column when the table has a header |
| `--macs-bed-out PATH` | Single-sample BED output |
| `--macs-bed-out-dir DIR` | Per-group BED output directory |

## Diagnostics

These write intermediate tracks. They are for inspection and comparison
against MACS3's own subcommands, not part of normal peak calling.

| Flag | Corresponds to |
|---|---|
| `--pileup-bdg PATH` | Binned Tn5-extended *cut* pileup. Differs from `macs3 pileup -f FRAG`. |
| `--frag-span-pileup-bdg PATH` | Fragment-span pileup. Matches `macs3 pileup -f FRAG` when `--frag-pileup-macs3-uint8-counts` is set. |
| `--frag-lambda-bdg PATH` | MACS3 `control_lambda` |
| `--frag-score-ppois-bdg PATH` | `macs3 bdgcmp -m ppois` |
| `--frag-score-fe-bdg PATH` | `macs3 bdgcmp -m FE` |
| `--frag-score-pseudocount X` | `macs3 bdgcmp -p` |
| `--bdgpeakcall-input PATH` / `--bdgpeakcall-output PATH` | `macs3 bdgpeakcall` on an existing score track |

## Exit status

`0` on success. Non-zero on a usage error, a missing or malformed input, or an
I/O failure. Parse errors name the offending file and line.
