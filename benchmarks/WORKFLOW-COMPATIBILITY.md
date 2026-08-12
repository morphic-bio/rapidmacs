# Real-world workflow compatibility audit

Audit date: 2026-08-11. RapidMACS base revision:
`d3b4c7d8f04e8555e9d7bbb7c88012a75f74c799` plus the working-tree changes
described here. Reference caller: MACS3 3.0.3. Compiler: GCC 12.3.0.

This is a correctness and compatibility audit, not a timing benchmark. No
runtime or speed claim is made from these runs. Comparisons use complete files;
"pass" means both narrowPeak and summits are byte-identical.

## Inputs

| Role | Path | SHA-256 |
|---|---|---|
| PBMC 100K fragments | `fragments.tsv` | `ecc6825c65900fa85a19390ca5e9c333c8afff4aae84a809481583036243c119` |
| ChIP treatment | `ENCFF686KKV.bam` | `37ee55258385f908e139e934dcb542360e91ce0dc4dca4867cc35f5507fbe4ea` |
| ChIP control | `ENCFF181CXT.bam` | `400347e0ebf66079f47bb06593129d0237bf1ceac7579dd9442ace0d2a57ad1f` |
| CUT&RUN treatment | `SRR14255054.bam` | `28f4f45a31b4a90da5a68c44b093b166f0c8bde009b04e00d549636d5fcae95a` |
| CUT&RUN control | `SRR14255053.bam` | `33f74032ce65e141bb137582201162f4dddee3950e1983e8963097dc9f57149d` |

The BEDPE matrix inputs were derived losslessly from the listed BAMPE files
with `macs3 randsample -f BAMPE -p 100`. The exact conversion commands are in
the Phase 2 generated command log.

## Phase 1: Signac-style BED input

Pinned source: Signac `12d1ff82774890b08e55b3101f36fab83e674b88`,
[`R/peaks.R`](https://github.com/stuart-lab/signac/blob/12d1ff82774890b08e55b3101f36fab83e674b88/R/peaks.R).
The effective default invocation is BED, human genome size 2.7e9, q=0.05,
`--nomodel --extsize 200 --shift -100 --llocal 10000 --keep-dup 1`.

| Output | Rows | MACS3 MD5 | RapidMACS MD5 | Result |
|---|---:|---|---|---|
| narrowPeak (all 10 columns) | 3,367 | `860ca3753fdd5a79e7d7cf558cbb42d6` | `860ca3753fdd5a79e7d7cf558cbb42d6` | pass |
| summits BED | 3,367 | `0d21daa7cc5cca9c119383ad12176ec1` | `0d21daa7cc5cca9c119383ad12176ec1` | pass |

Recorded outputs are under
`build/workflow_compatibility/20260811T152047Z_phase1_signac/`. The commands
are preserved in both `macs3/time.txt` and `rapidmacs/time.txt`; the reusable
harness is `benchmarks/run_workflow_signac_parity.sh`.

## Phase 2: paired-end duplicate retention

RapidMACS accepts `--keep-dup all` or a positive integer for BAMPE and BEDPE.
The default remains `1`. Duplicate caps are applied independently to sorted
identical `(chrom, start, end)` runs in treatment and control. FRAG inputs are
already count-collapsed, so any explicit non-`all` value receives the MACS3
warning and is overridden to `all`.

| Assay | Format | `--keep-dup` | Peaks | narrowPeak MD5 | summits MD5 | Full-file result |
|---|---|---:|---:|---|---|---|
| ChIP-seq | BAMPE | 1 | 33,985 | `24740aa002242a59fa7afe9bc6679c9c` | `0268822f5398f1d6505eb8a1cd27b0e6` | pass |
| ChIP-seq | BAMPE | 2 | 33,985 | `aa32b2cc6c0ccbb3097cb113d7c57297` | `3a89710e35a898824b8c33b19cc2ccb3` | pass |
| ChIP-seq | BAMPE | all | 33,985 | `400b4eeda438608a35d70fd87b4f60c1` | `e8fa3bfec15cdaab955b63ca499763e4` | pass |
| ChIP-seq | BEDPE | 1 | 33,985 | `07cc2669084768eb7556f6426b9c4264` | `e91b02ab32b9a7c4304671144a466159` | pass |
| ChIP-seq | BEDPE | 2 | 33,985 | `b73e51e2ecf1035801d2a15d363d709e` | `5181dcfe06084371c9aef7b38717a554` | pass |
| ChIP-seq | BEDPE | all | 33,985 | `063a31ac3ab724f06b166524dc719575` | `e03aa840816f0c4aa1241b98be28a6cd` | pass |
| CUT&RUN | BAMPE | 1 | 9,097 | `bfadfec5434ebd63192a87557e6ac747` | `af17ce7847c6d74317d141ae75c2f2a5` | pass |
| CUT&RUN | BAMPE | 2 | 18,532 | `e877aeea24e8941f7bb3afedf983b9fa` | `af9e00bbb894a35e0972f3c7d091063c` | pass |
| CUT&RUN | BAMPE | all | 29,176 | `da06dbcdbe99a51a04730f7aeb04d1e1` | `7558c421b80ad3db662e9be15bf8760c` | pass |
| CUT&RUN | BEDPE | 1 | 9,097 | `1b1f3ce00f02bb4507b7896208c99555` | `a387944a8e6229580e3d2c9526811a0b` | pass |
| CUT&RUN | BEDPE | 2 | 18,532 | `9d3cb36e63962f60f0e974d4fc5b395c` | `691a88c8649841b8af87fbdfd9ad4d2d` | pass |
| CUT&RUN | BEDPE | all | 29,176 | `664df8bbd7ced22fd6f429d4b51b9fe2` | `958c1bb6dcc190877b07b7764309660a` | pass |

Harness: `benchmarks/run_workflow_keepdup_parity.sh`. Final results and exact
shell-escaped commands are under
`build/workflow_compatibility/20260811_phase2_keepdup/` in `results.tsv` and
`commands.sh`.

## Phase 3: upstream workflow defaults and blockers

The exact upstream source revisions and effective commands were inspected,
then the MACS3 command and closest supported RapidMACS command were run on the
same full PBMC input.

In the command table, `$PBMC` is the full PBMC path in the Inputs table,
`$MACS3` is
`macs3`, and
`$RAPIDMACS` is `rapidmacs`. `$OUT` is the tool's
subdirectory under the recorded phase output directory. The generated command
logs contain the same commands with every path fully expanded.

| Tool | Exact MACS3 command | Closest RapidMACS command | Byte-identical | Blocking behavior |
|---|---|---|---|---|
| Signac | `$MACS3 callpeak -t "$PBMC" -f BED -g 2700000000 -n signac_atac --outdir "$OUT/macs3" -q 0.05 --nomodel --extsize 200 --shift -100 --llocal 10000 --keep-dup 1` | `$RAPIDMACS --macs-profile signac-atac --macs-input-format BED --macs-bed-input "$PBMC" --macs-bed-peaks-out "$OUT/rapidmacs/signac_atac_peaks.narrowPeak" --macs-bed-summits-out "$OUT/rapidmacs/signac_atac_summits.bed" --macs-bed-summary "$OUT/rapidmacs/summary.tsv" --outdir "$OUT/rapidmacs" -n signac_atac` | yes | none |
| ArchR | `$MACS3 callpeak -g 2700000000 --name archr --treatment "$PBMC" --outdir "$OUT/macs3" --format BED --call-summits --keep-dup all --nomodel --nolambda --shift -75 --extsize 150 -q 0.1` | `$RAPIDMACS --macs-input-format BED --macs-bed-input "$PBMC" --outdir "$OUT/rapidmacs" --macs-bed-peaks-out "$OUT/rapidmacs/archr_peaks.narrowPeak" --macs-bed-summits-out "$OUT/rapidmacs/archr_summits.bed" --effective-genome-size 2700000000 --qvalue 0.1 --nomodel --shift -75 --extsize 150 --keep-dup all --call-summits --llocal 10000 -n archr` | no | BED `--nolambda` is unsupported |
| pycisTopic | `$MACS3 callpeak --treatment "$PBMC" --name pycistopic --outdir "$OUT/macs3" --format BEDPE --gsize hs --qvalue 0.05 --nomodel --shift 73 --extsize 146 --keep-dup all --call-summits --nolambda` | `$RAPIDMACS --preset cutrun -f BEDPE -t "$PBMC" -g hs -q 0.05 --keep-dup all -n pycistopic --out-prefix "$OUT/rapidmacs/pycistopic"` | no | paired `--nolambda` and paired `--call-summits` are unsupported |

| Tool | Pinned source behavior | MACS3 peaks / MD5 | Closest RapidMACS peaks / MD5 | BED3 | Summits | Blocker |
|---|---|---|---|---|---|---|
| Signac | BED; g=2.7e9; q=.05; nomodel; ext=200; shift=-100; keep=1 | 3,367 / `860ca3753fdd5a79e7d7cf558cbb42d6` | 3,367 / same | equal | MD5 equal | none |
| ArchR | BED; g=2.7e9; q=.1; nomodel; nolambda; ext=150; shift=-75; keep=all; call-summits | 12,923 / `eb6af8927c8f7422f7fa40940dafe119` | 5,520 / `5b4828c5aa1abd74c158aa8580fa892e` | different | different | BED `--nolambda` is unsupported |
| pycisTopic | BEDPE; caller-supplied g=`hs`; q=.05; nomodel; nolambda; keep=all; call-summits (shift/ext are neutral in paired mode) | 9,181 / `0fb01536d6ba4e6c73811718d7ad3275` | 581 / `18202be591f81ac23045fc5583975b4e` | different | different | paired `--nolambda` and paired `--call-summits` are unsupported |

ArchR source: `6feec354ad6c8052ddbc4626a2ca2d858ed465bf`,
[`R/ReproduciblePeakSet.R`](https://github.com/GreenleafLab/ArchR/blob/6feec354ad6c8052ddbc4626a2ca2d858ed465bf/R/ReproduciblePeakSet.R).
pycisTopic source: `787ce422a37f5975b0ebb9e7b19eeaed44847501`,
[`pseudobulk_peak_calling.py`](https://github.com/aertslab/pycisTopic/blob/787ce422a37f5975b0ebb9e7b19eeaed44847501/src/pycisTopic/pseudobulk_peak_calling.py).

Summits MD5s for the two blocked comparisons are:

| Tool | MACS3 summits MD5 | RapidMACS summits MD5 |
|---|---|---|
| ArchR | `9ea1ee52175d7d5930d8c3c9b7dd9dd4` | `b681f8239bdec744f91c918a81a71246` |
| pycisTopic | `090a329741c87fb894f446f5de279376` | `7b0efb2a7635f45796d62c1b00178bcc` |

Harness: `benchmarks/run_workflow_tool_audit.sh`. Outputs, full MD5 inputs,
and exact commands are under
`build/workflow_compatibility/20260811_phase3_tool_audit/`.

At the end of Phase 3, the next compatibility feature identified was
`--nolambda`, because it would unblock ArchR while removing local-window work.
Paired `--call-summits` was the second required gap for exact pycisTopic
compatibility. Phase 4 below closes both gaps.

## Phase 4: `--nolambda` and paired `--call-summits`

Both Phase 3 blockers are implemented as opt-in flags; existing defaults are
unchanged. `--nolambda` is shared by BED, FRAG, BAMPE, and BEDPE. It preserves
MACS3's genome-wide lambda and treatment/control depth-scaling arithmetic while
skipping treatment-derived and control-derived local-window construction.
`--call-summits` is likewise shared and implements MACS3's smoothed derivative,
peak-shape filter, fallback summit, and grouped subpeak naming.

The Phase 3 MACS3 3.0.3 files were reused as the immutable reference; MACS3 was
not rerun. RapidMACS was run serially for ArchR and then pycisTopic with the
now-exact commands recorded in the generated `commands.sh`.

| Tool | RapidMACS compatibility command difference from Phase 3 | Rows | narrowPeak MD5 | summits MD5 | Result |
|---|---|---:|---|---|---|
| ArchR | added `--nolambda`; existing BED `--call-summits` now performs subpeak calling | 12,923 | `eb6af8927c8f7422f7fa40940dafe119` | `9ea1ee52175d7d5930d8c3c9b7dd9dd4` | full-file pass |
| pycisTopic | uses paired callpeak outputs with `--nolambda --call-summits` | 9,181 | `0fb01536d6ba4e6c73811718d7ad3275` | `090a329741c87fb894f446f5de279376` | full-file pass |

For both workflows the row count, full 10-column narrowPeak MD5, BED3, and
full summits BED MD5 equal MACS3. Results, the reference pointer, stderr, and
fully expanded commands are under
`build/workflow_compatibility/20260811_phase4_flags/`. The reusable harness is
`benchmarks/run_workflow_tool_audit.sh`; set `MACS_REFERENCE_ROOT` to a prior
completed audit to reuse its unchanged MACS3 files.

An additional controlled BAMPE check used the committed CUT&RUN chr21 fixtures
with `--keep-dup 1 --nolambda --call-summits`. All 182 rows were byte-identical
(`narrowPeak` MD5 `a126a4655e086fbab5a68c381a990793`; summits MD5
`a8fd175775c34bbe6e96127eb5ad983c`). Outputs and logs are under
`build/workflow_compatibility/20260811_phase4_control_check/`.

## Verification

`make` succeeds. `make test` passes all 10 unit tests, the `--keep-dup` CLI
regression, the legacy and MACS3 100K FRAG smokes, and the ChIP-seq and
CUT&RUN paired-control smokes. The regression explicitly verifies the
unchanged default (`--keep-dup 1`), values `2` and `all`, invalid paired-end
values, and FRAG warn-and-override behavior. The focused compatibility unit
test verifies that `--nolambda` allocates no local-window events, constructs
the exact global background/position array, selects MACS3's smoothed summit
instead of the old plateau midpoint, and writes grouped subpeak names.
