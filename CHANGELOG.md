# Changelog

## 1.0.2 (2026-10-12)

- Fix the BED single-end background total under numeric duplicate removal.
  In v1.0.1, a chromosome strand containing exactly one tag contributed to
  the filtered total. MACS3 3.0.3 retains that tag in the pileup but excludes
  it from the total used to scale genome-wide background. This could change
  reported fold enrichment, P and q scores, and threshold-dependent results.
  `--keep-dup all` continues to count every tag. BAMPE/BEDPE treatment and
  control singleton fragments already follow MACS3's paired-end rule, which
  includes them in the total.
- Fix BED summit column 5 under `--pvalue`: v1.0.1 wrote minus-log10 q scores
  instead of MACS3's minus-log10 P scores. The shared writer now accepts the
  summit metric explicitly; q-value mode continues to write q scores.
- Fix FRAG summit endpoint multiplicities when `frag_pileup_max_count` (CLI
  `--frag-pileup-max-count`, corresponding to MACS3 `--max-count`) is set.
  In v1.0.1, pileups and background totals used the capped count, while
  endpoints used for summit placement were repeated by the uncapped count.
  This could change summits and their associated fold enrichment, P, q,
  and integer scores. Event, sorted-batch, and sweep workspaces now use the
  same effective count throughout. Uncapped reference outputs are unchanged.

Validation includes original MACS3 3.0.3 narrowPeak and summit goldens for
singleton strands in P and q modes with duplicate limits 1, 2, and `all`.
On a 5,955,396-row BED input, regenerated MACS3 references and RapidMACS are
byte-identical for all 23,932 peaks/summits in P mode and all 30,677 in q mode.
The same five-column input, parsed as FRAG with `--max-count 1`, matches all
19,487 peaks and summits. A 120-row synthetic FRAG fixture matches original
MACS3 for caps 1, 2, 5, and no cap through the CLI and four library entry points.
The cap-1 and cap-2 goldens fail on v1.0.1. Existing 100K legacy/FRAG,
ChIP-seq, and CUT&RUN goldens remain byte-identical.

Consumers must update their RapidMACS source revision to `v1.0.2` and rebuild
the CLI or rebuild and relink `librapidmacs`.
