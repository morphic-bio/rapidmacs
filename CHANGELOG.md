# Changelog

## 1.0.2 (proposed, unreleased)

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

Validation includes original MACS3 3.0.3 narrowPeak and summit goldens for
singleton strands in P and q modes with duplicate limits 1, 2, and `all`.
On a 5,955,396-row BED input, regenerated MACS3 references and RapidMACS are
byte-identical for all 23,932 peaks/summits in P mode and all 30,677 in q mode.

No release or tag has been created. Consumers must update their RapidMACS
source revision and rebuild the CLI or rebuild and relink `librapidmacs`.
