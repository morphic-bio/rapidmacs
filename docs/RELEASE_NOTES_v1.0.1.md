# RapidMACS v1.0.1 Release Notes

RapidMACS v1.0.1 is a correctness release for exact MACS3 v3.0.3
`callpeak` compatibility.

## Correctness fixes

- Match `CallPeakUnit`'s strict `score > cutoff` region predicate in both
  p-value and q-value modes. A score exactly equal to the requested cutoff is
  no longer included in a peak boundary.
- Apply MACS3 `callpeak`'s default fold-enrichment cutoff of 1.0 after peak
  construction, excluding statistically significant depletion regions from
  narrow-peak output.

Both fixes apply to the standalone MACS3-compatible FRAG/BAMPE/BEDPE/BED
callers and to applications embedding `librapidmacs`. The historical custom
streaming caller is unchanged. The separate `bdgpeakcall` diagnostic retains
its inclusive threshold, matching MACS3's `bdgpeakcall` command.

## Validation

- Exact-cutoff regressions cover p-value and optimized q-value thresholding.
- The full RapidMACS unit suite and the available ChIP-seq and CUT&RUN parity
  smokes pass.
- A full-scale ATAC comparison revalidated complete narrowPeak and summit
  outputs byte-for-byte against MACS3 v3.0.3 after the fixes.

There are no public API or command-line compatibility breaks in this release.
