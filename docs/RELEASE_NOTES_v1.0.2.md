# RapidMACS v1.0.2

This patch release fixes three cases where v1.0.1 differed from original
MACS3 3.0.3 `callpeak` output.

- **Single-end BED background totals:** with numeric `--keep-dup`, MACS3
  retains chromosome strands containing exactly one tag in the pileup but
  excludes them from the filtered total used for genome-wide background.
  RapidMACS now follows that rule. `--keep-dup all` still counts every tag.
- **BED summit scores under `--pvalue`:** summit column 5 now contains
  minus-log10 P rather than minus-log10 q. Q-value mode is unchanged.
- **Capped FRAG counts:** `--frag-pileup-max-count N` (MACS3 `--max-count N`)
  now applies to summit endpoint multiplicities as well as pileups and
  background totals. Previously, uncapped endpoint counts could change
  summit positions and their associated scores.

## Validation

All comparisons below require byte-identical complete narrowPeak and summit
files against newly generated original MACS3 3.0.3 references:

- An 85-row BED fixture passes P/q modes with `--keep-dup 1`, `2`, and `all`.
- A 120-row FRAG fixture passes caps 1, 2, 5, and uncapped mode through the
  CLI and four library entry points. The new regressions fail on v1.0.1.
- A 5,955,396-row BED input matches all 23,932 peaks/summits in P mode and
  all 30,677 in q mode. The same five-column input, parsed as FRAG with cap 1,
  matches all 19,487 peaks/summits.
- All 11 unit tests and all test targets pass, including both 100K smokes.
  Existing 100K, ChIP-seq, and CUT&RUN goldens are unchanged. Fresh FRAG,
  BAMPE, and BEDPE parity checks also pass.

BAMPE/BEDPE treatment and control singleton accounting already matched MACS3
and remains unchanged.

## Installation

The container supports Linux amd64 and arm64:

```sh
docker pull ghcr.io/morphic-bio/rapidmacs:1.0.2
```

Release tarballs contain the executable, static library, headers, changelog,
and documentation. Verify the matching `.sha256` file before unpacking:

```sh
sha256sum -c rapidmacs-v1.0.2-linux-amd64.tar.gz.sha256
tar xzf rapidmacs-v1.0.2-linux-amd64.tar.gz
```

Use the `arm64` tarball on Linux arm64. Source consumers should update their
pinned revision to `v1.0.2` and rebuild. Applications embedding `librapidmacs`
must rebuild and relink; the shared BED writer has an additional optional
argument. Existing source calls and CLI flags remain supported.
