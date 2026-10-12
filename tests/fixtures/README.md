# Paired-end ChIP-seq and CUT&RUN fixtures

The binary fixtures in this directory are small, single-locus subsets of public
paired-end datasets.  The full files live outside the repository under
`$DATA_DIR` and are used for the paper benchmarks.
Both MACS3 and RapidMACS consume the same subset BAM, so alignment is not part of
the peak-caller timing or parity comparison.

## ChIP-seq: ENCODE CTCF with matched input

- treatment experiment: ENCSR634OAQ, CTCF ChIP-seq in NCI-H929
- treatment alignment: ENCFF686KKV, biological replicate 1, GRCh38,
  paired-ended, md5 `bedc2707c0be9f6d61c09edf811d655e`
- declared control experiment: ENCSR409FGC, input-library control in NCI-H929
- control alignment: ENCFF181CXT, biological replicate 1, GRCh38,
  paired-ended, md5 `ffec8961e15783bafb161fb7ab534568`
- test interval: `chr21:20000000-30000000`

The accessions and treatment/control relationship come from the ENCODE
experiment metadata.  Recreate the full inputs and fixture with:

```bash
data=$DATA_DIR/chipseq/ENCSR634OAQ_ENCSR409FGC
mkdir -p "$data"
curl -L 'https://encode-public.s3.amazonaws.com/2021/02/02/e22433a5-541c-4892-9c00-890b62bda1a3/ENCFF686KKV.bam' \
  -o "$data/ENCFF686KKV.bam"
curl -L 'https://encode-public.s3.amazonaws.com/2020/10/21/9aee9c46-df07-4bc1-93b7-9f91bbc0614b/ENCFF181CXT.bam' \
  -o "$data/ENCFF181CXT.bam"
printf '%s  %s\n' \
  bedc2707c0be9f6d61c09edf811d655e "$data/ENCFF686KKV.bam" \
  ffec8961e15783bafb161fb7ab534568 "$data/ENCFF181CXT.bam" | md5sum -c -
samtools index -@ 8 "$data/ENCFF686KKV.bam"
samtools index -@ 8 "$data/ENCFF181CXT.bam"
samtools view -bh "$data/ENCFF686KKV.bam" chr21:20000000-30000000 \
  -o chipseq_ENCFF686KKV.chr21_20m_30m.bam
samtools view -bh "$data/ENCFF181CXT.bam" chr21:20000000-30000000 \
  -o chipseq_ENCFF181CXT.chr21_20m_30m.bam
```

## CUT&RUN: GEO K562 CTCF with matched IgG

- series: GSE172236
- treatment: GSM5244531 / SRX10617507 / SRR14255054, K562 sgLUC CTCF
  CUT&RUN, 150-bp paired-end
- control: GSM5244530 / SRX10617506 / SRR14255053, K562 IgG control,
  150-bp paired-end
- ENA FASTQ md5 (treatment R1/R2): `c3a10f1dbadba9e21bd2d6006d968bb2`,
  `110e63a1298591f5cad84b942bcebb0a`
- ENA FASTQ md5 (control R1/R2): `ce1b6be737b06fc0fc54ad253cddf507`,
  `9b38a4c0a41d6ec7dcfb90ab83127ee5`
- test interval: `chr21:20000000-30000000`

The GEO family record identifies the IgG/CTCF pair.  ENA supplies the original
FASTQs.  The commands below align both samples identically with Chromap's
ChIP-seq mapping parameters, without its preset's upstream PCR-deduplication;
MACS3's default coordinate deduplication is therefore tested directly.

```bash
data=$DATA_DIR/cutrun/GSE172236
ref=$DATA_DIR/refdata-gex-GRCh38-2024-A/fasta/genome.fa
index=$DATA_DIR/pbmc_unsorted_3k_100k/chromap_index/genome.index
mkdir -p "$data"
for run in SRR14255053 SRR14255054; do
  for mate in 1 2; do
    curl -L "https://ftp.sra.ebi.ac.uk/vol1/fastq/${run:0:6}/${run: -3}/${run}/${run}_${mate}.fastq.gz" \
      -o "$data/${run}_${mate}.fastq.gz"
  done
  chromap -x "$index" -r "$ref" \
    -1 "$data/${run}_1.fastq.gz" -2 "$data/${run}_2.fastq.gz" \
    -e 8 -s 2 -f 500,1000 -l 2000 -q 30 --min-read-length 30 \
    --SAM -t 24 -o "$data/${run}.sam"
  samtools sort -@ 8 -o "$data/${run}.bam" "$data/${run}.sam"
  samtools index -@ 8 "$data/${run}.bam"
done
samtools view -bh "$data/SRR14255054.bam" chr21:20000000-30000000 \
  -o cutrun_SRR14255054.chr21_20m_30m.bam
samtools view -bh "$data/SRR14255053.bam" chr21:20000000-30000000 \
  -o cutrun_SRR14255053.chr21_20m_30m.bam
```

## Golden generation

Goldens are generated with MACS3 3.0.3, without `--call-summits` (the default
one-summit-per-peak path implemented by RapidMACS):

```bash
macs3 callpeak \
  -t $DATA_DIR/pbmc_unsorted_3k_100k/fragments.tsv \
  -f FRAG -g hs -p 1e-5 --min-length 200 --max-gap 30 -n atac_full \
  --outdir tests/golden/100k_macs3_frag
macs3 callpeak -t tests/fixtures/chipseq_ENCFF686KKV.chr21_20m_30m.bam \
  -c tests/fixtures/chipseq_ENCFF181CXT.chr21_20m_30m.bam \
  -f BAMPE -g hs -q 0.05 -n chip_chr21 \
  --outdir tests/golden/chipseq_ENCFF686KKV
macs3 callpeak -t tests/fixtures/cutrun_SRR14255054.chr21_20m_30m.bam \
  -c tests/fixtures/cutrun_SRR14255053.chr21_20m_30m.bam \
  -f BAMPE -g hs -q 0.05 -n cutrun_chr21 \
  --outdir tests/golden/cutrun_SRR14255054
```

The smoke tests require byte-identical narrowPeak and summit files.  The full
benchmark harness applies the same full-file MD5 parity gate to every measured
run before it writes the timing summary.

## Single-end singleton strands

`singleton_strands.bed` is a synthetic 85-row BED6 fixture: three dense peaks,
singletons on both strands (including a singleton opposite a dense strand),
and a two-copy duplicate run on its own chromosome. MACS3 3.0.3's filtered
background totals are 73 with `--keep-dup 1`, 77 with `2`, and 85 with `all`.
Singleton tags remain in pileups; a duplicate run reduced to one still counts.

`tests/run_singleton_parity.sh` checks all columns of narrowPeak and summits
against `tests/golden/singleton_strands/`. It uses `--nomodel --shift -100
--extsize 200 -g 10000`, P cutoff `1e-5` or q cutoff `0.05`, and name `singleton`.
The small effective genome makes the filtered total affect peak scores.
To revalidate with an original MACS3 3.0.3 executable (or a container wrapper):

```sh
MACS3=/path/to/macs3 tests/run_singleton_parity.sh
MACS3=/path/to/macs3 UPDATE_GOLDENS=1 tests/run_singleton_parity.sh
```

Reference generation checks the MACS3 version and copies only its output;
RapidMACS output is never used to generate goldens. `RAPIDMACS` can select a
different candidate executable, for example to demonstrate failure on v1.0.1.

## Fragment count caps

`capped_counts.fragments.tsv` contains 120 synthetic fragments on two
chromosomes, with counts from 1 through 15 and three separated peak regions.
Overlapping fragments create tied summit plateaus whose selection depends on
endpoint multiplicity. At caps 1 and 2, v1.0.1 selects a different summit for
one of the three peaks. The uncapped case protects existing behavior.

Original MACS3 3.0.3 generated `tests/golden/capped_counts/{0,1,2,5}/` with
`callpeak -f FRAG -g 2913022398 -p 1e-5 --min-length 200 --max-gap 30
--max-count N -n capped`. Here 0 disables the cap. Revalidate or regenerate:

```sh
MACS3=/path/to/macs3 tests/run_frag_max_count_parity.sh
MACS3=/path/to/macs3 UPDATE_GOLDENS=1 tests/run_frag_max_count_parity.sh
```

The CLI smoke compares both complete files. `test_frag_max_count.cc` also
compares those goldens through the in-memory fragment pipeline, event iterator,
sorted iterator/batch merge, and direct sweep workspace, using two threads.
