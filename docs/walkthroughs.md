# Walkthroughs

Five worked examples, from a first run to embedding the library. Each states
what you should see, so you can tell whether it worked.

- [1. Call peaks on a 10x fragments file](#1-call-peaks-on-a-10x-fragments-file)
- [2. ChIP-seq with a matched input control](#2-chip-seq-with-a-matched-input-control)
- [3. Prove byte-identity yourself](#3-prove-byte-identity-yourself)
- [4. Drop RapidMACS into an existing MACS3 workflow](#4-drop-rapidmacs-into-an-existing-macs3-workflow)
- [5. Call peaks from your own C++ program](#5-call-peaks-from-your-own-c-program)

---

## 1. Call peaks on a 10x fragments file

The common single-cell ATAC case: Cell Ranger ARC gave you a fragments file and
you want narrow peaks instead of ARC's wider ones.

```sh
make -j

curl -fLO https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_unsorted_3k/pbmc_unsorted_3k_atac_fragments.tsv.gz

bin/rapidmacs --preset shiftTAG -f FRAG \
  -t pbmc_unsorted_3k_atac_fragments.tsv.gz \
  -g hs -n atac -p 1e-5 \
  --bdgpeakcall-min-len 200 --bdgpeakcall-max-gap 30 \
  --peak-caller-threads 8 \
  --macs3-frag-narrowpeak atac.narrowPeak \
  --macs3-frag-summits atac.summits.bed
```

You should get **50,003 peaks**:

```sh
$ wc -l < atac.narrowPeak
50003
$ md5sum atac.narrowPeak atac.summits.bed
2c0de7c506514b0afef1e4339c386c9f  atac.narrowPeak
c3805d1d3744ed78e3d589a523ad0b09  atac.summits.bed
```

Those are the same MD5 sums MACS3 v3.0.3 produces from the same input, and the
same ones the paper reports. If yours differ, something is wrong — say so in an
issue rather than working around it.

The equivalent MACS3 command, for reference:

```sh
macs3 callpeak -t pbmc_unsorted_3k_atac_fragments.tsv.gz -f FRAG \
  -g hs -n atac -p 1e-5 --min-length 200 --max-gap 30 --outdir macs3_out
```

On the paper's host that takes 747.55 s. RapidMACS at 8 threads takes roughly
20 s.

### Why `--preset shiftTAG`

It sets the score cutoff MACS3 uses for `-p 0.01`. Here `-p 1e-5` follows it
and wins, because presets are applied first and explicit flags override them.
Putting `--preset` after `-p` would silently discard your threshold.

---

## 2. ChIP-seq with a matched input control

```sh
bin/rapidmacs --preset chip \
  -t ENCFF686KKV.bam -c ENCFF181CXT.bam \
  -g hs -n ctcf \
  --peak-caller-threads 8 \
  --out-prefix out/ctcf
```

`--preset chip` selects BAMPE input, `-q 0.05`, and MACS3's dynamic paired-end
peak geometry. It requires `-c`; without a control it exits with an error
rather than quietly calling peaks against a genome-wide background.

For CUT&RUN use `--preset cutrun`, which is the same but treats the control as
optional:

```sh
bin/rapidmacs --preset cutrun \
  -t SRR14255054.bam -c SRR14255053.bam \
  -g hs -n ctcf_cutrun --out-prefix out/ctcf_cutrun
```

Both write `<prefix>.narrowPeak` and companion files, matching

```sh
macs3 callpeak -t treat.bam -c control.bam -f BAMPE -g hs -n NAME -q 0.05 --outdir OUT
```

BAM input must be coordinate-sorted. Duplicate handling follows MACS3:
`--keep-dup 1` by default, `--keep-dup all` to retain everything.

### One thing to know about MACS3 here

Released MACS3 v3.0.3 writes its internal control pileup to a file named
`tmp_b_c.txt` in the working directory on every run that has a control. On the
ENCODE ChIP-seq dataset that write costs 207.6 s of a 616.90 s run. It is
harmless but it is why the ChIP-seq speedup in the paper is large. See
[`paper/README.md`](../paper/README.md).

---

## 3. Prove byte-identity yourself

The claim is that the files are the same bytes, so check it directly rather
than taking our word for it.

```sh
paper/scripts/setup_macs3_reference.sh VARIANT=released
MACS3=.macs3-reference/macs3-3.0.3-venv/bin/macs3

# same input, both callers
$MACS3 callpeak -t fragments.tsv.gz -f FRAG -g hs -n cmp \
  -p 1e-5 --min-length 200 --max-gap 30 --outdir macs3_out

bin/rapidmacs --preset shiftTAG -f FRAG -t fragments.tsv.gz -g hs -n cmp \
  -p 1e-5 --bdgpeakcall-min-len 200 --bdgpeakcall-max-gap 30 \
  --macs3-frag-narrowpeak rapid.narrowPeak \
  --macs3-frag-summits rapid.summits.bed

cmp macs3_out/cmp_peaks.narrowPeak rapid.narrowPeak && echo "narrowPeak identical"
cmp macs3_out/cmp_summits.bed      rapid.summits.bed && echo "summits identical"
```

`cmp` is the right tool: it compares bytes, not peak sets within a tolerance.

For the full published comparison across all six configurations, including
timings, use [`paper/scripts/run_10x_pbmc_example.sh`](../paper/scripts/run_10x_pbmc_example.sh),
which does the download, the checksum verification, both callers, and the
comparison, and exits non-zero if anything disagrees.

---

## 4. Drop RapidMACS into an existing MACS3 workflow

Signac, ArchR, and SCENIC+/pycisTopic all shell out to a peak caller. Because
the output files are identical, you can substitute RapidMACS part-way through
an analysis without invalidating comparisons against work already done.

ArchR and pycisTopic both call with `--nolambda --call-summits`. The BED form:

```sh
bin/rapidmacs --macs-profile signac-atac \
  --macs-bed-input tags.bed --macs-input-format BED \
  --nolambda --call-summits \
  -n pseudobulk \
  --macs-bed-peaks-out pseudobulk.narrowPeak \
  --macs-bed-summits-out pseudobulk.summits.bed
```

`--macs-profile signac-atac` resolves the whole Signac parameter set —
`--nomodel --extsize 200 --shift -100 -g 2700000000 -q 0.05 --llocal 10000` —
so you do not have to keep those constants in sync by hand.

If you start from fragments rather than BED, project them first, optionally
splitting by a barcode-to-group table so each pseudobulk gets its own file:

```sh
bin/rapidmacs --input-fragments fragments.tsv.gz \
  --barcode-group-tsv cell_groups.tsv --group-column cluster \
  --macs-bed-out-dir beds/
```

Verified command lines and checksums for all three tools are in
[`paper/WORKFLOW-COMPATIBILITY.md`](../paper/WORKFLOW-COMPATIBILITY.md). Note
that all three default to `macs2`; the parity claim is against `macs3` v3.0.3.

---

## 5. Call peaks from your own C++ program

The reason the library exists: an aligner that has just produced fragments can
call peaks on them in-process, instead of writing a temporary file and starting
a Python process to read it back. On the single-cell dataset above that avoids
moving 2.3 GB through the filesystem.

```cpp
// peaks_from_memory.cc
#include <cstdio>
#include <string>
#include <vector>

#include "rapidmacs/fragments.h"
#include "rapidmacs/macs3_frag_peak_pipeline.h"

int main() {
  // Fragments your program already holds. chrom_id indexes chrom_names.
  // FragmentRecord has default member initialisers, so it is not a C++11
  // aggregate: fill the fields rather than brace-initialising it.
  std::vector<macs3::FragmentRecord> fragments;
  auto add = [&fragments](int32_t chrom_id, int32_t start, int32_t end,
                          uint32_t count) {
    macs3::FragmentRecord f;
    f.chrom_id = chrom_id;
    f.start = start;
    f.end = end;
    f.count = count;
    fragments.push_back(f);
  };
  add(0, 1000, 1200, 1);
  add(0, 1050, 1300, 2);
  add(1, 5000, 5250, 1);

  std::vector<std::string> chrom_names;
  chrom_names.push_back("chr1");
  chrom_names.push_back("chr2");

  // Sorts on construction, so the sorted-iterator contract holds.
  auto iter = macs3::WrapVectorFragmentIterator(std::move(fragments),
                                                std::move(chrom_names));

  chromap::peaks::Macs3FragPeakPipelineParams params;
  params.name_prefix = "inmem";
  params.peak_caller_threads = 4;

  std::string work_dir, error;
  const bool ok = chromap::peaks::RunMacs3FragPeakPipelineFromSortedIterator(
      *iter, params, chromap::peaks::Macs3FragPeakPipelinePaths{},
      "inmem.narrowPeak", "inmem.summits.bed",
      /*keep_intermediates_dir=*/"", /*work_dir_parent=*/"",
      &work_dir, &error);
  if (!ok) {
    std::fprintf(stderr, "rapidmacs: %s\n", error.c_str());
    return 1;
  }
  return 0;
}
```

```sh
g++ -std=c++11 -Iinclude peaks_from_memory.cc lib/librapidmacs.a \
  -lhts -lm -lz -lpthread -fopenmp -o peaks_from_memory
```

To stream instead of buffering — the case that actually matters at scale —
implement `macs3::FragmentIterator` over your own data and hand it to the same
function. Emit fragments chromosome-grouped and start-sorted and peak memory
stays in the tens of kilobytes regardless of input size. If your order is
wrong, the call returns `false` and says so rather than producing wrong peaks.

Full interface and parameter reference: [api.md](api.md).
