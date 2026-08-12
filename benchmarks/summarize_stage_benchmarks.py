#!/usr/bin/env python3
"""Summarize paired CLI-ingestion and in-memory library-pipeline timings."""

import csv
import statistics
import sys
from collections import defaultdict


METRICS = (
    "load_s",
    "library_pipeline_s",
    "build_workspace_s",
    "treat_s",
    "lambda_s",
    "ppois_regions_s",
    "output_s",
    "profiled_total_s",
)


def main(raw_path: str, output_path: str) -> None:
    groups = defaultdict(list)
    with open(raw_path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            groups[(row["assay"], int(row["threads"]))].append(row)

    fields = ["assay", "threads", "n"]
    for metric in METRICS:
        stem = metric[:-2] if metric.endswith("_s") else metric
        fields.extend((f"{stem}_median_s", f"{stem}_min_s", f"{stem}_max_s"))

    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for (assay, threads), rows in sorted(groups.items()):
            out = {"assay": assay, "threads": threads, "n": len(rows)}
            for metric in METRICS:
                values = [float(row[metric]) for row in rows]
                stem = metric[:-2] if metric.endswith("_s") else metric
                out[f"{stem}_median_s"] = f"{statistics.median(values):.3f}"
                out[f"{stem}_min_s"] = f"{min(values):.3f}"
                out[f"{stem}_max_s"] = f"{max(values):.3f}"
            writer.writerow(out)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(f"usage: {sys.argv[0]} RAW.tsv SUMMARY.tsv")
    main(sys.argv[1], sys.argv[2])
