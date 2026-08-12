#!/usr/bin/env python3
"""Summarize repeated RapidMACS/MACS3 runs as median and sample range."""

import csv
import statistics
import sys
from collections import defaultdict


def main(raw_path: str, output_path: str) -> None:
    groups = defaultdict(list)
    with open(raw_path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            groups[(row["assay"], row["tool"], int(row["threads"]))].append(row)

    fields = [
        "assay", "tool", "threads", "n", "wall_median_s", "wall_min_s",
        "wall_max_s", "rss_median_kb", "rss_min_kb", "rss_max_kb", "peaks",
        "narrowpeak_md5", "bed3_md5", "summits_md5", "summits_bed3_md5",
    ]
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for (assay, tool, threads), rows in sorted(groups.items()):
            walls = [float(r["wall_s"]) for r in rows]
            rss = [int(r["rss_kb"]) for r in rows]
            for key in (
                "peaks", "narrowpeak_md5", "bed3_md5", "summits_md5",
                "summits_bed3_md5",
            ):
                if len({r[key] for r in rows}) != 1:
                    raise SystemExit(f"non-reproducible {key}: {assay} {tool} T={threads}")
            writer.writerow({
                "assay": assay,
                "tool": tool,
                "threads": threads,
                "n": len(rows),
                "wall_median_s": f"{statistics.median(walls):.3f}",
                "wall_min_s": f"{min(walls):.3f}",
                "wall_max_s": f"{max(walls):.3f}",
                "rss_median_kb": f"{statistics.median(rss):.0f}",
                "rss_min_kb": min(rss),
                "rss_max_kb": max(rss),
                "peaks": rows[0]["peaks"],
                "narrowpeak_md5": rows[0]["narrowpeak_md5"],
                "bed3_md5": rows[0]["bed3_md5"],
                "summits_md5": rows[0]["summits_md5"],
                "summits_bed3_md5": rows[0]["summits_bed3_md5"],
            })


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(f"usage: {sys.argv[0]} RAW.tsv SUMMARY.tsv")
    main(sys.argv[1], sys.argv[2])
