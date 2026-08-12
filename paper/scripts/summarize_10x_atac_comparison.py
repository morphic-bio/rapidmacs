#!/usr/bin/env python3
"""Build the supplied-10x ATAC MACS3/RapidMACS comparison table."""

import argparse
import csv
from pathlib import Path


PARITY_FIELDS = (
    "peaks",
    "narrowpeak_md5",
    "bed3_md5",
    "summits_md5",
    "summits_bed3_md5",
)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit(f"empty summary: {path}")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--macs3", required=True, type=Path)
    parser.add_argument("--rapid", required=True, action="append", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    macs_rows = read_rows(args.macs3)
    if len(macs_rows) != 1 or macs_rows[0]["tool"] != "MACS3":
        raise SystemExit("the MACS3 summary must contain exactly one MACS3 row")
    macs = macs_rows[0]
    if int(macs["n"]) != 3:
        raise SystemExit(f"MACS3 reference must have n=3, found n={macs['n']}")

    rapid_rows: list[dict[str, str]] = []
    for path in args.rapid:
        rapid_rows.extend(read_rows(path))
    if {int(row["threads"]) for row in rapid_rows} != {1, 4, 24}:
        raise SystemExit("RapidMACS summaries must contain threads 1, 4, and 24")
    for row in rapid_rows:
        if row["tool"] != "RapidMACS" or int(row["n"]) != 3:
            raise SystemExit("each RapidMACS configuration must have n=3")
        for field in PARITY_FIELDS:
            if row[field] != macs[field]:
                raise SystemExit(f"RapidMACS parity failure at T={row['threads']}: {field}")

    output_rows = [macs, *sorted(rapid_rows, key=lambda row: int(row["threads"]))]
    fields = [
        "dataset",
        "input_format",
        "tool",
        "threads",
        "n",
        "wall_median_s",
        "wall_min_s",
        "wall_max_s",
        "rss_median_kb",
        "peaks",
        "speedup_vs_macs3",
        "exact_peak_and_summit_parity",
        "narrowpeak_md5",
        "summits_md5",
    ]
    macs_median = float(macs["wall_median_s"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in output_rows:
            writer.writerow(
                {
                    "dataset": "10x PBMC 3k supplied fragments",
                    "input_format": "FRAG (5-column)",
                    "tool": row["tool"],
                    "threads": row["threads"],
                    "n": row["n"],
                    "wall_median_s": row["wall_median_s"],
                    "wall_min_s": row["wall_min_s"],
                    "wall_max_s": row["wall_max_s"],
                    "rss_median_kb": row["rss_median_kb"],
                    "peaks": row["peaks"],
                    "speedup_vs_macs3": (
                        f"{macs_median / float(row['wall_median_s']):.2f}"
                    ),
                    "exact_peak_and_summit_parity": "yes",
                    "narrowpeak_md5": row["narrowpeak_md5"],
                    "summits_md5": row["summits_md5"],
                }
            )


if __name__ == "__main__":
    main()
