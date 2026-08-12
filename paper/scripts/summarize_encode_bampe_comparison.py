#!/usr/bin/env python3
"""Build the matched ENCSR095QNB BAMPE comparison table."""

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


def read_summary(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit(f"empty summary: {path}")
    return rows


def add_comparison(
    output_rows: list[dict[str, str]],
    input_label: str,
    macs_path: Path,
    rapid_paths: list[Path],
) -> None:
    macs_rows = read_summary(macs_path)
    if len(macs_rows) != 1 or macs_rows[0]["tool"] != "MACS3":
        raise SystemExit(f"expected one MACS3 row: {macs_path}")
    macs = macs_rows[0]
    rows = [macs]
    for path in rapid_paths:
        rows.extend(read_summary(path))

    seen_threads: set[int] = set()
    for row in rows[1:]:
        if row["tool"] != "RapidMACS":
            raise SystemExit(f"expected RapidMACS row in comparison for {input_label}")
        threads = int(row["threads"])
        if threads in seen_threads:
            raise SystemExit(f"duplicate RapidMACS T={threads} for {input_label}")
        seen_threads.add(threads)
        for field in PARITY_FIELDS:
            if row[field] != macs[field]:
                raise SystemExit(
                    f"parity failure for {input_label} T={threads}: {field}"
                )

    if seen_threads != {1, 4, 24}:
        raise SystemExit(
            f"expected RapidMACS threads 1, 4, 24 for {input_label}; "
            f"found {sorted(seen_threads)}"
        )

    macs_median = float(macs["wall_median_s"])
    for row in sorted(rows, key=lambda item: (item["tool"] != "MACS3", int(item["threads"]))):
        output_rows.append(
            {
                "input": input_label,
                "tool": row["tool"],
                "threads": row["threads"],
                "n": row["n"],
                "wall_median_s": row["wall_median_s"],
                "wall_min_s": row["wall_min_s"],
                "wall_max_s": row["wall_max_s"],
                "rss_median_kb": row["rss_median_kb"],
                "peaks": row["peaks"],
                "speedup_vs_macs3": f"{macs_median / float(row['wall_median_s']):.2f}",
                "exact_peak_and_summit_parity": "yes",
                "narrowpeak_md5": row["narrowpeak_md5"],
                "summits_md5": row["summits_md5"],
            }
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--original-macs3", required=True, type=Path)
    parser.add_argument("--original-rapid", required=True, action="append", type=Path)
    parser.add_argument("--chromap-macs3", required=True, type=Path)
    parser.add_argument("--chromap-rapid", required=True, action="append", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    output_rows: list[dict[str, str]] = []
    add_comparison(
        output_rows,
        "original_ENCODE_ENCFF646NWY",
        args.original_macs3,
        args.original_rapid,
    )
    add_comparison(
        output_rows,
        "Chromap_FASTQ_realignment",
        args.chromap_macs3,
        args.chromap_rapid,
    )

    fields = list(output_rows[0])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)


if __name__ == "__main__":
    main()
