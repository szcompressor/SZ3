#!/usr/bin/env python3
"""Merge the metrics CSVs a workflow run produced into the compression baseline.

The integration jobs each write an integration_metrics.csv and keep it as an artifact. Download
those, point this at them, and commit the result. Do that when a change is meant to move the
compressed size; the driver fails a case whose size has drifted more than the tolerance it
carries, which is what makes an unintended move visible.

Usage:
    python3 update_compression_baseline.py <metrics.csv> [<metrics.csv> ...]
"""

import argparse
import csv
import json
import os

BASELINE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "compression_baseline.json")
NOTE = ("Compressed size in bytes per dataset/field|algo|error bound, from a run on master. "
        "Regenerate with update_compression_baseline.py when a change is meant to move it.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_files", nargs="+")
    parser.add_argument("--baseline", default=BASELINE)
    parser.add_argument("--prune", action="store_true",
                        help="drop entries the given CSVs do not mention")
    args = parser.parse_args()

    try:
        with open(args.baseline) as handle:
            baseline = json.load(handle)
    except FileNotFoundError:
        baseline = {}

    seen, changed, added = set(), 0, 0
    for path in args.csv_files:
        with open(path, newline="") as handle:
            for row in csv.DictReader(handle):
                if not row.get("bytes"):
                    continue
                key = f"{row['dataset']}/{row['field']}|{row['algo']}|{row['error_bound']}"
                value = int(row["bytes"])
                seen.add(key)
                if key not in baseline:
                    added += 1
                elif baseline[key] != value:
                    print(f"{key}: {baseline[key]} -> {value}")
                    changed += 1
                baseline[key] = value

    if args.prune:
        for key in [k for k in baseline if k not in seen and not k.startswith("_")]:
            print(f"dropping {key}")
            del baseline[key]

    baseline["_comment"] = NOTE
    ordered = {"_comment": baseline.pop("_comment")}
    ordered.update(dict(sorted(baseline.items())))
    with open(args.baseline, "w") as handle:
        json.dump(ordered, handle, indent=2)
        handle.write("\n")
    print(f"{len(ordered) - 1} entries in {args.baseline}: {added} added, {changed} changed")


if __name__ == "__main__":
    main()
