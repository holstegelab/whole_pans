#!/usr/bin/env python3
"""Create one sample/reference PAV analysis configuration."""

import argparse
import csv
import json
import re
from pathlib import Path


SAFE_NAME_RE = re.compile(r"^[A-Za-z0-9_-]+$")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--hap1", required=True)
    parser.add_argument("--hap2", required=True)
    parser.add_argument("--config-output", required=True)
    parser.add_argument("--assemblies-output", required=True)
    parser.add_argument(
        "--extra-config-json",
        default="{}",
        help="JSON object merged into config.json after the reference path",
    )
    args = parser.parse_args()

    if not SAFE_NAME_RE.fullmatch(args.sample):
        raise ValueError(
            "PAV sample names may contain only letters, numbers, underscores, and dashes: "
            + args.sample
        )
    extra = json.loads(args.extra_config_json)
    if not isinstance(extra, dict):
        raise ValueError("--extra-config-json must decode to a JSON object")

    config_path = Path(args.config_output)
    assemblies_path = Path(args.assemblies_output)
    config_path.parent.mkdir(parents=True, exist_ok=True)
    assemblies_path.parent.mkdir(parents=True, exist_ok=True)

    config = {"reference": str(Path(args.reference).resolve())}
    config.update(extra)
    config_path.write_text(json.dumps(config, indent=2, sort_keys=True) + "\n")
    with assemblies_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["NAME", "HAP_h1", "HAP_h2"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow(
            {
                "NAME": args.sample,
                "HAP_h1": str(Path(args.hap1).resolve()),
                "HAP_h2": str(Path(args.hap2).resolve()),
            }
        )


if __name__ == "__main__":
    main()
