#!/usr/bin/env python3
"""Resolve fixed TMEM pairs and select phenotype-blind random samples."""

from __future__ import annotations

import argparse
import csv
import random
from collections import Counter, defaultdict
from pathlib import Path

import yaml


FIELDS = [
    "sample_id", "role", "pair_id", "expected_allele",
    "expectation_basis", "batch", "r1", "r2", "coverage_proxy",
    "random_seed", "selection_reason",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def load_manifest(path: Path, columns: dict[str, int]) -> dict[str, dict[str, str]]:
    result = {}
    with path.open(newline="") as handle:
        for line_number, fields in enumerate(csv.reader(handle, delimiter="\t"), start=1):
            if not fields or all(not value.strip() for value in fields):
                continue
            def value(name: str) -> str:
                index = int(columns[name]) - 1
                return fields[index].strip() if index < len(fields) else ""

            sample = value("sample")
            data_type = value("data_type")
            if not sample or data_type != "fastq_paired":
                continue
            if sample in result:
                raise ValueError(f"duplicate sample {sample} across manifests (latest {path}:{line_number})")
            result[sample] = {
                "batch": value("batch"),
                "r1": value("r1"),
                "r2": value("r2"),
            }
    return result


def read_coverage(path: Path, column: str) -> dict[str, float]:
    if not path.is_file():
        return {}
    result = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sample = row.get("sample") or row.get("sample_id")
            try:
                result[str(sample)] = float(row.get(column, ""))
            except (TypeError, ValueError):
                pass
    return result


def allocate_quotas(batch_counts: Counter, total: int) -> dict[str, int]:
    batches = sorted(batch_counts)
    if not batches or total <= 0:
        return {}
    weight = sum(batch_counts.values())
    exact = {batch: total * batch_counts[batch] / weight for batch in batches}
    quota = {batch: int(exact[batch]) for batch in batches}
    remaining = total - sum(quota.values())
    for batch in sorted(batches, key=lambda key: (-(exact[key] - quota[key]), key))[:remaining]:
        quota[batch] += 1
    return quota


def main() -> None:
    args = parse_args()
    with args.config.open() as handle:
        config = yaml.safe_load(handle)
    pilot = config["tmem_pilot"]
    columns = pilot["manifest_columns"]
    available = {}
    for manifest in pilot["fastq_manifests"]:
        available.update(load_manifest(Path(manifest), columns))

    coverage = read_coverage(
        Path(pilot.get("coverage_proxy_table", "")),
        pilot.get("coverage_proxy_column", "mapped_records"),
    )
    fixed_rows = []
    fixed_samples = set()
    batch_counts: Counter = Counter()
    for pair_number, pair in enumerate(pilot["pairs"], start=1):
        if len(pair) != 2:
            raise ValueError(f"TMEM pair {pair_number} must contain candidate and control")
        for sample, role in zip(pair, ("tmem_candidate", "tmem_matched_control")):
            if sample not in available:
                raise ValueError(f"required TMEM pilot sample is absent from FASTQ manifests: {sample}")
            fixed_samples.add(sample)
            batch = available[sample]["batch"] or sample.split("_", 3)[2]
            batch_counts[batch] += 1
            fixed_rows.append({
                "sample_id": sample,
                "role": role,
                "pair_id": f"pair_{pair_number:02d}",
                "expected_allele": (
                    "insertion_present" if role == "tmem_candidate" else "insertion_absent"
                ),
                "expectation_basis": "prior_strict_junction_screen_not_independent_truth",
                "batch": batch,
                "r1": available[sample]["r1"],
                "r2": available[sample]["r2"],
                "coverage_proxy": coverage.get(sample, ""),
                "random_seed": "",
                "selection_reason": "predefined_candidate_control_pair",
            })

    seed = int(pilot.get("random_seed", 1))
    rng = random.Random(seed)
    random_count = int(pilot.get("random_blinded_count", 4))
    quotas = allocate_quotas(batch_counts, random_count)
    by_batch: dict[str, list[str]] = defaultdict(list)
    for sample, row in available.items():
        if sample not in fixed_samples:
            batch = row["batch"] or sample.split("_", 3)[2]
            by_batch[batch].append(sample)

    fixed_coverage_by_batch: dict[str, list[float]] = defaultdict(list)
    for row in fixed_rows:
        value = coverage.get(row["sample_id"])
        if value is not None:
            fixed_coverage_by_batch[row["batch"]].append(value)

    random_rows = []
    selected = set()
    for batch, quota in quotas.items():
        values = fixed_coverage_by_batch.get(batch, [])
        target = sorted(values)[len(values) // 2] if values else None
        candidates = list(by_batch.get(batch, []))
        rng.shuffle(candidates)
        candidates.sort(key=lambda sample: (
            abs(coverage[sample] - target) if target is not None and sample in coverage else float("inf"),
        ))
        for sample in candidates[:quota]:
            selected.add(sample)
            random_rows.append({
                "sample_id": sample,
                "role": "tmem_random_blinded",
                "pair_id": "",
                "expected_allele": "blinded",
                "expectation_basis": "phenotype_blind_no_expected_allele",
                "batch": batch,
                "r1": available[sample]["r1"],
                "r2": available[sample]["r2"],
                "coverage_proxy": coverage.get(sample, ""),
                "random_seed": seed,
                "selection_reason": "phenotype_blind_batch_stratified_coverage_matched_when_available",
            })

    if len(random_rows) < random_count:
        remaining = [sample for sample in available if sample not in fixed_samples | selected]
        rng.shuffle(remaining)
        target_values = [value for value in coverage.values()]
        target = sorted(target_values)[len(target_values) // 2] if target_values else None
        remaining.sort(key=lambda sample: (
            abs(coverage[sample] - target) if target is not None and sample in coverage else float("inf"),
        ))
        for sample in remaining[: random_count - len(random_rows)]:
            batch = available[sample]["batch"] or sample.split("_", 3)[2]
            random_rows.append({
                "sample_id": sample,
                "role": "tmem_random_blinded",
                "pair_id": "",
                "expected_allele": "blinded",
                "expectation_basis": "phenotype_blind_no_expected_allele",
                "batch": batch,
                "r1": available[sample]["r1"],
                "r2": available[sample]["r2"],
                "coverage_proxy": coverage.get(sample, ""),
                "random_seed": seed,
                "selection_reason": "phenotype_blind_fallback_coverage_matched_when_available",
            })
    if len(random_rows) != random_count:
        raise ValueError(f"could select only {len(random_rows)} of {random_count} blinded samples")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(args.output) + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(fixed_rows + sorted(random_rows, key=lambda row: (row["batch"], row["sample_id"])))
    temporary.replace(args.output)


if __name__ == "__main__":
    main()
