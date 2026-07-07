#!/usr/bin/env python3
"""Combine assembly metrics and classify assemblies for graph inclusion."""

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import yaml


ASSEMBLY_RE = re.compile(
    r"^(?P<sample>.+?)\.hifi\.hifiasm\.bp\.(?P<haplotype>hap[12])\.p_ctg(?:\..*)?$"
)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--assembly-output", required=True)
    parser.add_argument("--sample-output", required=True)
    parser.add_argument("--included-output", required=True)
    parser.add_argument("--excluded-output", required=True)
    parser.add_argument(
        "--allow-missing-mates",
        action="store_true",
        help=(
            "Keep a passing haplotype when its mate is absent from the input manifest. "
            "Use this for post-decontamination QC of assemblies already selected by QC."
        ),
    )
    return parser.parse_args()


def read_one_row(path):
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"Expected one data row in {path}, found {len(rows)}")
    return rows[0]


def number(value):
    return float(str(value).replace(",", ""))


def parse_seqkit(path):
    row = read_one_row(path)

    def get(*names):
        for name in names:
            if name in row:
                return number(row[name])
        raise KeyError(f"None of {names} found in {path}; columns are {sorted(row)}")

    total = get("sum_len")
    n_bases = get("sum_n")
    return {
        "contig_count": int(get("num_seqs")),
        "total_length_bp": int(total),
        "contig_n50_bp": int(get("N50")),
        "largest_contig_bp": int(get("max_len")),
        "gc_percent": get("GC(%)", "GC"),
        "n_bases": int(n_bases),
        "n_percent": 100.0 * n_bases / total if total else 100.0,
    }


def parse_compleasm(path):
    values = {}
    with open(path) as handle:
        for line in handle:
            match = re.match(r"^\s*([SDFIMN]):\s*([0-9.]+)%?", line)
            if match:
                values[match.group(1)] = float(match.group(2))
    required = {"S", "D", "F", "I", "M"}
    if not required.issubset(values):
        raise ValueError(f"Could not parse compleasm summary {path}; found categories {sorted(values)}")
    return {
        "compleasm_single_percent": values["S"],
        "compleasm_duplicated_percent": values["D"],
        "compleasm_complete_percent": values["S"] + values["D"],
        "compleasm_fragmented_percent": values["F"] + values["I"],
        "compleasm_missing_percent": values["M"],
    }


def parse_alignment(path, prefix):
    row = read_one_row(path)
    return {
        f"{prefix}_query_aligned_percent": number(row["query_aligned_percent"]),
        f"{prefix}_reference_covered_percent": number(row["reference_covered_percent"]),
        f"{prefix}_alignment_identity_percent": number(row["alignment_identity_percent"]),
        f"{prefix}_alignment_count": int(number(row["alignment_count"])),
    }


def reason(metric, operator, threshold, value):
    return f"{metric}{operator}{threshold:g} (observed {value:.3f})"


def classify(row, thresholds):
    fail = thresholds["fail"]
    warn = thresholds["warn"]
    failures = []
    warnings = []

    checks = [
        ("total_length_bp", ">=", "min_total_length_bp"),
        ("total_length_bp", "<=", "max_total_length_bp"),
        ("contig_n50_bp", ">=", "min_contig_n50_bp"),
        ("n_percent", "<=", "max_n_percent"),
        ("compleasm_complete_percent", ">=", "min_compleasm_complete_percent"),
        ("compleasm_duplicated_percent", "<=", "max_compleasm_duplicated_percent"),
        ("best_query_aligned_percent", ">=", "min_best_query_aligned_percent"),
        ("best_reference_covered_percent", ">=", "min_best_reference_covered_percent"),
        ("best_alignment_identity_percent", ">=", "min_best_alignment_identity_percent"),
    ]

    for metric, operator, threshold_name in checks:
        value = row[metric]
        hard_limit = float(fail[threshold_name])
        warn_limit = float(warn[threshold_name])
        hard_bad = value < hard_limit if operator == ">=" else value > hard_limit
        warn_bad = value < warn_limit if operator == ">=" else value > warn_limit
        if hard_bad:
            failures.append(reason(metric, operator, hard_limit, value))
        elif warn_bad:
            warnings.append(reason(metric, operator, warn_limit, value))

    if failures:
        return "FAIL", "; ".join(failures), "; ".join(warnings)
    if warnings:
        return "WARN", "", "; ".join(warnings)
    return "PASS", "", ""


def write_tsv(path, rows, fields):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    with open(args.config) as handle:
        config = yaml.safe_load(handle)
    thresholds = config["thresholds"]
    results = Path(args.results_dir)

    with open(args.manifest, newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))

    assembly_rows = []
    by_sample = defaultdict(list)
    for entry in manifest:
        assembly_id = entry["assembly_id"]
        match = ASSEMBLY_RE.match(assembly_id)
        if not match:
            raise ValueError(
                f"Assembly ID does not match expected hifiasm naming pattern: {assembly_id}"
            )

        row = {
            "assembly_id": assembly_id,
            "sample": match.group("sample"),
            "haplotype": match.group("haplotype"),
            "path": entry["path"],
        }
        row.update(parse_seqkit(results / "stats" / f"{assembly_id}.seqkit.tsv"))
        row.update(parse_compleasm(results / "compleasm" / assembly_id / "summary.txt"))
        row.update(
            parse_alignment(
                results / "alignment_metrics" / "CHM13" / f"{assembly_id}.tsv", "chm13"
            )
        )
        row.update(
            parse_alignment(
                results / "alignment_metrics" / "hg38" / f"{assembly_id}.tsv", "hg38"
            )
        )
        row["best_query_aligned_percent"] = max(
            row["chm13_query_aligned_percent"], row["hg38_query_aligned_percent"]
        )
        row["best_reference_covered_percent"] = max(
            row["chm13_reference_covered_percent"], row["hg38_reference_covered_percent"]
        )
        row["best_alignment_identity_percent"] = max(
            row["chm13_alignment_identity_percent"], row["hg38_alignment_identity_percent"]
        )
        row["assembly_status"], row["fail_reasons"], row["warning_reasons"] = classify(
            row, thresholds
        )
        assembly_rows.append(row)
        by_sample[row["sample"]].append(row)

    sample_rows = []
    included_paths = []
    excluded_rows = []
    for sample in sorted(by_sample):
        haplotypes = by_sample[sample]
        hap_by_name = {row["haplotype"]: row for row in haplotypes}
        pair_errors = []
        hap_counts = defaultdict(int)
        for row in haplotypes:
            hap_counts[row["haplotype"]] += 1
        missing = sorted({"hap1", "hap2"} - set(hap_by_name))
        if missing and not args.allow_missing_mates:
            pair_errors.append("missing " + ",".join(missing))
        duplicated = sorted(hap for hap, count in hap_counts.items() if count > 1)
        if duplicated:
            pair_errors.append("duplicate assembly for " + ",".join(duplicated))

        failed = [row for row in haplotypes if row["assembly_status"] == "FAIL"]
        warned = [row for row in haplotypes if row["assembly_status"] == "WARN"]
        sample_reasons = list(pair_errors)
        if missing and args.allow_missing_mates:
            sample_reasons.append(
                "mate absent from post-decontamination input: " + ",".join(missing)
            )
        for row in failed:
            sample_reasons.append(f"{row['haplotype']}: {row['fail_reasons']}")

        if pair_errors or len(failed) == len(haplotypes):
            sample_status = "FAIL"
        elif failed or missing:
            sample_status = "PARTIAL"
        elif warned:
            sample_status = "WARN"
        else:
            sample_status = "PASS"

        sample_rows.append(
            {
                "sample": sample,
                "sample_status": sample_status,
                "hap1_status": hap_by_name.get("hap1", {}).get("assembly_status", "MISSING"),
                "hap2_status": hap_by_name.get("hap2", {}).get("assembly_status", "MISSING"),
                "reasons": "; ".join(sample_reasons),
            }
        )

        if pair_errors:
            excluded_haplotypes = haplotypes
            included_haplotypes = []
        else:
            excluded_haplotypes = failed
            included_haplotypes = [
                row for row in haplotypes if row["assembly_status"] != "FAIL"
            ]

        for row in sorted(excluded_haplotypes, key=lambda item: item["haplotype"]):
            excluded_rows.append(
                {
                    "sample": sample,
                    "assembly_id": row["assembly_id"],
                    "path": row["path"],
                    "sample_status": sample_status,
                    "reasons": "; ".join(sample_reasons),
                }
            )
        included_paths.extend(
            row["path"]
            for row in sorted(included_haplotypes, key=lambda item: item["haplotype"])
        )

    assembly_fields = list(assembly_rows[0])
    write_tsv(args.assembly_output, assembly_rows, assembly_fields)
    write_tsv(
        args.sample_output,
        sample_rows,
        ["sample", "sample_status", "hap1_status", "hap2_status", "reasons"],
    )
    write_tsv(
        args.excluded_output,
        excluded_rows,
        ["sample", "assembly_id", "path", "sample_status", "reasons"],
    )
    Path(args.included_output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.included_output, "w") as handle:
        for path in included_paths:
            handle.write(path + "\n")

    counts = defaultdict(int)
    for row in sample_rows:
        counts[row["sample_status"]] += 1
    print(
        f"Assemblies: {len(assembly_rows)}; samples: {len(sample_rows)}; "
        f"PASS={counts['PASS']}, WARN={counts['WARN']}, "
        f"PARTIAL={counts['PARTIAL']}, FAIL={counts['FAIL']}"
    )


if __name__ == "__main__":
    main()
