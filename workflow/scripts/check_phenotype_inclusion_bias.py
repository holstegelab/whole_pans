#!/usr/bin/env python3
"""Check whether graph inclusion differs between two phenotype groups."""

import argparse
import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path


ASSEMBLY_RE = re.compile(
    r"^(?P<sample>.+?)\.hifi\.hifiasm\.bp\."
    r"(?P<haplotype>hap[12])\.p_ctg(?:\..*)?$"
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compare graph inclusion rates between AD and centenarian samples. "
            "A sample is included when at least one haplotype occurs in the graph "
            "inclusion list."
        )
    )
    parser.add_argument(
        "--ad",
        required=True,
        type=Path,
        help="Text file containing one AD sample ID per line",
    )
    parser.add_argument(
        "--centenarians",
        required=True,
        type=Path,
        help="Text file containing one centenarian sample ID per line",
    )
    parser.add_argument(
        "--included",
        required=True,
        type=Path,
        help=(
            "Graph inclusion list containing one assembly FASTA path or assembly "
            "ID per line"
        ),
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        help="Optional TSV output for group-level inclusion statistics",
    )
    parser.add_argument(
        "--sample-output",
        type=Path,
        help="Optional TSV output with one row per phenotype sample",
    )
    return parser.parse_args()


def read_sample_ids(path):
    sample_ids = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            value = line.strip()
            if not value or value.startswith("#"):
                continue
            if any(character.isspace() for character in value):
                raise ValueError(
                    f"Expected one sample ID at {path}:{line_number}, found {value!r}"
                )
            sample_ids.append(value)

    if not sample_ids:
        raise ValueError(f"No sample IDs found in {path}")
    duplicates = sorted(
        sample for sample in set(sample_ids) if sample_ids.count(sample) > 1
    )
    if duplicates:
        raise ValueError(
            f"Duplicate sample IDs in {path}: {', '.join(duplicates[:10])}"
        )
    return set(sample_ids)


def read_included_haplotypes(path):
    included = defaultdict(set)
    seen_assemblies = set()
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            value = line.strip()
            if not value or value.startswith("#"):
                continue
            name = Path(value).name
            match = ASSEMBLY_RE.match(name)
            if not match:
                raise ValueError(
                    f"Cannot parse hifiasm assembly name at {path}:{line_number}: "
                    f"{value!r}"
                )
            assembly = (match.group("sample"), match.group("haplotype"))
            if assembly in seen_assemblies:
                raise ValueError(
                    f"Duplicate included assembly at {path}:{line_number}: "
                    f"{assembly[0]} {assembly[1]}"
                )
            seen_assemblies.add(assembly)
            included[assembly[0]].add(assembly[1])
    return included


def group_statistics(samples, included_haplotypes):
    counts = {0: 0, 1: 0, 2: 0}
    included_haplotype_count = 0
    for sample in samples:
        haplotypes = included_haplotypes.get(sample, set())
        count = len(haplotypes)
        if count > 2:
            raise ValueError(f"More than two haplotypes found for sample {sample}")
        counts[count] += 1
        included_haplotype_count += count

    total = len(samples)
    included_samples = counts[1] + counts[2]
    possible_haplotypes = 2 * total
    return {
        "total_samples": total,
        "included_samples": included_samples,
        "sample_inclusion_rate": included_samples / total,
        "fully_included_samples": counts[2],
        "partially_included_samples": counts[1],
        "excluded_samples": counts[0],
        "included_haplotypes": included_haplotype_count,
        "possible_haplotypes": possible_haplotypes,
        "haplotype_inclusion_rate": (
            included_haplotype_count / possible_haplotypes
        ),
    }


def hypergeometric_probability(x, row_ad, row_cent, included_total):
    return (
        math.comb(row_ad, x)
        * math.comb(row_cent, included_total - x)
        / math.comb(row_ad + row_cent, included_total)
    )


def fisher_exact_two_sided(ad_included, ad_excluded, cent_included, cent_excluded):
    row_ad = ad_included + ad_excluded
    row_cent = cent_included + cent_excluded
    included_total = ad_included + cent_included
    lower = max(0, included_total - row_cent)
    upper = min(row_ad, included_total)
    observed = hypergeometric_probability(
        ad_included, row_ad, row_cent, included_total
    )
    return min(
        1.0,
        sum(
            hypergeometric_probability(x, row_ad, row_cent, included_total)
            for x in range(lower, upper + 1)
            if hypergeometric_probability(x, row_ad, row_cent, included_total)
            <= observed + 1e-12
        ),
    )


def odds_ratio_and_ci(a, b, c, d):
    corrected = any(value == 0 for value in (a, b, c, d))
    if corrected:
        a, b, c, d = (value + 0.5 for value in (a, b, c, d))
    odds_ratio = a * d / (b * c)
    standard_error = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    lower = math.exp(math.log(odds_ratio) - 1.96 * standard_error)
    upper = math.exp(math.log(odds_ratio) + 1.96 * standard_error)
    return odds_ratio, lower, upper, corrected


def write_summary(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "phenotype",
        "total_samples",
        "included_samples",
        "sample_inclusion_rate_percent",
        "fully_included_samples",
        "partially_included_samples",
        "excluded_samples",
        "included_haplotypes",
        "possible_haplotypes",
        "haplotype_inclusion_rate_percent",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def write_samples(path, groups, included_haplotypes):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "sample",
        "phenotype",
        "included_haplotypes",
        "included_haplotype_count",
        "inclusion_status",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for phenotype, samples in groups.items():
            for sample in sorted(samples):
                haplotypes = sorted(included_haplotypes.get(sample, set()))
                count = len(haplotypes)
                status = {0: "excluded", 1: "partial", 2: "full"}[count]
                writer.writerow(
                    {
                        "sample": sample,
                        "phenotype": phenotype,
                        "included_haplotypes": ",".join(haplotypes),
                        "included_haplotype_count": count,
                        "inclusion_status": status,
                    }
                )


def percent(value):
    return f"{100 * value:.2f}%"


def main():
    args = parse_args()
    groups = {
        "AD": read_sample_ids(args.ad),
        "centenarian": read_sample_ids(args.centenarians),
    }
    overlap = sorted(groups["AD"] & groups["centenarian"])
    if overlap:
        raise ValueError(
            "Samples occur in both phenotype lists: " + ", ".join(overlap[:10])
        )

    included_haplotypes = read_included_haplotypes(args.included)
    phenotype_samples = groups["AD"] | groups["centenarian"]
    unclassified = sorted(set(included_haplotypes) - phenotype_samples)
    missing_from_included = phenotype_samples - set(included_haplotypes)

    statistics = {
        phenotype: group_statistics(samples, included_haplotypes)
        for phenotype, samples in groups.items()
    }
    ad = statistics["AD"]
    cent = statistics["centenarian"]

    a = ad["included_samples"]
    b = ad["excluded_samples"]
    c = cent["included_samples"]
    d = cent["excluded_samples"]
    risk_difference = ad["sample_inclusion_rate"] - cent["sample_inclusion_rate"]
    risk_ratio = (
        ad["sample_inclusion_rate"] / cent["sample_inclusion_rate"]
        if cent["sample_inclusion_rate"]
        else math.inf
    )
    odds_ratio, odds_lower, odds_upper, corrected = odds_ratio_and_ci(a, b, c, d)
    fisher_p = fisher_exact_two_sided(a, b, c, d)

    source_ad_share = len(groups["AD"]) / len(phenotype_samples)
    included_known = a + c
    included_ad_share = a / included_known if included_known else math.nan

    print("Phenotype representation among graph-included samples")
    print(f"Included list: {args.included}")
    print()
    print(
        "phenotype\ttotal\tincluded_any\trate\tfull_pair\tpartial_pair\t"
        "excluded\tincluded_haps\thaplotype_rate"
    )
    for phenotype in ("AD", "centenarian"):
        values = statistics[phenotype]
        print(
            f"{phenotype}\t{values['total_samples']}\t"
            f"{values['included_samples']}\t"
            f"{percent(values['sample_inclusion_rate'])}\t"
            f"{values['fully_included_samples']}\t"
            f"{values['partially_included_samples']}\t"
            f"{values['excluded_samples']}\t"
            f"{values['included_haplotypes']}/{values['possible_haplotypes']}\t"
            f"{percent(values['haplotype_inclusion_rate'])}"
        )

    print()
    print(f"AD minus centenarian inclusion rate: {100 * risk_difference:+.2f} percentage points")
    print(f"AD/centenarian risk ratio: {risk_ratio:.3f}")
    correction_note = " (0.5 correction for a zero cell)" if corrected else ""
    print(
        f"Odds ratio: {odds_ratio:.3f} "
        f"(95% CI {odds_lower:.3f}-{odds_upper:.3f}){correction_note}"
    )
    print(f"Fisher exact two-sided p-value: {fisher_p:.6g}")
    print(f"AD share of phenotype cohort: {percent(source_ad_share)}")
    print(f"AD share of included phenotype samples: {percent(included_ad_share)}")
    print(f"Included samples without AD/centenarian phenotype: {len(unclassified)}")
    print(f"Phenotype samples with no included haplotype: {len(missing_from_included)}")

    if unclassified:
        print(
            "Warning: included samples absent from both phenotype lists: "
            + ", ".join(unclassified[:10]),
            file=sys.stderr,
        )

    if args.summary_output:
        rows = []
        for phenotype in ("AD", "centenarian"):
            values = statistics[phenotype]
            rows.append(
                {
                    "phenotype": phenotype,
                    "total_samples": values["total_samples"],
                    "included_samples": values["included_samples"],
                    "sample_inclusion_rate_percent": (
                        f"{100 * values['sample_inclusion_rate']:.6f}"
                    ),
                    "fully_included_samples": values["fully_included_samples"],
                    "partially_included_samples": values[
                        "partially_included_samples"
                    ],
                    "excluded_samples": values["excluded_samples"],
                    "included_haplotypes": values["included_haplotypes"],
                    "possible_haplotypes": values["possible_haplotypes"],
                    "haplotype_inclusion_rate_percent": (
                        f"{100 * values['haplotype_inclusion_rate']:.6f}"
                    ),
                }
            )
        write_summary(args.summary_output, rows)

    if args.sample_output:
        write_samples(args.sample_output, groups, included_haplotypes)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        raise SystemExit(f"error: {error}") from error
