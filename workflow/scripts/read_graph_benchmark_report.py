#!/usr/bin/env python3
"""Render the compact benchmark report without making truth-free accuracy claims."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", required=True, type=Path)
    parser.add_argument("--leakage", required=True, type=Path)
    parser.add_argument("--metrics", required=True, type=Path)
    parser.add_argument("--contrasts", required=True, type=Path)
    parser.add_argument("--simulation", required=True, type=Path)
    parser.add_argument("--rescue", required=True, type=Path)
    parser.add_argument("--tmem", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    args = parse_args()
    preflight = rows(args.preflight)
    leakage = rows(args.leakage)
    metrics = rows(args.metrics)
    contrasts = rows(args.contrasts)
    simulation = rows(args.simulation)
    rescue = rows(args.rescue)
    tmem = rows(args.tmem)
    lines = [
        "# Read-to-graph benchmark report", "",
        "Real-read results below are mapping quality/control metrics, not mapping accuracy. Accuracy is reported only for truth-labelled simulated reads.", "",
        "## Graph availability and compatibility", "",
        "| Graph | Available | Type | Backbone compatible | Long mapper | Short mapper | QC status |",
        "|---|---:|---|---:|---|---|---|",
    ]
    for row in preflight:
        lines.append(
            f"| `{row['graph_id']}` | {row['available']} | {row['graph_type']} | "
            f"{row['backbone_compatible']} | {row['long_read_mapper']} | "
            f"{row['short_read_mapper']} | {row['qc_status']} |"
        )
    unavailable = [row for row in preflight if row.get("available") != "true"]
    if unavailable:
        lines.extend(["", "Unavailable/blocked artifacts:", ""])
        lines.extend(f"- `{row['graph_id']}`: {row.get('problems') or 'unavailable'}" for row in unavailable)

    lines.extend(["", "## Membership leakage", ""])
    failed = [row for row in leakage if row.get("audit_status") == "FAIL"]
    if failed:
        lines.append("**FAIL: held-out leakage or unauditable membership was detected. Headline comparison is blocked.**")
        lines.extend(f"- {row['sample_id']} / {row['graph_id']}: {row['reason']}" for row in failed)
    elif any(row.get("audit_status") == "PASS" for row in leakage):
        heldout = {row["sample_id"] for row in leakage if row.get("audit_status") == "PASS"}
        lines.append(f"PASS for {len(heldout)} held-out sample(s) across all declared headline graphs.")
    else:
        lines.append("No `held_out_primary` samples are configured; headline generalization results are therefore unavailable.")

    lines.extend(["", "## Real-read mapping QC by sample role", ""])
    role_counts = Counter(row.get("sample_role", "") for row in metrics)
    if role_counts:
        lines.extend(f"- `{role}`: {count} sample/graph/mapper result rows" for role, count in sorted(role_counts.items()))
        lines.append("MAPQ is interpreted only within a mapper; Minigraph and minimap2 MAPQ values are not directly compared.")
    else:
        lines.append("No real-read mapping metrics were produced.")

    lines.extend(["", "## Held-out paired contrasts", ""])
    heldout_contrasts = [row for row in contrasts if row.get("headline_eligible") == "true" and row.get("status") == "PASS"]
    if heldout_contrasts:
        grouped = defaultdict(int)
        for row in heldout_contrasts:
            grouped[row["contrast_id"]] += 1
        lines.extend(f"- `{contrast}`: {count} per-sample metric deltas" for contrast, count in sorted(grouped.items()))
        lines.append("Deltas are right minus left and remain per-sample; no read-level p-values are reported.")
    else:
        lines.append("No eligible held-out paired contrasts are available. Control groups are not substituted into the headline estimate.")

    lines.extend(["", "## Graph-membership control contrasts", ""])
    control_contrasts = [
        row for row in contrasts
        if row.get("headline_eligible") != "true" and row.get("status") == "PASS"
    ]
    if control_contrasts:
        transitions = Counter(
            (row.get("sample_role", ""), row.get("membership_transition", ""))
            for row in control_contrasts
        )
        lines.extend(
            f"- `{role}` / `{transition}`: {count} metric deltas"
            for (role, transition), count in sorted(transitions.items())
        )
        lines.append("These are graph-membership controls, not held-out generalization estimates.")
    else:
        lines.append("No graph-membership control contrasts were produced.")

    lines.extend(["", "## Simulated mapping accuracy", ""])
    if simulation:
        lines.append(f"Truth-labelled simulation rows: {len(simulation)}. Simulation results remain separate from real-read QC.")
    else:
        lines.append("Simulation mode was disabled or no phased truth assemblies were selected for simulation.")

    lines.extend(["", "## V2 rescue capture", ""])
    if rescue:
        counts = Counter(row.get("capture_status", "") for row in rescue)
        lines.extend(f"- `{status}`: {count}" for status, count in sorted(counts.items()))
        independently_tested = sum(row.get("independent_carrier_tested") == "true" for row in rescue)
        lines.append(f"Rows with an independent carrier configured: {independently_tested}.")
    else:
        lines.append("No rescue-capture results were produced.")
    lines.append("V2 is experimental: structural integrity passed, but overall cohort-accounting QC is FAIL. That pre-existing QC failure is provenance, not a mapping outcome.")

    lines.extend(["", "## TMEM candidate/control pilot", ""])
    if tmem:
        by_role = Counter(row.get("sample_role", "") for row in tmem)
        by_call = Counter(row.get("evidence_call", "") for row in tmem)
        lines.append("Roles: " + ", ".join(f"{key}={value}" for key, value in sorted(by_role.items())) + ".")
        lines.append("Calls: " + ", ".join(f"{key}={value}" for key, value in sorted(by_call.items())) + ".")
        candidate_calls = Counter(row.get("evidence_call", "") for row in tmem if row.get("sample_role") == "tmem_candidate")
        control_calls = Counter(row.get("evidence_call", "") for row in tmem if row.get("sample_role") == "tmem_matched_control")
        candidate_concordance = Counter(row.get("evidence_concordance", "") for row in tmem if row.get("sample_role") == "tmem_candidate")
        control_concordance = Counter(row.get("evidence_concordance", "") for row in tmem if row.get("sample_role") == "tmem_matched_control")
        lines.append("Candidate calls: " + ", ".join(f"{key}={value}" for key, value in sorted(candidate_calls.items())) + ".")
        lines.append("Matched-control calls: " + ", ".join(f"{key}={value}" for key, value in sorted(control_calls.items())) + ".")
        lines.append("Candidate expectation concordance: " + ", ".join(f"{key}={value}" for key, value in sorted(candidate_concordance.items())) + ".")
        lines.append("Matched-control expectation concordance: " + ", ".join(f"{key}={value}" for key, value in sorted(control_concordance.items())) + ".")
    else:
        lines.append("No TMEM pilot evidence was produced.")
    lines.append("The ~18.7 kb insertion call requires mapped support at both junctions; allele-edge k-mers and the reference path are reported as corroborating evidence. Fixed expectations come from the prior strict screen and are not independent truth. Internal repetitive sequence and the non-discriminating AluYb8 result are excluded as carrier evidence.")

    lines.extend([
        "", "## Limitations and unsupported claims", "",
        "- Graph-member, rescue-source and held-out samples are never pooled into one accuracy estimate.",
        "- Real reads without truth support mapping QC comparisons only, not correct-locus accuracy claims.",
        "- Minigraph SV-level rGFA and Minigraph-Cactus GBZ are different representations; unsupported conversions are not performed.",
        "- Short reads cannot span the full TMEM insertion; evidence is junction/allele based.",
        "- No inferential p-values treat reads as biological replicates.",
        "- Unified SV counts refer to normalized events/clusters, never raw VCF record counts.",
    ])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(args.output) + ".tmp")
    temporary.write_text("\n".join(lines) + "\n")
    temporary.replace(args.output)


if __name__ == "__main__":
    main()
