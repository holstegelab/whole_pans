#!/usr/bin/env python3
"""Render the preliminary SV-level rGFA short-read mapping QC report."""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter
from pathlib import Path


LEAD = (
    "This is a preliminary comparison of short-read mapping to SV-level rGFAs indexed "
    "with artificial path-cover haplotypes. It is not a base-resolution Minigraph-Cactus "
    "benchmark and does not estimate correct-locus mapping accuracy without independent truth."
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("smoke", "full"), required=True)
    parser.add_argument("--index-preflight", required=True, type=Path)
    parser.add_argument("--selected", required=True, type=Path)
    parser.add_argument("--membership", required=True, type=Path)
    parser.add_argument("--index-manifest", action="append", default=[])
    parser.add_argument(
        "--index-benchmark", action="append", default=[], metavar="GRAPH=PATH"
    )
    parser.add_argument("--metrics", required=True, type=Path)
    parser.add_argument("--contrasts", required=True, type=Path)
    parser.add_argument("--target-presence", required=True, type=Path)
    parser.add_argument("--evidence", required=True, type=Path)
    parser.add_argument("--tools", required=True, type=Path)
    parser.add_argument("--provenance", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def cell(value: object) -> str:
    return str(value if value not in (None, "") else "NA").replace("|", "\\|").replace("\n", " ")


def number(value: object, digits: int = 4) -> str:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return cell(value)
    if not math.isfinite(parsed):
        return "NA"
    return f"{parsed:.{digits}f}"


def scaled(value: object, divisor: float, digits: int = 2) -> str:
    try:
        return number(float(value) / divisor, digits)
    except (TypeError, ValueError):
        return "NA"


def markdown_table(headers: list[str], table_rows: list[list[object]]) -> list[str]:
    output = [
        "| " + " | ".join(headers) + " |",
        "|" + "|".join("---" for _ in headers) + "|",
    ]
    output.extend("| " + " | ".join(cell(value) for value in row) + " |" for row in table_rows)
    return output


def main() -> None:
    args = parse_args()
    preflight = rows(args.index_preflight)
    selected = rows(args.selected)
    membership = rows(args.membership)
    manifests = [row for value in args.index_manifest for row in rows(Path(value))]
    index_benchmarks = {}
    for value in args.index_benchmark:
        graph, separator, path = value.partition("=")
        if not separator:
            raise ValueError(f"--index-benchmark must be GRAPH=PATH: {value}")
        benchmark_rows = rows(Path(path))
        index_benchmarks[graph] = benchmark_rows[-1] if benchmark_rows else {}
    metrics = rows(args.metrics)
    contrasts = rows(args.contrasts)
    targets = rows(args.target_presence)
    evidence = rows(args.evidence)
    tools = rows(args.tools)
    provenance = rows(args.provenance)

    lines = [LEAD, "", "# SV-graph short-read mapping QC", ""]
    if args.mode == "smoke":
        lines.extend([
            "**Smoke-mode result:** this deterministic subset is for index and genome-wide mapping QC. "
            "Absence of TMEM junction reads at this depth is not a genotype result.", "",
        ])

    lines.extend(["## Graph and index preflight", ""])
    lines.extend(markdown_table(
        ["Graph", "Source GiB", "P lines", "W lines", "Indexing mode", "Graph QC"],
        [[
            row.get("graph_id"),
            scaled(row.get("source_size_bytes"), 1024**3, 2),
            row.get("p_path_lines"), row.get("w_walk_lines"),
            row.get("indexing_mode"), row.get("graph_qc_status"),
        ] for row in preflight],
    ))
    lines.extend(["", (
        "Every derived index uses `rgfa_artificial_path_cover`. P/W counts describe source GFA "
        "lines but do not turn these SV-level rGFAs into real haplotype-embedded base-resolution graphs."
    ), ""])

    if manifests:
        lines.extend(markdown_table(
            [
                "Graph", "VG", "GBZ GiB", "All indexes GiB", "Index wall h",
                "Peak RSS MB", "Artificial cover", "Cover log seen", "Validation",
            ],
            [[
                row.get("graph_id"), row.get("vg_version"),
                scaled(row.get("gbz_size_bytes"), 1024**3, 2),
                scaled(row.get("total_index_size_bytes"), 1024**3, 2),
                scaled(index_benchmarks.get(row.get("graph_id", ""), {}).get("s"), 3600, 2),
                number(index_benchmarks.get(row.get("graph_id", ""), {}).get("max_rss"), 1),
                row.get("artificial_path_cover_created"),
                row.get("artificial_path_cover_log_seen"), row.get("validation"),
            ] for row in sorted(manifests, key=lambda value: value.get("graph_id", ""))],
        ))
        lines.append("")
        unconfirmed_cover_logs = [
            row.get("graph_id", "") for row in manifests
            if row.get("artificial_path_cover_log_seen") != "true"
        ]
        if unconfirmed_cover_logs:
            lines.extend([
                "**Index warning:** the completion manifest records artificial path-cover indexing, "
                "but no explicit path-cover message was detected in the VG log for: "
                + ", ".join(f"`{graph}`" for graph in unconfirmed_cover_logs)
                + ". Inspect those autoindex logs before interpreting mappings.",
                "",
            ])

    lines.extend(["## Selected samples", ""])
    role_counts = Counter(row.get("role", "") for row in selected)
    lines.append(
        f"Frozen manifest selection: {len(selected)} samples; "
        + ", ".join(f"`{role}`={count}" for role, count in sorted(role_counts.items()))
        + f". Expected mapping rows: {len(selected) * len(preflight)}."
    )
    lines.extend(["", *markdown_table(
        ["Sample", "Role", "Batch", "R1 lanes", "R2 lanes"],
        [[
            row.get("sample_id"), row.get("role"), row.get("batch"),
            len(row.get("r1", "").split(",")), len(row.get("r2", "").split(",")),
        ] for row in selected],
    ), ""])

    lines.extend(["## Graph membership/leakage audit", ""])
    audit_counts = Counter(row.get("audit_status", "") for row in membership)
    lines.append(
        "Membership rows: " + ", ".join(f"{status}={count}" for status, count in sorted(audit_counts.items())) + "."
    )
    uncertain = [row for row in membership if row.get("audit_status") == "UNCERTAIN"]
    failed = [row for row in membership if row.get("audit_status") == "FAIL"]
    if failed:
        lines.append("**FAIL: at least one selected sample is recorded in a target graph membership manifest.**")
    if uncertain:
        lines.append(
            "Membership uncertainty remains and is not converted into an absence claim: "
            + "; ".join(
                f"`{row['sample_id']}`/`{row['graph_id']}` ({row.get('reason')})" for row in uncertain[:12]
            )
            + ("; ..." if len(uncertain) > 12 else "") + "."
        )
    lines.append("")

    lines.extend(["## Per-sample genome-wide mapping QC", ""])
    if not metrics:
        lines.extend(["No mapping metric rows were produced.", ""])
    else:
        lines.extend(markdown_table(
            ["Sample", "Role", "Graph", "Mapped", "MAPQ≥20", "Median MAPQ", "Identity", "Aligned", "Clip bp/read", "Fallback"],
            [[
                row.get("sample_id"), row.get("sample_role"), row.get("graph_id"),
                number(row.get("mapped_read_fraction")), number(row.get("mapq_ge_20_fraction")),
                number(row.get("median_mapq"), 2), number(row.get("mean_alignment_identity")),
                number(row.get("mean_query_aligned_fraction")), number(row.get("mean_clipped_bases_per_read"), 2),
                row.get("single_end_fallback"),
            ] for row in metrics],
        ))
        lines.extend(["", (
            "Mapped-read fractions use independently counted FASTQ reads as the denominator. MAPQ "
            "fractions use primary mapped reads. Proper-pair percentage is unavailable from GAF and is not inferred."
        ), ""])

    lines.extend(["## Paired graph differences from HPRC", ""])
    summaries = [
        row for row in contrasts
        if row.get("row_type") in {"summary_overall", "summary_role"}
        and row.get("metric") in {
            "mapped_read_fraction", "mapq_ge_20_fraction", "mean_alignment_identity",
            "mean_query_aligned_fraction", "mean_clipped_bases_per_read",
        }
    ]
    if summaries:
        lines.extend(markdown_table(
            ["Contrast", "Group", "Metric", "N", "Median Δ", "Range"],
            [[
                row.get("contrast_id"), row.get("summary_group"), row.get("metric"),
                row.get("n_samples"), number(row.get("median_delta")),
                f"{number(row.get('min_delta'))} to {number(row.get('max_delta'))}",
            ] for row in summaries],
        ))
        lines.extend(["", (
            "The sample is the comparison unit. Deltas are right graph minus left graph; no read-level "
            "inferential p-values are reported for this small selected cohort."
        ), ""])
    else:
        lines.extend(["No complete paired contrast summaries were available.", ""])

    lines.extend(["## Pairing behavior and resource use", ""])
    fallback = [row for row in metrics if row.get("single_end_fallback") == "true"]
    learning_fail = [row for row in metrics if row.get("paired_fragment_learning_succeeded") != "true"]
    if fallback:
        lines.append("**Single-end fallback warning:** " + "; ".join(
            f"`{row['sample_id']}`/`{row['graph_id']}` ({row.get('single_end_fallback_read_count')} reads)"
            for row in fallback
        ) + ".")
    else:
        lines.append("No explicit Giraffe single-end fallback was detected in completed mapping logs.")
    if learning_fail:
        lines.append("Fragment learning was not confirmed for: " + "; ".join(
            f"`{row['sample_id']}`/`{row['graph_id']}`" for row in learning_fail
        ) + ".")
    if metrics:
        lines.extend(["", *markdown_table(
            ["Sample", "Graph", "Wall min", "Peak RSS MB", "Reads/s", "Fragment mean±SD"],
            [[
                row.get("sample_id"), row.get("graph_id"),
                scaled(row.get("mapping_wall_time_seconds"), 60, 2),
                number(row.get("peak_rss_mb"), 1), number(row.get("reads_per_second"), 1),
                f"{cell(row.get('learned_fragment_mean'))}±{cell(row.get('learned_fragment_sd'))}",
            ] for row in metrics],
        )])
    lines.append("")

    lines.extend(["## TMEM graph-content presence and callability", ""])
    if targets:
        lines.extend(markdown_table(
            ["Graph", "Present", "Callable", "ALT query", "Coverage", "Identity", "Failure/ambiguity"],
            [[
                row.get("graph_id"), row.get("target_present"), row.get("callable_status"),
                row.get("alternate_query_name"), number(row.get("alternate_query_coverage")),
                number(row.get("alternate_query_identity")),
                row.get("failure_reason") or row.get("ambiguity"),
            ] for row in targets],
        ))
        lines.extend(["", (
            "Targets were derived independently for each graph from reference and alternate queries with "
            "at least 1 kb flanks, assembly-to-graph GAF mappings, base-level CIGAR, and orientation-aware traversals."
        ), ""])
    else:
        lines.extend(["TMEM graph-content audit was not produced; junction analysis is blocked.", ""])

    lines.extend(["## TMEM mapped-read evidence", ""])
    if evidence:
        show = [row for row in evidence if row.get("sample_role") in {"tmem_candidate", "tmem_random_blinded", "tmem_matched_control"}]
        lines.extend(markdown_table(
            ["Sample", "Role", "Graph", "Callability", "Left", "Right", "Min", "Reference", "Both pass"],
            [[
                row.get("sample_id"), row.get("sample_role"), row.get("graph_id"),
                row.get("evidence_status"), row.get("left_alt_junction_distinct_primary_reads"),
                row.get("right_alt_junction_distinct_primary_reads"), row.get("minimum_junction_support"),
                row.get("reference_traversal_distinct_primary_reads"), row.get("both_junctions_meet_threshold"),
            ] for row in show],
        ))
        lines.extend(["", (
            "A candidate is supported only by distinct primary read names crossing both graph-specific "
            "alternate junctions at the configured MAPQ/coverage threshold. Random samples remain blinded/unknown "
            "and are not labelled concordant or discordant."
        ), ""])
    else:
        lines.extend(["No TMEM read-evidence rows were produced.", ""])

    lines.extend(["## Cautious recommendation", ""])
    if not metrics:
        lines.append("No graph can yet be recommended from this design because mapping results are unavailable.")
    elif args.mode == "smoke":
        lines.append(
            "Use this run only to decide whether index construction, paired mapping, denominators, and warnings "
            "are technically acceptable. Inspect every fallback and clipping/MAPQ contrast before enabling full-depth mapping."
        )
    else:
        lines.append(
            "Advance a local or HPRC+local graph only if its within-sample mapped fraction is accompanied by "
            "stable or improved MAPQ, clipping, identity, and paired behavior. Treat improved TMEM junction support "
            "as a graph-content result only where the graph-specific target audit is callable. Independent truth-labelled "
            "reads or orthogonal genotypes are required before making accuracy, sensitivity, or specificity claims."
        )
    lines.append("")

    lines.extend([
        "## Explicit limitations", "",
        "- These are SV-level rGFAs indexed with artificial path-cover haplotypes, not haplotype-embedded base-resolution GBZ graphs.",
        "- Results are “SV-graph short-read mapping QC”; real-read mapping rates are not mapping accuracy.",
        "- This benchmark does not predict performance of a future Minigraph-Cactus graph.",
        "- Sensitivity and specificity cannot be estimated without independent truth.",
        "- TMEM-positive candidates were selected using previous local-graph short-read evidence, favoring the local representation; this is not independent genotype truth.",
        "- Random samples are blinded/unknown, not confirmed noncarriers.",
        "- If TMEM is absent from HPRC, missing junction support is a graph-content result, not mapper failure; absent/ambiguous targets are not reported as zero-support negative genotypes.",
        "- Local V2 remains experimental because its cohort-accounting QC is currently FAIL.",
        "- Better mapping rate alone is insufficient if MAPQ, clipping, identity, or paired behavior worsens.",
        "- Internal repetitive insertion sequence alone does not count as carrier evidence.",
        "- Smoke-mode TMEM absence is never interpreted as a genotype result.",
        "- The selected cohort is small and enriched by prior evidence; reads are not independent biological replicates.",
        "",
        "## Tool/provenance status", "",
    ])
    lines.extend(markdown_table(
        ["Tool", "Version", "Status", "Required check"],
        [[row.get("tool"), row.get("version"), row.get("status"), row.get("required_check")] for row in tools],
    ))
    lines.extend(["", f"Provenance keys recorded: {len(provenance)}.", ""])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(args.output) + ".tmp")
    temporary.write_text("\n".join(lines) + "\n")
    temporary.replace(args.output)


if __name__ == "__main__":
    main()
