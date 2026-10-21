#!/usr/bin/env python3
"""Summarize DRAGMAP flat-reference alignments and compare them with graph mapping."""

from __future__ import annotations

import argparse
import math
import statistics
from collections import Counter
from pathlib import Path

from personalized_sv_statistics import (
    benchmark_values,
    histogram_median,
    parse_coverage,
    parse_fastq_stats,
    parse_flagstat,
    ratio,
    read_tsv,
    write_tsv,
)


FLAT_FIELDS = [
    "sample_id", "sample_group", "batch", "target", "reference_fasta",
    "reference_sha256", "reference_size_bytes", "hash_size_bytes",
    "mapper_version", "input_reads", "input_read_pairs", "input_bases",
    "primary_mapped_reads", "mapped_read_fraction", "mapq_0_fraction",
    "mapq_ge_5_fraction", "mapq_ge_20_fraction", "mapq_ge_40_fraction",
    "mean_mapq", "median_mapq", "mean_alignment_identity",
    "fraction_read_bases_aligned", "clipped_bases_per_read",
    "secondary_multimapping_fraction", "proper_pair_fraction",
    "reference_mapped_fraction", "reference_breadth_fraction",
    "reference_mean_depth", "mapping_mode", "index_wall_time_seconds",
    "index_peak_rss_mb", "mapping_wall_time_seconds",
    "mapping_cpu_time_seconds", "mapping_peak_rss_mb", "bam_size_bytes",
]

PAIRED_METRICS = [
    "mapped_read_fraction", "mapq_0_fraction", "mapq_ge_5_fraction",
    "mapq_ge_20_fraction", "mapq_ge_40_fraction", "mean_mapq",
    "median_mapq", "mean_alignment_identity", "fraction_read_bases_aligned",
    "clipped_bases_per_read", "secondary_multimapping_fraction",
    "proper_pair_fraction", "mapping_wall_time_seconds",
    "mapping_cpu_time_seconds", "mapping_peak_rss_mb",
]

PAIR_FIELDS = [
    "sample_id", "sample_group", "batch", "graph_target", "flat_target",
    "metric", "graph_value", "flat_value", "graph_minus_flat",
    "input_reads", "input_bases",
]
SUMMARY_FIELDS = [
    "graph_target", "flat_target", "metric", "n_samples",
    "mean_graph_value", "mean_flat_value", "mean_difference",
    "median_difference",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    metrics = commands.add_parser("metrics", help="summarize one flat-reference BAM")
    for option in (
        "bam", "fastq-stats", "coverage", "flagstat", "reference-info",
        "mapping-benchmark", "index-benchmark", "output",
    ):
        metrics.add_argument(f"--{option}", required=True, type=Path)
    for option in ("sample", "sample-group", "batch", "target", "mapper-version"):
        metrics.add_argument(f"--{option}", required=True)

    compare = commands.add_parser("compare", help="pair four graph maps with two flat references")
    compare.add_argument("--samples", required=True, type=Path)
    compare.add_argument("--flat", action="append", required=True, type=Path)
    compare.add_argument("--graph", action="append", required=True, type=Path)
    for option in ("flat-output", "pairs-output", "summary-output", "report-output"):
        compare.add_argument(f"--{option}", required=True, type=Path)
    return parser.parse_args()


def bam_mapping_metrics(path: Path, input_reads: int, input_bases: int) -> dict[str, object]:
    """Count primary read alignments while retaining secondary record frequency."""
    import pysam

    primary = secondary = aligned_bases = clipped_bases = 0
    identity_sum = mapq_sum = 0.0
    mapq_histogram: Counter[int] = Counter()
    with pysam.AlignmentFile(str(path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary:
                secondary += 1
                continue
            if read.is_supplementary or read.is_unmapped:
                continue
            primary += 1
            mapq = read.mapping_quality
            if 0 <= mapq < 255:
                mapq_histogram[mapq] += 1
                mapq_sum += mapq

            cigar = read.cigartuples
            if not cigar:
                raise ValueError(f"mapped BAM record has no CIGAR: {read.query_name}")
            query_aligned = sum(length for op, length in cigar if op in (0, 1, 7, 8))
            block_length = sum(length for op, length in cigar if op in (0, 1, 2, 7, 8))
            clipped = sum(length for op, length in cigar if op in (4, 5))
            if block_length == 0:
                raise ValueError(f"mapped BAM record has no aligned bases: {read.query_name}")
            if not read.has_tag("NM"):
                raise ValueError(f"mapped BAM record has no NM tag: {read.query_name}")
            nm = read.get_tag("NM")
            if not 0 <= nm <= block_length:
                raise ValueError(f"invalid NM tag for BAM record: {read.query_name}")
            identity_sum += (block_length - nm) / block_length
            aligned_bases += query_aligned
            clipped_bases += clipped

    if primary > input_reads:
        raise ValueError(f"{path}: {primary} primary mapped records exceed {input_reads} FASTQ reads")
    return {
        "primary_mapped_reads": primary,
        "mapped_read_fraction": ratio(primary, input_reads),
        "mapq_0_fraction": ratio(mapq_histogram[0], input_reads),
        "mapq_ge_5_fraction": ratio(sum(n for q, n in mapq_histogram.items() if q >= 5), input_reads),
        "mapq_ge_20_fraction": ratio(sum(n for q, n in mapq_histogram.items() if q >= 20), input_reads),
        "mapq_ge_40_fraction": ratio(sum(n for q, n in mapq_histogram.items() if q >= 40), input_reads),
        "mean_mapq": ratio(mapq_sum, sum(mapq_histogram.values())),
        "median_mapq": histogram_median(mapq_histogram),
        "mean_alignment_identity": ratio(identity_sum, primary),
        "fraction_read_bases_aligned": ratio(aligned_bases, input_bases),
        "clipped_bases_per_read": ratio(clipped_bases, primary),
        "secondary_multimapping_fraction": ratio(secondary, input_reads),
    }


def reference_info(path: Path, target: str) -> dict[str, str]:
    required = (
        "target", "reference_fasta", "reference_sha256",
        "reference_size_bytes", "hash_size_bytes",
    )
    rows = read_tsv(path)
    matches = [row for row in rows if row.get("target") == target]
    if len(matches) != 1 or any(not matches[0].get(field) for field in required):
        raise ValueError(f"{path}: expected one complete reference-info row for {target}")
    return matches[0]


def write_metrics(args: argparse.Namespace) -> None:
    reads, bases = parse_fastq_stats(args.fastq_stats)
    if not reads or not bases:
        raise ValueError(f"{args.fastq_stats}: input FASTQ read/base counts must be positive")
    info = reference_info(args.reference_info, args.target)
    breadth, depth = parse_coverage(args.coverage)
    reference_mapped, proper_pair = parse_flagstat(args.flagstat)
    wall, rss, cpu = benchmark_values(args.mapping_benchmark)
    index_wall, index_rss, _ = benchmark_values(args.index_benchmark)
    row = {
        "sample_id": args.sample,
        "sample_group": args.sample_group,
        "batch": args.batch,
        "target": args.target,
        "reference_fasta": info["reference_fasta"],
        "reference_sha256": info["reference_sha256"],
        "reference_size_bytes": info["reference_size_bytes"],
        "hash_size_bytes": info["hash_size_bytes"],
        "mapper_version": args.mapper_version,
        "input_reads": reads,
        "input_read_pairs": reads // 2,
        "input_bases": bases,
        **bam_mapping_metrics(args.bam, reads, bases),
        "proper_pair_fraction": proper_pair,
        "reference_mapped_fraction": reference_mapped,
        "reference_breadth_fraction": breadth,
        "reference_mean_depth": depth,
        "mapping_mode": "flat_reference",
        "index_wall_time_seconds": index_wall,
        "index_peak_rss_mb": index_rss,
        "mapping_wall_time_seconds": wall,
        "mapping_cpu_time_seconds": cpu,
        "mapping_peak_rss_mb": rss,
        "bam_size_bytes": args.bam.stat().st_size,
    }
    write_tsv(args.output, [row], FLAT_FIELDS)


def numeric(value: str) -> float | None:
    if value == "" or value is None:
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"non-finite metric value: {value}")
    return result


def formatted(value: float | None) -> str:
    return "" if value is None else f"{value:.10g}"


def read_rows(paths: list[Path]) -> list[dict[str, str]]:
    rows = []
    for path in paths:
        new_rows = read_tsv(path)
        if not new_rows:
            raise ValueError(f"{path}: no data rows")
        rows.extend(new_rows)
    return rows


def indexed_rows(
    rows: list[dict[str, str]], samples: dict[str, dict[str, str]],
    expected_targets: int, label: str,
) -> tuple[dict[tuple[str, str], dict[str, str]], list[str]]:
    targets = sorted({row.get("target", "") for row in rows})
    if "" in targets or len(targets) != expected_targets:
        raise ValueError(f"expected {expected_targets} distinct {label} targets, found {targets}")
    indexed = {}
    for row in rows:
        sample = row.get("sample_id", "")
        target = row["target"]
        key = sample, target
        if sample not in samples:
            raise ValueError(f"unexpected {label} sample: {sample}")
        if key in indexed:
            raise ValueError(f"duplicate {label} sample/target: {sample}/{target}")
        for field in ("sample_group", "batch"):
            if row.get(field) != samples[sample][field]:
                raise ValueError(f"{label} {sample}/{target}: {field} disagrees with manifest")
        for field in ("input_reads", "input_bases"):
            if not row.get(field) or int(row[field]) <= 0:
                raise ValueError(f"{label} {sample}/{target}: missing or invalid {field}")
        indexed[key] = row
    expected = {(sample, target) for sample in samples for target in targets}
    if set(indexed) != expected:
        missing = sorted(expected - set(indexed))
        raise ValueError(f"incomplete {label} sample/target rows; first missing: {missing[:5]}")
    return indexed, targets


def mean_column(rows: list[dict[str, str]], field: str) -> str:
    values = [numeric(row.get(field, "")) for row in rows]
    present = [value for value in values if value is not None]
    return formatted(statistics.mean(present)) if present else ""


def markdown_number(value: str, scale: float = 1) -> str:
    return "N/A" if value == "" else f"{float(value) * scale:+.2f}"


def render_report(
    path: Path, samples_path: Path, flat_targets: list[str], graph_targets: list[str],
    flat: dict[tuple[str, str], dict[str, str]],
    graph: dict[tuple[str, str], dict[str, str]],
    summary: list[dict[str, object]], sample_count: int,
) -> None:
    summary_by_key = {(row["graph_target"], row["flat_target"], row["metric"]): row for row in summary}
    lines = [
        "# Flat reference versus pangenome graph read mapping", "",
        f"The same {sample_count} pilot samples were aligned to {len(flat_targets)} flat references "
        f"and {len(graph_targets)} graphs. All contrasts pair the same sample and FASTQ read/base counts.",
        f"Sample manifest: `{samples_path}`.", "",
        "## Paired mapping differences", "",
        "Values are mean graph minus flat differences across samples. Percentage-point columns "
        "are differences in proportions multiplied by 100. Positive values favor the graph only "
        "for metrics where larger values are desirable.", "",
        "| Graph | Flat reference | n | Mapped (pp) | MAPQ ≥20 (pp) | Identity (pp) | "
        "Aligned bases (pp) | Proper pairs (pp) | Wall time (s) | Peak RSS (MB) |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    display = (
        ("mapped_read_fraction", 100), ("mapq_ge_20_fraction", 100),
        ("mean_alignment_identity", 100), ("fraction_read_bases_aligned", 100),
        ("proper_pair_fraction", 100), ("mapping_wall_time_seconds", 1),
        ("mapping_peak_rss_mb", 1),
    )
    for graph_target in graph_targets:
        for flat_target in flat_targets:
            entries = [summary_by_key[graph_target, flat_target, metric] for metric, _ in display]
            n = min(int(entry["n_samples"]) for entry in entries)
            cells = [markdown_number(str(entry["mean_difference"]), scale) for entry, (_, scale) in zip(entries, display)]
            lines.append(f"| {graph_target} | {flat_target} | {n} | " + " | ".join(cells) + " |")

    lines.extend([
        "", "The complete per-sample differences and mean/median summaries are in the companion TSV files.",
        "", "## Reference coverage (descriptive)", "",
        "Coverage is calculated from each alignment's BAM. Graph BAMs were surjected to their "
        "CHM13 backbone; flat BAMs use the indicated linear reference.", "",
        "| Alignment | Target | Samples | Mean breadth | Mean depth |",
        "|---|---|---:|---:|---:|",
    ])
    for label, targets, indexed, breadth_field, depth_field in (
        ("Graph backbone", graph_targets, graph, "common_backbone_breadth_fraction", "common_backbone_mean_depth"),
        ("Flat reference", flat_targets, flat, "reference_breadth_fraction", "reference_mean_depth"),
    ):
        for target in targets:
            rows = [row for (sample, row_target), row in indexed.items() if row_target == target]
            breadth = mean_column(rows, breadth_field)
            depth = mean_column(rows, depth_field)
            lines.append(f"| {label} | {target} | {len(rows)} | {breadth or 'N/A'} | {depth or 'N/A'} |")

    lines.extend([
        "", "## Interpretation", "",
        "- CHM13 flat coverage and graph CHM13 backbone coverage may share a reference sequence, "
        "but graph surjection and alignment placement can differ. hg38 has a different sequence "
        "and denominator. Coverage is therefore shown descriptively, without paired differences.",
        "- MAPQ scales differ between DRAGMAP and Giraffe. Threshold contrasts describe their "
        "reported confidence, not independently measured mapping accuracy.",
        "- Flat identity uses BAM NM and CIGAR; graph identity uses GAF matches and block length. "
        "The estimates are close analogues, not identical definitions.",
        "- Mapping benchmarks include the steps in each Snakemake mapping rule, including output "
        "processing; index construction is recorded separately.",
        "- When a sample has multiple FASTQ pairs, DRAGMAP maps each pair separately and merges "
        "the sorted BAMs; the graph workflow streams its mate files together. Runtime includes "
        "those respective steps.",
    ])
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(path) + ".tmp")
    temporary.write_text("\n".join(lines) + "\n")
    temporary.replace(path)


def compare(args: argparse.Namespace) -> None:
    sample_rows = read_tsv(args.samples)
    if not sample_rows:
        raise ValueError(f"{args.samples}: sample manifest has no data rows")
    samples = {}
    for row in sample_rows:
        sample = row.get("sample_id", "")
        if not sample or sample in samples or not row.get("sample_group") or not row.get("batch"):
            raise ValueError(f"invalid or duplicate sample manifest row: {sample}")
        samples[sample] = row

    flat_rows = read_rows(args.flat)
    graph_rows = read_rows(args.graph)
    flat, flat_targets = indexed_rows(flat_rows, samples, 2, "flat")
    graph, graph_targets = indexed_rows(graph_rows, samples, 4, "graph")
    for sample in samples:
        rows = [flat[sample, target] for target in flat_targets]
        rows += [graph[sample, target] for target in graph_targets]
        for field in ("input_reads", "input_bases"):
            counts = {int(row[field]) for row in rows}
            if len(counts) != 1:
                raise ValueError(f"{sample}: {field} differs across graph and flat alignments")

    ordered_samples = list(samples)
    write_tsv(
        args.flat_output,
        [flat[sample, target] for sample in ordered_samples for target in flat_targets],
        FLAT_FIELDS,
    )

    pairs = []
    summary = []
    for graph_target in graph_targets:
        for flat_target in flat_targets:
            for metric in PAIRED_METRICS:
                graph_values = []
                flat_values = []
                differences = []
                for sample in ordered_samples:
                    graph_row = graph[sample, graph_target]
                    flat_row = flat[sample, flat_target]
                    graph_value = numeric(graph_row.get(metric, ""))
                    flat_value = numeric(flat_row.get(metric, ""))
                    difference = (
                        graph_value - flat_value
                        if graph_value is not None and flat_value is not None else None
                    )
                    if difference is not None:
                        graph_values.append(graph_value)
                        flat_values.append(flat_value)
                        differences.append(difference)
                    pairs.append({
                        "sample_id": sample,
                        "sample_group": samples[sample]["sample_group"],
                        "batch": samples[sample]["batch"],
                        "graph_target": graph_target,
                        "flat_target": flat_target,
                        "metric": metric,
                        "graph_value": formatted(graph_value),
                        "flat_value": formatted(flat_value),
                        "graph_minus_flat": formatted(difference),
                        "input_reads": flat_row["input_reads"],
                        "input_bases": flat_row["input_bases"],
                    })
                summary.append({
                    "graph_target": graph_target,
                    "flat_target": flat_target,
                    "metric": metric,
                    "n_samples": len(differences),
                    "mean_graph_value": formatted(statistics.mean(graph_values)) if differences else "",
                    "mean_flat_value": formatted(statistics.mean(flat_values)) if differences else "",
                    "mean_difference": formatted(statistics.mean(differences)) if differences else "",
                    "median_difference": formatted(statistics.median(differences)) if differences else "",
                })
    write_tsv(args.pairs_output, pairs, PAIR_FIELDS)
    write_tsv(args.summary_output, summary, SUMMARY_FIELDS)
    render_report(
        args.report_output, args.samples, flat_targets, graph_targets,
        flat, graph, summary, len(samples),
    )


def main() -> None:
    args = parse_args()
    if args.command == "metrics":
        write_metrics(args)
    else:
        compare(args)


if __name__ == "__main__":
    main()
