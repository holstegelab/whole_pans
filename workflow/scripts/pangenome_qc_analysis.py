#!/usr/bin/env python3
"""QC analysis for the initial Minigraph SV pangenome rGFA.

This implements the analysis plan in pangenome/pangenome_qc_plan.md. The GFA is
large, so the parser streams it once and stores a compact segment checkpoint.
"""

import argparse
import csv
import datetime as dt
import gzip
import os
import re
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd


ASSEMBLY_RE = re.compile(
    r"^(?P<sample>.+?)\.hifi\.hifiasm\.bp\.(?P<hap>hap[12])\.p_ctg(?:\.clean)?$"
)
CANONICAL_CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y", "M", "MT"]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a reproducible QC report for a Minigraph rGFA pangenome."
    )
    parser.add_argument("--gfa", required=True)
    parser.add_argument("--graph-summary", required=True)
    parser.add_argument("--ordered-assemblies", required=True)
    parser.add_argument("--mash-distances", required=True)
    parser.add_argument("--tool-versions", required=True)
    parser.add_argument("--build-log", required=True)
    parser.add_argument("--build-benchmark", required=True)
    parser.add_argument("--post-qc-assembly", required=True)
    parser.add_argument("--post-qc-sample", required=True)
    parser.add_argument("--post-qc-included", required=True)
    parser.add_argument("--contamination-summary", required=True)
    parser.add_argument("--vg", default="")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--ln-sample-limit", type=int, default=100000)
    parser.add_argument("--mad-k", type=float, default=3.5)
    parser.add_argument("--low-matching-fraction", type=float, default=0.90)
    return parser.parse_args()


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def ensure_dirs(output_dir):
    output = Path(output_dir)
    for name in ["data", "figures", "tables", "report", "logs", "metadata"]:
        (output / name).mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(output / "metadata" / "mplconfig"))
    os.environ.setdefault("XDG_CACHE_HOME", str(output / "metadata" / "cache"))
    (output / "metadata" / "mplconfig").mkdir(parents=True, exist_ok=True)
    (output / "metadata" / "cache").mkdir(parents=True, exist_ok=True)
    return output


def parse_optional_tags(fields):
    tags = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def parse_int(value):
    if value in {None, ""}:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def normalize_chrom(value):
    if value is None:
        return ""
    chrom = str(value).strip()
    for prefix in ("CHM13.", "GRCh38.", "hg38."):
        if chrom.startswith(prefix):
            chrom = chrom[len(prefix) :]
            break
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    if chrom in {"m", "M"}:
        return "MT"
    return chrom.upper() if chrom in {"x", "y", "mt"} else chrom


def fasta_id(path_or_name):
    name = Path(str(path_or_name)).name
    for suffix in (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name


def clean_base_id(assembly_id):
    return str(assembly_id).removesuffix(".clean")


def parse_sample_hap(assembly_id):
    match = ASSEMBLY_RE.match(str(assembly_id))
    if not match:
        return "", ""
    return match.group("sample"), match.group("hap")


def read_metric_tsv(path):
    metrics = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            metrics[row["metric"]] = row["value"]
    return metrics


def n50(lengths):
    values = sorted((int(value) for value in lengths if int(value) > 0), reverse=True)
    if not values:
        return 0
    half = sum(values) / 2
    running = 0
    for value in values:
        running += value
        if running >= half:
            return value
    return 0


def numeric_quantiles(series, prefix):
    values = pd.Series(series).dropna()
    rows = []
    if values.empty:
        for q in [0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0]:
            rows.append((f"{prefix}_q{int(q * 100):02d}", 0))
        return rows
    for q in [0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0]:
        rows.append((f"{prefix}_q{int(q * 100):02d}", float(values.quantile(q))))
    return rows


def stable_name_conflict_stats(seg_df):
    """Summarize rGFA stable sequence names associated with multiple ranks."""
    empty = {
        "stable_names_total": 0,
        "stable_names_multiple_ranks": 0,
        "stable_name_conflict_segments": 0,
        "stable_name_conflict_bp": 0,
    }
    if seg_df.empty or not {"sn", "sr", "ln"}.issubset(seg_df.columns):
        return empty

    valid = seg_df[(seg_df["sn"].astype(str) != "") & (seg_df["sr"] >= 0)]
    if valid.empty:
        return empty
    ranks_per_name = valid.groupby("sn")["sr"].nunique()
    conflict_names = set(ranks_per_name[ranks_per_name > 1].index)
    affected = valid[valid["sn"].isin(conflict_names)]
    return {
        "stable_names_total": int(len(ranks_per_name)),
        "stable_names_multiple_ranks": int(len(conflict_names)),
        "stable_name_conflict_segments": int(len(affected)),
        "stable_name_conflict_bp": int(affected["ln"].sum()),
    }


def write_csv(path, rows, fieldnames):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def assign_reference_chromosome(seg_df, link_pairs):
    """Assign every segment to a reference chromosome via graph connectivity.

    In a Minigraph rGFA the ``SN`` tag is the reference chromosome only for
    rank-0 (``SR=0``) segments; for ``SR>0`` segments ``SN`` is the origin
    CONTIG name from the donor assembly (e.g. ``h1tg000048l``). Grouping
    non-reference bp directly by ``SN`` therefore produces a per-contig
    breakdown, not a per-chromosome one.

    Each chromosome forms its own connected component anchored by its rank-0
    backbone, so we label every component by the reference chromosome that
    dominates its rank-0 bp and propagate that label to the non-reference
    segments in the same component. Returns (assigned_chrom Series aligned to
    seg_df.index, component Series, n_components, n_reference_components).
    """
    seg_ids = seg_df["segment_id"].to_numpy()
    index = {sid: i for i, sid in enumerate(seg_ids)}
    n = len(seg_ids)
    parent = np.arange(n, dtype=np.int64)

    def find(x):
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:
            parent[x], x = root, parent[x]
        return root

    for left, right in link_pairs:
        li = index.get(left)
        ri = index.get(right)
        if li is None or ri is None:
            continue
        rl, rr = find(li), find(ri)
        if rl != rr:
            parent[rl] = rr

    roots = np.fromiter((find(i) for i in range(n)), dtype=np.int64, count=n)
    comp = pd.Series(roots, index=seg_df.index)

    ref_mask = (seg_df["sr"].to_numpy() == 0)
    ref = pd.DataFrame(
        {
            "comp": roots[ref_mask],
            "chrom": seg_df.loc[ref_mask, "chrom"].to_numpy(),
            "ln": seg_df.loc[ref_mask, "ln"].to_numpy(),
        }
    )
    if ref.empty:
        dominant = pd.Series(dtype=str)
    else:
        dominant = (
            ref.groupby(["comp", "chrom"])["ln"].sum()
            .reset_index()
            .sort_values("ln", ascending=False)
            .drop_duplicates("comp")
            .set_index("comp")["chrom"]
        )
    assigned = comp.map(dominant)
    n_components = int(np.unique(roots).size)
    n_reference_components = int(ref["comp"].nunique()) if not ref.empty else 0
    return assigned, comp, n_components, n_reference_components


def parse_gfa(gfa_path, output_dir, ln_sample_limit):
    segment_rows = []
    segment_ids = set()
    degree = Counter()
    segment_rank = Counter()
    segment_rank_bp = Counter()
    link_rank = Counter()
    link_pairs = []
    record_counts = Counter()
    missing_tags = Counter()
    ln_sampled = 0
    ln_mismatches = 0
    duplicate_segments = 0

    with open_text(gfa_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            record_type = fields[0]
            record_counts[record_type] += 1

            if record_type == "S":
                if len(fields) < 3:
                    missing_tags["malformed_S"] += 1
                    continue
                segment_id = fields[1]
                sequence = fields[2]
                tags = parse_optional_tags(fields[3:])
                for tag in ["SN", "SO", "SR", "LN"]:
                    if tag not in tags:
                        missing_tags[tag] += 1

                sr = parse_int(tags.get("SR"))
                ln = parse_int(tags.get("LN"))
                so = parse_int(tags.get("SO"))
                sn = tags.get("SN", "")
                if sr is None:
                    sr = -1
                if ln is None:
                    ln = len(sequence) if sequence != "*" else 0

                if sequence != "*" and tags.get("LN") is not None and ln_sampled < ln_sample_limit:
                    ln_sampled += 1
                    if ln != len(sequence):
                        ln_mismatches += 1

                if segment_id in segment_ids:
                    duplicate_segments += 1
                segment_ids.add(segment_id)
                segment_rank[sr] += 1
                segment_rank_bp[sr] += ln
                segment_rows.append(
                    {
                        "segment_id": segment_id,
                        "sr": sr,
                        "ln": ln,
                        "sn": sn,
                        "chrom": normalize_chrom(sn),
                        "so": -1 if so is None else so,
                        "is_ref": sr == 0,
                    }
                )

            elif record_type == "L":
                if len(fields) < 5:
                    missing_tags["malformed_L"] += 1
                    continue
                left = fields[1]
                right = fields[3]
                tags = parse_optional_tags(fields[5:])
                sr = parse_int(tags.get("SR"))
                if sr is None:
                    sr = -1
                link_rank[sr] += 1
                degree[left] += 1
                degree[right] += 1
                link_pairs.append((left, right))

    dangling_links = sum(
        1 for left, right in link_pairs if left not in segment_ids or right not in segment_ids
    )
    dangling_endpoints = sum(
        count for segment_id, count in degree.items() if segment_id not in segment_ids
    )

    seg_df = pd.DataFrame(segment_rows)

    # Assign every segment to a reference chromosome via graph connectivity.
    # `chrom` (from the SN tag) is only correct for rank-0 segments; downstream
    # per-chromosome aggregation must use `chrom_assigned`.
    n_components = n_reference_components = 0
    if not seg_df.empty:
        assigned, comp, n_components, n_reference_components = assign_reference_chromosome(
            seg_df, link_pairs
        )
        seg_df["chrom_assigned"] = assigned.to_numpy()
        seg_df["component"] = comp.to_numpy()
    else:
        seg_df["chrom_assigned"] = pd.Series(dtype=str)
        seg_df["component"] = pd.Series(dtype="int64")

    seg_path = output_dir / "data" / "seg_records.parquet"
    seg_df.to_parquet(seg_path, index=False)

    ranks = sorted(set(segment_rank) | set(link_rank))
    rank_rows = [
        {
            "sr": rank,
            "segment_count": segment_rank[rank],
            "segment_bp": segment_rank_bp[rank],
            "link_count": link_rank[rank],
        }
        for rank in ranks
    ]
    rank_tally = pd.DataFrame(rank_rows)
    rank_tally.to_csv(output_dir / "data" / "rank_tally.tsv", sep="\t", index=False)

    degree_rows = [
        {"segment_id": segment_id, "degree": degree[segment_id]}
        for segment_id in sorted(segment_ids)
    ]
    degree_df = pd.DataFrame(degree_rows)
    degree_df.to_parquet(output_dir / "data" / "degree_records.parquet", index=False)

    parse_stats = {
        "segments": len(seg_df),
        "links": record_counts["L"],
        "paths": record_counts["P"] + record_counts["W"],
        "segment_bp": int(seg_df["ln"].sum()) if not seg_df.empty else 0,
        "max_sr": int(seg_df["sr"].max()) if not seg_df.empty else -1,
        "record_counts": record_counts,
        "missing_tags": missing_tags,
        "ln_sampled": ln_sampled,
        "ln_mismatches": ln_mismatches,
        "duplicate_segments": duplicate_segments,
        "dangling_links": dangling_links,
        "dangling_endpoints": dangling_endpoints,
        "n_components": n_components,
        "n_reference_components": n_reference_components,
        **stable_name_conflict_stats(seg_df),
    }
    return seg_df, rank_tally, degree_df, parse_stats


def parse_build_log(path):
    assembly_count = ""
    warnings = 0
    warning_categories = Counter()
    rgfa_consistency_warnings = 0
    real_time = ""
    cpu_time = ""
    peak_rss = ""
    pattern = re.compile(r"Assemblies to add after CHM13 and hg38:\s+(\d+)")
    resource_pattern = re.compile(
        r"Real time:\s+([^;]+)\s+sec;\s+CPU:\s+([^;]+)\s+sec;\s+Peak RSS:\s+(.+)"
    )
    warning_pattern = re.compile(r"^\[W(?:::(?P<category>[^\]]+))?\]")
    with open(path, errors="replace") as handle:
        for line in handle:
            warning_match = warning_pattern.match(line)
            if warning_match:
                warnings += 1
                warning_categories[warning_match.group("category") or "general"] += 1
                if (
                    "inconsistent rGFA" in line
                    or "associated with different ranks" in line
                ):
                    rgfa_consistency_warnings += 1
            match = pattern.search(line)
            if match:
                assembly_count = match.group(1)
            match = resource_pattern.search(line)
            if match:
                real_time, cpu_time, peak_rss = match.groups()
    return {
        "assembly_count": assembly_count,
        "warning_lines": warnings,
        "warning_categories": warning_categories,
        "rgfa_consistency_warnings": rgfa_consistency_warnings,
        "real_time_sec_from_log": real_time,
        "cpu_time_sec_from_log": cpu_time,
        "peak_rss_from_log": peak_rss,
    }


def read_optional_tsv(path):
    if not path or not Path(path).exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def integrity_checks(
    seg_df,
    rank_tally,
    parse_stats,
    args,
    ordered,
    mash,
    post_sample,
    post_included,
    build_log_info,
):
    rows = []
    summary = read_metric_tsv(args.graph_summary)

    def add(check, expected, observed, status, details=""):
        rows.append(
            {
                "check": check,
                "expected": str(expected),
                "observed": str(observed),
                "status": status,
                "details": details,
            }
        )

    for metric in ["segments", "links", "paths", "segment_bp"]:
        expected = summary.get(metric, "")
        observed = parse_stats[metric]
        add(
            f"parser_matches_graph_summary_{metric}",
            expected,
            observed,
            "PASS" if str(observed) == str(expected) else "FAIL",
        )

    expected_max_rank = len(ordered) - 1
    max_rank = parse_stats["max_sr"]
    add(
        "max_sr_rank_matches_ordered_sources",
        expected_max_rank,
        max_rank,
        "PASS" if max_rank == expected_max_rank else "FAIL",
    )

    observed_ranks = set(int(value) for value in rank_tally["sr"] if int(value) >= 0)
    missing_ranks = sorted(set(range(expected_max_rank + 1)) - observed_ranks)
    add(
        "all_expected_ranks_contribute_segments",
        "no missing ranks",
        ",".join(map(str, missing_ranks)) if missing_ranks else "none",
        "PASS" if not missing_ranks else "WARN",
    )

    add(
        "link_endpoint_referential_integrity",
        0,
        parse_stats["dangling_links"],
        "PASS" if parse_stats["dangling_links"] == 0 else "FAIL",
        f"dangling endpoint references: {parse_stats['dangling_endpoints']}",
    )

    for tag in ["SN", "SO", "SR", "LN"]:
        missing = parse_stats["missing_tags"][tag]
        add(
            f"segment_tag_{tag}_complete",
            0,
            missing,
            "PASS" if missing == 0 else "FAIL",
        )

    add(
        "sampled_ln_matches_sequence_length",
        0,
        parse_stats["ln_mismatches"],
        "PASS" if parse_stats["ln_mismatches"] == 0 else "FAIL",
        f"sampled S records with sequence present: {parse_stats['ln_sampled']}",
    )
    add(
        "duplicate_segment_ids",
        0,
        parse_stats["duplicate_segments"],
        "PASS" if parse_stats["duplicate_segments"] == 0 else "FAIL",
    )
    stable_name_conflicts = parse_stats.get("stable_names_multiple_ranks", 0)
    add(
        "stable_sequence_names_have_single_source_rank",
        0,
        stable_name_conflicts,
        "PASS" if stable_name_conflicts == 0 else "FAIL",
        (
            f"affected segments: {parse_stats.get('stable_name_conflict_segments', 0)}; "
            f"affected bp: {parse_stats.get('stable_name_conflict_bp', 0)}"
        ),
    )
    rgfa_log_warnings = build_log_info.get("rgfa_consistency_warnings", 0)
    add(
        "build_log_rgfa_consistency_warnings",
        0,
        rgfa_log_warnings,
        "PASS" if rgfa_log_warnings == 0 else "FAIL",
        "Warnings containing inconsistent rGFA or source-rank conflicts.",
    )
    warning_lines = build_log_info.get("warning_lines", 0)
    add(
        "build_log_warning_lines",
        0,
        warning_lines,
        "PASS" if warning_lines == 0 else "WARN",
        "; ".join(
            f"{name}={count}"
            for name, count in sorted(
                build_log_info.get("warning_categories", {}).items()
            )
        ),
    )

    assembly_rows = ordered[ordered["role"] == "assembly"].copy()
    expected_assembly_count = len(post_included)
    add(
        "ordered_assembly_count",
        expected_assembly_count,
        len(assembly_rows),
        "PASS" if len(assembly_rows) == expected_assembly_count else "FAIL",
    )
    build_assembly_count = build_log_info.get("assembly_count", "")
    add(
        "build_log_assembly_count",
        expected_assembly_count,
        build_assembly_count,
        "PASS" if build_assembly_count == str(expected_assembly_count) else "FAIL",
    )

    ordered_ids = set(assembly_rows["assembly_id"].map(str))
    included_ids = set(fasta_id(path) for path in post_included)
    missing_from_graph = sorted(included_ids - ordered_ids)
    extra_in_graph = sorted(ordered_ids - included_ids)
    add(
        "post_qc_included_set_matches_graph_assemblies",
        f"same {expected_assembly_count} assembly IDs",
        f"missing={len(missing_from_graph)}; extra={len(extra_in_graph)}",
        "PASS" if not missing_from_graph and not extra_in_graph else "FAIL",
        (
            f"missing: {','.join(missing_from_graph[:10])}; "
            f"extra: {','.join(extra_in_graph[:10])}"
        ),
    )

    ordered_mash = pd.to_numeric(
        assembly_rows["mash_distance"].replace("NA", np.nan), errors="coerce"
    )
    inversions = int((ordered_mash.diff().dropna() < 0).sum())
    add(
        "ordered_assemblies_mash_distance_monotone",
        0,
        inversions,
        "PASS" if inversions == 0 else "FAIL",
    )

    n_comp = parse_stats.get("n_components", 0)
    n_ref_comp = parse_stats.get("n_reference_components", 0)
    add(
        "connected_components_equal_reference_chromosomes",
        n_ref_comp,
        n_comp,
        "PASS" if n_comp == n_ref_comp and n_comp > 0 else "WARN",
        (
            "Each chromosome should form one component anchored by its rank-0 "
            "backbone; equality validates the connectivity-based chromosome "
            "assignment used for per-chromosome non-reference content."
        ),
    )

    if not post_sample.empty and "sample_status" in post_sample:
        counts = post_sample["sample_status"].value_counts().to_dict()
        for status in ["PASS", "WARN", "PARTIAL", "FAIL"]:
            add(
                f"post_qc_sample_status_{status}",
                "reported",
                counts.get(status, 0),
                "INFO",
            )

    checks = pd.DataFrame(rows)
    checks.to_csv(Path(args.output_dir) / "tables" / "integrity_checks.csv", index=False)
    return checks


def graph_overview(seg_df, degree_df, parse_stats, output_dir):
    total_bp = int(seg_df["ln"].sum())
    ref = seg_df[seg_df["sr"] == 0]
    nonref = seg_df[seg_df["sr"] > 0]
    degree_values = degree_df["degree"] if not degree_df.empty else pd.Series(dtype=int)
    rows = []

    def add(metric, value):
        rows.append({"metric": metric, "value": value})

    add("segments", len(seg_df))
    add("links", parse_stats["links"])
    add("paths_or_walks", parse_stats["paths"])
    add("ranks_contributing_nonref_segments", int(nonref["sr"].nunique()))
    add("total_segment_bp", total_bp)
    add("reference_backbone_bp", int(ref["ln"].sum()))
    add("non_reference_bp", int(nonref["ln"].sum()))
    add(
        "non_reference_bp_percent",
        round(100 * int(nonref["ln"].sum()) / total_bp, 6) if total_bp else 0,
    )
    add("mean_segment_length_bp", round(float(seg_df["ln"].mean()), 3))
    add("median_segment_length_bp", round(float(seg_df["ln"].median()), 3))
    add("segment_n50_bp", n50(seg_df["ln"]))
    add("reference_segment_n50_bp", n50(ref["ln"]))
    add("nonref_segment_n50_bp", n50(nonref["ln"]) if not nonref.empty else 0)
    add("largest_segment_bp", int(seg_df["ln"].max()) if not seg_df.empty else 0)
    add("largest_nonref_segment_bp", int(nonref["ln"].max()) if not nonref.empty else 0)
    add("connected_components", parse_stats.get("n_components", 0))
    add("reference_chromosome_components", parse_stats.get("n_reference_components", 0))
    add("link_segment_ratio", round(parse_stats["links"] / len(seg_df), 6))
    add("mean_degree", round(float(degree_values.mean()), 6) if len(degree_values) else 0)
    add("median_degree", round(float(degree_values.median()), 6) if len(degree_values) else 0)
    add("max_degree", int(degree_values.max()) if len(degree_values) else 0)
    add("hub_nodes_degree_ge_4", int((degree_values >= 4).sum()))

    for metric, value in numeric_quantiles(seg_df["ln"], "segment_length_bp"):
        add(metric, round(value, 3))
    for metric, value in numeric_quantiles(degree_values, "degree"):
        add(metric, round(value, 3))

    bins = [
        ("lt_50bp", 0, 50),
        ("50bp_1kb", 50, 1000),
        ("1kb_10kb", 1000, 10000),
        ("gt_10kb", 10000, None),
        ("gt_100kb", 100000, None),
    ]
    for label, low, high in bins:
        for group_name, frame in [("ref", ref), ("nonref", nonref)]:
            if high is None:
                count = int((frame["ln"] > low).sum())
            else:
                count = int(((frame["ln"] >= low) & (frame["ln"] < high)).sum())
            add(f"{group_name}_segments_{label}", count)

    overview = pd.DataFrame(rows)
    overview.to_csv(output_dir / "tables" / "graph_overview_stats.csv", index=False)
    return overview


def per_source_contribution(rank_tally, ordered, output_dir):
    rank = rank_tally.rename(columns={"sr": "rank"}).copy()
    ordered = ordered.copy()
    ordered["rank"] = pd.to_numeric(ordered["order"], errors="coerce").astype("Int64") - 1
    merged = ordered.merge(rank, on="rank", how="left")
    for column in ["segment_count", "segment_bp", "link_count"]:
        merged[column] = pd.to_numeric(merged[column], errors="coerce").fillna(0).astype(int)
    merged["novel_segment_count"] = np.where(merged["rank"] > 0, merged["segment_count"], 0)
    merged["novel_bp"] = np.where(merged["rank"] > 0, merged["segment_bp"], 0)
    merged["cumulative_nonref_bp"] = merged["novel_bp"].cumsum()
    no_hg38 = merged["novel_bp"].where(merged["rank"] != 1, 0)
    merged["cumulative_nonref_bp_without_hg38"] = no_hg38.cumsum()
    merged["mash_distance"] = pd.to_numeric(
        merged["mash_distance"].replace("NA", np.nan), errors="coerce"
    )
    keep = [
        "rank",
        "order",
        "role",
        "assembly_id",
        "path",
        "mash_distance",
        "segment_count",
        "segment_bp",
        "link_count",
        "novel_segment_count",
        "novel_bp",
        "cumulative_nonref_bp",
        "cumulative_nonref_bp_without_hg38",
    ]
    merged[keep].to_csv(output_dir / "tables" / "per_source_contribution.csv", index=False)

    fit = {"heaps_a": np.nan, "heaps_b": np.nan, "last_decile_mean_novel_bp": np.nan}
    fit_frame = merged[(merged["rank"] > 1) & (merged["cumulative_nonref_bp_without_hg38"] > 0)]
    if len(fit_frame) >= 3:
        x = np.arange(1, len(fit_frame) + 1, dtype=float)
        y = fit_frame["cumulative_nonref_bp_without_hg38"].to_numpy(dtype=float)
        slope, intercept = np.polyfit(np.log(x), np.log(y), 1)
        fit["heaps_a"] = float(np.exp(intercept))
        fit["heaps_b"] = float(slope)
        last_decile = max(1, int(np.ceil(len(fit_frame) * 0.1)))
        fit["last_decile_mean_novel_bp"] = float(fit_frame["novel_bp"].tail(last_decile).mean())
    return merged, fit


def per_chromosome_nonref(seg_df, output_dir):
    """Per-chromosome reference and non-reference content.

    Non-reference segments are grouped by ``chrom_assigned`` (the reference
    chromosome inferred from graph connectivity), NOT by the raw ``SN`` tag,
    which for donor segments is a contig name. Reference length per chromosome
    is the summed rank-0 backbone bp (each chromosome is one connected
    component, so this equals the CHM13 chromosome length in the graph).
    """
    group_col = "chrom_assigned" if "chrom_assigned" in seg_df.columns else "chrom"
    ref = seg_df[seg_df["sr"] == 0].copy()
    nonref = seg_df[seg_df["sr"] > 0].copy()

    if ref.empty:
        ref_lengths = pd.DataFrame(columns=["chrom", "ref_bp"])
    else:
        ref_lengths = ref.groupby(group_col, as_index=False)["ln"].sum()
        ref_lengths = ref_lengths.rename(columns={group_col: "chrom", "ln": "ref_bp"})

    if nonref.empty:
        nonref_summary = pd.DataFrame(columns=["chrom", "nonref_bp", "nonref_segments"])
    else:
        nonref_summary = (
            nonref.groupby(group_col, as_index=False)
            .agg(nonref_bp=("ln", "sum"), nonref_segments=("segment_id", "count"))
            .rename(columns={group_col: "chrom"})
        )

    chrom = ref_lengths.merge(nonref_summary, on="chrom", how="outer").fillna(0)
    chrom = chrom[chrom["chrom"].astype(str) != ""]
    chrom["ref_bp"] = chrom["ref_bp"].astype(int)
    chrom["nonref_bp"] = chrom["nonref_bp"].astype(int)
    chrom["nonref_segments"] = chrom["nonref_segments"].astype(int)
    chrom["nonref_bp_per_mb"] = np.where(
        chrom["ref_bp"] > 0, chrom["nonref_bp"] / (chrom["ref_bp"] / 1_000_000), np.nan
    )

    def sort_key(name):
        stripped = normalize_chrom(name)
        try:
            return CANONICAL_CHROMS.index(stripped)
        except ValueError:
            return 10_000

    chrom["sort_key"] = chrom["chrom"].map(sort_key)
    chrom = chrom.sort_values(["sort_key", "chrom"]).drop(columns=["sort_key"])
    chrom.to_csv(output_dir / "tables" / "per_chromosome_nonref.csv", index=False)
    return chrom


def mash_outliers(args, mash, assembly_qc, contamination, output_dir):
    if mash.empty:
        out = pd.DataFrame()
        out.to_csv(output_dir / "tables" / "mash_outliers.csv", index=False)
        return out, pd.DataFrame()

    mash = mash.copy()
    mash["mash_distance"] = pd.to_numeric(mash["mash_distance"], errors="coerce")
    hashes = mash["matching_hashes"].str.extract(r"(?P<matching>\d+)/(?P<total>\d+)")
    mash["matching_hashes_n"] = pd.to_numeric(hashes["matching"], errors="coerce")
    mash["sketch_hashes_n"] = pd.to_numeric(hashes["total"], errors="coerce")
    mash["matching_hash_fraction"] = mash["matching_hashes_n"] / mash["sketch_hashes_n"]
    parsed = mash["assembly_id"].map(parse_sample_hap)
    mash["sample"] = [item[0] for item in parsed]
    mash["haplotype"] = [item[1] for item in parsed]
    median = float(mash["mash_distance"].median())
    mad = float((mash["mash_distance"] - median).abs().median())
    if mad == 0:
        mash["mash_mad_score"] = 0.0
        high_threshold = median
    else:
        mash["mash_mad_score"] = (mash["mash_distance"] - median) / mad
        high_threshold = median + args.mad_k * mad
    mash["high_mash_distance"] = mash["mash_distance"] > high_threshold
    mash["low_matching_hash_fraction"] = mash["matching_hash_fraction"] < args.low_matching_fraction

    if not assembly_qc.empty:
        qc = assembly_qc.copy()
        qc["assembly_id_clean"] = qc["assembly_id"].map(clean_base_id)
        mash["assembly_id_clean"] = mash["assembly_id"].map(clean_base_id)
        mash = mash.merge(
            qc[
                [
                    "assembly_id_clean",
                    "assembly_status",
                    "warning_reasons",
                    "fail_reasons",
                    "contig_n50_bp",
                    "best_reference_covered_percent",
                ]
            ],
            on="assembly_id_clean",
            how="left",
        )

    if not contamination.empty:
        contam = contamination.copy()
        contam["assembly_id_clean"] = contam["assembly_id"].map(clean_base_id)
        if "assembly_id_clean" not in mash:
            mash["assembly_id_clean"] = mash["assembly_id"].map(clean_base_id)
        keep_cols = [
            col
            for col in [
                "assembly_id_clean",
                "removed_contigs",
                "removed_whole_contig_bp",
                "review_contigs",
                "review_bp",
            ]
            if col in contam.columns
        ]
        mash = mash.merge(contam[keep_cols], on="assembly_id_clean", how="left")

    paired = (
        mash.pivot_table(index="sample", columns="haplotype", values="mash_distance", aggfunc="first")
        .reset_index()
        .rename_axis(None, axis=1)
    )
    if {"hap1", "hap2"}.issubset(paired.columns):
        paired["hap_abs_difference"] = (paired["hap1"] - paired["hap2"]).abs()

    flagged = mash[mash["high_mash_distance"] | mash["low_matching_hash_fraction"]].copy()
    flagged.to_csv(output_dir / "tables" / "mash_outliers.csv", index=False)
    return flagged, paired


def setup_plotting(output_dir):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme(style="whitegrid", context="talk")
    plt.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
        }
    )
    return plt, sns


def plot_segment_lengths(seg_df, output_dir):
    plt, sns = setup_plotting(output_dir)
    frame = seg_df[seg_df["ln"] > 0].copy()
    ref_ln = frame.loc[frame["sr"] == 0, "ln"].to_numpy()
    nonref_ln = frame.loc[frame["sr"] > 0, "ln"].to_numpy()
    # Log-spaced bins: seaborn's linear default (bins=80) collapses a 1 bp -
    # 32 Mb range into one visible block. Build the bins explicitly on a log
    # scale so the distribution shape is visible.
    max_ln = int(frame["ln"].max()) if not frame.empty else 1
    bins = np.logspace(0, np.log10(max_ln) + 0.01, 60)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist(ref_ln, bins=bins, histtype="step", linewidth=1.6, color="#4C72B0",
            label=f"reference SR=0 (n={len(ref_ln):,})")
    ax.hist(nonref_ln, bins=bins, histtype="step", linewidth=1.6, color="#DD8452",
            label=f"non-reference SR>0 (n={len(nonref_ln):,})")
    ax.set_xscale("log")
    ax.set_yscale("log")
    for arr, color in [(ref_ln, "#4C72B0"), (nonref_ln, "#DD8452")]:
        if len(arr):
            ax.axvline(np.median(arr), color=color, ls=":", lw=1)
    ax.set_xlabel("Segment length (bp, log scale)")
    ax.set_ylabel("Segment count (log scale)")
    ax.set_title("Segment length distribution")
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(output_dir / "figures" / "segment_length_distribution.png")
    plt.close(fig)


def plot_growth(per_source, output_dir):
    plt, _sns = setup_plotting(output_dir)
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    x = per_source["rank"].astype(int)
    axes[0].plot(x, per_source["cumulative_nonref_bp"] / 1e6, label="including hg38")
    axes[0].plot(
        x,
        per_source["cumulative_nonref_bp_without_hg38"] / 1e6,
        label="hg38 contribution removed",
    )
    axes[0].set_ylabel("Cumulative non-reference Mb")
    axes[0].set_title("Pangenome growth curve")
    axes[0].legend()
    colors = np.where(per_source["rank"] == 1, "#d95f02", "#1b9e77")
    axes[1].bar(x, per_source["novel_bp"] / 1e6, color=colors, width=0.85)
    axes[1].set_xlabel("SR rank / source addition order")
    axes[1].set_ylabel("Novel Mb first introduced")
    axes[1].set_title("Per-source contribution")
    fig.tight_layout()
    fig.savefig(output_dir / "figures" / "pangenome_growth_curve.png")
    plt.close(fig)


def plot_chromosome(chrom, output_dir):
    plt, _sns = setup_plotting(output_dir)
    plot_df = chrom[(chrom["chrom"].astype(str) != "") & (chrom["nonref_bp"] > 0)].copy()
    # Canonical chromosome order (chr1..22, X, Y), most-recent at top of a
    # horizontal bar reads top-to-bottom, so reverse for display.
    order = {name: i for i, name in enumerate(CANONICAL_CHROMS)}
    plot_df["__o"] = plot_df["chrom"].map(lambda c: order.get(normalize_chrom(c), 10_000))
    plot_df = plot_df.sort_values("__o", ascending=False).drop(columns="__o")
    labels = ["chr" + str(c) if not str(c).startswith("chr") else str(c)
              for c in plot_df["chrom"]]
    fig, axes = plt.subplots(1, 2, figsize=(12, max(5, 0.32 * len(plot_df))), sharey=True)
    axes[0].barh(labels, plot_df["nonref_bp"] / 1e6, color="#4C72B0",
                 edgecolor="white", linewidth=0.4)
    axes[0].set_xlabel("Non-reference sequence (Mb)")
    axes[0].set_ylabel("Chromosome")
    axes[0].set_title("Total novel sequence per chromosome")
    axes[1].barh(labels, plot_df["nonref_bp_per_mb"] / 1e3, color="#55A868",
                 edgecolor="white", linewidth=0.4)
    axes[1].set_xlabel("Non-reference density (kb per reference Mb)")
    axes[1].set_title("Novel sequence density")
    fig.suptitle("Non-reference content assigned to chromosome via graph connectivity",
                 y=1.0, fontsize=11)
    fig.tight_layout()
    fig.savefig(output_dir / "figures" / "per_chromosome_nonref.png")
    plt.close(fig)


def plot_mash(mash, paired, output_dir):
    plt, sns = setup_plotting(output_dir)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    sns.histplot(data=mash, x="mash_distance", hue="haplotype", bins=35, ax=axes[0])
    axes[0].set_title("Mash distance to CHM13")
    axes[0].set_xlabel("Mash distance")
    if {"hap1", "hap2"}.issubset(paired.columns):
        axes[1].scatter(paired["hap1"], paired["hap2"], s=28, alpha=0.8)
        low = min(paired["hap1"].min(), paired["hap2"].min())
        high = max(paired["hap1"].max(), paired["hap2"].max())
        axes[1].plot([low, high], [low, high], color="black", linewidth=1)
        axes[1].set_xlabel("hap1 Mash distance")
        axes[1].set_ylabel("hap2 Mash distance")
    axes[1].set_title("Haplotype-pair concordance")
    fig.tight_layout()
    fig.savefig(output_dir / "figures" / "mash_distance_diversity.png")
    plt.close(fig)


def record_tool_versions(args, output_dir):
    rows = []

    def add(tool, version, command=""):
        rows.append({"tool": tool, "version": version, "command": command})

    add("python", sys.version.split()[0], "python --version")
    for module_name, module in [("pandas", pd), ("numpy", np)]:
        add(module_name, getattr(module, "__version__", "unknown"))
    try:
        import matplotlib
        import seaborn
        import pyarrow

        add("matplotlib", matplotlib.__version__)
        add("seaborn", seaborn.__version__)
        add("pyarrow", pyarrow.__version__)
    except Exception as exc:
        add("plot/parquet_dependency_error", str(exc))

    for tool, command in [
        ("gfatools", ["gfatools"]),
        ("vg", [args.vg, "version"] if args.vg else ["vg", "version"]),
    ]:
        executable = command[0]
        if not executable or not (Path(executable).exists() or shutil.which(executable)):
            add(tool, "not found", " ".join(command))
            continue
        try:
            completed = subprocess.run(
                command,
                text=True,
                capture_output=True,
                check=False,
                timeout=20,
            )
            text = (completed.stdout or completed.stderr).strip().splitlines()
            add(tool, text[0] if text else "available", " ".join(command))
        except Exception as exc:
            add(tool, f"version check failed: {exc}", " ".join(command))

    tool_versions = pd.read_csv(args.tool_versions, sep="\t")
    for _, row in tool_versions.iterrows():
        add(f"build_{row['tool']}", row["version"], "recorded during graph build")

    pd.DataFrame(rows).to_csv(output_dir / "metadata" / "qc_tool_versions.tsv", sep="\t", index=False)
    return pd.DataFrame(rows)


def write_manifest(output_dir):
    rows = []
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and "mplconfig" not in path.parts:
            rows.append(
                {
                    "path": str(path.relative_to(output_dir)),
                    "bytes": path.stat().st_size,
                    "updated": dt.datetime.fromtimestamp(path.stat().st_mtime).isoformat(),
                }
            )
    pd.DataFrame(rows).to_csv(output_dir / "metadata" / "run_manifest.tsv", sep="\t", index=False)


def write_report(
    output_dir,
    checks,
    overview,
    per_source,
    heaps_fit,
    chrom,
    mash_out,
    paired,
    build_log_info,
    benchmark,
):
    status_counts = checks["status"].value_counts().to_dict()
    metric = dict(zip(overview["metric"], overview["value"]))
    total_nonref = metric.get("non_reference_bp", 0)
    total_bp = metric.get("total_segment_bp", 0)
    top_chrom = chrom.sort_values("nonref_bp", ascending=False).head(5)
    last_novel = heaps_fit.get("last_decile_mean_novel_bp", np.nan)
    if pd.notna(last_novel) and last_novel > 1_000_000:
        openness = "still adding megabase-scale sequence in the last decile"
    else:
        openness = "showing partial saturation by this simple last-decile metric"
    if status_counts.get("FAIL", 0):
        integrity_verdict = "FAIL: do not use this graph downstream until failures are resolved"
    elif status_counts.get("WARN", 0):
        integrity_verdict = "WARN: structurally usable only after warning review"
    else:
        integrity_verdict = "PASS: no integrity failures or warnings detected"
    warning_categories = build_log_info.get("warning_categories", {})
    warning_breakdown = "; ".join(
        f"{name}={count}" for name, count in sorted(warning_categories.items())
    ) or "none"

    lines = [
        "# Pangenome QC Report",
        "",
        "## Build Provenance",
        "",
        "- Graph type: Minigraph rGFA, SV-level `-cxggs` graph.",
        f"- Build log assemblies after CHM13/hg38: {build_log_info.get('assembly_count', 'not found')}.",
        f"- Build log warning lines: {build_log_info.get('warning_lines', 0)} ({warning_breakdown}).",
        f"- rGFA consistency warnings: {build_log_info.get('rgfa_consistency_warnings', 0)}.",
    ]
    if not benchmark.empty:
        row = benchmark.iloc[0].to_dict()
        lines.extend(
            [
                f"- Benchmark wall time: {row.get('h:m:s', row.get('s', 'NA'))}.",
                f"- Benchmark max RSS MB: {row.get('max_rss', 'NA')}.",
                f"- Benchmark CPU time seconds: {row.get('cpu_time', 'NA')}.",
            ]
        )

    lines.extend(
        [
            "",
            "## Integrity Verdict",
            "",
            f"- Overall verdict: **{integrity_verdict}**.",
            (
                f"- Checks: PASS={status_counts.get('PASS', 0)}, "
                f"WARN={status_counts.get('WARN', 0)}, "
                f"FAIL={status_counts.get('FAIL', 0)}, INFO={status_counts.get('INFO', 0)}."
            ),
            "- See [integrity_checks.csv](../tables/integrity_checks.csv).",
            "",
            "## Topology And Content",
            "",
            f"- Segments: {metric.get('segments', 'NA')}; links: {metric.get('links', 'NA')}.",
            f"- Total segment bp: {total_bp}; non-reference bp: {total_nonref}.",
            f"- Non-reference fraction: {metric.get('non_reference_bp_percent', 'NA')}%.",
            f"- Segment N50: {metric.get('segment_n50_bp', 'NA')} bp.",
            f"- Hub nodes with degree >= 4: {metric.get('hub_nodes_degree_ge_4', 'NA')}.",
            "",
            "![Segment length distribution](../figures/segment_length_distribution.png)",
            "",
            "## Non-reference Content And Openness",
            "",
            (
                "- Novel sequence is assigned to the rank that first introduced it. "
                "It is not per-sample uniqueness because the rGFA has no paths/walks."
            ),
            f"- Heaps-style exponent without hg38: {heaps_fit.get('heaps_b', np.nan):.4g}.",
            f"- Mean novel bp/source in the last decile: {last_novel:.3f}.",
            f"- Simple interpretation: the graph is {openness}.",
            "",
            "![Pangenome growth curve](../figures/pangenome_growth_curve.png)",
            "",
            "## Chromosomal Distribution",
            "",
            (
                "- Non-reference segments are assigned to a reference chromosome "
                "via graph connectivity (each chromosome is one connected "
                "component anchored by its rank-0 backbone), not by the raw `SN` "
                "tag, which for donor segments is a contig name."
            ),
            "- Top chromosomes by non-reference bp:",
        ]
    )
    for _, row in top_chrom.iterrows():
        chrom_name = str(row["chrom"])
        if not chrom_name.startswith("chr"):
            chrom_name = "chr" + chrom_name
        lines.append(
            f"  - {chrom_name}: {int(row['nonref_bp'])} bp "
            f"({row['nonref_bp_per_mb']:.0f} bp per reference Mb)"
        )
    lines.extend(
        [
            "",
            "![Per-chromosome non-reference content](../figures/per_chromosome_nonref.png)",
            "",
            "## Input Diversity And Outliers",
            "",
            f"- Mash outlier rows flagged: {len(mash_out)}.",
            f"- Paired samples with hap1/hap2 Mash distances: {len(paired)}.",
            "- See [mash_outliers.csv](../tables/mash_outliers.csv).",
            "",
            "![Mash distance diversity](../figures/mash_distance_diversity.png)",
            "",
            "## Limitations",
            "",
            "- Minigraph rGFA has no per-sample `P`/`W` traversals here.",
            "- Per-sample graph usage requires `minigraph --call`, `vg giraffe`, or another mapping step.",
            "- `-cxggs` is SV-level; SNPs and small indels are out of scope.",
            "- Original `/gpfs/...` paths in metadata are treated as metadata on this mount.",
            "",
            "## Recommended Next Steps",
            "",
            "- Review any `FAIL` rows in `integrity_checks.csv` before downstream use.",
            "- Review `mash_outliers.csv` jointly with assembly QC and decontamination summaries.",
            "- Run optional `gfatools bubble` if variant-site inventory is needed.",
        ]
    )
    report = output_dir / "report" / "pangenome_qc_report.md"
    report.write_text("\n".join(lines) + "\n")


def main():
    args = parse_args()
    output_dir = ensure_dirs(args.output_dir)

    seg_df, rank_tally, degree_df, parse_stats = parse_gfa(
        args.gfa, output_dir, args.ln_sample_limit
    )
    ordered = pd.read_csv(args.ordered_assemblies, sep="\t", dtype=str).fillna("")
    mash = pd.read_csv(args.mash_distances, sep="\t", dtype=str).fillna("")
    post_sample = read_optional_tsv(args.post_qc_sample)
    assembly_qc = read_optional_tsv(args.post_qc_assembly)
    contamination = read_optional_tsv(args.contamination_summary)
    with open(args.post_qc_included) as handle:
        post_included = [line.strip() for line in handle if line.strip()]
    build_log_info = parse_build_log(args.build_log)

    checks = integrity_checks(
        seg_df,
        rank_tally,
        parse_stats,
        args,
        ordered,
        mash,
        post_sample,
        post_included,
        build_log_info,
    )
    overview = graph_overview(seg_df, degree_df, parse_stats, output_dir)
    per_source, heaps_fit = per_source_contribution(rank_tally, ordered, output_dir)
    chrom = per_chromosome_nonref(seg_df, output_dir)
    mash_out, paired = mash_outliers(args, mash, assembly_qc, contamination, output_dir)
    if Path(args.build_benchmark).exists():
        benchmark = pd.read_csv(args.build_benchmark, sep="\t")
    else:
        benchmark = pd.DataFrame()

    plot_segment_lengths(seg_df, output_dir)
    plot_growth(per_source, output_dir)
    plot_chromosome(chrom, output_dir)
    parsed = mash["assembly_id"].map(parse_sample_hap)
    mash_plot = mash.copy()
    mash_plot["haplotype"] = [item[1] for item in parsed]
    mash_plot["mash_distance"] = pd.to_numeric(mash_plot["mash_distance"], errors="coerce")
    plot_mash(mash_plot, paired, output_dir)

    record_tool_versions(args, output_dir)
    write_report(
        output_dir,
        checks,
        overview,
        per_source,
        heaps_fit,
        chrom,
        mash_out,
        paired,
        build_log_info,
        benchmark,
    )
    write_manifest(output_dir)
    print(f"Wrote pangenome QC analysis to {output_dir}")


if __name__ == "__main__":
    main()
