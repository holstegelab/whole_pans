#!/usr/bin/env python3
"""Compare two completed seeded Minigraph-Cactus builds without loading graphs."""

import argparse
import csv
import gzip
import io
import json
import sys
from pathlib import Path

from minigraph_cactus_common import atomic_bytes, read_json


def vcf_samples(path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#CHROM\t"):
                fields = line.rstrip("\n").split("\t")
                return fields[9:] if len(fields) > 9 else []
    raise ValueError(f"VCF has no #CHROM header: {path}")


def graph_stats(path):
    rows = {}
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"#ref_chrom", "nodes", "edges", "length"}.issubset(reader.fieldnames or []):
            raise ValueError(f"unexpected Cactus graph-stats header: {path}")
        for row in reader:
            chrom = row["#ref_chrom"]
            if not chrom or chrom in rows:
                raise ValueError(f"empty or duplicate reference chromosome in {path}")
            rows[chrom] = {key: (None if row[key] == "NA" else int(row[key]))
                           for key in ("nodes", "edges", "length")}
    if not rows:
        raise ValueError(f"Cactus graph-stats has no chromosome rows: {path}")
    return rows


def build_metrics(settings):
    root = Path(settings["results_root"]) / settings["build_id"]
    preflight_path = root / "metadata" / "preflight.json"
    preflight = read_json(preflight_path)
    seed = read_json(root / "stages" / "minigraph" / "seed-normalization.json")
    graph = root / "graph"
    prefix = settings["out_name"]
    stats_path = graph / f"{prefix}.stats" / "graph-stats.tsv"
    chrom_stats = graph_stats(stats_path)
    files = {"clipped_gfa_bytes": graph / f"{prefix}.gfa.gz",
             "gbz_bytes": graph / f"{prefix}.gbz",
             "vcf_bytes": graph / f"{prefix}.vcf.gz",
             "raw_vcf_bytes": graph / f"{prefix}.raw.vcf.gz",
             "hal_bytes": graph / f"{prefix}.full.hal",
             "giraffe_gbz_bytes": graph / f"{prefix}.d2.gbz"}
    metrics = {"haplotypes": preflight["haplotype_count"],
               "local_haplotypes": preflight["local_haplotype_count"],
               "input_bases": preflight["total_input_bases"],
               "seed_segments": seed["counts"]["segments"],
               "seed_links": seed["counts"]["links"],
               "seed_segment_bp": seed["counts"]["segment_bp"],
               "clipped_nodes": sum(row["nodes"] or 0 for row in chrom_stats.values()),
               "clipped_edges": sum(row["edges"] or 0 for row in chrom_stats.values()),
               "clipped_graph_bp": sum(row["length"] or 0 for row in chrom_stats.values()),
               "vcf_samples": len(vcf_samples(graph / f"{prefix}.vcf.gz"))}
    metrics.update({name: path.stat().st_size for name, path in files.items()})
    return metrics, preflight, seed, stats_path, chrom_stats


def render(local_settings, hprc_settings, out_prefix):
    local, local_preflight, local_seed, local_stats, local_chroms = build_metrics(local_settings)
    hprc, hprc_preflight, hprc_seed, hprc_stats, hprc_chroms = build_metrics(hprc_settings)
    if local["local_haplotypes"] != hprc["local_haplotypes"]:
        raise ValueError("builds contain different numbers of local V1 haplotypes")
    local_members = {row["name"]: (row["path"], row["sha256"])
                     for row in local_preflight["inputs"] if row["name"] not in {"CHM13", "GRCh38"}}
    combined_members = {row["name"]: (row["path"], row["sha256"])
                        for row in hprc_preflight["inputs"] if row["name"] in local_members}
    if local_members != combined_members:
        raise ValueError("HPRC+local build does not contain exactly the local V1 FASTA content")
    rows = []
    for name in local:
        left, right = local[name], hprc[name]
        rows.append((name, left, right, right - left,
                     f"{right / left:.3f}" if left else "NA"))
    table = io.StringIO()
    writer = csv.writer(table, delimiter="\t", lineterminator="\n")
    writer.writerow(["metric", "local_v1", "hprc_plus_local_v1", "difference", "ratio"])
    writer.writerows(rows)
    out = Path(out_prefix)
    atomic_bytes(str(out) + ".tsv", table.getvalue().encode())
    chrom_table = io.StringIO()
    chrom_writer = csv.writer(chrom_table, delimiter="\t", lineterminator="\n")
    chrom_writer.writerow(["reference_chromosome", "local_nodes", "hprc_plus_local_nodes",
                           "local_edges", "hprc_plus_local_edges", "local_graph_bp",
                           "hprc_plus_local_graph_bp"])
    for chrom in sorted(set(local_chroms) | set(hprc_chroms)):
        left, right = local_chroms.get(chrom, {}), hprc_chroms.get(chrom, {})
        chrom_writer.writerow([chrom, left.get("nodes") if left.get("nodes") is not None else "NA",
                               right.get("nodes") if right.get("nodes") is not None else "NA",
                               left.get("edges") if left.get("edges") is not None else "NA",
                               right.get("edges") if right.get("edges") is not None else "NA",
                               left.get("length") if left.get("length") is not None else "NA",
                               right.get("length") if right.get("length") is not None else "NA"])
    atomic_bytes(str(out) + ".chromosomes.tsv", chrom_table.getvalue().encode())
    md = ["# Minigraph-Cactus: local V1 versus HPRC+local V1", "",
          "Both builds reused their frozen SV graph topology and recomputed mapping, split, alignment and join.",
          "The HPRC+local build contains additional HPRC haplotypes. Its HPRC CHM13/GRCh38 FASTAs also differ",
          "from the local build's reference FASTAs, so differences are not an isolated causal effect of HPRC addition.",
          "", "| Metric | Local V1 | HPRC+local V1 | Difference | Ratio |",
          "| --- | ---: | ---: | ---: | ---: |"]
    for name, left, right, difference, ratio in rows:
        md.append(f"| {name} | {left:,} | {right:,} | {difference:+,} | {ratio} |")
    md += ["", "## Provenance", "",
           f"- Local source GFA SHA-256: `{local_seed['source_sha256']}`",
           f"- HPRC+local source GFA SHA-256: `{hprc_seed['source_sha256']}`",
           f"- Local normalized GFA SHA-256: `{local_seed['normalized_sha256']}`",
           f"- HPRC+local normalized GFA SHA-256: `{hprc_seed['normalized_sha256']}`",
           f"- Local graph statistics: `{local_stats}`",
           f"- HPRC+local graph statistics: `{hprc_stats}`", "",
           f"Per-chromosome graph complexity is in `{out}.chromosomes.tsv`.", "",
           "File sizes and VCF sample counts are descriptive QC, not biological performance measures.",
           "Review split logs, graph statistics, reference paths and paired read-mapping accuracy before choosing a graph.", ""]
    atomic_bytes(str(out) + ".md", "\n".join(md).encode())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("local-settings", "hprc-settings", "out-prefix"):
        parser.add_argument(f"--{name}", required=True)
    args = parser.parse_args()
    render(json.loads(args.local_settings), json.loads(args.hprc_settings), args.out_prefix)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError, json.JSONDecodeError) as exc:
        sys.exit(f"Minigraph-Cactus comparison failed: {exc}")
