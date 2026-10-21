#!/usr/bin/env python3
"""Summarize junction and allele-edge evidence for the TMEM 18.7 kb pilot."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import Counter
from pathlib import Path


NODE_RE = re.compile(r"[<>]([0-9]+)")
COMPLEMENT = str.maketrans("ACGT", "TGCA")


FIELDS = [
    "sample_id", "sample_role", "pair_id", "expected_allele",
    "expectation_basis", "graph_id", "strict_left_junction_support",
    "strict_right_junction_support", "alt_left_edge_kmer_support",
    "alt_right_edge_kmer_support", "reference_path_support", "allele_balance",
    "evidence_call", "evidence_concordance", "ambiguous_status",
    "mapping_targets_compared", "target_agreement", "evidence_scope",
    "limitations",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    summarize = sub.add_parser("summarize")
    summarize.add_argument("--gaf", required=True, type=Path)
    summarize.add_argument("--targets", required=True, type=Path)
    summarize.add_argument("--vcf", required=True, type=Path)
    summarize.add_argument("--sample-manifest", required=True, type=Path)
    summarize.add_argument("--sample", required=True)
    summarize.add_argument("--graph", required=True)
    summarize.add_argument("--kmer-size", type=int, default=31)
    summarize.add_argument("--edge-window", type=int, default=250)
    summarize.add_argument("--min-mapq", type=int, default=5)
    summarize.add_argument("--min-junction", type=int, default=3)
    summarize.add_argument("--min-reference", type=int, default=3)
    summarize.add_argument("--output", required=True, type=Path)
    merge = sub.add_parser("merge")
    merge.add_argument("--input", action="append", default=[])
    merge.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def path_nodes(value: str) -> list[str]:
    return NODE_RE.findall(value)


def adjacent_pairs(nodes: list[str]) -> set[frozenset[str]]:
    return {frozenset((left, right)) for left, right in zip(nodes, nodes[1:]) if left != right}


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def kmers(sequence: str, k: int) -> set[str]:
    sequence = sequence.upper()
    return {
        sequence[index:index + k]
        for index in range(len(sequence) - k + 1)
        if set(sequence[index:index + k]) <= {"A", "C", "G", "T"}
    }


def target_definition(path: Path) -> dict[str, object]:
    rows = read_rows(path)
    if not rows:
        raise ValueError(f"no TMEM target insertions in {path}")
    insertion_nodes = set()
    left_pairs = set()
    right_pairs = set()
    anchor_pairs = set()
    for row in rows:
        nodes = {node for node in row.get("insertion_nodes", "").split(",") if node}
        traversal = path_nodes(row.get("alt_traversal", ""))
        left = row.get("left_anchor_node", "")
        right = row.get("right_anchor_node", "")
        insertion_nodes.update(nodes)
        anchor_pairs.add((left, right))
        for first, second in zip(traversal, traversal[1:]):
            if first == left and second in nodes or second == left and first in nodes:
                left_pairs.add(frozenset((first, second)))
            if first == right and second in nodes or second == right and first in nodes:
                right_pairs.add(frozenset((first, second)))
    return {
        "rows": rows,
        "insertion_nodes": insertion_nodes,
        "left_pairs": left_pairs,
        "right_pairs": right_pairs,
        "anchor_pairs": anchor_pairs,
    }


def read_vcf_alleles(vcf: Path, targets: list[dict[str, str]]) -> tuple[set[str], set[str]]:
    wanted = {}
    for row in targets:
        key = (row["chrom"], row["pos"], row["variant_id"])
        wanted.setdefault(key, set()).add(int(row["alt_index"].removeprefix("ALT")))
    refs, alts = set(), set()
    opener = gzip.open if str(vcf).endswith(".gz") else open
    with opener(vcf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            key = tuple(fields[:3])
            if key not in wanted:
                continue
            refs.add(fields[3].upper())
            available_alts = fields[4].split(",")
            for index in wanted[key]:
                alts.add(available_alts[index - 1].upper())
    if not refs or not alts:
        raise ValueError("TMEM target alleles were not found in the configured graph VCF")
    return refs, alts


def kmer_panels(refs: set[str], alts: set[str], k: int, edge_window: int) -> dict[str, set[str]]:
    ref_kmers = set().union(*(kmers(sequence, k) for sequence in refs))
    alt_kmers = set().union(*(kmers(sequence, k) for sequence in alts))
    left = set().union(*(kmers(sequence[:edge_window], k) for sequence in alts)) - ref_kmers
    right = set().union(*(kmers(sequence[-edge_window:], k) for sequence in alts)) - ref_kmers
    return {"left": left, "right": right, "ref": ref_kmers - alt_kmers}


def fastq_sequences(paths: list[Path]):
    for path in paths:
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as handle:
            while True:
                name = handle.readline()
                if not name:
                    break
                sequence = handle.readline().strip().upper()
                plus = handle.readline()
                quality = handle.readline()
                if not quality:
                    raise ValueError(f"truncated FASTQ near {name.strip()} in {path}")
                yield sequence


def count_kmers(paths: list[Path], panels: dict[str, set[str]], k: int) -> Counter:
    lookup = {}
    for group, values in panels.items():
        for value in values:
            lookup.setdefault(value, set()).add(group)
            lookup.setdefault(reverse_complement(value), set()).add(group)
    counts = Counter()
    for sequence in fastq_sequences(paths):
        for index in range(len(sequence) - k + 1):
            for group in lookup.get(sequence[index:index + k], ()):
                counts[group] += 1
    return counts


def concordance(expected: str, call: str) -> str:
    if expected == "blinded":
        return "blinded"
    if call == "no_call":
        return "no_call"
    alt_supported = call in {"ALT_supported", "ALT_and_reference_supported"}
    if expected == "insertion_present":
        return "concordant" if alt_supported else "discordant"
    if expected == "insertion_absent":
        return "discordant" if alt_supported else "concordant"
    return "unknown_expectation"


def summarize(args: argparse.Namespace) -> dict[str, object]:
    manifest_rows = {row["sample_id"]: row for row in read_rows(args.sample_manifest)}
    if args.sample not in manifest_rows:
        raise ValueError(f"{args.sample} absent from {args.sample_manifest}")
    sample = manifest_rows[args.sample]
    target = target_definition(args.targets)
    left = right = reference = 0
    with args.gaf.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12 or int(fields[11]) < args.min_mapq:
                continue
            nodes = path_nodes(fields[5])
            pairs = adjacent_pairs(nodes)
            node_set = set(nodes)
            left_hit = bool(pairs & target["left_pairs"])
            right_hit = bool(pairs & target["right_pairs"])
            left += left_hit
            right += right_hit
            if not node_set & target["insertion_nodes"]:
                for left_anchor, right_anchor in target["anchor_pairs"]:
                    if left_anchor in nodes and right_anchor in nodes and nodes.index(left_anchor) < nodes.index(right_anchor):
                        reference += 1
                        break
    refs, alts = read_vcf_alleles(args.vcf, target["rows"])
    panels = kmer_panels(refs, alts, args.kmer_size, args.edge_window)
    paths = [Path(value) for key in ("r1", "r2") for value in sample[key].split(",") if value]
    counts = count_kmers(paths, panels, args.kmer_size)
    alt_total = left + right
    allele_balance = alt_total / (alt_total + 2 * reference) if alt_total + reference else ""
    two_junction = left >= args.min_junction and right >= args.min_junction
    ref_supported = reference >= args.min_reference
    if two_junction and ref_supported:
        call, ambiguous = "ALT_and_reference_supported", "false"
    elif two_junction:
        call, ambiguous = "ALT_supported", "false"
    elif ref_supported and left == 0 and right == 0:
        call, ambiguous = "reference_supported", "false"
    else:
        call, ambiguous = "no_call", "true"
    return {
        "sample_id": args.sample,
        "sample_role": sample["role"],
        "pair_id": sample.get("pair_id", ""),
        "expected_allele": sample.get("expected_allele", ""),
        "expectation_basis": sample.get("expectation_basis", ""),
        "graph_id": args.graph,
        "strict_left_junction_support": left,
        "strict_right_junction_support": right,
        "alt_left_edge_kmer_support": counts["left"],
        "alt_right_edge_kmer_support": counts["right"],
        "reference_path_support": reference,
        "allele_balance": allele_balance,
        "evidence_call": call,
        "evidence_concordance": concordance(sample.get("expected_allele", ""), call),
        "ambiguous_status": ambiguous,
        "mapping_targets_compared": 1,
        "target_agreement": "single_supported_target",
        "evidence_scope": "junction_and_allele_edge_only",
        "limitations": "short_reads_do_not_span_18.7kb;internal_repetitive_sequence_and_AluYb8_not_used_as_carrier_evidence",
    }


def write(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(path) + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def merge(args: argparse.Namespace) -> None:
    rows = [row for path in args.input for row in read_rows(Path(path))]
    by_sample = {}
    for row in rows:
        by_sample.setdefault(row["sample_id"], []).append(row)
    for sample_rows in by_sample.values():
        calls = {row["evidence_call"] for row in sample_rows}
        agreement = "agree" if len(calls) == 1 else "disagree"
        if len(sample_rows) == 1:
            agreement = "single_supported_target"
        for row in sample_rows:
            row["mapping_targets_compared"] = len(sample_rows)
            row["target_agreement"] = agreement
    rows.sort(key=lambda row: (row.get("pair_id", ""), row["sample_role"], row["sample_id"], row["graph_id"]))
    write(args.output, rows)


def main() -> None:
    args = parse_args()
    if args.command == "summarize":
        write(args.output, [summarize(args)])
    else:
        merge(args)


if __name__ == "__main__":
    main()
