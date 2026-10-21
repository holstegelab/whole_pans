#!/usr/bin/env python3
"""Build TMEM allele queries and derive orientation-aware targets per whole-genome graph."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path


TARGET_FIELDS = [
    "graph_id", "target_present", "callable_status", "edge_namespace", "node_translation",
    "reference_traversal",
    "alternate_traversal", "left_junction", "right_junction",
    "reference_defining_edges", "insertion_nodes", "insertion_edges",
    "reference_query_coverage", "reference_query_identity", "reference_query_mapq",
    "alternate_query_name", "alternate_query_coverage", "alternate_query_identity",
    "alternate_query_mapq", "qualified_alternate_queries", "mapping_coverage",
    "mapping_identity", "base_level_cigar_present", "both_flanks_covered",
    "insertion_interval_covered", "derivation_source", "ambiguity", "failure_reason",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    queries = sub.add_parser("queries")
    queries.add_argument("--vcf", required=True, type=Path)
    queries.add_argument("--targets", required=True, type=Path)
    queries.add_argument("--reference", required=True, type=Path)
    queries.add_argument("--vcf-position-offset", required=True, type=int)
    queries.add_argument("--flank-bp", type=int, default=1000)
    queries.add_argument("--output-fasta", required=True, type=Path)
    queries.add_argument("--output-manifest", required=True, type=Path)

    audit = sub.add_parser("audit")
    audit.add_argument("--graph", required=True)
    audit.add_argument("--gaf", required=True, type=Path)
    audit.add_argument("--query-manifest", required=True, type=Path)
    audit.add_argument("--translation", type=Path,
                       help="Segment-to-node table exported from the mapping GBZ with vg gbwt -Z --translation")
    audit.add_argument("--min-coverage", type=float, default=0.95)
    audit.add_argument("--min-identity", type=float, default=0.90)
    audit.add_argument("--min-mapq", type=int, default=20)
    audit.add_argument("--min-flank-aligned", type=int, default=900)
    audit.add_argument("--output", required=True, type=Path)

    merge = sub.add_parser("merge")
    merge.add_argument("--input", action="append", default=[])
    merge.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(path) + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def read_single_fasta(path: Path) -> tuple[str, str]:
    name = ""
    sequence = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name:
                    raise ValueError(f"regional reference must have one sequence: {path}")
                name = line[1:].strip().split()[0]
            else:
                sequence.append(line.strip())
    if not name or not sequence:
        raise ValueError(f"empty regional reference FASTA: {path}")
    return name, "".join(sequence).upper()


def find_vcf_record(vcf: Path, target_rows: list[dict[str, str]]) -> tuple[str, int, str, list[str]]:
    keys = {(row["chrom"], row["pos"], row["variant_id"]) for row in target_rows}
    if len(keys) != 1:
        raise ValueError("target_insertions.tsv must describe one TMEM VCF locus")
    wanted = next(iter(keys))
    with open_text(vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 5 and tuple(fields[:3]) == wanted:
                return fields[0], int(fields[1]), fields[3].upper(), [value.upper() for value in fields[4].split(",")]
    raise ValueError(f"TMEM target {wanted} was not found in {vcf}")


def build_queries(args: argparse.Namespace) -> None:
    targets = read_tsv(args.targets)
    if not targets:
        raise ValueError(f"no target insertion rows in {args.targets}")
    chrom, pos, ref, alts = find_vcf_record(args.vcf, targets)
    reference_name, reference_sequence = read_single_fasta(args.reference)
    start = pos - args.vcf_position_offset - 1
    if start < args.flank_bp or start + len(ref) + args.flank_bp > len(reference_sequence):
        raise ValueError("regional reference does not contain the requested TMEM flanks")
    observed_ref = reference_sequence[start:start + len(ref)]
    if observed_ref != ref:
        raise ValueError(
            f"VCF REF does not match regional FASTA at derived offset: {ref[:40]} != {observed_ref[:40]}"
        )
    left = reference_sequence[start - args.flank_bp:start]
    right = reference_sequence[start + len(ref):start + len(ref) + args.flank_bp]
    wanted_alt_indices = sorted({int(row["alt_index"].removeprefix("ALT")) for row in targets})
    records = [("TMEM106B|REF", "reference", "REF", ref)]
    for index in wanted_alt_indices:
        if index < 1 or index > len(alts):
            raise ValueError(f"target ALT{index} is absent from the VCF record")
        records.append((f"TMEM106B|ALT{index}", "alternate", f"ALT{index}", alts[index - 1]))

    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    fasta_tmp = Path(str(args.output_fasta) + ".tmp")
    manifest_rows = []
    with fasta_tmp.open("w") as fasta:
        for query_name, allele_type, allele_index, allele in records:
            query = left + allele + right
            fasta.write(f">{query_name}\n")
            for offset in range(0, len(query), 80):
                fasta.write(query[offset:offset + 80] + "\n")
            manifest_rows.append({
                "query_name": query_name,
                "allele_type": allele_type,
                "allele_index": allele_index,
                "query_length": len(query),
                "left_flank_bp": len(left),
                "right_flank_bp": len(right),
                "allele_start": len(left),
                "allele_end": len(left) + len(allele),
                "allele_length": len(allele),
                "vcf_chrom": chrom,
                "vcf_pos": pos,
                "regional_reference_name": reference_name,
                "vcf_position_offset": args.vcf_position_offset,
                "derivation_source": f"{args.vcf};{args.targets};{args.reference}",
            })
    fasta_tmp.replace(args.output_fasta)
    write_tsv(args.output_manifest, manifest_rows, list(manifest_rows[0]))


def tags(fields: list[str]) -> dict[str, str]:
    result = {}
    for value in fields[12:]:
        parts = value.split(":", 2)
        if len(parts) == 3:
            result[parts[0]] = parts[2]
    return result


def path_tokens(value: str) -> list[str]:
    if not value.startswith((">", "<")):
        raise ValueError("Target GAF uses stable coordinates; regenerate it with minigraph --vc")
    tokens = [orientation + name for orientation, name in re.findall(r"([<>])([^<>]+)", value)]
    if "".join(tokens) != value or any(re.search(r":\d+-\d+$", token) for token in tokens):
        raise ValueError("Target GAF must contain source segment IDs, not stable intervals; use minigraph --vc")
    return tokens


def flip(token: str) -> str:
    return ("<" if token.startswith(">") else ">") + token[1:]


def canonical_traversal(path: str, strand: str) -> list[str]:
    tokens = path_tokens(path)
    return [flip(value) for value in reversed(tokens)] if strand == "-" else tokens


def parse_alignments(path: Path) -> dict[str, list[dict[str, object]]]:
    result: dict[str, list[dict[str, object]]] = defaultdict(list)
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            try:
                query_length = int(fields[1])
                query_start = int(fields[2])
                query_end = int(fields[3])
                matches = int(fields[9])
                block_length = int(fields[10])
                mapq = int(fields[11])
            except ValueError:
                continue
            if fields[5] == "*" or block_length == 0:
                continue
            record_tags = tags(fields)
            result[fields[0]].append({
                "query_length": query_length,
                "query_start": query_start,
                "query_end": query_end,
                "coverage": (query_end - query_start) / query_length if query_length else 0,
                "identity": matches / block_length if block_length else 0,
                "mapq": mapq,
                "traversal": canonical_traversal(fields[5], fields[4]),
                "cigar": record_tags.get("cg", ""),
                "primary": record_tags.get("tp", "P") != "S",
                "strand": fields[4],
                "line_number": line_number,
            })
    return result


def flank_aligned_bases(record: dict[str, object], start: int, end: int) -> int:
    """Count query bases aligned to graph sequence, excluding inserted/clipped bases."""
    cigar = str(record["cigar"])
    operations = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    if "".join(length + op for length, op in operations) != cigar:
        return 0
    if record.get("strand") == "-":
        operations.reverse()
    position = int(record["query_start"])
    aligned = 0
    for length_text, operation in operations:
        length = int(length_text)
        if operation in "M=X":
            aligned += max(0, min(end, position + length) - max(start, position))
        if operation in "M=XIS":
            position += length
    return aligned


def qualify(
    record: dict[str, object], query: dict[str, str], min_coverage: float,
    min_identity: float, min_mapq: int, min_flank_aligned: int,
) -> tuple[bool, dict[str, bool]]:
    qstart = int(record["query_start"])
    qend = int(record["query_end"])
    left_flank = int(query["left_flank_bp"])
    right_start = int(query["query_length"]) - int(query["right_flank_bp"])
    checks = {
        "query_length": int(record["query_length"]) == int(query["query_length"]),
        "coverage": float(record["coverage"]) >= min_coverage,
        "identity": float(record["identity"]) >= min_identity,
        "mapq": min_mapq <= int(record["mapq"]) < 255,
        "cigar": bool(record["cigar"]),
        "left_flank": flank_aligned_bases(record, 0, left_flank) >= min_flank_aligned,
        "right_flank": flank_aligned_bases(record, right_start, int(query["query_length"])) >= min_flank_aligned,
        "allele": qstart <= int(query["allele_start"]) and qend >= int(query["allele_end"]),
        "traversal": len(record["traversal"]) >= 2,
    }
    return all(checks.values()), checks


def common_ends(reference: list[str], alternate: list[str]) -> tuple[int, int]:
    prefix = 0
    while prefix < min(len(reference), len(alternate)) and reference[prefix] == alternate[prefix]:
        prefix += 1
    suffix = 0
    while (
        suffix < min(len(reference) - prefix, len(alternate) - prefix)
        and reference[-suffix - 1] == alternate[-suffix - 1]
    ):
        suffix += 1
    return prefix, suffix


def edge_text(edges: set[tuple[str, str]]) -> str:
    return ";".join(f"{left},{right}" for left, right in sorted(edges))


def derive_target(
    graph: str, query_rows: list[dict[str, str]], alignments: dict[str, list[dict[str, object]]],
    min_coverage: float, min_identity: float, min_mapq: int, min_flank_aligned: int,
) -> dict[str, object]:
    ref_query = next(row for row in query_rows if row["allele_type"] == "reference")

    def qualified_records(query: dict[str, str]) -> list[dict[str, object]]:
        result = []
        for record in alignments.get(query["query_name"], []):
            passed, checks = qualify(
                record, query, min_coverage, min_identity, min_mapq, min_flank_aligned
            )
            record = dict(record)
            record["checks"] = checks
            if passed:
                result.append(record)
        return sorted(
            result,
            key=lambda value: (value["mapq"], value["coverage"], value["identity"]),
            reverse=True,
        )

    refs = qualified_records(ref_query)
    best_ref = refs[0] if len(refs) == 1 else None
    qualified_alts = []
    ambiguous_alt_queries = []
    for query in query_rows:
        if query["allele_type"] != "alternate":
            continue
        records = qualified_records(query)
        if len(records) == 1:
            qualified_alts.append((query, records[0]))
        elif len(records) > 1:
            ambiguous_alt_queries.append(query["query_name"])
    qualified_alts.sort(
        key=lambda value: (value[1]["mapq"], value[1]["coverage"], value[1]["identity"]),
        reverse=True,
    )

    failure = []
    if not refs:
        failure.append("reference_query_not_uniquely_complete_high_identity_mapping")
    elif len(refs) > 1:
        failure.append("reference_query_has_multiple_complete_high_identity_mappings")
    if not qualified_alts:
        if ambiguous_alt_queries:
            failure.append("alternate_queries_have_multiple_complete_high_identity_mappings")
        else:
            failure.append("no_alternate_query_complete_high_identity_both-flank_mapping")

    left_edges: set[tuple[str, str]] = set()
    right_edges: set[tuple[str, str]] = set()
    reference_edges: set[tuple[str, str]] = set()
    insertion_edges: set[tuple[str, str]] = set()
    insertion_nodes: set[str] = set()
    valid_alts = []
    if best_ref is not None:
        reference = list(best_ref["traversal"])
        for query, record in qualified_alts:
            alternate = list(record["traversal"])
            prefix, suffix = common_ends(reference, alternate)
            if prefix < 1 or suffix < 1 or prefix + suffix >= len(alternate):
                continue
            left_edges.add((alternate[prefix - 1], alternate[prefix]))
            right_edges.add((alternate[-suffix - 1], alternate[-suffix]))
            alternate_internal = alternate[prefix:len(alternate) - suffix]
            insertion_nodes.update(
                {token[1:] for token in alternate_internal} - {token[1:] for token in reference}
            )
            insertion_edges.update(zip(alternate_internal, alternate_internal[1:]))
            ref_window = reference[prefix - 1:len(reference) - suffix + 1]
            alt_window = alternate[prefix - 1:len(alternate) - suffix + 1]
            reference_edges.update(set(zip(ref_window, ref_window[1:])) - set(zip(alt_window, alt_window[1:])))
            valid_alts.append((query, record, alternate))
    # An edge used by any ALT traversal cannot distinguish REF from ALT.
    all_alt_edges = set()
    for _, _, alternate in valid_alts:
        for left, right in zip(alternate, alternate[1:]):
            all_alt_edges.update(((left, right), (flip(right), flip(left))))
    reference_edges.difference_update(all_alt_edges)
    if best_ref is not None and qualified_alts and not valid_alts:
        failure.append("alternate_and_reference_traversals_lack_unambiguous_shared_flank_boundaries")

    all_ref_edges = set()
    if best_ref is not None:
        for left, right in zip(best_ref["traversal"], best_ref["traversal"][1:]):
            all_ref_edges.update(((left, right), (flip(right), flip(left))))
    ambiguous_boundaries = (
        len(left_edges) != 1 or len(right_edges) != 1 or not reference_edges
        or bool((left_edges | right_edges) & all_ref_edges)
    )
    if best_ref is not None and valid_alts and not ambiguous_boundaries:
        callable_status = "callable_complete_graph_specific_target"
        present = "true"
        ambiguity = "none"
    elif best_ref is not None and not qualified_alts and not ambiguous_alt_queries:
        callable_status = "not_callable_no_qualified_alternate"
        present = "unknown"
        ambiguity = "failed_query_mapping_does_not_prove_allele_absence"
    else:
        callable_status = "not_callable_target_ambiguous"
        present = "true" if best_ref is not None and valid_alts else "unknown"
        ambiguity = "nonunique_query_mapping_or_graph_traversal_boundary"
        if best_ref is not None and valid_alts and ambiguous_boundaries:
            failure.append("alternate_junction_or_reference_traversal_boundary_is_not_unique")

    best_alt_query, best_alt, _ = valid_alts[0] if valid_alts else (
        qualified_alts[0][0], qualified_alts[0][1], qualified_alts[0][1]["traversal"]
    ) if qualified_alts else ({}, None, [])
    derivation_sources = sorted({row["derivation_source"] for row in query_rows})
    return {
        "graph_id": graph,
        "target_present": present,
        "callable_status": callable_status,
        "edge_namespace": "source_gfa_segments",
        "node_translation": "",
        "reference_traversal": "" if best_ref is None else "".join(best_ref["traversal"]),
        "alternate_traversal": ";;".join("".join(value[2]) for value in valid_alts),
        "left_junction": edge_text(left_edges),
        "right_junction": edge_text(right_edges),
        "reference_defining_edges": edge_text(reference_edges),
        "insertion_nodes": ",".join(sorted(insertion_nodes)),
        "insertion_edges": edge_text(insertion_edges),
        "reference_query_coverage": "" if best_ref is None else best_ref["coverage"],
        "reference_query_identity": "" if best_ref is None else best_ref["identity"],
        "reference_query_mapq": "" if best_ref is None else best_ref["mapq"],
        "alternate_query_name": best_alt_query.get("query_name", ""),
        "alternate_query_coverage": "" if best_alt is None else best_alt["coverage"],
        "alternate_query_identity": "" if best_alt is None else best_alt["identity"],
        "alternate_query_mapq": "" if best_alt is None else best_alt["mapq"],
        "qualified_alternate_queries": ",".join(value[0]["query_name"] for value in valid_alts),
        "mapping_coverage": "" if best_alt is None else best_alt["coverage"],
        "mapping_identity": "" if best_alt is None else best_alt["identity"],
        "base_level_cigar_present": str(bool(best_ref and best_ref["cigar"] and best_alt and best_alt["cigar"])).lower(),
        "both_flanks_covered": str(bool(best_ref and valid_alts)).lower(),
        "insertion_interval_covered": str(bool(valid_alts)).lower(),
        "derivation_source": ";".join(derivation_sources) + ";minigraph_-c_-xasm_--vc_graph_specific_GAF",
        "ambiguity": ambiguity,
        "failure_reason": ";".join(failure),
    }


def translate_target(result: dict, translation_path: Path) -> dict:
    """Translate only allele-defining segments; never assume how vg chopped nodes."""
    edge_fields = ("left_junction", "right_junction", "reference_defining_edges", "insertion_edges")
    source_edges = {
        field: [tuple(edge.split(",")) for edge in result[field].split(";") if edge]
        for field in edge_fields
    }
    wanted = set(filter(None, result["insertion_nodes"].split(",")))
    for edges in source_edges.values():
        wanted.update(token[1:] for edge in edges for token in edge)
    translation = {}
    with translation_path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 3 or fields[0] != "T" or fields[1] not in wanted:
                continue
            nodes = fields[2].split(",")
            if not nodes or any(not node.isdigit() or int(node) <= 0 for node in nodes):
                raise ValueError(f"Invalid vg node IDs in translation for {fields[1]}")
            if fields[1] in translation:
                raise ValueError(f"Duplicate segment in vg translation: {fields[1]}")
            translation[fields[1]] = nodes
    missing = wanted - translation.keys()
    if missing:
        raise ValueError(
            "Mapping GBZ translation lacks target segments: " + ", ".join(sorted(missing)[:10])
            + ". Check that the GFA and existing read-mapping GBZ belong to the same graph."
        )

    def oriented_nodes(token: str) -> list[str]:
        nodes = translation[token[1:]]
        return [token[0] + node for node in (nodes if token[0] == ">" else reversed(nodes))]

    for field, edges in source_edges.items():
        result[field] = edge_text({
            (oriented_nodes(left)[-1], oriented_nodes(right)[0]) for left, right in edges
        })
    result["insertion_nodes"] = ",".join(sorted({
        node for segment in filter(None, result["insertion_nodes"].split(","))
        for node in translation[segment]
    }, key=int))
    result["edge_namespace"] = "vg_node_ids"
    result["node_translation"] = str(translation_path)
    return result


def audit(args: argparse.Namespace) -> None:
    queries = read_tsv(args.query_manifest)
    result = derive_target(
        args.graph, queries, parse_alignments(args.gaf), args.min_coverage,
        args.min_identity, args.min_mapq, args.min_flank_aligned,
    )
    if args.translation is not None:
        result = translate_target(result, args.translation)
    write_tsv(args.output, [result], TARGET_FIELDS)


def merge(args: argparse.Namespace) -> None:
    rows = [row for value in args.input for row in read_tsv(Path(value))]
    rows.sort(key=lambda row: row.get("graph_id", ""))
    write_tsv(args.output, rows, TARGET_FIELDS)


def main() -> None:
    args = parse_args()
    if args.command == "queries":
        build_queries(args)
    elif args.command == "audit":
        audit(args)
    else:
        merge(args)


if __name__ == "__main__":
    main()
