#!/usr/bin/env python3
"""Prepare provenance-checked SV follow-up inputs; never infer validated genotypes."""

import argparse
import csv
import gzip
import hashlib
import json
import shutil
from collections import Counter
from pathlib import Path


def open_text(path):
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def rows(path):
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_table(path, records, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    records = list(records)
    fields = fields or list(records[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def write_fasta(path, sequences):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for name, sequence in sequences.items():
            if not sequence or set(sequence) - set("ACGTN"):
                raise ValueError(f"Missing or unsupported allele sequence: {name}")
            handle.write(f">{name}\n")
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset:offset + 80] + "\n")


def read_fasta(path):
    result = {}
    name = None
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in result:
                    raise ValueError(f"Duplicate FASTA ID: {name}")
                result[name] = []
            elif line.strip():
                if name is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                result[name].append(line.strip().upper())
    return {name: "".join(parts) for name, parts in result.items()}


def correct_provenance(args):
    manifest = list(rows(args.manifest))
    if not manifest or len({r["assembly_id"] for r in manifest}) != len(manifest):
        raise ValueError("Build manifest is empty or has duplicate assemblies")
    # Reproduce hprc_cohort_samples() and its per-sample awk selection exactly.
    samples = list(dict.fromkeys(r["sample"] for r in manifest))
    ordered = [r for sample in samples for r in manifest if r["sample"] == sample]
    build_order = [dict(source_rank=args.first_rank + i, **r) for i, r in enumerate(ordered)]
    by_rank = {r["source_rank"]: r for r in build_order}
    tables = {Path(path).name: list(rows(path)) for path in args.tables}
    candidates = {}
    for table in tables.values():
        for row in table:
            key = row["candidate_id"]
            if key in candidates and any(candidates[key][field] != row[field] for field in
                                        ("inserted_node_id", "source_rank", "source_name")):
                raise ValueError(f"Conflicting candidate provenance: {key}")
            candidates[key] = row
    if not candidates:
        raise ValueError("No HPRC candidates supplied")
    wanted = {r["inserted_node_id"] for r in candidates.values()}
    nodes = {}
    with open_text(args.graph) as handle:
        for line in handle:
            if not line.startswith("S\t"):
                continue
            fields = line.rstrip().split("\t")
            if fields[1] not in wanted:
                continue
            tags = {part.split(":", 2)[0]: part.split(":", 2)[2] for part in fields[3:]}
            if fields[1] in nodes:
                raise ValueError(f"Duplicate graph node: {fields[1]}")
            nodes[fields[1]] = (fields[2].upper(), tags)
    corrected = {}
    for key, row in candidates.items():
        sequence, tags = nodes[row["inserted_node_id"]]
        rank = int(tags["SR"])
        donor = by_rank.get(rank)
        if donor is None:
            raise ValueError(f"{key}: SR {rank} is outside the declared local build order")
        prefix = donor["assembly_id"].removesuffix(".clean") + "."
        if not tags["SN"].startswith(prefix):
            raise ValueError(f"{key}: grouped build order disagrees with source contig {tags['SN']}; "
                             "check first_local_source_rank and the actual build manifest")
        if (str(rank) != row["source_rank"] or tags["SN"] != row["source_name"]
                or tags["SO"] != row["source_offset_0based"]):
            raise ValueError(f"{key}: saved candidate tags disagree with the supplied graph")
        if len(sequence) != int(row["insertion_size_bp"]) or set(sequence) - set("ACGTN"):
            raise ValueError(f"{key}: inserted node sequence/length is invalid")
        corrected[key] = dict(row, previous_introducing_sample=row["introducing_sample"],
                              previous_introducing_assembly=row["introducing_assembly"],
                              introducing_sample=donor["sample"],
                              introducing_assembly=donor["assembly_id"], source_assembly=donor["path"],
                              provenance_method="sample_grouped_build_order_and_source_contig",
                              provenance_status="PASS", local_provenance="true")

    # Bubble sequences are in reference traversal orientation; a GFA node can
    # be reverse-complemented relative to that traversal.
    by_node = {r["inserted_node_id"]: key for key, r in corrected.items()}
    alleles = {}
    with open_text(args.bubbles) as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 14:
                raise ValueError("Expected 14-column gfatools bubble output")
            segment_ids = fields[11].split(",")
            if len(segment_ids) != 3 or segment_ids[1] not in by_node:
                continue
            key = by_node[segment_ids[1]]
            row = corrected[key]
            if fields[0] != row["reference_path"] or int(fields[1]) != int(row["position_0based"]):
                continue
            if (int(fields[2]) != int(fields[1]) or list(map(int, fields[3:7])) != [3, 2, 0, 0]
                    or int(fields[7]) != int(row["insertion_size_bp"])):
                raise ValueError(f"{key}: candidate is no longer a simple insertion bubble")
            sequence = fields[13].upper()
            node_sequence = nodes[segment_ids[1]][0]
            if sequence not in (node_sequence, reverse_complement(node_sequence)):
                raise ValueError(f"{key}: bubble allele differs from its GFA node")
            if key in alleles:
                raise ValueError(f"{key}: multiple matching bubbles")
            alleles[key] = sequence
            row["sequence_sha256"] = hashlib.sha256(sequence.encode()).hexdigest()
            row["sequence_orientation"] = "reference_bubble_traversal"
    if set(alleles) != set(corrected):
        raise ValueError(f"Missing sequence-resolved bubbles: {sorted(set(corrected) - set(alleles))}")
    out = Path(args.output)
    for name, table in tables.items():
        write_table(out / name, [corrected[r["candidate_id"]] for r in table],
                    list(next(iter(corrected.values()))))
    write_table(out / "build_order.tsv", build_order)
    write_table(out / "all_candidates.tsv", [corrected[key] for key in sorted(corrected)])
    write_fasta(out / "insertion_alleles.fa", dict(sorted(alleles.items())))


def fetch_windows(path, windows):
    """Stream plain/gzip reference once; retain only requested small windows."""
    chunks = {key: [] for key in windows}
    by_chrom = {}
    for key, (chrom, start, end) in windows.items():
        if start < 0 or end <= start:
            raise ValueError(f"Invalid reference window for {key}")
        by_chrom.setdefault(chrom, []).append((key, start, end))
    chrom = None
    offset = 0
    seen = set()
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                chrom = line[1:].split()[0]
                if chrom in seen:
                    raise ValueError(f"Duplicate reference contig: {chrom}")
                seen.add(chrom)
                offset = 0
                continue
            sequence = line.strip().upper()
            end_offset = offset + len(sequence)
            for key, start, end in by_chrom.get(chrom, []):
                if start < end_offset and end > offset:
                    chunks[key].append(sequence[max(0, start - offset):min(len(sequence), end - offset)])
            offset = end_offset
    result = {key: "".join(parts) for key, parts in chunks.items()}
    for key, sequence in result.items():
        if len(sequence) != windows[key][2] - windows[key][1]:
            raise ValueError(f"Reference window missing/truncated for {key}; check build and contig names")
    return result


def bed_context(paths, candidates, padding):
    labels = {r["candidate_id"]: set() for r in candidates}
    if not paths:
        return {key: "NOT_ASSESSED" for key in labels}
    seen_chroms = set()
    for path in paths:
        with open_text(path) as handle:
            for line in handle:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                fields = line.rstrip().split("\t")
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                if start < 0 or end <= start:
                    raise ValueError(f"Invalid BED interval in {path}: {line[:100]}")
                seen_chroms.add(chrom)
                for row in candidates:
                    pos = int(row["position_0based"])
                    if chrom == row["chromosome"] and start < pos + padding + 1 and end > max(0, pos - padding):
                        labels[row["candidate_id"]].add(Path(path).name)
    return {r["candidate_id"]: (";".join(sorted(labels[r["candidate_id"]])) or
                               ("NO_OVERLAP" if r["chromosome"] in seen_chroms else "CONTIG_NOT_IN_BED"))
            for r in candidates}


def prepare_panel(args):
    if (not 2 <= args.max_loci <= 50 or args.flank < 100 or
            min(args.genic_loci, args.repeat_loci) < 0 or
            (not args.selection and args.genic_loci + args.repeat_loci >= args.max_loci)):
        raise ValueError("Use 2..50 loci, flank >=100, and quotas leaving space for TMEM")
    candidates = list(rows(args.candidates))
    alleles = read_fasta(args.alleles)
    by_id = {r["candidate_id"]: r for r in candidates}
    repeat = bed_context(args.repeat_bed, candidates, args.flank)
    segdup = bed_context(args.segdup_bed, candidates, args.flank)
    selected = {}
    if args.selection:
        for row in rows(args.selection):
            key, role = row["candidate_id"], row["panel_role"]
            if key not in by_id or key in selected or role not in {
                    "candidate", "genic_candidate", "repeat_review_control"}:
                raise ValueError(f"Invalid or duplicate panel selection: {row}")
            selected[key] = role
    else:
        eligible = [r for r in candidates if r["presentation_tier"] in {"A", "B"}]
        ranked = sorted(eligible, key=lambda r: (r["presentation_tier"], -int(r["insertion_size_bp"]), r["candidate_id"]))
        # Repeat-review controls are difficult candidates, not truth-negative alleles.
        difficult = sorted(eligible, key=lambda r: (
            repeat[r["candidate_id"]] in {"NOT_ASSESSED", "NO_OVERLAP", "CONTIG_NOT_IN_BED"},
            float(r["kmer5_diversity"]), r["candidate_id"]))
        for row in difficult[:args.repeat_loci]:
            selected[row["candidate_id"]] = "repeat_review_control"
        genic = [r for r in ranked if r["inside_annotated_gene"] == "true" and r["candidate_id"] not in selected]
        for row in genic[:args.genic_loci]:
            selected[row["candidate_id"]] = "genic_candidate"
        for row in ranked:
            if len(selected) >= args.max_loci - 1:
                break
            selected.setdefault(row["candidate_id"], "candidate")
    if not selected or len(selected) + 1 > args.max_loci:
        raise ValueError("Panel selection must contain 1..max_loci-1 HPRC loci (TMEM uses one slot)")
    chosen = [by_id[key] for key in selected]
    windows = {r["candidate_id"]: (r["chromosome"], int(r["position_0based"]) - args.flank,
                                   int(r["position_0based"]) + args.flank) for r in chosen}
    reference = fetch_windows(args.reference, windows)
    out = Path(args.output)
    sequences = {}
    queries = []
    panel = []
    for row in chosen:
        key = row["candidate_id"]
        sequence = alleles[key]
        if hashlib.sha256(sequence.encode()).hexdigest() != row["sequence_sha256"]:
            raise ValueError(f"Allele checksum mismatch: {key}")
        flank_sequence = reference[key]
        for allele_type, allele in (("REF", ""), ("ALT", sequence)):
            name = key + "|" + allele_type
            sequences[name] = flank_sequence[:args.flank] + allele + flank_sequence[args.flank:]
            queries.append(dict(query_name=name, candidate_id=key, allele_type=allele_type,
                                reference_build="CHM13", chromosome=row["chromosome"],
                                position_0=row["position_0based"], end_0=row["position_0based"],
                                allele_start=args.flank, allele_end=args.flank + len(allele),
                                query_length=len(sequences[name])))
        panel.append(dict(candidate_id=key, panel_role=selected[key], reference_build="CHM13",
                          chromosome=row["chromosome"], position_0=row["position_0based"],
                          end_0=row["position_0based"], svtype="INS", svlen=len(sequence),
                          donor_sample=row["introducing_sample"], donor_assembly=row["introducing_assembly"],
                          source_assembly=row["source_assembly"], source_contig=row["source_name"],
                          source_offset_0=row["source_offset_0based"], sequence_sha256=row["sequence_sha256"],
                          gene=row["nearest_gene"], repeat_context=repeat[key], segdup_context=segdup[key],
                          validation_status="PENDING_READ_AND_SOURCE_ASSEMBLY_VALIDATION"))
    write_fasta(out / "CHM13.queries.fa", sequences)
    tmem_sequences = read_fasta(args.tmem_fasta)
    tmem_rows = list(rows(args.tmem_targets))
    if not tmem_rows:
        raise ValueError("TMEM target manifest is empty")
    if {r["query_name"] for r in tmem_rows} != set(tmem_sequences):
        raise ValueError("TMEM query FASTA and target manifest disagree")
    ref_rows = [r for r in tmem_rows if r["allele_type"] == "reference"]
    if len(ref_rows) != 1:
        raise ValueError("Expected exactly one reference TMEM query")
    ref_length = int(ref_rows[0]["allele_length"])
    for row in tmem_rows:
        pos = int(row["whole_genome_pos_1based"]) - 1
        queries.append(dict(query_name=row["query_name"], candidate_id="TMEM106B",
                            allele_type=row["allele_index"], reference_build="hg38",
                            chromosome=row["whole_genome_chrom"], position_0=pos, end_0=pos + ref_length,
                            allele_start=row["allele_start"], allele_end=row["allele_end"],
                            query_length=len(tmem_sequences[row["query_name"]])))
    tmem = tmem_rows[0]
    panel.append(dict(candidate_id="TMEM106B", panel_role="assay_control", reference_build="hg38",
                      chromosome=tmem["whole_genome_chrom"], position_0=int(tmem["whole_genome_pos_1based"]) - 1,
                      end_0=int(tmem["whole_genome_pos_1based"]) - 1 + ref_length,
                      svtype="MULTIALLELIC_INS", svlen="see_query_manifest", donor_sample="",
                      donor_assembly="", source_assembly="", source_contig="", source_offset_0="",
                      sequence_sha256="see_query_FASTA", gene="TMEM106B", repeat_context="NOT_ASSESSED",
                      segdup_context="NOT_ASSESSED", validation_status="ASSAY_CONTROL_NOT_TRUTH_GENOTYPES"))
    write_fasta(out / "hg38.queries.fa", tmem_sequences)
    write_table(out / "panel.tsv", panel)
    write_table(out / "query_manifest.tsv", queries)
    worksheet = [dict(candidate_id=r["candidate_id"], donor_sample=r["donor_sample"],
                      source_assembly=r["source_assembly"], sample_to_validate="", read_path="",
                      independent_of_source="unknown", source_assembly_checked="not_checked",
                      left_junction_support="", right_junction_support="", reference_support="",
                      genome_wide_placement_review="pending", repeat_review="pending",
                      callability="unknown", decision="PENDING", notes="") for r in panel]
    write_table(out / "validation_worksheet.tsv", worksheet)


def freeze_tmem(args):
    evidence = list(rows(args.evidence))
    inventory_rows = list(rows(args.inventory))
    if Counter(r["target"] for r in inventory_rows) != Counter(args.graphs):
        raise ValueError("Expected exactly one inventory record per configured graph")
    inventory = {r["target"]: r for r in inventory_rows}
    expected = {(sample, graph) for sample in args.samples for graph in args.graphs}
    observed = Counter((r["sample_id"], r["target"]) for r in evidence)
    if set(observed) != expected or any(n != 1 for n in observed.values()):
        raise ValueError("Cannot freeze incomplete/duplicated TMEM evidence: expected one row per sample/graph")
    target_rows = list(rows(args.targets))
    if Counter(r["graph_id"] for r in target_rows) != Counter(args.graphs):
        raise ValueError("Expected exactly one TMEM target per configured graph")
    for row in evidence:
        if row["graph_checksum"] != inventory[row["target"]]["graph_checksum"]:
            raise ValueError("TMEM evidence and graph inventory checksums disagree")
        if row["counting_unit"] != "read_name_collapsed_fragments":
            raise ValueError("Cannot freeze old read-count-based TMEM evidence")
        if (row["target_callable_status"] != "callable_complete_graph_specific_target"
                or row["evidence_available"] != "true"):
            raise ValueError("Cannot freeze stale or unavailable TMEM evidence")
    for row in target_rows:
        if row["callable_status"] != "callable_complete_graph_specific_target":
            raise ValueError(f"Cannot freeze uncallable target: {row['graph_id']}")
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    artifacts = {"sv_evidence.tsv": args.evidence, "graph_inventory.tsv": args.inventory,
                 "tmem_graph_targets.tsv": args.targets, "report.md": args.report}
    for item in args.artifact:
        name, path = item.split("=", 1)
        if Path(name).name != name or name in artifacts:
            raise ValueError(f"Invalid/duplicate snapshot filename: {name}")
        artifacts[name] = path
    records = []
    for name, source in artifacts.items():
        destination = out / name
        shutil.copyfile(source, destination)
        records.append(dict(artifact=name, source=str(source), sha256=sha256(destination),
                            bytes=destination.stat().st_size))
    effective = out / "effective_config.json"
    effective.write_text(json.dumps(json.loads(args.effective_config), indent=2, sort_keys=True) + "\n")
    records.append(dict(artifact=effective.name, source="effective_snakemake_config",
                        sha256=sha256(effective), bytes=effective.stat().st_size))
    controls = [dict(sample_id=r["sample_id"], graph=r["target"], prior_group=r["sample_group"],
                     evidence_class=r["call_class"], reference_fragments=r["reference_path_reads"],
                     alt_left_fragments=r["left_alt_junction_reads"], alt_right_fragments=r["right_alt_junction_reads"],
                     status="PRIOR_SCREEN_CONTROL_NOT_INDEPENDENT_TRUTH")
                for r in evidence if r["sample_group"] == "high_insertion_absence_probability"]
    write_table(out / "negative_control_evidence.tsv", controls, ["sample_id", "graph", "prior_group",
                "evidence_class", "reference_fragments", "alt_left_fragments", "alt_right_fragments", "status"])
    control_path = out / "negative_control_evidence.tsv"
    records.append(dict(artifact=control_path.name, source="derived_from_frozen_sv_evidence",
                        sha256=sha256(control_path), bytes=control_path.stat().st_size))
    write_table(out / "checksums.tsv", records)


def report(args):
    panel = list(rows(args.panel))
    queries = {r["query_name"]: r for r in rows(args.queries)}
    alignments = {key: [] for key in queries}
    for specification in args.paf:
        build, path = specification.split("=", 1)
        with open_text(path) as handle:
            for line in handle:
                f = line.rstrip().split("\t")
                if len(f) < 12 or f[0] not in queries or queries[f[0]]["reference_build"] != build:
                    raise ValueError(f"Unexpected PAF record in {path}")
                alignments[f[0]].append(f)
    placement = []
    for key, query in queries.items():
        hits = alignments[key]
        near = [f for f in hits if f[5] == query["chromosome"] and
                int(f[7]) <= int(query["end_0"]) + args.flank and
                int(f[8]) >= max(0, int(query["position_0"]) - args.flank)]
        anchored = [f for f in near if int(f[2]) <= 100 and int(f[3]) >= int(f[1]) - 100
                    and 20 <= int(f[11]) < 255]
        placement.append(dict(query_name=key, candidate_id=query["candidate_id"],
                              reference_build=query["reference_build"], reported_alignments=len(hits),
                              expected_locus_alignments=len(near), two_flank_span_alignments=len(anchored),
                              placement_status="REVIEW_ALIGNMENT" if anchored else "UNRESOLVED_PLACEMENT",
                              validation_status="NOT_VALIDATION"))
    write_table(args.placement, placement)
    roles = Counter(r["panel_role"] for r in panel)
    text = ["# Review follow-up preparation", "", f"Prepared {len(panel)} loci: {dict(roles)}.", "",
            "The HPRC donor labels were checked against sample-grouped build order and graph source contigs. "
            "Insertion sequences were matched to simple bubbles and written in reference traversal orientation.", "",
            "Coordinates are zero-based. HPRC candidates use CHM13; TMEM uses GRCh38. "
            "No cross-reference harmonization or off-reference annotation was performed.", "",
            "The complete configured TMEM sample/graph panel is copied to tmem_frozen with checksums. "
            "Graph hashes are inherited from the checked inventory; freezing does not rehash the large graphs. "
            "ALT_SUPPORTED, HET_LIKE and REF_SUPPORTED remain evidence categories, not calibrated genotypes. "
            "Prior high-absence samples are controls selected by an earlier assay, not independent truth.", "",
            "Genome-wide PAFs and placement_review.tsv screen the allele queries for placement problems. "
            "An alignment spanning the query ends does not prove both junctions or the inserted sequence: "
            "inspect its CIGAR and alternative hits. Unmapped ALT queries are unresolved, not rejected alleles.", "",
            "Repeat-review controls include low-complexity candidates selected by 5-mer diversity when BEDs "
            "are absent. NOT_ASSESSED and CONTIG_NOT_IN_BED are not evidence of repeat-free sequence. "
            "BED annotations cover reference flanks and do not annotate the novel inserted sequence.", "",
            "Copy validation_worksheet.tsv before filling it in. Check original donor assemblies, "
            "both junctions, REF evidence, callability and repeat/paralog recruitment. Record whether "
            "reads are independent of the source assembly. None of these candidates is promoted to a "
            "validated allele by this preparation workflow.", "",
            "The corrected V2 QC report retains its existing-evidence accounting scope. "
            "A nearby changed bubble does not establish allele capture, and incomplete cohort accounting "
            "does not mean the small rescue pilot failed structurally.", "",
            "Older V2 qc tables and original sv_examples donor labels are superseded by this run's outputs. "
            "Original inputs and reports are retained. Cohort genotyping and association analysis are later stages.", ""]
    Path(args.output).write_text("\n".join(text))


def main():
    csv.field_size_limit(100_000_000)
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("provenance")
    for flag in ("manifest", "graph", "bubbles", "output"):
        p.add_argument("--" + flag, required=True)
    p.add_argument("--tables", nargs="+", required=True)
    p.add_argument("--first-rank", type=int, required=True)
    p.set_defaults(func=correct_provenance)
    p = sub.add_parser("panel")
    for flag in ("candidates", "alleles", "reference", "tmem-fasta", "tmem-targets", "output"):
        p.add_argument("--" + flag, required=True)
    p.add_argument("--selection")
    p.add_argument("--repeat-bed", action="append", default=[])
    p.add_argument("--segdup-bed", action="append", default=[])
    p.add_argument("--max-loci", type=int, default=30)
    p.add_argument("--genic-loci", type=int, default=15)
    p.add_argument("--repeat-loci", type=int, default=4)
    p.add_argument("--flank", type=int, default=1000)
    p.set_defaults(func=prepare_panel)
    p = sub.add_parser("freeze")
    for flag in ("evidence", "inventory", "targets", "report", "output", "effective-config"):
        p.add_argument("--" + flag, required=True)
    p.add_argument("--samples", nargs="+", required=True)
    p.add_argument("--graphs", nargs="+", required=True)
    p.add_argument("--artifact", action="append", default=[])
    p.set_defaults(func=freeze_tmem)
    p = sub.add_parser("report")
    for flag in ("panel", "queries", "placement", "output"):
        p.add_argument("--" + flag, required=True)
    p.add_argument("--paf", action="append", required=True)
    p.add_argument("--flank", type=int, default=1000)
    p.set_defaults(func=report)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
