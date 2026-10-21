"""Build-aware annotation joins for the bounded, technically passing panel.

Interval searches belong to bedtools rules. This script retains exact insertion
boundaries, transcript/evidence metadata and explicit missing-resource states.
"""

import csv
from collections import defaultdict
from pathlib import Path


FEATURE_FIELDS = ["chromosome", "start_0", "end_0", "feature_id", "feature_type",
                  "gene_id", "gene_name", "transcript_id", "strand", "evidence",
                  "cell_type", "build", "release"]
FEATURE_TYPES = {"gene", "tss", "exon", "CDS", "splice", "UTR", "intron", "repeat", "segdup", "cCRE"}


def read(path):
    if not path:
        return []
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, records, fields):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    records = list(records)
    fields = list(dict.fromkeys(fields + [key for row in records for key in row]))
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def joined(values):
    return ";".join(sorted({str(value) for value in values if value and value != "."}))


def true(value):
    return str(value).lower() == "true"


def prepare(sm):
    cfg = sm.params.settings
    build, release = cfg["build"], cfg["release"]
    features = read(sm.input.features[0]) if sm.input.features else []
    if features and (not build or not release):
        raise ValueError("Freeze annotation build and release before annotating")
    seen = set()
    for row in features:
        if set(FEATURE_FIELDS) - row.keys() or row["feature_id"] in seen:
            raise ValueError("Feature schema incomplete or feature_id duplicated")
        seen.add(row["feature_id"])
        start, end = int(row["start_0"]), int(row["end_0"])
        if start < 0 or end <= start or row["feature_type"] not in FEATURE_TYPES:
            raise ValueError("Invalid 0-based half-open feature interval/type")
        if row["build"] != build or not row["release"].strip("."):
            raise ValueError("Feature build differs from configuration or source release is missing")
        if row["feature_type"] in {"gene", "tss", "exon", "CDS", "splice", "UTR", "intron"}:
            if not row["gene_id"] or row["gene_id"] == ".":
                raise ValueError("Gene/transcript feature lacks gene_id")
        if row["feature_type"] in {"tss", "exon", "CDS", "splice", "UTR", "intron"}:
            if row["transcript_id"] in {"", "."}:
                raise ValueError("Transcript-specific feature lacks transcript_id")
        if row["feature_type"] == "tss" and (end-start != 1 or row["strand"] not in {"+", "-"}):
            raise ValueError("Each TSS must be one base with a known transcript strand")
        if row["feature_type"] == "cCRE" and not row["evidence"].strip("."):
            raise ValueError("cCRE rows require explicit downloaded source/release evidence")
    lock = {}
    for row in read(sm.input.lock):
        if row["graph_id"] != sm.params.main_graph:
            continue
        key = row["allele_id"]
        if key in lock:
            raise ValueError("Measurement lock must have one row per allele")
        lock[key] = row
    reviews = {}
    for row in (read(sm.input.review[0]) if sm.input.review else []):
        key = row["allele_id"], row["build"]
        if key in reviews:
            raise ValueError("Duplicate annotation-coordinate review")
        reviews[key] = row
    panel = read(sm.input.panel)
    passing = [row for row in panel if lock.get(row["allele_id"], {}).get("status") == "LOCKED"
               and lock.get(row["allele_id"], {}).get("technical_status") == "PASS"]
    if not passing:
        raise ValueError("Annotation requires at least one calibrated LOCKED allele for the main graph")
    sequences, name = {}, ""
    with open(sm.input.fasta) as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in sequences:
                    raise ValueError("Duplicate catalogue FASTA identifier")
                sequences[name] = []
            elif line.strip():
                sequences[name].append(line.strip())
    loci, queries = [], []
    feature_contigs = {r["chromosome"] for r in features}
    flank = int(cfg.get("flank_bp", 1000))
    if flank < 1:
        raise ValueError("Annotation flank_bp must be positive")
    with open(sm.output.inserted, "w") as inserted:
        for row in passing:
            key = row["allele_id"]
            if row["status"] not in {"PASS", "PASS_CONTROL"}:
                raise ValueError("Technically passing allele has unresolved source placement")
            if lock[key].get("catalogue_graph_sha256") != row["graph_sha256"]:
                raise ValueError("Lock catalogue_graph_sha256 must identify the panel's source GFA")
            if lock[key].get("locked_before_phenotypes") != "yes":
                raise ValueError("Annotation requires phenotype-blind calibrated measurement rules")
            for field in ("calibration_evidence", "threshold_sensitivity_evidence", "depth_batch_review", "technical_review"):
                if not lock[key].get(field) or not Path(lock[key][field]).is_file():
                    raise ValueError(f"Locked annotation allele lacks existing {field}")
            same_build = row["reference_build"] == build
            coordinate = reviews.get((key, build), {})
            reviewed = coordinate and all(true(coordinate.get(field)) for field in (
                "left_flank_verified", "right_flank_verified", "left_flank_unique", "right_flank_unique",
                "orientation_verified", "synteny_verified"))
            reviewed = reviewed and Path(coordinate.get("evidence", "")).is_file()
            status = "READY" if same_build or reviewed else "UNRESOLVED_COORDINATES"
            if same_build and not coordinate:
                chrom, start, end = row["chromosome"], int(row["reference_start_0"]), int(row["reference_end_0"])
                evidence = "CATALOGUE_RESOLVED_NATIVE_BUILD"
            elif reviewed:
                chrom, start, end = coordinate["chromosome"], int(coordinate["start_0"]), int(coordinate["end_0"])
                evidence = coordinate["evidence"]
            else:
                chrom, start, end, evidence = "", "", "", "TWO_SIDED_CORRESPONDENCE_NOT_ESTABLISHED"
                status = "UNRESOLVED_COORDINATES"
            if status == "READY":
                if start < 0 or end < start or (row["sv_type"] == "INS" and end != start):
                    raise ValueError("Invalid affected-reference interval; an insertion is a zero-length boundary")
                if not features:
                    status = "NOT_ASSESSED"
                elif chrom not in feature_contigs:
                    status = "ANNOTATION_CONTIG_UNAVAILABLE"
                else:
                    # This interval is only a search proxy. Final context
                    # and distances below use the original zero-length boundary.
                    queries.append((chrom, max(0, start-1) if start == end else start,
                                    max(start+1, end), key+"|locus"))
                    # Searching each outer edge also preserves exact-distance
                    # ties at abutting features despite BED's base distances.
                    if start:
                        queries.append((chrom, start-1, start, key+"|edge_left"))
                    queries.append((chrom, end, end+1, key+"|edge_right"))
                    if start:
                        queries.append((chrom, max(0, start-flank), start, key+"|left_flank"))
                    queries.append((chrom, end, end+flank, key+"|right_flank"))
            record = dict(row, annotation_build=build, annotation_release=release,
                          annotation_chromosome=chrom, annotation_start_0=start, annotation_end_0=end,
                          status=status, reason=evidence)
            sequence = "".join(sequences.get(key+"|ALT", []))
            left, right = int(row["allele_start"]), int(row["allele_end"])
            if row["sv_type"] in {"INS", "MULTIALLELIC_INS"} and sequence and 0 <= left < right <= len(sequence):
                inserted.write(f">{key}\n{sequence[left:right]}\n")
                record["inserted_sequence_length"] = right-left
            else:
                record["inserted_sequence_length"] = ""
            loci.append(record)
    write(sm.output.loci, loci, ["locus_id", "allele_id", "graph_version", "status", "reason"])
    write(sm.output.metadata, features, FEATURE_FIELDS)
    for path, kind in ((sm.output.features, None), (sm.output.genes, "gene"), (sm.output.tss, "tss")):
        with open(path, "w") as handle:
            for row in features:
                if kind is None or row["feature_type"] == kind:
                    handle.write("\t".join(row[field] for field in FEATURE_FIELDS[:4])+"\n")
    with open(sm.output.queries, "w") as handle:
        for row in queries:
            handle.write("\t".join(map(str, row))+"\n")


def hits(path, features):
    """Join native BED4 query + BED4 feature results, retaining all ties."""
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 8 and fields[7] in features:
                allele, query_type = fields[3].rsplit("|", 1)
                yield allele, query_type, features[fields[7]]


def distance(locus, feature):
    start, end = int(locus["annotation_start_0"]), int(locus["annotation_end_0"])
    left, right = int(feature["start_0"]), int(feature["end_0"])
    if feature["feature_type"] == "tss":
        return max(start-left, left-end, 0)
    return max(start-right, left-end, 0)


def summarize(sm):
    cfg = sm.params.settings
    loci = {r["allele_id"]: r for r in read(sm.input.loci)}
    features = {r["feature_id"]: r for r in read(sm.input.metadata)}
    overlaps, nearby, nearest = defaultdict(list), defaultdict(list), defaultdict(list)
    evidence = []
    for allele, query_type, feature in hits(sm.input.overlaps, features):
        if query_type in {"edge_left", "edge_right"}:
            continue
        relation = "REFERENCE_FLANK_OVERLAP"
        if query_type == "locus":
            start, end = int(loci[allele]["annotation_start_0"]), int(loci[allele]["annotation_end_0"])
            left, right = int(feature["start_0"]), int(feature["end_0"])
            relation = ("INSERTION_INSIDE_FEATURE" if left < start < right else "INSERTION_AT_FEATURE_BOUNDARY") if start == end else "AFFECTED_REFERENCE_OVERLAP"
        overlaps[allele].append((query_type, feature, relation))
        evidence.append(dict(allele_id=allele, locus_id=loci[allele]["locus_id"],
                             graph_version=loci[allele]["graph_version"], query_type=query_type,
                             relation=relation, status="CONTEXT_ONLY", reason="NO_FUNCTIONAL_EFFECT_INFERRED", **feature))
    for allele, query_type, feature in hits(sm.input.nearby, features):
        if query_type == "locus" and distance(loci[allele], feature) <= 100000:
            nearby[allele].append(feature)
    for path in (sm.input.nearest_gene, sm.input.nearest_tss):
        for allele, query_type, feature in hits(path, features):
            if query_type in {"locus", "edge_left", "edge_right"}:
                nearest[allele, feature["feature_type"]].append(feature)
    links = defaultdict(list)
    for row in (read(sm.input.links[0]) if sm.input.links else []):
        required = {"feature_id", "target_gene_id", "target_gene_name", "evidence_type", "source", "cell_type"}
        if required - row.keys() or any(not row[field].strip(".") for field in ("target_gene_id", "evidence_type", "source")):
            raise ValueError("Regulatory links require target, evidence type and explicit source")
        if row["feature_id"] not in features or features[row["feature_id"]]["feature_type"] != "cCRE":
            raise ValueError("Regulatory link does not refer to a frozen cCRE feature")
        links[row["feature_id"]].append(row)
    repeat_hits = defaultdict(list)
    if sm.input.repeatmasker:
        with open(sm.input.repeatmasker[0]) as handle:
            for line in handle:
                fields = line.split()
                if not fields or not fields[0].isdigit():
                    continue
                if len(fields) < 11 or fields[4] not in loci:
                    raise ValueError("Unexpected RepeatMasker query/output format")
                start, end = int(fields[5])-1, int(fields[6])
                if not 0 <= start < end <= int(loci[fields[4]]["inserted_sequence_length"]):
                    raise ValueError("RepeatMasker insertion coordinates are outside the submitted sequence")
                repeat_hits[fields[4]].append(f"{fields[9]}:{fields[10]}:{start}-{end}")
    assessed = set(cfg.get("assessed_features", []))
    if assessed - FEATURE_TYPES:
        raise ValueError("Unknown assessed_features category")
    result = []
    for allele, locus in loci.items():
        row = dict(locus)
        entries = overlaps[allele]
        available = locus["status"] == "READY"
        for kind in sorted(FEATURE_TYPES - {"tss"}):
            matching = [(feature, relation) for query, feature, relation in entries if query == "locus" and feature["feature_type"] == kind]
            row[kind+"_features"] = joined(f"{f['feature_id']}:{f['gene_id']}:{f['transcript_id']}:{relation}" for f, relation in matching)
            if not row[kind+"_features"]:
                row[kind+"_features"] = "NO_OVERLAP" if available and kind in assessed else "NOT_ASSESSED"
        row["gene_context"] = joined(feature["feature_type"] for query, feature, relation in entries
                                     if query == "locus" and feature["feature_type"] in {"exon", "CDS", "splice", "UTR", "intron"})
        if not row["gene_context"]:
            row["gene_context"] = "INTERGENIC" if row["gene_features"] == "NO_OVERLAP" else "TRANSCRIPT_CONTEXT_NOT_ASSESSED"
        for kind in ("gene", "tss"):
            candidates = nearest[allele, kind]
            best = min((distance(locus, f) for f in candidates), default=None)
            chosen = [f for f in candidates if distance(locus, f) == best]
            row["nearest_"+kind] = joined(f"{f['gene_id']}:{f['gene_name']}:{f['transcript_id']}" for f in chosen) or "NOT_ASSESSED"
            row["nearest_"+kind+"_distance_bp"] = best if best is not None else ""
        for window in (2000, 10000, 100000):
            genes = {f["gene_id"] for f in nearby[allele] if distance(locus, f) <= window}
            row[f"tss_genes_within_{window}bp"] = joined(genes) if genes else ("NONE" if available and "tss" in assessed else "NOT_ASSESSED")
            row[f"tss_gene_count_{window}bp"] = len(genes) if available and "tss" in assessed else ""
        row["all_tss_within_100kb"] = joined(f"{f['gene_id']}:{f['gene_name']}:{f['transcript_id']}:{f['strand']}:{f['start_0']}:{distance(locus, f)}bp" for f in nearby[allele])
        for side in ("left_flank", "right_flank"):
            for kind in ("repeat", "segdup"):
                row[side+"_"+kind] = joined(f["feature_id"] for query, f, rel in entries
                                             if query == side and f["feature_type"] == kind) or (
                                                 "NO_OVERLAP" if available and kind in assessed else "NOT_ASSESSED")
        target_links = [link for query, feature, relation in entries if query == "locus" and feature["feature_type"] == "cCRE" for link in links[feature["feature_id"]]]
        row["regulatory_targets"] = joined(f"{r['target_gene_id']}:{r['target_gene_name']}:{r['evidence_type']}:{r['cell_type']}:{r['source']}" for r in target_links) or "NOT_ASSESSED"
        row["inserted_sequence_repeats"] = joined(repeat_hits[allele]) or (
            "NO_HIT_IN_SUPPLIED_LIBRARY" if sm.input.repeatmasker and locus["inserted_sequence_length"] else "NOT_ASSESSED")
        row["repeat_library_release"] = cfg.get("repeatmasker", {}).get("library_release", "") if sm.input.repeatmasker else ""
        row["status"] = "ANNOTATED_CONTEXT" if available else locus["status"]
        result.append(row)
    write(sm.output.context, result, ["locus_id", "allele_id", "graph_version", "annotation_build", "annotation_release", "status", "reason"])
    write(sm.output.evidence, evidence, ["locus_id", "allele_id", "graph_version", "query_type", "relation", "status", "reason"]+FEATURE_FIELDS)
    counts = []
    for window in (2000, 10000, 100000):
        assessed_rows = [r for r in result if r[f"tss_gene_count_{window}bp"] != ""]
        near = [r for r in assessed_rows if int(r[f"tss_gene_count_{window}bp"]) > 0]
        counts.append(dict(window_bp=window, allele_denominator=len(assessed_rows), allele_count=len(near),
                           locus_denominator=len({r["locus_id"] for r in assessed_rows}),
                           locus_count=len({r["locus_id"] for r in near}), build=cfg["build"], release=cfg["release"],
                           status="DESCRIPTIVE_ONLY", reason="SELECTED_PANEL_NO_ENRICHMENT_TEST"))
    write(sm.output.counts, counts, ["window_bp", "allele_denominator", "allele_count", "locus_denominator", "locus_count", "build", "release", "status", "reason"])


with open(snakemake.log[0], "w") as annotation_log:
    try:
        if snakemake.params.task == "prepare":
            prepare(snakemake)
        else:
            summarize(snakemake)
        annotation_log.write("Completed build-aware annotation table joins.\n")
    except Exception as error:
        annotation_log.write(f"{type(error).__name__}: {error}\n")
        raise
