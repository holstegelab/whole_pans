#!/usr/bin/env python3
"""Insertion-first catalogue, conservative resource review, and explicit pilot panel.

Uses the existing follow-up reader/reconstruction helpers. No unmapped query is
called absent; representation decisions require a recorded locus-level review.
"""

import argparse
import hashlib
from collections import Counter, defaultdict
from pathlib import Path

from review_followup import open_text, read_fasta, reverse_complement, rows, write_fasta, write_table


RESOURCES = ("HPRC", "GRCh38", "CHM13")
STATUSES = {"REPRESENTED", "PARTIAL_OR_DIFFERENT_ALLELE", "NOT_REPRESENTED_AT_RESOLVED_LOCUS",
            "UNRESOLVED", "NOT_ASSESSED"}


def write(path, records, fields):
    """Keep headers even when no candidates pass; retain original provenance columns."""
    records = list(records)
    names = list(dict.fromkeys(fields + [key for row in records for key in row]))
    write_table(path, [{name: row.get(name, "") for name in names} for row in records], names)


def digest(*values):
    return hashlib.sha256("\0".join(map(str, values)).encode()).hexdigest()


def yes(value):
    return str(value).lower() == "true"


def windows(path, wanted):
    """One streaming FASTA scan; missing/truncated windows remain explicit failures."""
    by_contig = defaultdict(list)
    result = {key: [] for key in wanted}
    for key, (chrom, start, end) in wanted.items():
        if start >= 0 and end > start:
            by_contig[chrom].append((key, start, end))
    chrom, offset = None, 0
    seen = set()
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                chrom, offset = line[1:].split()[0], 0
                if chrom in seen:
                    raise ValueError(f"Duplicate source/reference contig in {path}: {chrom}")
                seen.add(chrom)
                continue
            sequence = line.strip().upper()
            for key, start, end in by_contig[chrom]:
                if start < offset + len(sequence) and end > offset:
                    result[key].append(sequence[max(0, start-offset):min(len(sequence), end-offset)])
            offset += len(sequence)
    return {key: "".join(chunks) for key, chunks in result.items()
            if len("".join(chunks)) == wanted[key][2] - wanted[key][1]}


def catalogue(args):
    if (args.flank < 100 or not args.graph_version or len(args.graph_sha256) != 64
            or set(args.graph_sha256.lower()) - set("0123456789abcdef")):
        raise ValueError("Require flank >=100 and frozen graph version/SHA256")
    original = list(rows(args.candidates))
    if not original or len({r["candidate_id"] for r in original}) != len(original):
        raise ValueError("Candidate catalogue is empty or contains duplicate candidate IDs")
    insertions = read_fasta(args.alleles)
    review = {}
    for row in rows(args.source_review) if args.source_review else []:
        if row["candidate_id"] in review or row["candidate_id"] not in {r["candidate_id"] for r in original}:
            raise ValueError("Duplicate or unknown candidate in source review")
        review[row["candidate_id"]] = row
    flank, out = args.flank, Path(args.output)
    rejected, attempts, records, crosswalk, queries = [], [], [], [], []
    sequences, anchored, unique = {}, {}, {}
    usable = []
    for row in original:
        key = row["candidate_id"]
        reason = ""
        if row.get("sv_type") != "INS" or int(row.get("reference_length_bp", -1)) != 0:
            reason = "UNSUPPORTED_NON_INSERTION_REQUIRES_EXACT_CONFIGURATION"
        elif row.get("provenance_status") != "PASS":
            reason = "CORRECTED_PROVENANCE_NOT_PASS"
        elif key not in insertions or digest(insertions[key]) != row.get("sequence_sha256"):
            reason = "MISSING_INSERTION_OR_CHECKSUM_DISAGREEMENT"
        if reason:
            rejected.append(dict(row, status="REJECTED", reason=reason))
            attempts.append(dict(candidate_id=key, status="REJECTED", reason=reason))
        else:
            usable.append(row)
    reference = windows(args.reference, {r["candidate_id"]: (
        r["chromosome"], int(r["position_0based"])-flank, int(r["position_0based"])+flank)
        for r in usable})
    by_source = defaultdict(list)
    for row in usable:
        by_source[row["source_assembly"]].append(row)
    source = {}
    for path, group in by_source.items():
        if Path(path).is_file():
            source.update(windows(path, {r["candidate_id"]: (r["source_name"],
                int(r["source_offset_0based"])-flank,
                int(r["source_offset_0based"])+len(insertions[r["candidate_id"]])+flank) for r in group}))
    for row in usable:
        key = row["candidate_id"]
        inserted, ref = insertions[key], reference.get(key, "")
        observed = source.get(key, "")
        reason, orientation = "", ""
        if not ref or set(ref) - set("ACGTN"):
            reason = "REFERENCE_FLANKS_MISSING_TRUNCATED_OR_INVALID"
        elif not observed or set(observed) - set("ACGTN"):
            reason = "SOURCE_FASTA_CONTIG_OR_COMPLETE_FLANKS_UNAVAILABLE"
        else:
            centre = observed[flank:-flank]
            if centre == inserted and reverse_complement(centre) != inserted:
                orientation = "+"
            elif reverse_complement(centre) == inserted and centre != inserted:
                orientation, observed = "-", reverse_complement(observed)
            elif centre == inserted:
                decision = review.get(key, {})
                orientation = decision.get("source_orientation", "")
                if orientation == "-":
                    observed = reverse_complement(observed)
                elif orientation != "+":
                    reason = "PALINDROMIC_INSERTION_ORIENTATION_REQUIRES_FLANK_REVIEW"
            else:
                reason = "SOURCE_OFFSET_INSERTION_SEQUENCE_DISAGREEMENT"
        if reason:
            rejected.append(dict(row, status="REJECTED", reason=reason))
            attempts.append(dict(candidate_id=key, status="REJECTED", reason=reason))
            continue
        if review.get(key, {}).get("source_orientation") not in {None, "", orientation}:
            raise ValueError(f"Source review contradicts exact source orientation: {key}")
        exact_anchors = observed[:flank] == ref[:flank] and observed[-flank:] == ref[-flank:]
        decision = review.get(key, {})
        reviewed = all(yes(decision.get(field)) for field in
                       ("source_locus_resolved", "left_anchor_verified", "right_anchor_verified",
                        "orientation_verified")) and Path(decision.get("evidence", "")).is_file()
        status = "PASS" if exact_anchors or reviewed else "SOURCE_INSERT_CONFIRMED_ANCHORS_PENDING"
        reason = "EXACT_SOURCE_AND_REFERENCE_ANCHORS" if exact_anchors else (
            decision.get("reason", "RECORDED_SOURCE_ANCHOR_REVIEW") if reviewed else "REVIEW_SOURCE_FLANK_PLACEMENT")
        locus = "L" + digest(args.graph_version, args.graph_sha256, row["reference_build"],
                             row["chromosome"], row["position_0based"], ref)[:20]
        allele = locus + ".A" + digest(inserted)[:16]
        record = dict(row, locus_id=locus, allele_id=allele, graph_version=args.graph_version,
                      graph_sha256=args.graph_sha256, status=status, reason=reason,
                      source_orientation=orientation, source_sequence_status="EXACT_INSERTION_VERIFIED",
                      source_flank_review_evidence=decision.get("evidence", ""),
                      reference_start_0=row["position_0based"], reference_end_0=row["position_0based"],
                      inserted_source_start_0=row["source_offset_0based"],
                      inserted_source_end_0=int(row["source_offset_0based"])+len(inserted),
                      allele_start=flank, allele_end=flank+len(inserted),
                      left_flank_bp=flank, right_flank_bp=flank,
                      repeat_status=row.get("excluded_region") or "NOT_ASSESSED",
                      segdup_status="NOT_ASSESSED", placement_status="UNIQUE_LOCUS_PLACEMENT_NOT_YET_REVIEWED")
        attempts.append(dict(candidate_id=key, locus_id=locus, allele_id=allele, status=status, reason=reason))
        crosswalk.append(dict(candidate_id=key, locus_id=locus, allele_id=allele,
                              reference_build=row["reference_build"], chromosome=row["chromosome"],
                              start_0=row["position_0based"], end_0=row["position_0based"],
                              vcf_anchor_pos_1based=row["position_0based"],
                              source_assembly=row["source_assembly"], source_contig=row["source_name"],
                              source_start_0=row["source_offset_0based"],
                              source_end_0=int(row["source_offset_0based"])+len(inserted),
                              source_orientation=orientation, graph_nodes=row.get("segment_ids", ""),
                              graph_version=args.graph_version, status=status, reason=reason))
        if allele in unique:
            unique[allele]["candidate_id"] += ";" + key
            continue
        unique[allele] = record
        records.append(record)
        for kind, sequence in (("REF", ref), ("ALT", observed),
                               ("LEFT", observed[:flank]), ("RIGHT", observed[-flank:])):
            name = allele + "|" + kind
            sequences[name] = sequence
            if kind in {"REF", "ALT"}:
                anchored[name] = sequence
            queries.append(dict(query_name=name, candidate_id=key, locus_id=locus, allele_id=allele,
                                query_type=kind, query_length=len(sequence),
                                allele_start=flank if kind in {"REF", "ALT"} else "",
                                allele_end=flank+len(inserted) if kind == "ALT" else (
                                    flank if kind == "REF" else ""),
                                reference_build=row["reference_build"], chromosome=row["chromosome"],
                                start_0=row["position_0based"], end_0=row["position_0based"],
                                graph_version=args.graph_version, status=status, reason=reason))
    if args.tmem_fasta or args.tmem_targets:
        if not (args.tmem_fasta and args.tmem_targets):
            raise ValueError("Supply both TMEM FASTA and target table")
        add_tmem(args, records, queries, sequences, anchored, crosswalk)
    if args.controls or args.control_fasta:
        if not (args.controls and args.control_fasta):
            raise ValueError("Supply both known-control table and query FASTA")
        add_controls(args, records, queries, sequences, anchored, crosswalk)
    # Native coordinates are not a liftover. Two independently placed flanks
    # and synteny must establish the corresponding interval in the other build.
    for native in list(crosswalk):
        native["correspondence_status"] = "NATIVE_COORDINATE"
        counterpart = "CHM13" if native["reference_build"] == "GRCh38" else "GRCh38"
        crosswalk.append(dict(candidate_id=native["candidate_id"], locus_id=native["locus_id"],
                              allele_id=native["allele_id"], reference_build=counterpart,
                              graph_version=args.graph_version, correspondence_status="UNRESOLVED",
                              status="UNRESOLVED", reason="REQUIRES_INDEPENDENT_TWO_FLANK_SYNTENY_REVIEW"))
    write(out / "alleles.tsv", records, ["candidate_id", "locus_id", "allele_id", "status", "reason"])
    write(out / "query_manifest.tsv", queries, ["query_name", "locus_id", "allele_id", "query_type", "query_length"])
    write(out / "coordinate_crosswalk.tsv", crosswalk, ["candidate_id", "locus_id", "allele_id", "status", "reason"])
    write(out / "rejections.tsv", rejected, ["candidate_id", "status", "reason"])
    write(out / "candidate_flow.tsv", attempts, ["candidate_id", "locus_id", "allele_id", "status", "reason"])
    descriptions = Counter()
    outcomes = {row["candidate_id"]: row["status"] for row in attempts}
    for row in original:
        size = int(row.get("insertion_size_bp", 0))
        size_bin = "<1kb" if size < 1000 else ("1-10kb" if size < 10000 else ">=10kb")
        descriptions[row.get("sv_type", "UNKNOWN"), size_bin,
                     row.get("presentation_tier", "NOT_ASSESSED"), outcomes[row["candidate_id"]]] += 1
    write(out / "descriptive_totals.tsv", [dict(sv_type=kind, size_bin=size_bin,
          confidence_tier=tier, status=status, unit="candidate_records", count=count)
          for (kind, size_bin, tier, status), count in descriptions.items()],
          ["sv_type", "size_bin", "confidence_tier", "status", "unit", "count"])
    write_fasta(out / "alleles.fa", anchored)
    write_fasta(out / "locus_queries.fa", sequences)


def add_tmem(args, records, queries, sequences, anchored, crosswalk):
    """Keep preexisting TMEM query boundaries/build and every ALT separately."""
    fasta, targets = read_fasta(args.tmem_fasta), list(rows(args.tmem_targets))
    if set(fasta) != {r["query_name"] for r in targets}:
        raise ValueError("TMEM target manifest and FASTA IDs disagree")
    refs = [r for r in targets if r["allele_type"] == "reference"]
    if len(refs) != 1:
        raise ValueError("TMEM requires one exact reference query")
    ref = refs[0]
    locus = "L" + digest(args.graph_version, args.graph_sha256, "GRCh38", ref["whole_genome_chrom"],
                         ref["whole_genome_pos_1based"], fasta[ref["query_name"]])[:20]
    for row in targets:
        if row is ref:
            continue
        full = fasta[row["query_name"]]
        start, end = int(row["allele_start"]), int(row["allele_end"])
        if not 0 < start <= end < len(full):
            raise ValueError("Invalid TMEM allele/flank boundaries")
        allele = locus + ".A" + digest(full)[:16]
        position = int(row["whole_genome_pos_1based"])-1
        records.append(dict(candidate_id="TMEM106B", locus_id=locus, allele_id=allele,
                            reference_build="GRCh38", chromosome=row["whole_genome_chrom"],
                            reference_start_0=position, reference_end_0=position+int(ref["allele_length"]),
                            sv_type="MULTIALLELIC_INS", insertion_size_bp=end-start,
                            graph_version=args.graph_version, graph_sha256=args.graph_sha256,
                            nearest_gene="TMEM106B", sequence_sha256=digest(full),
                            allele_start=start, allele_end=end, left_flank_bp=start,
                            right_flank_bp=len(full)-end, source_sequence_status="PREEXISTING_CONTROL",
                            status="PASS_CONTROL", reason="FROZEN_PREEXISTING_QUERY_NOT_TRUTH_GENOTYPE"))
        crosswalk.append(dict(candidate_id="TMEM106B", locus_id=locus, allele_id=allele,
                              reference_build="GRCh38", chromosome=row["whole_genome_chrom"],
                              start_0=position, end_0=position+int(ref["allele_length"]),
                              graph_version=args.graph_version, status="PASS_CONTROL",
                              reason="PREEXISTING_CONTROL_COORDINATES"))
        for kind, sequence, astart, aend in (
                ("REF", fasta[ref["query_name"]], int(ref["allele_start"]), int(ref["allele_end"])),
                ("ALT", full, start, end), ("LEFT", full[:start], "", ""),
                ("RIGHT", full[end:], "", "")):
            name = allele + "|" + kind
            sequences[name] = sequence
            if kind in {"REF", "ALT"}:
                anchored[name] = sequence
            queries.append(dict(query_name=name, candidate_id="TMEM106B", locus_id=locus, allele_id=allele,
                                query_type=kind, query_length=len(sequence), allele_start=astart,
                                allele_end=aend, reference_build="GRCh38", chromosome=row["whole_genome_chrom"],
                                start_0=position, end_0=position+int(ref["allele_length"]),
                                graph_version=args.graph_version, status="PASS_CONTROL",
                                reason="PREEXISTING_CONTROL_COORDINATES"))


def add_controls(args, records, queries, sequences, anchored, crosswalk):
    """Import sequence-resolved controls outside the local-addition catalogue."""
    fasta = read_fasta(args.control_fasta)
    existing = {r["allele_id"] for r in records}
    seen = set()
    for row in rows(args.controls):
        key = row["control_id"]
        if key in seen:
            raise ValueError(f"Duplicate control ID: {key}")
        seen.add(key)
        if not Path(row["evidence"]).is_file():
            raise ValueError(f"Missing known-control sequence evidence: {key}")
        ref, alt = (fasta[key + "|" + kind] for kind in ("REF", "ALT"))
        left, right = int(row["allele_start"]), int(row["allele_end"])
        ref_left, ref_right = int(row["ref_allele_start"]), int(row["ref_allele_end"])
        start, end = int(row["start_0"]), int(row["end_0"])
        if (row["reference_build"] not in {"CHM13", "GRCh38"} or not 0 <= start <= end
                or not 0 < left <= right < len(alt) or not 0 < ref_left <= ref_right < len(ref)
                or ref_right-ref_left != end-start
                or alt[:left] != ref[:ref_left] or alt[right:] != ref[ref_right:]):
            raise ValueError(f"Inconsistent control coordinates/flanks: {key}")
        locus = "L" + digest(args.graph_version, args.graph_sha256, row["reference_build"],
                             row["chromosome"], start, ref)[:20]
        allele = locus + ".A" + digest(alt)[:16]
        if allele in existing:
            raise ValueError(f"Control duplicates an existing allele: {key}")
        existing.add(allele)
        common = dict(candidate_id=key, locus_id=locus, allele_id=allele,
                      reference_build=row["reference_build"], chromosome=row["chromosome"],
                      graph_version=args.graph_version, graph_sha256=args.graph_sha256,
                      status="PASS_CONTROL", reason="PREPARED_CONTROL_REQUIRES_RESOURCE_REVIEW")
        records.append(dict(common, reference_start_0=start, reference_end_0=end,
                            sv_type=row["sv_type"], insertion_size_bp=right-left,
                            allele_start=left, allele_end=right, sequence_sha256=digest(alt),
                            left_flank_bp=left, right_flank_bp=len(alt)-right,
                            source_sequence_status="PREPARED_CONTROL", control_evidence=row["evidence"]))
        crosswalk.append(dict(common, start_0=start, end_0=end,
                              correspondence_status="NATIVE_COORDINATE"))
        for kind, sequence, astart, aend in (
                ("REF", ref, ref_left, ref_right), ("ALT", alt, left, right),
                ("LEFT", alt[:left], "", ""), ("RIGHT", alt[right:], "", "")):
            name = allele + "|" + kind
            sequences[name] = sequence
            if kind in {"REF", "ALT"}:
                anchored[name] = sequence
            queries.append(dict(common, query_name=name, query_type=kind,
                                query_length=len(sequence), allele_start=astart, allele_end=aend,
                                start_0=start, end_0=end))


def representation(args):
    catalogue_rows = list(rows(args.catalogue))
    query_rows, search_rows = list(rows(args.queries)), list(rows(args.search_manifest))
    by_query = {r["query_name"]: r for r in query_rows}
    searches = {r["resource"]: r for r in search_rows}
    if len(by_query) != len(query_rows) or len(searches) != len(search_rows):
        raise ValueError("Duplicate query or resource in search manifests")
    if set(searches) - set(RESOURCES):
        raise ValueError("Resource names must be HPRC, GRCh38, CHM13")
    reviews = {}
    for row in rows(args.review):
        key = row["allele_id"], row["resource"]
        if (key in reviews or row["resource"] not in RESOURCES
                or row["allele_id"] not in {r["allele_id"] for r in catalogue_rows}):
            raise ValueError(f"Duplicate or unknown resource review: {key}")
        reviews[key] = row
    evidence, hits = [], Counter()
    for resource, search in searches.items():
        path = search.get("alignment_path", "")
        if not path or not Path(path).is_file():
            continue
        with open_text(path) as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.rstrip().split("\t")
                if len(f) < 12 or f[0] not in by_query:
                    raise ValueError(f"Malformed alignment or unknown query in {path}")
                query = by_query[f[0]]
                if int(f[1]) != int(query["query_length"]):
                    raise ValueError(f"Query length disagreement: {f[0]}")
                hits[query["allele_id"], resource, query["query_type"]] += 1
                evidence.append(dict(allele_id=query["allele_id"], locus_id=query["locus_id"], resource=resource,
                    query_name=f[0], query_type=query["query_type"], query_length=f[1], query_start_0=f[2],
                    query_end_0=f[3], strand=f[4], target=f[5], target_length=f[6], target_start_0=f[7],
                    target_end_0=f[8], matching_bases=f[9], alignment_block=f[10], mapq=f[11],
                    alignment_tags=";".join(f[12:]), alignment_path=path,
                    graph_version=query["graph_version"], resource_version=search["resource_version"],
                    status="ALIGNMENT_SCREEN_ONLY", reason="REQUIRES_LOCUS_AND_STRUCTURAL_REVIEW"))
    matrix = []
    for allele in catalogue_rows:
        for resource in RESOURCES:
            search, review = searches.get(resource, {}), reviews.get((allele["allele_id"], resource), {})
            available = bool(search.get("alignment_path")) and Path(search["alignment_path"]).is_file()
            status = "UNRESOLVED" if available else "NOT_ASSESSED"
            reason = "NO_COMPLETED_LOCUS_REVIEW" if available else "RESOURCE_OR_SEARCH_UNAVAILABLE"
            if review:
                requested = review["status"]
                if requested not in STATUSES:
                    raise ValueError(f"Unknown representation status: {requested}")
                if requested in {"REPRESENTED", "PARTIAL_OR_DIFFERENT_ALLELE", "NOT_REPRESENTED_AT_RESOLVED_LOCUS"}:
                    checksum = search.get("resource_sha256", "")
                    if (not search.get("resource_version") or len(checksum) != 64
                            or set(checksum.lower()) - set("0123456789abcdef")):
                        raise ValueError(f"Freeze resource version and SHA256 before a definitive decision: {resource}")
                    required = all(yes(review.get(field)) for field in (
                        "locus_resolved", "search_adequate", "structural_configuration_reviewed",
                        "both_flanks_syntenic", "comparator_locus_covered"))
                    if (not available or not yes(search.get("search_complete")) or not required
                            or not Path(review.get("evidence", "")).is_file() or not review.get("reason")):
                        raise ValueError(f"Unsupported definitive resource decision: {allele['allele_id']} {resource}")
                    if allele["status"] not in {"PASS", "PASS_CONTROL"}:
                        raise ValueError("Source-anchor-pending alleles cannot receive definitive resource status")
                    if requested == "REPRESENTED" and not yes(review.get("whole_allele_and_junctions_supported")):
                        raise ValueError("REPRESENTED requires whole-allele and junction evidence")
                status, reason = requested, review.get("reason", "RECORDED_REVIEW")
            matrix.append(dict(allele_id=allele["allele_id"], locus_id=allele["locus_id"], resource=resource,
                graph_version=allele["graph_version"], resource_version=search.get("resource_version", ""),
                resource_sha256=search.get("resource_sha256", ""), status=status, reason=reason,
                alt_alignment_count=hits[allele["allele_id"], resource, "ALT"],
                left_alignment_count=hits[allele["allele_id"], resource, "LEFT"],
                right_alignment_count=hits[allele["allele_id"], resource, "RIGHT"],
                locus_resolved=review.get("locus_resolved", "false"),
                search_adequate=review.get("search_adequate", "false"),
                structural_configuration_reviewed=review.get("structural_configuration_reviewed", "false"),
                both_flanks_syntenic=review.get("both_flanks_syntenic", "false"),
                comparator_locus_covered=review.get("comparator_locus_covered", "false"),
                whole_allele_and_junctions_supported=review.get("whole_allele_and_junctions_supported", "false"),
                resource_chromosome=review.get("resource_chromosome", ""),
                resource_start_0=review.get("resource_start_0", ""),
                resource_end_0=review.get("resource_end_0", ""),
                resource_strand=review.get("resource_strand", ""),
                evidence=review.get("evidence", ""), elsewhere_similarity=review.get("elsewhere_similarity", "NOT_REVIEWED"),
                threshold_sensitive=review.get("threshold_sensitive", "NOT_ASSESSED")))
    strict = {r["allele_id"] for r in catalogue_rows if all(
        m["status"] == "NOT_REPRESENTED_AT_RESOLVED_LOCUS" for m in matrix if m["allele_id"] == r["allele_id"])}
    for row in matrix:
        row["strict_three_resource"] = str(row["allele_id"] in strict).lower()
    out = Path(args.output)
    write(out / "resource_matrix.tsv", matrix, ["allele_id", "locus_id", "resource", "status", "reason"])
    write(out / "evidence.tsv", evidence, ["allele_id", "locus_id", "resource", "query_name", "status", "reason"])
    counts = [dict(resource=resource, status=status, allele_count=count) for (resource, status), count in
              Counter((r["resource"], r["status"]) for r in matrix).items()]
    write(out / "resource_counts.tsv", counts, ["resource", "status", "allele_count"])
    flow_path = Path(args.catalogue).parent / "candidate_flow.tsv"
    attempted = list(rows(flow_path))
    local = {r["allele_id"] for r in attempted if r.get("allele_id")}
    resolved = {allele for allele in local if all(m["status"] not in {"UNRESOLVED", "NOT_ASSESSED"}
                for m in matrix if m["allele_id"] == allele)}
    write(out / "flow.tsv", [dict(metric=name, count=value) for name, value in (
        ("all_attempted_candidate_records", len(attempted)),
        ("rejected_candidate_records", sum(r["status"] == "REJECTED" for r in attempted)),
        ("source_anchor_pending_candidate_records", sum(r["status"] == "SOURCE_INSERT_CONFIRMED_ANCHORS_PENDING" for r in attempted)),
        ("reconstructed_distinct_local_alt_alleles", len(local)),
        ("three_resource_resolved_distinct_local_alt_alleles", len(resolved)),
        ("strict_additional_distinct_local_alt_alleles", len(strict & local)),
        ("strict_additional_distinct_local_loci", len({r["locus_id"] for r in catalogue_rows if r["allele_id"] in strict & local})),
        ("strict_additional_attempted_candidate_records", sum(r.get("allele_id") in strict & local for r in attempted)))],
        ["metric", "count"])


def panel(args):
    catalogue_rows, selection = list(rows(args.catalogue)), list(rows(args.selection))
    roles = Counter(r["panel_role"] for r in selection)
    if roles != Counter(candidate_addition=20, represented_control=5, technical_control=5):
        raise ValueError("Freeze exactly 30 loci: 20 candidate additions, 5 represented and 5 technical controls")
    if len({r["locus_id"] for r in selection}) != 30 or any(not r["selection_reason"] for r in selection):
        raise ValueError("Require 30 distinct locus IDs and phenotype-blind selection reasons")
    selected = {r["locus_id"]: r for r in selection}
    known = {r["locus_id"] for r in catalogue_rows}
    if set(selected) - known:
        raise ValueError(f"Unknown selected loci: {sorted(set(selected) - known)}")
    chosen = [dict(row, panel_role=selected[row["locus_id"]]["panel_role"],
                   selection_reason=selected[row["locus_id"]]["selection_reason"])
              for row in catalogue_rows if row["locus_id"] in selected]
    if any(r["status"] not in {"PASS", "PASS_CONTROL"} for r in chosen):
        raise ValueError("Resolve source/anchor QC before freezing selected alleles")
    if any(r["panel_role"] == "candidate_addition" and r["status"] != "PASS" for r in chosen):
        raise ValueError("Candidate additions must come from the corrected local catalogue, not imported controls")
    tmem = [r for r in chosen if r["candidate_id"] == "TMEM106B"]
    if not tmem or any(r["panel_role"] != "represented_control" for r in tmem):
        raise ValueError("TMEM106B must be selected among the five represented controls")
    ids = {r["allele_id"] for r in chosen}
    matrix = list(rows(args.resource_matrix))
    for row in chosen:
        if row["panel_role"] == "represented_control" and not any(
                m["allele_id"] == row["allele_id"] and m["status"] == "REPRESENTED" for m in matrix):
            raise ValueError("Every represented-control ALT requires recorded representation in at least one resource")
    queries = [r for r in rows(args.queries) if r["allele_id"] in ids]
    fasta = read_fasta(args.fasta)
    sequences = {r["query_name"]: fasta[r["query_name"]] for r in queries}
    out = Path(args.output)
    write(out / "panel.tsv", chosen, ["locus_id", "allele_id", "panel_role", "selection_reason", "status", "reason"])
    write(out / "query_manifest.tsv", queries, ["query_name", "locus_id", "allele_id", "query_type", "query_length"])
    write_fasta(out / "locus_queries.fa", sequences)
    write(out / "counts.tsv", [dict(panel_role=role, loci=number,
          alt_alleles=sum(r["panel_role"] == role for r in chosen)) for role, number in sorted(roles.items())],
          ["panel_role", "loci", "alt_alleles"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("catalogue")
    for flag in ("candidates", "alleles", "reference", "graph-version", "graph-sha256", "output"):
        p.add_argument("--" + flag, required=True)
    p.add_argument("--source-review")
    p.add_argument("--tmem-fasta")
    p.add_argument("--tmem-targets")
    p.add_argument("--controls")
    p.add_argument("--control-fasta")
    p.add_argument("--flank", type=int, default=1000)
    p.set_defaults(func=catalogue)
    p = sub.add_parser("representation")
    for flag in ("catalogue", "queries", "search-manifest", "review", "output"):
        p.add_argument("--" + flag, required=True)
    p.set_defaults(func=representation)
    p = sub.add_parser("panel")
    for flag in ("catalogue", "queries", "fasta", "selection", "resource-matrix", "output"):
        p.add_argument("--" + flag, required=True)
    p.set_defaults(func=panel)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
