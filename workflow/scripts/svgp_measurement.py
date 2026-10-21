#!/usr/bin/env python3
"""Graph molecule evidence and conservative, phenotype-blind measurement tables.

Native vg/minimap2/Sniffles commands remain in Snakemake. Python is needed here
for graph alignment edits, molecule joins and explicit missing-value accounting.
"""

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path


MISSING = {"", ".", "NA", "unknown"}
EVIDENCE_FIELDS = [
    "person_id", "analysis_group", "locus_id", "allele_id", "graph_id", "graph_sha256",
    "catalogue_graph_sha256", "node_space_id", "reference_build", "ref_support", "alt_left_junction_support",
    "alt_right_junction_support", "alt_internal_depth", "flank_depth",
    "allele_balance", "mapping_ambiguity", "duplicate_handling",
    "unassigned_fragments", "carrier_status", "genotype", "dosage",
    "measurement_type", "confidence", "callable", "status", "reason",
]


def open_text(path):
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def rows(path):
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_table(path, records, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def unique(records, keys):
    result = {}
    for row in records:
        key = tuple(row[k] for k in keys)
        if key in result:
            raise ValueError(f"Duplicate {keys}: {key}")
        result[key] = row
    return result


def require_file(value):
    if value in MISSING or not Path(value).is_file() or not Path(value).stat().st_size:
        raise ValueError(f"Missing review/evidence file: {value}")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def audit_bundle(args):
    graph = json.loads(args.graph_config)
    manifest = unique(rows(args.manifest), ["component"])
    checked = []
    for component in ("gbz", "dist", "minimizer", "zipcodes"):
        row = manifest[(component,)]
        if (row["graph_id"] != args.graph or row["node_space_id"] != graph["node_space_id"]
                or row["path"] != graph[component] or row["status"] != "VERIFIED"):
            raise ValueError(f"{component}: bundle provenance/path mismatch")
        digest = sha256(row["path"])
        if digest != row["sha256"]:
            raise ValueError(f"{component}: checksum mismatch")
        if component == "gbz" and digest != graph["graph_sha256"]:
            raise ValueError("graph_sha256 must identify the exact mapped GBZ")
        checked.append(dict(row, observed_sha256=digest))
    reviews = unique(rows(args.reference_review), ["path_name"]) if graph.get("reference_paths") else {}
    for name in graph.get("reference_paths", []):
        row = reviews[(name,)]
        if (row["status"] != "VERIFIED" or row["path_role"] != "REFERENCE_FASTA"
                or row["node_space_id"] != graph["node_space_id"]
                or row["graph_id"] != args.graph
                or not re.fullmatch(r"[0-9a-f]{64}", row["path_sequence_sha256"])
                or row["path_sequence_sha256"] != row["reference_sequence_sha256"]):
            raise ValueError(f"{name}: exact reference spelling not verified")
        require_file(row["evidence"])
    write_table(args.output, checked, list(checked[0]))


def flip(handle):
    if not re.fullmatch(r"[<>][0-9]+", handle):
        raise ValueError(f"Expected oriented numeric node, got {handle}")
    return ("<" if handle[0] == ">" else ">") + handle[1:]


def fragment_name(name):
    # Only the explicit standard mate suffix is normalized, never sample IDs.
    return re.sub(r"/[12]$", "", name.split()[0])


def merge_intervals(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


def merge_length(intervals):
    return sum(end-start for start, end in merge_intervals(intervals))


def alignment_features(alignment, edge_lookup, node_lookup, anchor_matches):
    """Match oriented junctions with consecutive exact matches on both sides.

    A mismatch, insertion or deletion interrupts an anchoring run. Edges in a
    clipped/unmatched part of a graph path cannot contribute support.
    """
    mappings = alignment.get("path", {}).get("mapping", [])
    tokens, edits, boundaries, intervals = [], [], [], []
    cursor = 0
    for mapping in mappings:
        pos = mapping.get("position", {})
        node = str(pos.get("node_id", ""))
        reverse = pos.get("is_reverse", False)
        token = ("<" if reverse else ">") + node
        tokens.append(token)
        node_cursor = int(pos.get("offset", 0))
        for edit in mapping.get("edit", []):
            source, target = int(edit.get("from_length", 0)), int(edit.get("to_length", 0))
            exact = source == target and source > 0 and not edit.get("sequence")
            edits.append((cursor, cursor + target, exact))
            if exact and node in node_lookup:
                for key, feature, start, end, length in node_lookup[node]:
                    a, b = node_cursor, node_cursor + source
                    if reverse:
                        a, b = length - b, length - a
                    if max(a, start) < min(b, end):
                        intervals.append((key, feature, node, max(a, start), min(b, end)))
            cursor += target
            node_cursor += source
        boundaries.append(cursor)
    # Runs are joined only across exact edits; a deletion has target length 0
    # and must still break a run.
    runs, start, end = [], None, None
    for a, b, exact in edits:
        if exact:
            if start is None:
                start = a
            end = b
        elif start is not None:
            runs.append((start, end))
            start, end = None, None
    if start is not None:
        runs.append((start, end))
    hits = set()
    for index, edge in enumerate(zip(tokens, tokens[1:])):
        boundary = boundaries[index]
        if any(a <= boundary - anchor_matches and b >= boundary + anchor_matches
               for a, b in runs):
            hits.update(edge_lookup.get(edge, ()))
    return hits, intervals


def extract(args):
    graph = json.loads(args.graph_config)
    sample = unique(rows(args.samples), ["person_id"])[(args.person,)]
    panel = unique(rows(args.panel), ["locus_id", "allele_id"])
    targets = list(rows(args.targets))
    edge_lookup, node_lookup = defaultdict(list), defaultdict(list)
    available, denominator = defaultdict(set), Counter()
    target_intervals, node_lengths = defaultdict(list), {}
    for row in targets:
        key = (row["locus_id"], row["allele_id"])
        if key not in panel:
            raise ValueError(f"Target outside frozen panel: {key}")
        if (row["graph_id"] != args.graph or row["graph_sha256"] != graph["graph_sha256"]
                or row["node_space_id"] != graph["node_space_id"]
                or row["catalogue_graph_sha256"] != panel[key]["graph_sha256"]):
            raise ValueError(f"Assay node-space mismatch: {key}")
        if row["status"] != "VERIFIED":
            continue
        require_file(row["evidence"])
        feature = row["feature"]
        available[key].add(feature)
        if feature in {"ref", "left", "right"}:
            edge = (row["edge_from"], row["edge_to"])
            reverse = (flip(edge[1]), flip(edge[0]))
            edge_lookup[edge].append((key, feature))
            edge_lookup[reverse].append((key, feature))
        elif feature in {"flank", "internal"}:
            start, end, length = int(row["start0"]), int(row["end0"]), int(row["node_length"])
            if not 0 <= start < end <= length:
                raise ValueError(f"Invalid graph interval: {row}")
            node = row["node_id"]
            if node in node_lengths and node_lengths[node] != length:
                raise ValueError(f"Inconsistent target node length: {node}")
            node_lengths[node] = length
            target_intervals[key, feature, node].append((start, end))
        else:
            raise ValueError(f"Unknown target feature: {feature}")
    # The assay denominator and covered bases both use interval unions: repeated
    # or overlapping target records must not artificially lower depth.
    for (key, feature, node), intervals in target_intervals.items():
        for start, end in merge_intervals(intervals):
            node_lookup[node].append((key, feature, start, end, node_lengths[node]))
            denominator[key, feature] += end-start
    # Keep only panel-overlapping observations, then stream the potentially
    # whole-genome assignment table, retaining only those fragment names.
    observations, wanted = [], set()
    counts, depth = defaultdict(set), defaultdict(list)
    seen, ambiguous, unassigned = defaultdict(set), defaultdict(set), defaultdict(set)
    for library in rows(Path(args.mapping) / "libraries.tsv"):
        library_id = library["library_id"]
        process = subprocess.Popen(["vg", "view", "-aj", library["gam"]],
                                   stdout=subprocess.PIPE, text=True)
        try:
            for line in process.stdout:
                alignment = json.loads(line)
                hits, intervals = alignment_features(alignment, edge_lookup, node_lookup,
                                                     args.anchor_matches)
                touched = {key for key, feature in hits} | {x[0] for x in intervals}
                if not touched:
                    continue
                name = fragment_name(alignment["name"])
                raw_fragment = (library_id, name)
                wanted.add(raw_fragment)
                low_quality = (alignment.get("is_secondary", False)
                               or int(alignment.get("mapping_quality", 0)) < args.min_mapq)
                observations.append((raw_fragment, hits, intervals, touched, low_quality))
        finally:
            if process.stdout:
                process.stdout.close()
            if process.wait() != 0:
                raise RuntimeError(f"vg view failed for {library['gam']}")
    # These are reviewed physical-molecule assignments, not read-name dedup.
    molecules = {}
    for row in rows(args.molecules):
        key = (row["library_id"], fragment_name(row["read_name"]))
        if row["person_id"] != args.person or key not in wanted:
            continue
        if row["status"] != "VERIFIED" or not row["molecule_id"] or not row["method"]:
            continue
        require_file(row["evidence"])
        molecule = (row["library_id"], row["molecule_id"])
        if key in molecules and molecules[key] != molecule:
            raise ValueError(f"Conflicting molecule assignment: {key}")
        molecules[key] = molecule
    for raw_fragment, hits, intervals, touched, low_quality in observations:
        molecule = molecules.get(raw_fragment)
        for key in touched:
            seen[key].add(raw_fragment)
            if molecule is None:
                unassigned[key].add(raw_fragment)
        if low_quality:
            for key in touched:
                ambiguous[key].add(raw_fragment)
            continue
        if molecule is None:
            continue
        for key, feature in hits:
            counts[key, feature].add(molecule)
        for key, feature, node, a, b in intervals:
            depth[key, feature, node, molecule].append((a, b))
    conflicts = {}
    for key in panel:
        conflicts[key] = counts[key, "ref"] & (counts[key, "left"] | counts[key, "right"])
        for feature in ("ref", "left", "right"):
            counts[key, feature].difference_update(conflicts[key])
        ambiguous[key].update(fragment for fragment, molecule in molecules.items()
                              if molecule in conflicts[key])
    depth_bases = Counter()
    for (key, feature, node, molecule), intervals in depth.items():
        if molecule not in conflicts[key]:
            depth_bases[key, feature] += merge_length(intervals)
    output = []
    for key in panel:
        ref, left, right = (len(counts[key, feature]) for feature in ("ref", "left", "right"))
        alt = min(left, right)
        reason = ("ASSAY_TARGETS_UNVERIFIED" if not {"ref", "left", "right"} <= available[key]
                  else "DUPLICATE_ASSIGNMENTS_INCOMPLETE" if unassigned[key]
                  else "UNCALIBRATED")
        record = dict(zip(("locus_id", "allele_id"), key), person_id=args.person,
                      analysis_group=sample["analysis_group"],
                      graph_id=args.graph, graph_sha256=graph["graph_sha256"],
                      catalogue_graph_sha256=panel[key]["graph_sha256"],
                      node_space_id=graph["node_space_id"], reference_build=graph["reference_build"],
                      ref_support=ref, alt_left_junction_support=left,
                      alt_right_junction_support=right,
                      alt_internal_depth=(depth_bases[key, "internal"] / denominator[key, "internal"]
                                          if denominator[key, "internal"] else "."),
                      flank_depth=(depth_bases[key, "flank"] / denominator[key, "flank"]
                                   if denominator[key, "flank"] else "."),
                      allele_balance=alt / (ref + alt) if ref + alt else ".",
                      mapping_ambiguity=len(ambiguous[key]) / len(seen[key]) if seen[key] else ".",
                      duplicate_handling="reviewed_physical_molecule_ids_within_library",
                      unassigned_fragments=len(unassigned[key]), carrier_status="NO_CALL",
                      genotype="./.", dosage=".", measurement_type="uncalibrated_support",
                      confidence="NOT_CALIBRATED", callable="no", status="NO_CALL", reason=reason)
        output.append(record)
    write_table(args.output, output, EVIDENCE_FIELDS)


def apply_lock(record, lock, min_mapq, anchor_matches):
    record = dict(record)
    if record["reason"] != "UNCALIBRATED" or lock is None or lock["status"] != "LOCKED":
        return record
    if (lock["graph_sha256"] != record["graph_sha256"]
            or lock["catalogue_graph_sha256"] != record["catalogue_graph_sha256"]
            or lock["node_space_id"] != record["node_space_id"]
            or int(lock["min_mapq"]) != min_mapq
            or int(lock["anchor_matches"]) != anchor_matches
            or lock["locked_before_phenotypes"] != "yes"):
        raise ValueError("Measurement lock differs from extraction/graph or was not phenotype-blind")
    if lock["technical_status"] != "PASS":
        record["reason"] = "LOCUS_FAILED_TECHNICAL_REVIEW"
        return record
    for field in ("calibration_evidence", "threshold_sensitivity_evidence", "depth_batch_review", "technical_review"):
        require_file(lock[field])
    for field in ("min_ref_molecules", "min_alt_left_molecules", "min_alt_right_molecules"):
        if int(lock[field]) <= 0:
            raise ValueError(f"Locked {field} must be positive")
    if (not math.isfinite(float(lock["min_flank_depth"])) or float(lock["min_flank_depth"]) <= 0
            or not 0 <= float(lock["max_mapping_ambiguity"]) <= 1):
        raise ValueError("Invalid locked depth/ambiguity limits")
    if lock["measurement_type"] != "carrier":
        raise ValueError("This junction assay supports calibrated carrier status only; dosage needs a validated genotype parser")
    for field in ("flank_depth", "mapping_ambiguity"):
        if record[field] in MISSING:
            record["reason"] = "MISSING_COVERAGE_OR_AMBIGUITY"
            return record
    if (float(record["flank_depth"]) < float(lock["min_flank_depth"])
            or float(record["mapping_ambiguity"]) > float(lock["max_mapping_ambiguity"])):
        record["reason"] = "LOW_DEPTH_OR_AMBIGUOUS_MAPPING"
        return record
    ref = int(record["ref_support"])
    left, right = (int(record[f"alt_{side}_junction_support"]) for side in ("left", "right"))
    alt = left >= int(lock["min_alt_left_molecules"]) and right >= int(lock["min_alt_right_molecules"])
    if alt:
        carrier = "ALT_SUPPORTED"
    elif left == right == 0 and ref >= int(lock["min_ref_molecules"]):
        carrier = "REF_SUPPORTED"
    else:
        record["reason"] = "INSUFFICIENT_OR_ONE_SIDED_ALLELE_SUPPORT"
        return record
    record.update(carrier_status=carrier, measurement_type="carrier", confidence="LOCKED_CARRIER",
                  callable="yes", status="PASS", reason=".")
    # REF-supported carrier status is still not a calibrated diploid 0/0 call.
    return record


def summarize(args):
    locks = unique(rows(args.lock), ["graph_id", "locus_id", "allele_id"]) if args.lock else {}
    evidence = []
    for path in args.inputs:
        for row in rows(path):
            key = (row["graph_id"], row["locus_id"], row["allele_id"])
            evidence.append(apply_lock(row, locks.get(key), args.min_mapq, args.anchor_matches))
    unique(evidence, ["person_id", "graph_id", "locus_id", "allele_id"])
    if not evidence:
        raise ValueError("No attempted measurement cells; configure people and graphs")
    write_table(args.output, evidence, EVIDENCE_FIELDS)
    grouped = defaultdict(list)
    for row in evidence:
        grouped[row["graph_id"], row["analysis_group"], row["locus_id"], row["allele_id"]].append(row)
    summaries = []
    for (graph, analysis_group, locus, allele), cells in grouped.items():
        callable_n = sum(row["callable"] == "yes" for row in cells)
        carriers = sum(row["carrier_status"] == "ALT_SUPPORTED" for row in cells)
        summaries.append(dict(graph_id=graph, analysis_group=analysis_group, locus_id=locus, allele_id=allele,
                              attempted=len(cells), callable=callable_n, carriers=carriers,
                              no_call=len(cells) - callable_n, callability=callable_n / len(cells),
                              carrier_frequency=carriers / callable_n if callable_n else ".",
                              allele_count=".", callable_ploidy=".", allele_frequency=".",
                              measurement_type="carrier", status="PASS" if callable_n else "NO_CALL",
                              reason="carrier_only_no_allele_frequency"))
    write_table(args.summary, summaries, list(summaries[0]))
    if args.paired_summary:
        paired = []
        for person in sorted({r["person_id"] for r in evidence}):
            groups = defaultdict(dict)
            for row in evidence:
                if row["person_id"] == person:
                    groups[row["graph_id"]][row["locus_id"], row["allele_id"]] = row
            main, comparator = groups["hprc_plus_local"], groups["hprc_only"]
            shared = sorted(set(main) & set(comparator))
            both = [key for key in shared if main[key]["callable"] == comparator[key]["callable"] == "yes"]
            paired.append(dict(person_id=person, attempted_alleles=len(shared),
                               main_callable=sum(main[key]["callable"] == "yes" for key in shared),
                               comparator_callable=sum(comparator[key]["callable"] == "yes" for key in shared),
                               both_callable=len(both),
                               carrier_agree=sum(main[key]["carrier_status"] == comparator[key]["carrier_status"] for key in both),
                               status="DESCRIPTIVE", reason="paired_person_summary_not_accuracy_without_ONT"))
        write_table(args.paired_summary, paired, list(paired[0]))
    if args.matrix_dir:
        out = Path(args.matrix_dir)
        people = sorted({r["person_id"] for r in evidence})
        alleles = sorted({(r["locus_id"], r["allele_id"]) for r in evidence})
        cells = unique(evidence, ["person_id", "locus_id", "allele_id"])
        for name, field in (("carrier_matrix", "carrier_status"), ("dosage_matrix", "dosage")):
            table = []
            for locus, allele in alleles:
                record = dict(locus_id=locus, allele_id=allele,
                              measurement_type="carrier" if name == "carrier_matrix" else "NOT_CALIBRATED")
                for person in people:
                    value = cells[(person, locus, allele)][field]
                    record[person] = ("1" if value == "ALT_SUPPORTED" else "0" if value == "REF_SUPPORTED"
                                      else ".") if field == "carrier_status" else value
                table.append(record)
            write_table(out / f"{name}.tsv.gz", table,
                        ["locus_id", "allele_id", "measurement_type"] + people)
        write_table(out / "callability.tsv.gz", evidence, EVIDENCE_FIELDS)


def validate_ont(args):
    panel = unique(rows(args.panel), ["locus_id", "allele_id"])
    samples = unique(rows(args.samples), ["person_id"])
    paired = unique(rows(args.pairing), ["person_id"])
    adjudicated = unique(rows(args.adjudication), ["person_id", "locus_id", "allele_id"])
    illumina = unique((r for r in rows(args.illumina) if r["graph_id"] == args.graph),
                      ["person_id", "locus_id", "allele_id"])
    results = []
    for person in args.people:
        if (person,) not in samples:
            raise ValueError(f"Unknown ONT person: {person}")
        identity = paired.get((person,), {}).get("status", "NOT_PAIRED")
        for locus, allele in panel:
            key = (person, locus, allele)
            row = dict(adjudicated.get(key, {}))
            resolved = row.get("status") in {"RESOLVED_CARRIER", "RESOLVED_GENOTYPE"}
            if resolved:
                if (row["reference_build"] != args.reference_build
                        or row["identity_confirmed"] != "yes" or identity != "DECLARED_FINGERPRINT_CONFIRMED"
                        or row["development_or_held_out"] not in {"development", "held_out"}
                        or row["carrier_status"] not in {"ALT_SUPPORTED", "REF_SUPPORTED"}):
                    raise ValueError(f"Invalid resolved ONT identity/evidence cell: {key}")
                require_file(row["evidence"])
                for field in ("alt_left_molecules", "alt_right_molecules", "ref_molecules", "local_depth"):
                    if row[field] in MISSING or not math.isfinite(float(row[field])) or float(row[field]) < 0:
                        raise ValueError(f"Missing resolved ONT evidence: {key} {field}")
                if float(row["local_depth"]) <= 0 or (row["carrier_status"] == "ALT_SUPPORTED" and
                        min(float(row["alt_left_molecules"]), float(row["alt_right_molecules"])) <= 0):
                    raise ValueError(f"Resolved ONT ALT needs positive depth and both junctions: {key}")
                if row["carrier_status"] == "REF_SUPPORTED" and float(row["ref_molecules"]) <= 0:
                    raise ValueError(f"Resolved ONT REF needs reference support: {key}")
            short = illumina.get(key, {})
            evaluable = resolved and short.get("callable") == "yes"
            results.append(dict(row, person_id=person, locus_id=locus, allele_id=allele,
                                reference_build=args.reference_build,
                                status=row.get("status", "NOT_EVALUABLE"),
                                reason=row.get("reason", "ONT_READ_ADJUDICATION_REQUIRED"),
                                illumina_carrier=short.get("carrier_status", "NO_CALL"),
                                comparison_evaluable="yes" if evaluable else "no",
                                concordant=("yes" if short["carrier_status"] == row["carrier_status"] else "no")
                                           if evaluable else "."))
    if not results:
        raise ValueError("No configured ONT people; ONT validation is not completed")
    fields = list(dict.fromkeys(key for row in results for key in row))
    write_table(args.output, results, fields)
    summary = []
    for group in ("development", "held_out"):
        cells = [r for r in results if r.get("development_or_held_out") == group]
        comparable = [r for r in cells if r["comparison_evaluable"] == "yes"]
        nonref = [r for r in comparable if r.get("carrier_status") == "ALT_SUPPORTED"
                  or r["illumina_carrier"] == "ALT_SUPPORTED"]
        for scope, subset in (("all_resolved", comparable), ("nonreference_union", nonref)):
            summary.append(dict(validation_set=group, scope=scope, attempted=len(cells),
                                evaluable=len(subset), concordant=sum(r["concordant"] == "yes" for r in subset),
                                no_call=len(cells) - len(comparable), status="DESCRIPTIVE",
                                reason="correlated_cells_no_pooled_accuracy_interval"))
    write_table(args.summary, summary, list(summary[0]))


def release(args):
    reviews = [r for r in rows(args.review) if r["graph_id"] == args.graph]
    if len(reviews) != 1 or reviews[0]["status"] != "APPROVED":
        raise ValueError("A single reviewed cohort release is required")
    row = reviews[0]
    for field in ("pilot_evidence", "ont_evidence", "calibration_evidence", "compute_budget_evidence"):
        require_file(row[field])
    # Normalize only the path spelling; mounted and server paths must never be
    # resolved through symlinks or substituted for one another.
    for field, consumed in (("pilot_evidence", args.pilot), ("ont_evidence", args.ont)):
        if os.path.abspath(row[field]) != os.path.abspath(consumed):
            raise ValueError(f"Cohort release {field} must identify the actual consumed input")
    if not row["retention_policy"] or not row["reviewer"]:
        raise ValueError("Missing compute retention policy/reviewer")
    locks = [r for r in rows(args.lock) if r["graph_id"] == args.graph
             and r["status"] == "LOCKED" and r["technical_status"] == "PASS"]
    locked = {(r["locus_id"], r["allele_id"]) for r in locks}
    calibrated = {(r["locus_id"], r["allele_id"]) for r in rows(args.pilot)
                  if r["graph_id"] == args.graph and r["callable"] == "yes"} & locked
    if not calibrated:
        raise ValueError("No locked, callable pilot alleles")
    ont = list(rows(args.ont))
    if not any(r["status"] in {"RESOLVED_CARRIER", "RESOLVED_GENOTYPE"}
               and (r["locus_id"], r["allele_id"]) in calibrated for r in ont):
        raise ValueError("No resolved ONT cells for locked, callable main-graph alleles")
    write_table(args.output, [row], list(row))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    audit = commands.add_parser("audit-bundle")
    for name in ("manifest", "reference-review", "graph", "graph-config", "output"):
        audit.add_argument("--" + name, required=True)
    extraction = commands.add_parser("extract")
    for name in ("panel", "targets", "molecules", "mapping", "samples", "person", "graph", "graph-config", "output"):
        extraction.add_argument("--" + name, required=True)
    summary = commands.add_parser("summarize")
    summary.add_argument("--inputs", nargs="+", required=True)
    summary.add_argument("--lock")
    summary.add_argument("--matrix-dir")
    summary.add_argument("--paired-summary")
    for name in ("output", "summary"):
        summary.add_argument("--" + name, required=True)
    for command in (extraction, summary):
        command.add_argument("--min-mapq", type=int, required=True)
        command.add_argument("--anchor-matches", type=int, required=True)
    ont = commands.add_parser("validate-ont")
    for name in ("panel", "samples", "pairing", "adjudication", "illumina", "graph", "reference-build", "output", "summary"):
        ont.add_argument("--" + name, required=True)
    ont.add_argument("--people", nargs="+", required=True)
    gate = commands.add_parser("release")
    for name in ("review", "graph", "lock", "pilot", "ont", "output"):
        gate.add_argument("--" + name, required=True)
    args = parser.parse_args()
    if args.command in {"extract", "summarize"} and (args.min_mapq < 0 or args.anchor_matches <= 0):
        parser.error("Choose explicit valid extraction settings using the development controls")
    {"audit-bundle": audit_bundle, "extract": extract, "summarize": summarize,
     "validate-ont": validate_ont, "release": release}[args.command](args)


if __name__ == "__main__":
    main()
