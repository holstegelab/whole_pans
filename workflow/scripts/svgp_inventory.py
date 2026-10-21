#!/usr/bin/env python3
"""Join curated people, explicit read lanes and frozen resource metadata.

This inventory only stats raw files. It does not read FASTQ contents, hash large
resources, infer sample aliases, calculate relatedness or verify fingerprints.
Run from the server repository so relative data paths have their intended meaning.
"""

import argparse
import csv
import math
import os
from collections import Counter, defaultdict
from pathlib import Path


UNKNOWN = {"", ".", "na", "n/a", "unknown", "not_assessed", "none"}
SAMPLE_FIELDS = "person_id illumina_id fastq_manifest_id phenotype_id ont_id discovery_donor related_to_discovery_donor phenotype sex age_definition ancestry_covariates batch library_id read_length depth analysis_eligible exclusion_reason".split()
READ_FIELDS = "read_id person_id illumina_id fastq_manifest_id library_id lane batch r1 r2".split()
RESOURCE_FIELDS = "resource_id resource_type path build release scope sha256 checksum_source checksum_verified donor_manifest annotation_release software_version".split()
IDENTIFIERS = ("illumina_id", "fastq_manifest_id", "phenotype_id", "ont_id")


def protect_inputs(inputs, outputs):
    """Reject outputs aliasing inputs or each other without resolving mount paths."""
    previous = [Path(path) for path in inputs if path]
    for output in map(Path, outputs):
        for other in previous:
            if (os.path.abspath(output) == os.path.abspath(other)
                    or (output.exists() and other.exists() and output.samefile(other))):
                raise ValueError(f"Output path would overwrite an input or another output: {output}")
        previous.append(output)


def known(value):
    return str(value).strip().lower() not in UNKNOWN


def flag(value):
    value = str(value).strip().lower()
    if value in {"yes", "true", "1"}:
        return "yes"
    if value in {"no", "false", "0"}:
        return "no"
    if value in UNKNOWN:
        return "unknown"
    raise ValueError(f"Expected yes/no/unknown, received {value!r}")


def read_table(path, required):
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"Missing or duplicate TSV columns: {path}")
        absent = set(required) - set(reader.fieldnames)
        if absent:
            raise ValueError(f"{path}: missing columns {', '.join(sorted(absent))}")
        records = []
        for number, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{path}:{number}: malformed TSV row")
            records.append({key: value.strip() or "unknown" for key, value in row.items()})
    return records


def write_table(path, records, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for record in records:
            writer.writerow({key: record.get(key, "unknown") for key in fields})


def unique(records, field, label):
    result = {}
    for row in records:
        key = row[field]
        if not known(key) or key in result:
            raise ValueError(f"{label}: missing or duplicate {field}: {key!r}")
        result[key] = row
    return result


def file_state(path):
    if not known(path):
        return "NOT_PROVIDED"
    try:
        item = Path(path)
        if not item.is_file():
            return "UNAVAILABLE"
        return "PRESENT_NONEMPTY" if item.stat().st_size else "EMPTY"
    except OSError:
        return "UNAVAILABLE"


def inventory(args):
    samples = read_table(args.samples, SAMPLE_FIELDS)
    people = unique(samples, "person_id", "samples")
    reads = read_table(args.reads, READ_FIELDS)
    unique(reads, "read_id", "reads")
    resources = read_table(args.resources, RESOURCE_FIELDS)
    unique(resources, "resource_id", "resources")
    aliases = {}

    def add_alias(namespace, alias, person):
        if not known(alias):
            return
        if namespace not in IDENTIFIERS or person not in people:
            raise ValueError(f"Invalid crosswalk namespace/person: {namespace}, {person}")
        key = (namespace, alias)
        if key in aliases and aliases[key] != person:
            raise ValueError(f"Ambiguous identity: {namespace}={alias!r}")
        aliases[key] = person

    for person, sample in people.items():
        for namespace in IDENTIFIERS:
            add_alias(namespace, sample[namespace], person)
    if args.crosswalk:
        crosswalk = read_table(args.crosswalk, ["namespace", "alias", "person_id"])
        seen_aliases = set()
        for row in crosswalk:
            key = (row["namespace"], row["alias"])
            if not known(row["alias"]) or key in seen_aliases:
                raise ValueError(f"Missing or duplicate crosswalk key: {key}")
            seen_aliases.add(key)
            add_alias(*key, row["person_id"])

    missing, exclusions = [], []

    def issue(kind, identifier, reason, person="unknown", path="unknown"):
        missing.append(dict(kind=kind, identifier=identifier, person_id=person,
                            path=path, status="UNRESOLVED", reason=reason))

    for name, records, path in (("samples", samples, args.samples), ("reads", reads, args.reads),
                                ("resources", resources, args.resources)):
        if not records:
            issue("manifest", name, "header_only_manifest;cohort_or_resource_inventory_not_established", path=path)

    by_person = defaultdict(list)
    lane_keys, paths = set(), set()
    for row in reads:
        matches = set()
        if known(row["person_id"]):
            if row["person_id"] not in people:
                raise ValueError(f"Read {row['read_id']}: person_id absent from curated samples")
            matches.add(row["person_id"])
        unresolved = []
        for namespace in ("illumina_id", "fastq_manifest_id"):
            if known(row[namespace]):
                match = aliases.get((namespace, row[namespace]))
                if match:
                    matches.add(match)
                else:
                    unresolved.append(namespace)
        if len(matches) > 1:
            raise ValueError(f"Read {row['read_id']}: conflicting person/alias assignments")
        person = next(iter(matches), "unknown")
        row["person_id"] = person
        reasons = [f"unresolved_{field}" for field in unresolved]
        if not known(person):
            reasons.append("person_identity_unresolved")
        for field in ("library_id", "lane"):
            if not known(row[field]):
                reasons.append(f"missing_{field}")
        if known(person) and known(row["library_id"]) and known(row["lane"]):
            key = (person, row["library_id"], row["lane"])
            if key in lane_keys:
                raise ValueError(f"Duplicate person/library/lane key: {key}")
            lane_keys.add(key)
        for mate in ("r1", "r2"):
            row[f"{mate}_status"] = file_state(row[mate])
            if known(row[mate]):
                # Paths must be explicit and consistently spelled; no suffix/path inference.
                path = str(Path(row[mate]))
                if path in paths:
                    raise ValueError(f"FASTQ path reused in multiple mate/lane records: {path}")
                paths.add(path)
            if row[f"{mate}_status"] != "PRESENT_NONEMPTY":
                reasons.append(f"{mate}_{row[f'{mate}_status'].lower()}")
                issue("read", row["read_id"], reasons[-1], person, row[mate])
        row["status"] = "READY" if not reasons else "UNRESOLVED"
        row["reason"] = ";".join(reasons) or "explicit_lane_and_nonempty_mate_files;contents_not_scanned"
        if reasons:
            exclusions.append(dict(entity="read", identifier=row["read_id"], person_id=person,
                                   scope="measurement", status="EXCLUDED_OR_UNRESOLVED", reason=row["reason"]))
        if known(person):
            by_person[person].append(row)
        elif reasons:
            issue("identity", row["read_id"], row["reason"])

    pairing = []
    for person, sample in people.items():
        lanes = by_person[person]
        available = bool(lanes) and all(row["status"] == "READY" for row in lanes)
        phenotype = known(sample["phenotype_id"]) and known(sample["phenotype"])
        covariates = (flag(sample.get("covariates_available", "unknown")) == "yes"
                      and all(known(sample[field]) for field in ("sex", "age_definition", "ancestry_covariates")))
        snps = flag(sample.get("snp_data_available", "unknown"))
        manifest_complete = flag(sample.get("fastq_manifest_complete", "unknown"))
        unrelated = flag(sample.get("unrelated_primary", "unknown"))
        donor = flag(sample["discovery_donor"])
        relative = flag(sample["related_to_discovery_donor"])
        requested = flag(sample["analysis_eligible"])
        reasons = []
        if not available:
            reasons.append("illumina_lanes_missing_or_unresolved")
        if manifest_complete != "yes":
            reasons.append(f"fastq_manifest_completeness_{manifest_complete}")
        if not phenotype:
            reasons.append("phenotype_link_incomplete")
        if not covariates:
            reasons.append("covariates_unavailable_or_unknown")
        if snps != "yes":
            reasons.append(f"snp_data_{snps}")
        if donor != "no":
            reasons.append(f"discovery_donor_{donor}")
        if relative != "no":
            reasons.append(f"related_to_discovery_donor_{relative}")
        if unrelated != "yes":
            reasons.append(f"unrelated_primary_{unrelated}")
        if requested != "yes":
            reasons.append(f"primary_eligibility_not_approved_{requested}")
        if known(sample["exclusion_reason"]):
            reasons.append(sample["exclusion_reason"])
        try:
            if known(sample["depth"]) and (not math.isfinite(float(sample["depth"])) or float(sample["depth"]) < 0):
                raise ValueError("invalid depth")
        except ValueError:
            raise ValueError(f"{person}: depth must be nonnegative numeric or unknown") from None
        sample.update(requested_analysis_eligible=requested, analysis_eligible="yes" if not reasons else "no",
                      exclusion_reason=";".join(reasons) or "none", lane_count=len(lanes),
                      reads_available="yes" if available else "no", phenotype_linked="yes" if phenotype else "no",
                      covariates_available="yes" if covariates else "unknown_or_unavailable",
                      snp_data_available=snps, discovery_donor=donor, related_to_discovery_donor=relative,
                      fastq_manifest_complete=manifest_complete, unrelated_primary=unrelated,
                      analysis_group="independent_primary" if not reasons else
                      ("development_or_sensitivity" if donor == "yes" or relative == "yes" else "excluded_or_unresolved"),
                      status="READY" if not reasons else "EXCLUDED_OR_UNRESOLVED", reason=";".join(reasons) or "curated_primary_requirements_met")
        for reason in reasons:
            exclusions.append(dict(entity="person", identifier=person, person_id=person,
                                   scope="primary_phenotype", status="EXCLUDED_OR_UNRESOLVED", reason=reason))
        if not lanes:
            issue("read_manifest", person, "no_explicit_lane_records", person)
        for field in ("phenotype_id", "sex", "age_definition", "ancestry_covariates"):
            if not known(sample[field]):
                issue("metadata", field, f"{field}_unknown", person)
        for field in ("snp_data_available", "fastq_manifest_complete", "unrelated_primary",
                      "discovery_donor", "related_to_discovery_donor"):
            if sample[field] == "unknown":
                issue("metadata", field, f"{field}_unknown", person)
        fingerprint = sample.get("ont_fingerprint_status", "not_assessed").lower()
        evidence = sample.get("ont_fingerprint_evidence", "unknown")
        if fingerprint not in {"confirmed", "mismatch", "not_assessed", "unknown"}:
            raise ValueError(f"{person}: unsupported ont_fingerprint_status: {fingerprint}")
        ont_status = file_state(sample.get("ont_path", "unknown"))
        linked = known(sample["ont_id"]) and available
        confirmed = linked and fingerprint == "confirmed" and file_state(evidence) == "PRESENT_NONEMPTY"
        pair_status = "DECLARED_FINGERPRINT_CONFIRMED" if confirmed else ("MISMATCH" if fingerprint == "mismatch" else "ID_LINK_ONLY" if linked else "NOT_PAIRED")
        sample["ont_pairing_status"] = pair_status
        pairing.append(dict(person_id=person, illumina_id=sample["illumina_id"], ont_id=sample["ont_id"],
                            ont_path=sample.get("ont_path", "unknown"), ont_path_status=ont_status,
                            fingerprint_status=fingerprint, fingerprint_evidence=evidence, status=pair_status,
                            reason="curated_fingerprint_result;not_recomputed_by_inventory" if confirmed else
                            "ID_crosswalk_is_not_fingerprint_validation"))
        if known(sample["ont_id"]) and (ont_status != "PRESENT_NONEMPTY" or not confirmed):
            issue("ont", sample["ont_id"], "ONT_path_or_SNP_fingerprint_confirmation_missing", person)

    for row in resources:
        row["path_status"] = file_state(row["path"])
        checksum = row["sha256"].lower()
        checksum_valid = len(checksum) == 64 and all(char in "0123456789abcdef" for char in checksum)
        verified = flag(row["checksum_verified"])
        row["checksum_status"] = ("DECLARED_VERIFIED_NOT_REHASHED" if checksum_valid and verified == "yes" and known(row["checksum_source"])
                                  else "RECORDED_NOT_VERIFIED" if checksum_valid else "MISSING_OR_INVALID")
        reasons = []
        if row["path_status"] != "PRESENT_NONEMPTY":
            reasons.append(f"resource_{row['path_status'].lower()}")
        for field in ("build", "release", "scope"):
            if not known(row[field]):
                reasons.append(f"missing_{field}")
        if row["checksum_status"] != "DECLARED_VERIFIED_NOT_REHASHED":
            reasons.append("checksum_verification_not_documented")
        row["status"] = "RECORDED" if not reasons else "INCOMPLETE"
        row["reason"] = ";".join(reasons) or "frozen_metadata_and_existing_checksum_recorded;resource_not_rehashed"
        for reason in reasons:
            issue("resource", row["resource_id"], reason, path=row["path"])

    flow = []
    remaining = samples
    checks = [("reads_available", lambda r: r["reads_available"] == "yes"),
              ("phenotype_linked", lambda r: r["phenotype_linked"] == "yes"),
              ("covariates_available", lambda r: r["covariates_available"] == "yes"),
              ("snp_data_available", lambda r: r["snp_data_available"] == "yes"),
              ("independent_unrelated_non_donor", lambda r: r["discovery_donor"] == r["related_to_discovery_donor"] == "no" and r["unrelated_primary"] == "yes"),
              ("final_primary", lambda r: r["analysis_eligible"] == "yes")]
    flow.append(dict(step="curated_people", count=len(samples), denominator=len(samples),
                     status="INVENTORIED" if samples else "MISSING_INPUTS",
                     reason="one_row_per_explicit_person_id" if samples else "manifest_not_populated;not_an_observed_empty_cohort"))
    for name, check in checks:
        before = len(remaining)
        remaining = [row for row in remaining if check(row)]
        flow.append(dict(step=name, count=len(remaining), denominator=before,
                         status="CUMULATIVE" if samples else "MISSING_INPUTS",
                         reason="unknown_values_do_not_pass" if samples else "manifest_not_populated"))
    for name, check in [("graph_donors", lambda r: r["discovery_donor"] == "yes"),
                        ("discovery_relatives", lambda r: r["related_to_discovery_donor"] == "yes"),
                        ("ont_id_linked", lambda r: r["ont_pairing_status"] in {"ID_LINK_ONLY", "DECLARED_FINGERPRINT_CONFIRMED"}),
                        ("ont_fingerprint_confirmed", lambda r: r["ont_pairing_status"] == "DECLARED_FINGERPRINT_CONFIRMED")]:
        flow.append(dict(step=name, count=sum(check(row) for row in samples), denominator=len(samples),
                         status="SEPARATE_COUNT" if samples else "MISSING_INPUTS",
                         reason="not_a_cumulative_gate" if samples else "manifest_not_populated"))
    groups = Counter((r["phenotype"], r["batch"], r["depth"], r["analysis_group"]) for r in samples)
    summary = [dict(phenotype=key[0], batch=key[1], depth=key[2], analysis_group=key[3], n_people=value,
                    status="DESCRIPTIVE", reason="exact_curated_depth;unknown_retained") for key, value in sorted(groups.items())]
    out = Path(args.output) / "metadata"
    computed = "requested_analysis_eligible lane_count reads_available phenotype_linked covariates_available snp_data_available fastq_manifest_complete unrelated_primary analysis_group ont_pairing_status status reason".split()
    sample_fields = list(dict.fromkeys([*SAMPLE_FIELDS, *computed, *(key for row in samples for key in row)]))
    write_table(out / "samples.tsv", samples, sample_fields)
    write_table(out / "reads.tsv", reads, [*READ_FIELDS, "r1_status", "r2_status", "status", "reason"])
    write_table(out / "resources.tsv", resources, [*RESOURCE_FIELDS, "path_status", "checksum_status", "status", "reason"])
    write_table(out / "sample_flow.tsv", flow, ["step", "count", "denominator", "status", "reason"])
    write_table(out / "phenotype_by_batch_depth.tsv", summary, ["phenotype", "batch", "depth", "analysis_group", "n_people", "status", "reason"])
    write_table(out / "donor_overlap.tsv", samples, ["person_id", "discovery_donor", "related_to_discovery_donor", "analysis_group", "status", "reason"])
    write_table(out / "ont_pairing.tsv", pairing, ["person_id", "illumina_id", "ont_id", "ont_path", "ont_path_status", "fingerprint_status", "fingerprint_evidence", "status", "reason"])
    write_table(out / "missing_inputs.tsv", missing, ["kind", "identifier", "person_id", "path", "status", "reason"])
    write_table(out / "exclusions.tsv", exclusions, ["entity", "identifier", "person_id", "scope", "status", "reason"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("samples", "reads", "resources", "output"):
        parser.add_argument(f"--{name}", required=True)
    parser.add_argument("--crosswalk", help="Optional TSV: namespace, alias, person_id")
    args = parser.parse_args()
    try:
        inventory(args)
    except (OSError, ValueError) as error:
        parser.exit(1, f"Inventory error: {error}\n")


if __name__ == "__main__":
    main()
