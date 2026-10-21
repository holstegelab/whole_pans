#!/usr/bin/env python3
"""Import prespecified existing ExpansionHunter repeats; never infer SV dosage.

Prepare selected catalog metadata once, then let Snakemake stream one VCF per job.
Inputs are unchanged. Caller PASS remains unvalidated external-call evidence.
"""

import argparse
import gzip
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

from svgp_inventory import known, protect_inputs, read_table, write_table


CATALOG_FIELDS = "locus_id variant_id chromosome start0 end0 variant_type locus_structure catalog catalog_sha256 catalog_status".split()
EVIDENCE_FIELDS = "source_id source_vcf caller_sample_id person_id identity_status reference_build reference_fasta reference_sha256 caller_version caller_version_source locus_id variant_id chromosome start0 end0 variant_type locus_structure catalog catalog_sha256 catalog_status repeat_unit reference_repeat_units reference_length_bp vcf_pos1 vcf_end vcf_ref vcf_alt vcf_filter gt repcn repci so adsp adfl adir locus_coverage allele1_repeat_units allele1_ci_low allele1_ci_high allele2_repeat_units allele2_ci_low allele2_ci_high allele_count call_status review_status".split()


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else Path(path).open()


def prepare_catalog(args):
    protect_inputs([args.catalog, args.loci], [args.output])
    selected = read_table(args.loci, ["locus_id", "variant_id"])
    wanted = {}
    for row in selected:
        if not known(row["locus_id"]) or not known(row["variant_id"]) or row["variant_id"] in wanted:
            raise ValueError("Selected repeats require unique, explicit locus_id and variant_id")
        wanted[row["variant_id"]] = row
    if not wanted:
        raise ValueError("Choose a small prespecified repeat subset in --loci before extraction")
    found = {}
    with Path(args.catalog).open("rb") as handle:
        raw = handle.read()
    checksum = hashlib.sha256(raw).hexdigest()
    for entry in json.loads(raw):
        regions = entry["ReferenceRegion"]
        regions = regions if isinstance(regions, list) else [regions]
        ids = entry.get("VariantId")
        if ids is None:
            ids = [entry["LocusId"]] if len(regions) == 1 else [f"{entry['LocusId']}_{region}" for region in regions]
        ids = ids if isinstance(ids, list) else [ids]
        types = entry["VariantType"]
        types = types if isinstance(types, list) else [types]
        if len(ids) != len(regions) or len(types) != len(regions):
            raise ValueError(f"Inconsistent catalog component arrays: {entry['LocusId']}")
        for variant, region, kind in zip(ids, regions, types):
            if variant not in wanted:
                continue
            if variant in found or wanted[variant]["locus_id"] != entry["LocusId"]:
                raise ValueError(f"Duplicate variant or selected locus/catalog mismatch: {variant}")
            match = re.fullmatch(r"(.+):(\d+)-(\d+)", region)
            if not match or int(match[2]) >= int(match[3]):
                raise ValueError(f"Invalid catalog reference interval: {region}")
            if kind not in {"Repeat", "RareRepeat"}:
                raise ValueError(f"Selected variant is not a repeat: {variant}")
            found[variant] = dict(wanted[variant], chromosome=match[1], start0=match[2], end0=match[3],
                                  variant_type=kind, locus_structure=entry["LocusStructure"], catalog_status="FOUND")
    records = []
    for variant, row in wanted.items():
        records.append(dict(found.get(variant, dict(row, catalog_status="NOT_IN_CATALOG")),
                            catalog=args.catalog, catalog_sha256=checksum))
    write_table(Path(args.output), records, CATALOG_FIELDS)


def parse_call(fields, catalog):
    """Retain raw fields and partial calls; repeat counts are per allele, not dosage."""
    info = dict(part.split("=", 1) for part in fields[7].split(";") if "=" in part)
    keys, values = fields[8].split(":"), fields[9].split(":")
    if len(keys) != len(set(keys)) or len(values) > len(keys):
        raise ValueError("Malformed FORMAT fields")
    data = dict(zip(keys, values))
    result = dict(repeat_unit=info.get("RU", "unknown"), reference_repeat_units=info.get("REF", "unknown"),
                  reference_length_bp=info.get("RL", "unknown"), vcf_pos1=fields[1], vcf_end=info.get("END", "unknown"),
                  vcf_ref=fields[3], vcf_alt=fields[4], vcf_filter=fields[6],
                  gt=data.get("GT", "."), repcn=data.get("REPCN", "."), repci=data.get("REPCI", "."),
                  so=data.get("SO", "."), adsp=data.get("ADSP", "."),
                  adfl=data.get("ADFL", "."), adir=data.get("ADIR", "."),
                  locus_coverage=data.get("LC", "unknown"))
    # EH repeat POS is the 1-based preceding anchor: repeat start0 equals POS.
    if catalog["catalog_status"] != "FOUND":
        return dict(result, call_status="NOT_IN_CATALOG")
    if (fields[0], fields[1], info.get("END")) != (catalog["chromosome"], catalog["start0"], catalog["end0"]):
        return dict(result, call_status="CATALOG_COORDINATE_MISMATCH")
    gt = re.split(r"[/|]", result["gt"])
    counts = result["repcn"].split("/")
    cis = result["repci"].split("/")
    result["allele_count"] = len(gt)
    complete = len(gt) in {1, 2} and len(counts) == len(gt) and len(cis) == len(gt)
    any_called = False
    for index in range(min(2, len(gt))):
        count = counts[index] if index < len(counts) else "."
        ci = re.fullmatch(r"(\d+)-(\d+)", cis[index]) if index < len(cis) else None
        if count.isdigit():
            result[f"allele{index+1}_repeat_units"] = count
            any_called = True
        if ci:
            result[f"allele{index+1}_ci_low"], result[f"allele{index+1}_ci_high"] = ci.groups()
        if (not gt[index].isdigit() or not count.isdigit() or not ci
                or (ci and not int(ci[1]) <= int(count) <= int(ci[2]))):
            complete = False
    if "LowDepth" in fields[6].split(";"):
        status = "LOW_DEPTH"
    elif not complete:
        status = "PARTIAL_GENOTYPE" if any_called else "NO_CALL"
    elif fields[6] != "PASS":
        status = "FILTERED_OR_UNASSESSED"
    else:
        status = "CALLER_PASS"
    return dict(result, call_status=status)


def extract(args):
    protect_inputs([args.vcf, args.catalog, args.crosswalk, args.run_log],
                   [args.out_evidence, args.out_status])
    catalog_rows = read_table(args.catalog, CATALOG_FIELDS)
    selected = {row["variant_id"]: row for row in catalog_rows}
    if not selected or len(selected) != len(catalog_rows):
        raise ValueError("Prepared selected catalog must contain unique variants and at least one repeat")
    aliases = {}
    for row in read_table(args.crosswalk, ["source_id", "caller_sample_id", "person_id", "review_status"]):
        if row["source_id"] != args.source_id:
            continue
        alias = row["caller_sample_id"]
        if not known(alias) or alias in aliases:
            raise ValueError(f"Missing or duplicate ExpansionHunter crosswalk label: {alias}")
        aliases[alias] = row
    version, version_source = "unknown", "unknown"
    if args.run_log and Path(args.run_log).is_file():
        with Path(args.run_log).open() as handle:
            for _ in range(20):
                line = handle.readline()
                match = re.search(r"Starting ExpansionHunter v([^\]\s]+)", line)
                if match:
                    version, version_source = match[1], args.run_log
                    break
    base = dict(source_id=args.source_id, source_vcf=args.vcf, caller_sample_id="unknown",
                person_id="unknown", identity_status="UNRESOLVED", reference_build=args.reference_build,
                reference_fasta=args.reference_fasta, reference_sha256=args.reference_sha256,
                caller_version=version, caller_version_source=version_source, review_status="NOT_VALIDATED")
    calls, total, file_status = {}, 0, "VCF_UNAVAILABLE"
    if Path(args.vcf).is_file() and Path(args.vcf).stat().st_size:
        with open_text(args.vcf) as handle:
            for line in handle:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM\t"):
                    if file_status != "VCF_UNAVAILABLE":
                        raise ValueError("Repeated ExpansionHunter VCF sample header")
                    header = line.rstrip("\r\n").split("\t")
                    if (len(header) != 10 or not known(header[9])
                            or header[:9] != ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]):
                        raise ValueError("Expected one-sample ExpansionHunter VCF")
                    base["caller_sample_id"] = header[9]
                    alias = aliases.get(header[9], {})
                    if alias.get("review_status") == "REVIEWED" and known(alias.get("person_id", "unknown")):
                        base.update(person_id=alias["person_id"], identity_status="LINKED")
                    file_status = "HEADER_ONLY_OR_TRUNCATED"
                    continue
                if line.startswith("#") or not line.strip():
                    continue
                if file_status == "VCF_UNAVAILABLE":
                    raise ValueError("VCF body precedes sample header")
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != 10:
                    raise ValueError("Malformed ExpansionHunter VCF row")
                total += 1
                info = dict(part.split("=", 1) for part in fields[7].split(";") if "=" in part)
                variant = info.get("VARID", info.get("REPID", ""))
                if variant not in selected:
                    continue
                if variant in calls:
                    raise ValueError(f"Duplicate selected repeat in VCF: {variant}")
                calls[variant] = parse_call(fields, selected[variant])
        if total:
            file_status = "READ_TO_EOF_COMPLETENESS_UNVERIFIED"
        elif file_status == "VCF_UNAVAILABLE":
            file_status = "NO_SAMPLE_HEADER"
    records = []
    for variant, metadata in selected.items():
        missing = "NOT_IN_CATALOG" if metadata["catalog_status"] != "FOUND" else "RECORD_NOT_REPORTED"
        if file_status in {"VCF_UNAVAILABLE", "NO_SAMPLE_HEADER", "HEADER_ONLY_OR_TRUNCATED"}:
            missing = file_status
        records.append(dict(base, **metadata, **calls.get(variant, {"call_status": missing})))
    write_table(Path(args.out_evidence), records, EVIDENCE_FIELDS)
    counts = Counter(row["call_status"] for row in records)
    summary = dict(base, file_status=file_status, records_scanned=total, selected_variants=len(selected),
                   reported_selected_variants=len(calls), call_status_counts=json.dumps(dict(sorted(counts.items()))))
    write_table(Path(args.out_status), [summary], list(summary))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    catalog = commands.add_parser("catalog", help="Prepare a small selected catalog once")
    for name in ("catalog", "loci", "output"):
        catalog.add_argument(f"--{name}", required=True)
    catalog.set_defaults(func=prepare_catalog)
    call = commands.add_parser("extract", help="Stream one existing VCF; no new genotyping")
    for name in ("vcf", "catalog", "crosswalk", "out-evidence", "out-status"):
        call.add_argument(f"--{name}", required=True)
    call.add_argument("--source-id", default="expansionhunter")
    call.add_argument("--run-log", default="")
    for name in ("reference-build", "reference-fasta", "reference-sha256"):
        call.add_argument(f"--{name}", default="unknown")
    call.set_defaults(func=extract)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
