#!/usr/bin/env python3
"""Import one existing single-sample gCNV segment VCF without calling or filtering."""

import argparse
import csv
import gzip
import hashlib
import re
from collections import Counter
from pathlib import Path

from svgp_inventory import protect_inputs


EVIDENCE_FIELDS = "source_id caller_sample_id person_id identity_status source_path source_record_id record_number chromosome start0 end0 svtype reported_cn raw_gt np qa qs qss qse qual filter reported_state reference_build caller_version reference_path".split()
STATUS_FIELDS = "source_id caller_sample_id person_id identity_status source_path file_status import_status records reference_segments cnv_segments missing_cn_segments unsupported_segments sha256 header_sha256 reference_build caller_version reference_path autosomal_reference_cn allosomal_contigs source_role interpretation".split()


def writer(path, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = path.open("w", newline="")
    table = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    table.writeheader()
    return handle, table


def identity(path, source_id, sample):
    mappings = {}
    with Path(path).open(newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        required = {"source_id", "caller_sample_id", "person_id", "review_status"}
        if not rows.fieldnames or len(rows.fieldnames) != len(set(rows.fieldnames)) or not required <= set(rows.fieldnames):
            raise ValueError(f"{path}: expected unique columns {sorted(required)}")
        for number, row in enumerate(rows, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{path}:{number}: malformed crosswalk row")
            row = {key: value.strip() for key, value in row.items()}
            key = row["source_id"], row["caller_sample_id"]
            if not all(key) or key in mappings:
                raise ValueError(f"{path}:{number}: missing or duplicate source/sample key")
            mappings[key] = row
    row = mappings.get((source_id, sample))
    if row and row["review_status"] == "REVIEWED":
        if row["person_id"].lower() in {"", ".", "unknown", "na", "n/a", "none"}:
            raise ValueError(f"Reviewed crosswalk lacks person_id for {source_id}/{sample}")
        return row["person_id"], "REVIEWED"
    return "unknown", "UNRESOLVED"


def header_value(line, pattern):
    match = re.search(pattern, line)
    return match.group(1) if match else "unknown"


def import_vcf(args):
    source = Path(args.vcf)
    protect_inputs([source, args.crosswalk], [args.out_evidence, args.out_status])
    status = dict.fromkeys(STATUS_FIELDS, "unknown")
    status.update(source_id=args.source_id, source_path=str(source), reference_build=args.reference_build,
                  caller_sample_id=args.expected_sample or "unknown", identity_status="UNRESOLVED",
                  source_role="existing_linear_cnv", import_status="NO_RECORDS",
                  interpretation="raw_segments_only;no_graph_allele_equivalence;missing_is_not_reference")
    counts = Counter()
    evidence_handle, evidence = writer(args.out_evidence, EVIDENCE_FIELDS)
    try:
        if not source.is_file():
            status.update(file_status="UNAVAILABLE", import_status="NOT_IMPORTED")
        elif not source.stat().st_size:
            status.update(file_status="EMPTY", import_status="NOT_IMPORTED")
        else:
            status["file_status"] = "PRESENT_NONEMPTY"
            digest = hashlib.sha256()
            with source.open("rb") as compressed:
                for block in iter(lambda: compressed.read(1024 * 1024), b""):
                    digest.update(block)
            status["sha256"] = digest.hexdigest()
            header_digest = hashlib.sha256()
            opener = gzip.open if str(source).endswith(".gz") else open
            sample = None
            is_segments = False
            with opener(source, "rt") as handle:
                for number, line in enumerate(handle, 1):
                    if line.startswith("#"):
                        header_digest.update(line.encode("utf-8"))
                        if line.startswith("##FORMAT=<ID=NP,"):
                            is_segments = True
                        if line.startswith("##reference="):
                            status["reference_path"] = line.rstrip().split("=", 1)[1]
                        if line.startswith("##GATKCommandLine=<ID=PostprocessGermlineCNVCalls,"):
                            status["caller_version"] = header_value(line, r'Version="([^\"]+)"')
                            reference = header_value(line, r'--reference\s+([^\s"]+)')
                            if reference != "unknown":
                                status["reference_path"] = reference
                            status["autosomal_reference_cn"] = header_value(line, r'--autosomal-ref-copy-number\s+(\d+)')
                            allosomes = re.findall(r'--allosomal-contig\s+([^\s"]+)', line)
                            status["allosomal_contigs"] = ",".join(allosomes) or "not_declared"
                        if line.startswith("#CHROM\t"):
                            fields = line.rstrip("\r\n").split("\t")
                            if (sample is not None or len(fields) != 10 or not is_segments
                                    or fields[:9] != ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
                                    or not fields[9].strip()):
                                raise ValueError(f"{source}: expected one sample and segment FORMAT/NP")
                            sample = fields[9]
                            if args.expected_sample and sample != args.expected_sample:
                                raise ValueError(f"{source}: sample differs from expected header sample")
                            person, identity_status = identity(args.crosswalk, args.source_id, sample)
                            status.update(caller_sample_id=sample, person_id=person, identity_status=identity_status)
                        continue
                    if sample is None:
                        raise ValueError(f"{source}:{number}: missing #CHROM header")
                    fields = line.rstrip("\r\n").split("\t")
                    if len(fields) != 10:
                        raise ValueError(f"{source}:{number}: expected 10 VCF fields")
                    chrom, pos, record_id, ref, alt, qual, filt, raw_info, raw_format, call = fields
                    info = dict(item.split("=", 1) for item in raw_info.split(";") if "=" in item)
                    keys, values = raw_format.split(":"), call.split(":")
                    if len(keys) != len(set(keys)) or len(values) > len(keys) or not {"CN", "NP"} <= set(keys):
                        raise ValueError(f"{source}:{number}: invalid segment FORMAT")
                    fmt = dict(zip(keys, values))
                    start, end = int(pos) - 1, int(info["END"])
                    if start < 0 or end <= start:
                        raise ValueError(f"{source}:{number}: invalid POS/END")
                    cn, gt = fmt.get("CN", "."), fmt.get("GT", ".")
                    if cn != "." and (not cn.isdigit() or int(cn) < 0):
                        raise ValueError(f"{source}:{number}: CN must be a nonnegative integer or '.'")
                    # A duplication may have GT=./. while CN is still reported by gCNV.
                    svtype = {"<DEL>": "DEL", "<DUP>": "DUP"}.get(alt, "unknown")
                    if cn == ".":
                        state = "MISSING_CN"
                    elif svtype in {"DEL", "DUP"}:
                        state = "REPORTED_CNV"
                    elif alt == "." and gt in {"0", "0/0", "0|0"}:
                        svtype, state = "REF", "REPORTED_REFERENCE_SEGMENT"
                    else:
                        state = "UNSUPPORTED_ALLELE_REPRESENTATION"
                    counts["records"] += 1
                    counts[state] += 1
                    row = dict(status, source_record_id=record_id, record_number=counts["records"],
                               chromosome=chrom, start0=start, end0=end, svtype=svtype, reported_cn=cn,
                               raw_gt=gt, qual=qual, filter=filt, reported_state=state)
                    row.update({key.lower(): fmt.get(key, ".") for key in ("NP", "QA", "QS", "QSS", "QSE")})
                    evidence.writerow({key: row.get(key, "unknown") for key in EVIDENCE_FIELDS})
            if sample is None:
                raise ValueError(f"{source}: missing #CHROM header")
            status["header_sha256"] = header_digest.hexdigest()
            if counts["records"]:
                status["import_status"] = "IMPORTED_RAW_EVIDENCE"
    finally:
        evidence_handle.close()
    status.update(records=counts["records"], reference_segments=counts["REPORTED_REFERENCE_SEGMENT"],
                  cnv_segments=counts["REPORTED_CNV"], missing_cn_segments=counts["MISSING_CN"],
                  unsupported_segments=counts["UNSUPPORTED_ALLELE_REPRESENTATION"])
    status_handle, table = writer(args.out_status, STATUS_FIELDS)
    try:
        table.writerow(status)
    finally:
        status_handle.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--source-id", required=True)
    parser.add_argument("--crosswalk", required=True)
    parser.add_argument("--expected-sample", default="")
    parser.add_argument("--reference-build", default="unknown")
    parser.add_argument("--out-evidence", required=True)
    parser.add_argument("--out-status", required=True)
    args = parser.parse_args()
    try:
        import_vcf(args)
    except (OSError, ValueError, KeyError, EOFError) as error:
        parser.exit(1, f"gCNV import failed: {error}\n")


if __name__ == "__main__":
    main()
