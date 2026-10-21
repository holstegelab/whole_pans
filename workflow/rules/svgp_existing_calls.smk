"""Inspect existing callers; import one reviewed callset without rerunning callers."""

import gzip
import hashlib

SGP_EXISTING = SGP.get("existing_calls", {})
SGP_GCNV = SGP_EXISTING.get("gcnv", {})
SGP_EH = SGP_EXISTING.get("expansionhunter", {})
SGP_EXISTING_CROSSWALK = SGP_EXISTING.get("crosswalk", "config/svgp_existing_call_crosswalk.example.tsv")
SGP_EXISTING_OUT = SGP_OUT + "/existing_calls"
SGP_CNV_ROOT = Path(SGP_GCNV.get("root", "."))
SGP_EH_ROOT = Path(SGP_EH.get("root", "."))
SGP_CALL_FOLDERS = {
    "gcnv_standard": [SGP_CNV_ROOT / "GATK_gCNV" / ("GENOTYPED_CALLS_segments_" + str(i)) for i in (0, 1)],
    "gcnv_wgs_autosomes": [SGP_CNV_ROOT / "GATK_gCNV/wgs_autosomes_all_samples/GENOTYPED_CALLS_segments_0"],
    "expansionhunter": [SGP_EH_ROOT / "results"],
}


def sgp_existing_files(wildcards):
    return sorted({str(path) for folders in SGP_CALL_FOLDERS.values() for folder in folders
                   for pattern in ("*.vcf", "*.vcf.gz") for path in folder.glob(pattern)})


def sgp_existing_records():
    manifest = checkpoints.svgp_existing_call_inventory.get().output.files
    with open(manifest, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sgp_selected_calls(kind):
    source = SGP_GCNV.get("selected_callset", "") if kind == "gcnv" else "expansionhunter"
    if source not in SGP_CALL_FOLDERS or (kind == "gcnv" and not source.startswith("gcnv_")):
        raise ValueError("Choose one gcnv selected_callset: gcnv_standard or gcnv_wgs_autosomes; do not pool overlapping runs")
    records = [row for row in sgp_existing_records() if row["source_id"] == source]
    if not records:
        raise ValueError(f"No existing VCF files inventoried for {source}; inspect source_summary.tsv")
    return records


def sgp_existing_row(wildcards, kind):
    matches = [row for row in sgp_selected_calls(kind) if row["file_id"] == wildcards.file_id]
    if len(matches) != 1:
        raise ValueError(f"Unknown or unselected existing {kind} source file: {wildcards.file_id}")
    return matches[0]


def sgp_existing_parts(kind, part):
    return [f"{SGP_EXISTING_OUT}/{kind}/{row['file_id']}.{part}.tsv" for row in sgp_selected_calls(kind)]


rule existing_calls_inventory:
    input:
        SGP_EXISTING_OUT + "/source_files.tsv",
        SGP_EXISTING_OUT + "/source_summary.tsv"


checkpoint svgp_existing_call_inventory:
    input:
        vcfs=sgp_existing_files
    output:
        files=SGP_EXISTING_OUT + "/source_files.tsv",
        summary=SGP_EXISTING_OUT + "/source_summary.tsv"
    params:
        source_folders={source: list(map(str, folders)) for source, folders in SGP_CALL_FOLDERS.items()}
    threads: 1
    resources:
        mem_mb=1000,
        runtime_min=30
    log: SGP_OUT + "/logs/existing_calls_inventory.log"
    run:
        records, summaries = [], []
        for source, folder_names in params.source_folders.items():
            folders = list(map(Path, folder_names))
            paths = sorted({path for folder in folders for pattern in ("*.vcf", "*.vcf.gz") for path in folder.glob(pattern)})
            samples = []
            for path in paths:
                file_id = hashlib.sha256((source + "\0" + str(path)).encode()).hexdigest()[:20]
                record = dict(file_id=file_id, source_id=source, path=str(path), caller_sample_id="unknown",
                              status="HEADER_NOT_READ", reason="unknown")
                opener = gzip.open if str(path).endswith(".gz") else open
                try:
                    with opener(path, "rt") as handle:
                        for line in handle:
                            if line.startswith("#CHROM\t"):
                                fields = line.rstrip("\r\n").split("\t")
                                if (len(fields) != 10 or not fields[9].strip()
                                        or fields[:9] != ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]):
                                    raise ValueError("Expected a single-sample VCF")
                                record.update(caller_sample_id=fields[9], status="HEADER_AVAILABLE",
                                              reason="literal_VCF_sample_only;identity_and_call_QC_unreviewed")
                                samples.append(fields[9])
                                break
                            if not line.startswith("#"):
                                break
                    if record["status"] == "HEADER_NOT_READ":
                        record.update(status="HEADER_UNAVAILABLE", reason="no_single_sample_CHROM_header")
                except (OSError, ValueError, EOFError) as error:
                    record.update(status="HEADER_UNAVAILABLE", reason=str(error))
                records.append(record)
            summaries.append(dict(source_id=source, source_folders=";".join(map(str, folders)),
                                  files=len(paths), distinct_header_samples=len(set(samples)),
                                  repeated_sample_files=len(samples)-len(set(samples)),
                                  status="INVENTORIED_HEADERS_ONLY" if paths else "UNAVAILABLE_OR_EMPTY",
                                  reason="overlapping_callsets_are_not_independent;files_are_not_people"))
        for path, rows, fields in ((output.files, records, ["file_id", "source_id", "path", "caller_sample_id", "status", "reason"]),
                                   (output.summary, summaries, list(summaries[0]))):
            with open(path, "w", newline="") as handle:
                writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
        Path(log[0]).write_text("Read existing VCF headers only; did not rerun callers or read BAM/CRAM data.\n")


rule gcnv_import:
    input: SGP_EXISTING_OUT + "/gcnv/calls.tsv.gz", SGP_EXISTING_OUT + "/gcnv/import_status.tsv"


rule repeat_import:
    input: SGP_EXISTING_OUT + "/expansionhunter/calls.tsv.gz", SGP_EXISTING_OUT + "/expansionhunter/import_status.tsv"


rule existing_calls:
    input:
        SGP_EXISTING_OUT + "/gcnv/calls.tsv.gz",
        SGP_EXISTING_OUT + "/expansionhunter/calls.tsv.gz"


rule svgp_import_gcnv:
    input:
        script=str(SGP_ROOT / "workflow/scripts/svgp_gcnv.py"),
        common=str(SGP_ROOT / "workflow/scripts/svgp_inventory.py"),
        vcf=lambda wc: sgp_existing_row(wc, "gcnv")["path"],
        crosswalk=SGP_EXISTING_CROSSWALK
    output:
        evidence=SGP_EXISTING_OUT + "/gcnv/{file_id}.evidence.tsv",
        status=SGP_EXISTING_OUT + "/gcnv/{file_id}.status.tsv"
    params:
        source=lambda wc: sgp_existing_row(wc, "gcnv")["source_id"],
        sample=lambda wc: sgp_existing_row(wc, "gcnv")["caller_sample_id"],
        build=SGP_GCNV.get("reference_build", "unknown")
    threads: 1
    resources:
        mem_mb=1000,
        runtime_min=30,
        existing_call_jobs=1
    log: SGP_OUT + "/logs/existing_calls/gcnv/{file_id}.log"
    conda: SGP_TABLE_ENV
    shell:
        "python {input.script:q} --vcf {input.vcf:q} --source-id {params.source:q} "
        "--expected-sample {params.sample:q} --reference-build {params.build:q} "
        "--crosswalk {input.crosswalk:q} --out-evidence {output.evidence:q} "
        "--out-status {output.status:q} > {log:q} 2>&1"


rule svgp_prepare_repeat_catalog:
    input:
        script=str(SGP_ROOT / "workflow/scripts/svgp_expansionhunter.py"),
        common=str(SGP_ROOT / "workflow/scripts/svgp_inventory.py"),
        catalog=SGP_EH.get("catalog", str(SGP_EH_ROOT / "variant_catalog.masked_reference.json")),
        loci=SGP_EH.get("loci", "config/svgp_existing_repeat_loci.example.tsv")
    output: SGP_EXISTING_OUT + "/expansionhunter/selected_catalog.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime_min=15
    log: SGP_OUT + "/logs/existing_calls/repeat_catalog.log"
    conda: SGP_TABLE_ENV
    shell:
        "python {input.script:q} catalog --catalog {input.catalog:q} --loci {input.loci:q} "
        "--output {output:q} > {log:q} 2>&1"


rule svgp_import_expansionhunter:
    input:
        script=str(SGP_ROOT / "workflow/scripts/svgp_expansionhunter.py"),
        common=str(SGP_ROOT / "workflow/scripts/svgp_inventory.py"),
        vcf=lambda wc: sgp_existing_row(wc, "expansionhunter")["path"],
        catalog=SGP_EXISTING_OUT + "/expansionhunter/selected_catalog.tsv",
        crosswalk=SGP_EXISTING_CROSSWALK
    output:
        evidence=SGP_EXISTING_OUT + "/expansionhunter/{file_id}.evidence.tsv",
        status=SGP_EXISTING_OUT + "/expansionhunter/{file_id}.status.tsv"
    params:
        build=SGP_EH.get("reference_build", "unknown"),
        reference=SGP_EXISTING.get("reference_fasta", "unknown"),
        reference_sha256=SGP_EXISTING.get("reference_sha256") or "unknown"
    threads: 1
    resources:
        mem_mb=1000,
        runtime_min=30,
        existing_call_jobs=1
    log: SGP_OUT + "/logs/existing_calls/expansionhunter/{file_id}.log"
    conda: SGP_TABLE_ENV
    shell:
        "python {input.script:q} extract --vcf {input.vcf:q} --catalog {input.catalog:q} "
        "--crosswalk {input.crosswalk:q} --source-id expansionhunter --reference-build {params.build:q} "
        "--reference-fasta {params.reference:q} --reference-sha256 {params.reference_sha256:q} "
        "--out-evidence {output.evidence:q} --out-status {output.status:q} > {log:q} 2>&1"


rule svgp_collect_existing_calls:
    input:
        evidence=lambda wc: sgp_existing_parts(wc.caller, "evidence"),
        statuses=lambda wc: sgp_existing_parts(wc.caller, "status")
    output:
        calls=SGP_EXISTING_OUT + "/{caller}/calls.tsv.gz",
        status=SGP_EXISTING_OUT + "/{caller}/import_status.tsv"
    wildcard_constraints:
        caller="gcnv|expansionhunter"
    threads: 1
    resources:
        mem_mb=1000,
        runtime_min=60
    log: SGP_OUT + "/logs/existing_calls/collect_{caller}.log"
    run:
        for paths, output_path, compressed in ((input.evidence, output.calls, True), (input.statuses, output.status, False)):
            opener = gzip.open if compressed else open
            with opener(output_path, "wt", newline="") as destination:
                writer = None
                for path in paths:
                    with open(path, newline="") as source:
                        reader = csv.DictReader(source, delimiter="\t")
                        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
                            raise ValueError(f"Missing or duplicate existing-call columns: {path}")
                        if writer is None:
                            writer = csv.DictWriter(destination, reader.fieldnames, delimiter="\t", lineterminator="\n")
                            writer.writeheader()
                        if reader.fieldnames != writer.fieldnames:
                            raise ValueError("Existing-call schema mismatch across source files")
                        rows = 0
                        for row in reader:
                            if None in row or any(value is None for value in row.values()):
                                raise ValueError(f"Malformed existing-call row: {path}")
                            writer.writerow(row)
                            rows += 1
                        if not compressed and rows != 1:
                            raise ValueError(f"Expected one existing-call status row: {path}")
        Path(log[0]).write_text("Streamed existing evidence/status tables; no association or graph-allele equivalence inferred.\n")
