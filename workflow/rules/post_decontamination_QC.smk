POST_QC_CONFIG = workflow.source_path("../../config/config.yaml")
configfile: POST_QC_CONFIG

POST_QC_ENTRY_OUTDIR = config["results"]["post_decontamination_qc"]


rule post_decontamination_qc:
    input:
        f"{POST_QC_ENTRY_OUTDIR}/summary/assembly_qc.tsv",
        f"{POST_QC_ENTRY_OUTDIR}/summary/sample_qc.tsv",
        f"{POST_QC_ENTRY_OUTDIR}/summary/graph_included_assemblies.txt",
        f"{POST_QC_ENTRY_OUTDIR}/summary/graph_excluded_assemblies.tsv",
        complete=lambda wildcards: POST_QC_COMPLETE_MARKER,


include: "decontamination.smk"


import csv
from pathlib import Path


POST_QC_OUTDIR = config["results"]["post_decontamination_qc"]
POST_QC_COMPLETE_MARKER = (
    f"{POST_QC_OUTDIR}/summary/"
    f"all_cleaned_assemblies.{DECONTAM_ASSEMBLY_MANIFEST_SIGNATURE[:12]}.complete"
)


def cleaned_assembly_id(path):
    name = strip_fasta_suffix(path)
    if not name.endswith(".clean"):
        raise WorkflowError(
            f"Expected a cleaned FASTA name ending in '.clean': {path}"
        )
    return name[: -len(".clean")]


rule post_qc_assembly_manifest:
    input:
        graph_list=f"{DECONTAM_OUTDIR}/summary/graph_cleaned_assemblies.txt"
    output:
        f"{POST_QC_OUTDIR}/resources/all_cleaned_assemblies.tsv"
    run:
        expected_ids = ASSEMBLY_IDS
        cleaned = {}
        with open(input.graph_list) as handle:
            for line in handle:
                path = line.strip()
                if not path:
                    continue
                assembly_id = cleaned_assembly_id(path)
                if assembly_id in cleaned:
                    raise WorkflowError(
                        f"Duplicate assembly ID in cleaned assembly list: {assembly_id}"
                    )
                cleaned[assembly_id] = path

        missing = sorted(set(expected_ids) - set(cleaned))
        unexpected = sorted(set(cleaned) - set(expected_ids))
        if missing or unexpected:
            details = []
            if missing:
                details.append("missing cleaned assemblies: " + ", ".join(missing))
            if unexpected:
                details.append("unexpected cleaned assemblies: " + ", ".join(unexpected))
            raise WorkflowError("; ".join(details))

        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        with open(output[0], "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["assembly_id", "path"])
            for assembly_id in expected_ids:
                writer.writerow([assembly_id, cleaned[assembly_id]])


rule summarize_post_decontamination_qc:
    input:
        manifest=f"{POST_QC_OUTDIR}/resources/all_cleaned_assemblies.tsv",
        config=f"{QC_OUTDIR}/resources/resolved_qc_config.json",
        stats=all_cleaned_stats
    output:
        assembly=f"{POST_QC_OUTDIR}/summary/assembly_qc.tsv",
        sample=f"{POST_QC_OUTDIR}/summary/sample_qc.tsv",
        included=f"{POST_QC_OUTDIR}/summary/graph_included_assemblies.txt",
        excluded=f"{POST_QC_OUTDIR}/summary/graph_excluded_assemblies.tsv",
        complete=POST_QC_COMPLETE_MARKER
    params:
        results_dir=DECONTAM_OUTDIR,
        script=SUMMARIZE_QC_SCRIPT
    log:
        f"{POST_QC_OUTDIR}/logs/summarize_qc.log"
    conda:
        TOOLS_ENV
    threads: 1
    resources:
        mem_mb=8000,
        runtime_min=60
    shell:
        r"""
        python {params.script:q} \
          --manifest {input.manifest:q} \
          --config {input.config:q} \
          --results-dir {params.results_dir:q} \
          --seqkit-suffix .clean.seqkit.tsv \
          --sequence-only \
          --assembly-output {output.assembly:q} \
          --sample-output {output.sample:q} \
          --included-output {output.included:q} \
          --excluded-output {output.excluded:q} \
          --complete-marker {output.complete:q} > {log:q} 2>&1
        """
