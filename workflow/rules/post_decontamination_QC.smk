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


def post_qc_assembly_input(wildcards):
    return f"{DECONTAM_OUTDIR}/cleaned/{wildcards.assembly}.clean.fa.gz"


def all_post_qc_stats(wildcards):
    return expand(
        f"{POST_QC_OUTDIR}/stats/{{assembly}}.seqkit.tsv",
        assembly=ASSEMBLY_IDS,
    )


def all_original_compleasm(wildcards):
    return expand(
        f"{QC_OUTDIR}/compleasm/{{assembly}}/summary.txt",
        assembly=ASSEMBLY_IDS,
    )


def all_post_qc_alignment_metrics(reference):
    def paths(wildcards):
        return expand(
            f"{POST_QC_OUTDIR}/alignment_metrics/{reference}/{{assembly}}.tsv",
            assembly=ASSEMBLY_IDS,
        )

    return paths


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


rule post_qc_fasta_stats:
    input:
        assembly=post_qc_assembly_input
    output:
        temp(f"{POST_QC_OUTDIR}/stats/{{assembly}}.seqkit.tsv")
    log:
        f"{POST_QC_OUTDIR}/logs/seqkit/{{assembly}}.log"
    benchmark:
        f"{POST_QC_OUTDIR}/benchmarks/seqkit/{{assembly}}.tsv"
    conda:
        TOOLS_ENV
    threads: 1
    resources:
        mem_mb=4000,
        runtime_min=60
    shell:
        r"""
        seqkit stats --all --tabular {input.assembly:q} > {output:q} 2> {log:q}
        """


rule post_qc_align_to_reference:
    input:
        assembly=post_qc_assembly_input,
        index=f"{QC_OUTDIR}/resources/references/{{reference}}.mmi"
    output:
        temp(f"{POST_QC_OUTDIR}/alignments/{{reference}}/{{assembly}}.paf.gz")
    log:
        f"{POST_QC_OUTDIR}/logs/minimap2/{{reference}}.{{assembly}}.log"
    benchmark:
        f"{POST_QC_OUTDIR}/benchmarks/minimap2/{{reference}}.{{assembly}}.tsv"
    conda:
        TOOLS_ENV
    threads: config["resources"]["minimap2_threads"]
    resources:
        mem_mb=config["resources"]["minimap2_mem_mb"],
        runtime_min=config["resources"]["minimap2_runtime_min"]
    shell:
        r"""
        set -o pipefail
        minimap2 -x asm5 --secondary=no -c -t {threads} {input.index:q} {input.assembly:q} 2> {log:q} \
          | gzip -c > {output:q}
        """


rule post_qc_paf_metrics:
    input:
        paf=f"{POST_QC_OUTDIR}/alignments/{{reference}}/{{assembly}}.paf.gz",
        query_stats=f"{POST_QC_OUTDIR}/stats/{{assembly}}.seqkit.tsv",
        reference_stats=f"{QC_OUTDIR}/resources/references/{{reference}}.seqkit.tsv"
    output:
        temp(f"{POST_QC_OUTDIR}/alignment_metrics/{{reference}}/{{assembly}}.tsv")
    params:
        min_mapq=ALIGNMENT_MIN_MAPQ,
        script=PAF_METRICS_SCRIPT
    log:
        f"{POST_QC_OUTDIR}/logs/paf_metrics/{{reference}}.{{assembly}}.log"
    conda:
        TOOLS_ENV
    threads: 1
    resources:
        mem_mb=8000,
        runtime_min=60
    shell:
        r"""
        python {params.script:q} \
          --paf {input.paf:q} \
          --query-stats {input.query_stats:q} \
          --reference-stats {input.reference_stats:q} \
          --reference {wildcards.reference:q} \
          --min-mapq {params.min_mapq} \
          --output {output:q} > {log:q} 2>&1
        """


rule summarize_post_decontamination_qc:
    input:
        manifest=f"{POST_QC_OUTDIR}/resources/all_cleaned_assemblies.tsv",
        config=f"{QC_OUTDIR}/resources/resolved_qc_config.json",
        stats=all_post_qc_stats,
        original_compleasm=all_original_compleasm,
        chm13=all_post_qc_alignment_metrics("CHM13"),
        hg38=all_post_qc_alignment_metrics("hg38")
    output:
        assembly=f"{POST_QC_OUTDIR}/summary/assembly_qc.tsv",
        sample=f"{POST_QC_OUTDIR}/summary/sample_qc.tsv",
        included=f"{POST_QC_OUTDIR}/summary/graph_included_assemblies.txt",
        excluded=f"{POST_QC_OUTDIR}/summary/graph_excluded_assemblies.tsv",
        complete=POST_QC_COMPLETE_MARKER
    params:
        results_dir=POST_QC_OUTDIR,
        compleasm_results_dir=QC_OUTDIR,
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
          --compleasm-results-dir {params.compleasm_results_dir:q} \
          --assembly-output {output.assembly:q} \
          --sample-output {output.sample:q} \
          --included-output {output.included:q} \
          --excluded-output {output.excluded:q} \
          --complete-marker {output.complete:q} > {log:q} 2>&1
        """
