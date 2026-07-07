POST_QC_CONFIG = workflow.source_path("../../config/config.yaml")
configfile: POST_QC_CONFIG

POST_QC_ENTRY_OUTDIR = config["results"]["post_decontamination_qc"]


rule post_decontamination_qc:
    input:
        f"{POST_QC_ENTRY_OUTDIR}/summary/assembly_qc.tsv",
        f"{POST_QC_ENTRY_OUTDIR}/summary/sample_qc.tsv",
        f"{POST_QC_ENTRY_OUTDIR}/summary/graph_included_assemblies.txt",
        f"{POST_QC_ENTRY_OUTDIR}/summary/graph_excluded_assemblies.tsv",


include: "decontamination.smk"


import csv
from pathlib import Path


POST_QC_OUTDIR = config["results"]["post_decontamination_qc"]


def cleaned_assembly_id(path):
    name = strip_fasta_suffix(path)
    if not name.endswith(".clean"):
        raise WorkflowError(
            f"Expected a cleaned FASTA name ending in '.clean': {path}"
        )
    return name[: -len(".clean")]


def post_qc_assembly_input(wildcards):
    return f"{DECONTAM_OUTDIR}/cleaned/{wildcards.assembly}.clean.fa.gz"


def selected_post_qc_stats(wildcards):
    return expand(
        f"{POST_QC_OUTDIR}/stats/{{assembly}}.seqkit.tsv",
        assembly=selected_assembly_ids(wildcards),
    )


def selected_post_qc_compleasm(wildcards):
    return expand(
        f"{POST_QC_OUTDIR}/compleasm/{{assembly}}/summary.txt",
        assembly=selected_assembly_ids(wildcards),
    )


def selected_post_qc_alignment_metrics(reference):
    def paths(wildcards):
        return expand(
            f"{POST_QC_OUTDIR}/alignment_metrics/{reference}/{{assembly}}.tsv",
            assembly=selected_assembly_ids(wildcards),
        )

    return paths


rule post_qc_assembly_manifest:
    input:
        graph_list=f"{DECONTAM_OUTDIR}/summary/graph_cleaned_assemblies.txt"
    output:
        f"{POST_QC_OUTDIR}/resources/assemblies.tsv"
    run:
        expected_ids = selected_assembly_ids(wildcards)
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


rule post_qc_compleasm:
    input:
        assembly=post_qc_assembly_input,
        lineage=f"{QC_OUTDIR}/resources/compleasm/{config['compleasm']['lineage']}_{config['compleasm']['odb']}.ready"
    output:
        summary=temp(f"{POST_QC_OUTDIR}/compleasm/{{assembly}}/summary.txt")
    params:
        library=f"{QC_OUTDIR}/resources/compleasm/library",
        lineage=config["compleasm"]["lineage"],
        odb=config["compleasm"]["odb"]
    log:
        f"{POST_QC_OUTDIR}/logs/compleasm/{{assembly}}.log"
    benchmark:
        f"{POST_QC_OUTDIR}/benchmarks/compleasm/{{assembly}}.tsv"
    conda:
        COMPLEASM_ENV
    threads: config["resources"]["compleasm_threads"]
    resources:
        mem_mb=config["resources"]["compleasm_mem_mb"],
        runtime_min=config["resources"]["compleasm_runtime_min"]
    shell:
        r"""
        scratch_root="${{TMPDIR:-/tmp}}"
        if [ ! -d "$scratch_root" ] || [ ! -w "$scratch_root" ]; then
            echo "TMPDIR is not a writable directory: $scratch_root" >&2
            exit 1
        fi

        workdir=$(mktemp -d -p "$scratch_root" "compleasm.post-clean.{wildcards.assembly}.XXXXXX")
        trap 'rm -rf "$workdir"' EXIT
        printf 'Compleasm scratch directory: %s\n' "$workdir" > {log:q}
        compleasm run \
          -a {input.assembly:q} \
          -o "$workdir/result" \
          -t {threads} \
          -l {params.lineage:q} \
          -L {params.library:q} \
          --odb {params.odb:q} >> {log:q} 2>&1
        mkdir -p "$(dirname {output.summary:q})"
        cp "$workdir/result/summary.txt" {output.summary:q}
        test -s {output.summary:q}
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
        manifest=f"{POST_QC_OUTDIR}/resources/assemblies.tsv",
        config=f"{QC_OUTDIR}/resources/resolved_qc_config.json",
        stats=selected_post_qc_stats,
        compleasm=selected_post_qc_compleasm,
        chm13=selected_post_qc_alignment_metrics("CHM13"),
        hg38=selected_post_qc_alignment_metrics("hg38")
    output:
        assembly=f"{POST_QC_OUTDIR}/summary/assembly_qc.tsv",
        sample=f"{POST_QC_OUTDIR}/summary/sample_qc.tsv",
        included=f"{POST_QC_OUTDIR}/summary/graph_included_assemblies.txt",
        excluded=f"{POST_QC_OUTDIR}/summary/graph_excluded_assemblies.tsv"
    params:
        results_dir=POST_QC_OUTDIR,
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
          --assembly-output {output.assembly:q} \
          --sample-output {output.sample:q} \
          --included-output {output.included:q} \
          --excluded-output {output.excluded:q} \
          --allow-missing-mates > {log:q} 2>&1
        """
