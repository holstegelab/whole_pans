DECONTAM_CONFIG = workflow.source_path("../../config/config.yaml")
configfile: DECONTAM_CONFIG

DECONTAM_ENTRY_OUTDIR = config["results"]["decontamination"]


rule decontamination:
    input:
        f"{DECONTAM_ENTRY_OUTDIR}/summary/contamination_summary.tsv",
        f"{DECONTAM_ENTRY_OUTDIR}/summary/contig_actions.tsv",
        f"{DECONTAM_ENTRY_OUTDIR}/summary/review_candidates.tsv",
        f"{DECONTAM_ENTRY_OUTDIR}/summary/graph_cleaned_assemblies.txt",


include: "QC.smk"


import csv
import json
import os
import re
import shlex
from pathlib import Path


QC_RESULTS_DIR = config["results"]["qc"]
DECONTAM_OUTDIR = config["results"]["decontamination"]
DECONTAM_ENV = f"{WORKFLOW_ROOT}/workflow/envs/tools.yaml"
CLASSIFY_CONTAMINATION_SCRIPT = (
    f"{WORKFLOW_ROOT}/workflow/scripts/classify_contamination.py"
)
FILTER_FASTA_SCRIPT = f"{WORKFLOW_ROOT}/workflow/scripts/filter_fasta.py"
SUMMARIZE_DECONTAMINATION_SCRIPT = (
    f"{WORKFLOW_ROOT}/workflow/scripts/summarize_decontamination.py"
)
DATABASES = config["nonhuman_databases"]
DATABASE_NAMES = sorted(DATABASES)
EXISTING_QC_INCLUDED = config.get("existing_qc", {}).get("included_assemblies")


def strip_fasta_suffix(path):
    name = os.path.basename(path)
    for suffix in (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    raise ValueError(f"Unsupported FASTA suffix: {path}")


if not DATABASE_NAMES:
    raise WorkflowError("No non-human BLAST databases configured")


def decontamination_assembly_input(wildcards):
    return ASSEMBLIES[wildcards.assembly]


def qc_included_list(wildcards):
    if EXISTING_QC_INCLUDED:
        return EXISTING_QC_INCLUDED
    return checkpoints.summarize_qc.get().output.included


def qc_selection_dependency(wildcards):
    if EXISTING_QC_INCLUDED:
        # The completed selection is deliberately treated as external state.
        # Declaring it as an input would connect it to summarize_qc and rebuild
        # cleaned temporary files from the original QC run.
        return []
    return checkpoints.summarize_qc.get().output.included


def selected_assembly_ids(wildcards):
    included_list = qc_included_list(wildcards)
    selected = []
    seen = set()
    with open(included_list) as handle:
        for line in handle:
            path = line.strip()
            if not path:
                continue
            assembly_id = strip_fasta_suffix(path)
            if assembly_id not in ASSEMBLIES:
                raise WorkflowError(
                    f"QC selected an assembly not discovered at startup: {path}"
                )
            if assembly_id in seen:
                raise WorkflowError(
                    f"Duplicate assembly ID in QC inclusion list: {assembly_id}"
                )
            seen.add(assembly_id)
            selected.append(assembly_id)
    if not selected:
        raise WorkflowError(
            "QC did not select any assemblies for decontamination. "
            "Review the QC summary and thresholds."
        )
    return sorted(selected)


def selected_decisions(wildcards):
    return expand(
        f"{DECONTAM_OUTDIR}/decisions/{{assembly}}.contigs.tsv",
        assembly=selected_assembly_ids(wildcards),
    )


def selected_cleaned_assemblies(wildcards):
    return expand(
        f"{DECONTAM_OUTDIR}/cleaned/{{assembly}}.clean.fa.gz",
        assembly=selected_assembly_ids(wildcards),
    )


def selected_cleaned_stats(wildcards):
    return expand(
        f"{DECONTAM_OUTDIR}/stats/{{assembly}}.clean.seqkit.tsv",
        assembly=selected_assembly_ids(wildcards),
    )


def database_marker(wildcards):
    return DATABASES[wildcards.database]["marker"]


def database_prefix(wildcards):
    return DATABASES[wildcards.database]["prefix"]


def blast_arguments(wildcards):
    arguments = []
    for database in DATABASE_NAMES:
        path = f"{DECONTAM_OUTDIR}/blast/{database}/{wildcards.assembly}.tsv.gz"
        arguments.extend(["--blast", f"{database}={path}"])
    return " ".join(shlex.quote(value) for value in arguments)


wildcard_constraints:
    database="|".join(re.escape(name) for name in DATABASE_NAMES)


rule decontamination_assembly_manifest:
    input:
        selection=qc_selection_dependency
    output:
        f"{DECONTAM_OUTDIR}/resources/assemblies.tsv"
    run:
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        with open(output[0], "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["assembly_id", "source_path"])
            with open(qc_included_list(wildcards)) as included_handle:
                for line in included_handle:
                    path = line.strip()
                    if path:
                        writer.writerow([strip_fasta_suffix(path), path])


rule resolved_classification_config:
    output:
        f"{DECONTAM_OUTDIR}/resources/resolved_classification_config.json"
    params:
        content=json.dumps(
            {"blast": config["blast"], "classification": config["classification"]},
            sort_keys=True,
        )
    run:
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        with open(output[0], "w") as handle:
            json.dump(json.loads(params.content), handle, indent=2, sort_keys=True)
            handle.write("\n")


rule blast_nonhuman:
    input:
        assembly=decontamination_assembly_input,
        database=database_marker
    output:
        f"{DECONTAM_OUTDIR}/blast/{{database}}/{{assembly}}.tsv.gz"
    params:
        prefix=database_prefix,
        word_size=config["blast"]["word_size"],
        evalue=config["blast"]["evalue"],
        min_identity=config["blast"]["min_identity"],
        max_target_seqs=config["blast"]["max_target_seqs"]
    log:
        f"{DECONTAM_OUTDIR}/logs/blast/{{database}}.{{assembly}}.log"
    benchmark:
        f"{DECONTAM_OUTDIR}/benchmarks/blast/{{database}}.{{assembly}}.tsv"
    conda:
        DECONTAM_ENV
    threads: config["resources"]["blast_threads"]
    resources:
        mem_mb=config["resources"]["blast_mem_mb"],
        runtime_min=config["resources"]["blast_runtime_min"]
    shell:
        r"""
        set -o pipefail
        blastn \
          -query {input.assembly:q} \
          -db {params.prefix:q} \
          -task megablast \
          -word_size {params.word_size} \
          -best_hit_overhang 0.1 \
          -best_hit_score_edge 0.1 \
          -dust yes \
          -evalue {params.evalue} \
          -min_raw_gapped_score 100 \
          -penalty -5 \
          -perc_identity {params.min_identity} \
          -soft_masking true \
          -max_target_seqs {params.max_target_seqs} \
          -num_threads {threads} \
          -outfmt '6 qseqid qlen sseqid pident length qstart qend evalue bitscore' \
          2> {log:q} | gzip -c > {output:q}
        """


rule classify_contigs:
    input:
        assembly=decontamination_assembly_input,
        blast=[
            f"{DECONTAM_OUTDIR}/blast/{database}/{{assembly}}.tsv.gz"
            for database in DATABASE_NAMES
        ],
        chm13=f"{QC_RESULTS_DIR}/alignments/CHM13/{{assembly}}.paf.gz",
        hg38=f"{QC_RESULTS_DIR}/alignments/hg38/{{assembly}}.paf.gz",
        config=f"{DECONTAM_OUTDIR}/resources/resolved_classification_config.json"
    output:
        decisions=f"{DECONTAM_OUTDIR}/decisions/{{assembly}}.contigs.tsv",
        split_bed=f"{DECONTAM_OUTDIR}/decisions/{{assembly}}.split_nonhuman.bed",
        remove_list=f"{DECONTAM_OUTDIR}/decisions/{{assembly}}.remove_contigs.txt",
        review_list=f"{DECONTAM_OUTDIR}/decisions/{{assembly}}.review_contigs.txt"
    params:
        blast_args=blast_arguments,
        script=CLASSIFY_CONTAMINATION_SCRIPT
    log:
        f"{DECONTAM_OUTDIR}/logs/classify/{{assembly}}.log"
    conda:
        DECONTAM_ENV
    threads: 1
    resources:
        mem_mb=16000,
        runtime_min=120
    shell:
        r"""
        python {params.script:q} \
          --assembly {input.assembly:q} \
          --assembly-id {wildcards.assembly:q} \
          --human-paf {input.chm13:q} \
          --human-paf {input.hg38:q} \
          {params.blast_args} \
          --config {input.config:q} \
          --decisions {output.decisions:q} \
          --split-bed {output.split_bed:q} \
          --remove-list {output.remove_list:q} \
          --review-list {output.review_list:q} > {log:q} 2>&1
        """


rule clean_assembly:
    input:
        assembly=decontamination_assembly_input,
        decisions=f"{DECONTAM_OUTDIR}/decisions/{{assembly}}.contigs.tsv",
        split_bed=f"{DECONTAM_OUTDIR}/decisions/{{assembly}}.split_nonhuman.bed"
    output:
        cleaned=f"{DECONTAM_OUTDIR}/cleaned/{{assembly}}.clean.fa.gz",
        removed=f"{DECONTAM_OUTDIR}/removed/{{assembly}}.nonhuman.fa.gz",
        review=f"{DECONTAM_OUTDIR}/review/{{assembly}}.review.fa.gz",
        split_map=f"{DECONTAM_OUTDIR}/split_maps/{{assembly}}.split_map.tsv"
    params:
        script=FILTER_FASTA_SCRIPT
    log:
        f"{DECONTAM_OUTDIR}/logs/clean/{{assembly}}.log"
    benchmark:
        f"{DECONTAM_OUTDIR}/benchmarks/clean/{{assembly}}.tsv"
    conda:
        DECONTAM_ENV
    threads: 1
    resources:
        mem_mb=8000,
        runtime_min=120
    shell:
        r"""
        python {params.script:q} \
          --assembly {input.assembly:q} \
          --decisions {input.decisions:q} \
          --split-bed {input.split_bed:q} \
          --cleaned {output.cleaned:q} \
          --removed {output.removed:q} \
          --review {output.review:q} \
          --split-map {output.split_map:q} > {log:q} 2>&1
        """


rule cleaned_stats:
    input:
        f"{DECONTAM_OUTDIR}/cleaned/{{assembly}}.clean.fa.gz"
    output:
        f"{DECONTAM_OUTDIR}/stats/{{assembly}}.clean.seqkit.tsv"
    log:
        f"{DECONTAM_OUTDIR}/logs/seqkit/{{assembly}}.clean.log"
    conda:
        DECONTAM_ENV
    threads: 1
    resources:
        mem_mb=4000,
        runtime_min=60
    shell:
        r"""
        seqkit stats --all --tabular {input:q} > {output:q} 2> {log:q}
        """


rule summarize_decontamination:
    input:
        manifest=f"{DECONTAM_OUTDIR}/resources/assemblies.tsv",
        decisions=selected_decisions,
        cleaned=selected_cleaned_assemblies,
        stats=selected_cleaned_stats
    output:
        summary=f"{DECONTAM_OUTDIR}/summary/contamination_summary.tsv",
        actions=f"{DECONTAM_OUTDIR}/summary/contig_actions.tsv",
        review=f"{DECONTAM_OUTDIR}/summary/review_candidates.tsv",
        graph_list=f"{DECONTAM_OUTDIR}/summary/graph_cleaned_assemblies.txt"
    params:
        results_dir=DECONTAM_OUTDIR,
        script=SUMMARIZE_DECONTAMINATION_SCRIPT
    log:
        f"{DECONTAM_OUTDIR}/logs/summarize.log"
    conda:
        DECONTAM_ENV
    threads: 1
    resources:
        mem_mb=8000,
        runtime_min=60
    shell:
        r"""
        python {params.script:q} \
          --manifest {input.manifest:q} \
          --results-dir {params.results_dir:q} \
          --summary {output.summary:q} \
          --actions {output.actions:q} \
          --review {output.review:q} \
          --graph-list {output.graph_list:q} > {log:q} 2>&1
        """
