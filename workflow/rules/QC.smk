QC_CONFIG = workflow.source_path("../../config/config.yaml")
configfile: QC_CONFIG


import csv
import glob
import json
import os
from pathlib import Path


ASSEMBLY_DIR = config["assemblies"]["directory"]
QC_OUTDIR = config["results"]["qc"]
WORKFLOW_ROOT = config["workflow_root"]
TOOLS_ENV = f"{WORKFLOW_ROOT}/workflow/envs/tools.yaml"
COMPLEASM_ENV = f"{WORKFLOW_ROOT}/workflow/envs/compleasm.yaml"
PAF_METRICS_SCRIPT = f"{WORKFLOW_ROOT}/workflow/scripts/paf_metrics.py"
SUMMARIZE_QC_SCRIPT = f"{WORKFLOW_ROOT}/workflow/scripts/summarize_qc.py"
REFERENCES = config["references"]
ALIGNMENT_MIN_MAPQ = config["alignment_min_mapq"]

missing_references = {"CHM13", "hg38"} - set(REFERENCES)
if missing_references:
    raise WorkflowError(
        "Missing required reference keys: " + ", ".join(sorted(missing_references))
    )


def strip_fasta_suffix(path):
    name = os.path.basename(path)
    for suffix in (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    raise ValueError(f"Unsupported FASTA suffix: {path}")


assembly_files = []
for pattern in config["assemblies"]["patterns"]:
    assembly_files.extend(glob.glob(os.path.join(ASSEMBLY_DIR, pattern)))

ASSEMBLIES = {}
for path in sorted(set(assembly_files)):
    assembly_id = strip_fasta_suffix(path)
    if assembly_id in ASSEMBLIES:
        raise WorkflowError(f"Duplicate assembly ID after removing FASTA suffix: {assembly_id}")
    ASSEMBLIES[assembly_id] = path

ASSEMBLY_IDS = sorted(ASSEMBLIES)
if not ASSEMBLY_IDS:
    raise WorkflowError(
        f"No assemblies found in {ASSEMBLY_DIR}. Check assemblies.directory and "
        "assemblies.patterns in config/config.yaml."
    )


def assembly_input(wildcards):
    return ASSEMBLIES[wildcards.assembly]


rule qc:
    input:
        f"{QC_OUTDIR}/summary/assembly_qc.tsv",
        f"{QC_OUTDIR}/summary/sample_qc.tsv",
        f"{QC_OUTDIR}/summary/graph_included_assemblies.txt",
        f"{QC_OUTDIR}/summary/graph_excluded_assemblies.tsv",


rule qc_assembly_manifest:
    input:
        list(ASSEMBLIES.values())
    output:
        temp(f"{QC_OUTDIR}/resources/assemblies.tsv")
    run:
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        with open(output[0], "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["assembly_id", "path"])
            for assembly_id in ASSEMBLY_IDS:
                writer.writerow([assembly_id, ASSEMBLIES[assembly_id]])


rule resolved_qc_config:
    output:
        f"{QC_OUTDIR}/resources/resolved_qc_config.json"
    params:
        content=json.dumps({"thresholds": config["thresholds"]}, sort_keys=True)
    run:
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        with open(output[0], "w") as handle:
            json.dump(json.loads(params.content), handle, indent=2, sort_keys=True)
            handle.write("\n")


rule fasta_stats:
    input:
        assembly=assembly_input
    output:
        temp(f"{QC_OUTDIR}/stats/{{assembly}}.seqkit.tsv")
    log:
        f"{QC_OUTDIR}/logs/seqkit/{{assembly}}.log"
    benchmark:
        f"{QC_OUTDIR}/benchmarks/seqkit/{{assembly}}.tsv"
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


rule compleasm_download:
    output:
        touch(f"{QC_OUTDIR}/resources/compleasm/{config['compleasm']['lineage']}_{config['compleasm']['odb']}.ready")
    params:
        library=f"{QC_OUTDIR}/resources/compleasm/library",
        lineage=config["compleasm"]["lineage"],
        odb=config["compleasm"]["odb"]
    log:
        f"{QC_OUTDIR}/logs/compleasm/download.log"
    conda:
        COMPLEASM_ENV
    threads: 1
    resources:
        mem_mb=4000,
        runtime_min=120
    shell:
        r"""
        mkdir -p {params.library:q}
        compleasm download {params.lineage:q} -L {params.library:q} --odb {params.odb:q} > {log:q} 2>&1
        """


rule compleasm:
    input:
        assembly=assembly_input,
        lineage=f"{QC_OUTDIR}/resources/compleasm/{config['compleasm']['lineage']}_{config['compleasm']['odb']}.ready"
    output:
        # Post-decontamination QC reuses this original-assembly result.
        summary=f"{QC_OUTDIR}/compleasm/{{assembly}}/summary.txt"
    params:
        library=f"{QC_OUTDIR}/resources/compleasm/library",
        lineage=config["compleasm"]["lineage"],
        odb=config["compleasm"]["odb"]
    log:
        f"{QC_OUTDIR}/logs/compleasm/{{assembly}}.log"
    benchmark:
        f"{QC_OUTDIR}/benchmarks/compleasm/{{assembly}}.tsv"
    conda:
        COMPLEASM_ENV
    threads: config["resources"]["compleasm_threads"]
    resources:
        mem_mb=config["resources"]["compleasm_mem_mb"],
        runtime_min=config["resources"]["compleasm_runtime_min"]
    shell:
        r"""
        # Compleasm creates many temporary files. Set TMPDIR in the shell or
        # cluster profile to use node-local scratch; /tmp is the fallback.
        scratch_root="${{TMPDIR:-/tmp}}"
        if [ ! -d "$scratch_root" ] || [ ! -w "$scratch_root" ]; then
            echo "TMPDIR is not a writable directory: $scratch_root" >&2
            exit 1
        fi

        workdir=$(mktemp -d -p "$scratch_root" "compleasm.{wildcards.assembly}.XXXXXX")
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


rule reference_index:
    input:
        lambda wildcards: REFERENCES[wildcards.reference]
    output:
        temp(f"{QC_OUTDIR}/resources/references/{{reference}}.mmi")
    log:
        f"{QC_OUTDIR}/logs/minimap2/{{reference}}.index.log"
    benchmark:
        f"{QC_OUTDIR}/benchmarks/minimap2/{{reference}}.index.tsv"
    conda:
        TOOLS_ENV
    threads: config["resources"]["minimap2_threads"]
    resources:
        mem_mb=config["resources"]["minimap2_index_mem_mb"],
        runtime_min=config["resources"]["minimap2_index_runtime_min"]
    shell:
        r"""
        minimap2 -x asm5 -t {threads} -d {output:q} {input:q} 2> {log:q}
        """


rule reference_stats:
    input:
        lambda wildcards: REFERENCES[wildcards.reference]
    output:
        temp(f"{QC_OUTDIR}/resources/references/{{reference}}.seqkit.tsv")
    log:
        f"{QC_OUTDIR}/logs/seqkit/{{reference}}.reference.log"
    conda:
        TOOLS_ENV
    threads: 1
    resources:
        mem_mb=4000,
        runtime_min=60
    shell:
        r"""
        seqkit stats --all --tabular {input:q} > {output:q} 2> {log:q}
        """


rule align_to_reference:
    input:
        assembly=assembly_input,
        index=f"{QC_OUTDIR}/resources/references/{{reference}}.mmi"
    output:
        temp(f"{QC_OUTDIR}/alignments/{{reference}}/{{assembly}}.paf.gz")
    log:
        f"{QC_OUTDIR}/logs/minimap2/{{reference}}.{{assembly}}.log"
    benchmark:
        f"{QC_OUTDIR}/benchmarks/minimap2/{{reference}}.{{assembly}}.tsv"
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


rule paf_metrics:
    input:
        paf=f"{QC_OUTDIR}/alignments/{{reference}}/{{assembly}}.paf.gz",
        query_stats=f"{QC_OUTDIR}/stats/{{assembly}}.seqkit.tsv",
        reference_stats=f"{QC_OUTDIR}/resources/references/{{reference}}.seqkit.tsv"
    output:
        temp(f"{QC_OUTDIR}/alignment_metrics/{{reference}}/{{assembly}}.tsv")
    params:
        min_mapq=ALIGNMENT_MIN_MAPQ,
        script=PAF_METRICS_SCRIPT
    log:
        f"{QC_OUTDIR}/logs/paf_metrics/{{reference}}.{{assembly}}.log"
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


checkpoint summarize_qc:
    input:
        manifest=f"{QC_OUTDIR}/resources/assemblies.tsv",
        config=f"{QC_OUTDIR}/resources/resolved_qc_config.json",
        stats=expand(f"{QC_OUTDIR}/stats/{{assembly}}.seqkit.tsv", assembly=ASSEMBLY_IDS),
        compleasm=expand(f"{QC_OUTDIR}/compleasm/{{assembly}}/summary.txt", assembly=ASSEMBLY_IDS),
        chm13=expand(f"{QC_OUTDIR}/alignment_metrics/CHM13/{{assembly}}.tsv", assembly=ASSEMBLY_IDS),
        hg38=expand(f"{QC_OUTDIR}/alignment_metrics/hg38/{{assembly}}.tsv", assembly=ASSEMBLY_IDS)
    output:
        assembly=f"{QC_OUTDIR}/summary/assembly_qc.tsv",
        sample=f"{QC_OUTDIR}/summary/sample_qc.tsv",
        included=f"{QC_OUTDIR}/summary/graph_included_assemblies.txt",
        excluded=f"{QC_OUTDIR}/summary/graph_excluded_assemblies.tsv"
    params:
        results_dir=QC_OUTDIR,
        script=SUMMARIZE_QC_SCRIPT
    log:
        f"{QC_OUTDIR}/logs/summarize_qc.log"
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
          --excluded-output {output.excluded:q} > {log:q} 2>&1
        """
