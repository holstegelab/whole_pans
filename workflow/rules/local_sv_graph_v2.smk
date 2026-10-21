"""Build a local SV graph from HiFi-validated rescue sequences."""

configfile: workflow.source_path("../../config/config.yaml")


import csv
import gzip
import re
import shlex
from pathlib import Path


PROJECT_ROOT = Path(workflow.current_basedir).resolve().parents[1]
workdir: str(PROJECT_ROOT)

CFG = config["local_sv_graph_v2"]
HIFI = CFG["hifi_reads"]
PLAN_CFG = CFG["plan"]
RESCUE_CFG = CFG["rescue"]
GRAPH_CFG = CFG["graph"]
QC_CFG = CFG.get("qc", {})
DISCOVERY_SCREEN = config.get("novel_sv_discovery", {}).get("graph_screen", {})

OUT = CFG["results"]
FROZEN_GRAPH = CFG["frozen_graph"]
FROZEN_SUMMARY = CFG["frozen_graph_summary"]
FROZEN_INVENTORY = CFG["frozen_graph_inventory"]
FROZEN_BUBBLES = CFG["frozen_graph_bubbles"]
CATALOG = CFG["catalog"]
CATALOG_EVIDENCE = CFG["catalog_evidence"]
GRAPH_CANDIDATES = CFG["graph_candidates"]
ASSEMBLY_MANIFEST = CFG["assembly_manifest"]
CENTENARIAN_SAMPLES = HIFI["centenarian_samples"]
HEADER_TEMPLATE = HIFI["header_template"]

ARCHIVE_SCRIPT = str(workflow.source_path("../scripts/archive_hifi_reads.py"))
RESCUE_SCRIPT = str(workflow.source_path("../scripts/build_clean_sv_rescues.py"))
GRAPH_SCRIPT = str(workflow.source_path("../scripts/build_incremental_sv_graph.py"))
QC_SCRIPT = str(workflow.source_path("../scripts/qc_local_sv_graph_v2.py"))
RESCUE_ENV = str(workflow.source_path("../envs/sv_rescue.yaml"))
GRAPH_ENV = str(workflow.source_path("../envs/sv_graph.yaml"))
GRAPH_ENV_INPUT = "workflow/envs/sv_graph.yaml"

READ_MANIFEST = f"{OUT}/hifi/hifi_archive_reads.tsv"
STAGED = f"{OUT}/hifi/archive/staged"
STAGED_BAM = f"{OUT}/hifi/archive/compute_bams/{{sample}}.bam"
FASTQ = f"{OUT}/hifi/fastq/{{sample}}.hifi.fastq.gz"
HEADER_QC = f"{OUT}/hifi/header_qc/{{sample}}.header_qc.tsv"
PLAN = f"{OUT}/rescue/plan.tsv"
RESCUE_RESULTS = f"{OUT}/rescue/results"
RESCUE_RESULT = f"{RESCUE_RESULTS}/{{rescue}}"
RESCUE_SUMMARY = f"{OUT}/rescue/rescue_summary.tsv"
RESCUE_FASTA = f"{OUT}/rescue/clean_rescues.fa"
GRAPH_INPUTS = f"{OUT}/rescue/clean_rescue_graph_inputs.tsv"
GRAPH = f"{OUT}/graphs/sv_pangenome.local.v2.rescued.minigraph.gfa"
GRAPH_SUMMARY = f"{OUT}/metadata/sv_pangenome.local.v2.graph_summary.tsv"
GRAPH_VERSIONS = f"{OUT}/metadata/tool_versions.tsv"
GRAPH_PROVENANCE = f"{OUT}/metadata/sv_pangenome.local.v2.provenance.tsv"
GRAPH_BUBBLES = f"{OUT}/metadata/sv_pangenome.local.v2.bubbles.tsv"
QC_DIR = f"{OUT}/qc"
QC_INTEGRITY = f"{QC_DIR}/tables/graph_integrity.tsv"
QC_GRAPH_METRICS = f"{QC_DIR}/tables/graph_qc_metrics.tsv"
QC_DELTA_SEGMENTS = f"{QC_DIR}/tables/delta_segments.tsv"
QC_DELTA_LINKS = f"{QC_DIR}/tables/delta_links.tsv"
QC_BUBBLE_DELTA = f"{QC_DIR}/tables/bubble_delta.tsv"
QC_RESCUES = f"{QC_DIR}/tables/rescued_sv_graph_qc.tsv"
QC_COHORT = f"{QC_DIR}/tables/cohort_completeness.tsv"
QC_REPORT = f"{QC_DIR}/report/local_sv_graph_v2_qc_report.md"

SAFE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")


def read_table(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_ids(path):
    with open(path) as handle:
        return {
            line.split()[0]
            for raw in handle
            if (line := raw.strip()) and not line.startswith("#")
        }


CENTENARIANS = read_ids(CENTENARIAN_SAMPLES)


def cohort(sample):
    if not SAFE_ID.fullmatch(sample):
        raise WorkflowError(f"Unsafe sample ID: {sample}")
    return "centenarian" if sample in CENTENARIANS else "ad"


def archive_bam(sample):
    return HIFI["archive_bam_pattern"].format(
        archive_root=HIFI["archive_root"],
        cohort=cohort(sample),
        sample=sample,
    )


def cli_options(values, defaults):
    """Convert a config section to quoted command-line options."""
    options = []
    for name, default in defaults.items():
        value = values.get(name, default)
        flag = "--" + name.replace("_", "-")
        if isinstance(default, bool):
            if value:
                options.append(flag)
        else:
            options.extend((flag, shlex.quote(str(value))))
    return " ".join(options)


PLAN_OPTIONS = cli_options(
    PLAN_CFG,
    {
        "confidence": "UNCERTAIN",
        "validation_status": "PENDING_HIFI",
        "require_methods": "graph_residual;assembly_dipcall;assembly_svim_asm",
        "min_independent_samples": 2,
        "max_events": 100,
        "source_qc_tiers": "fragmented_rescue;not_recommended",
        "breakpoint_distance": 500,
        "length_similarity": 0.70,
        "flank_bp": 50000,
        "min_flank_bp": 20000,
        "allow_contig_end": False,
    },
)

RESCUE_OPTIONS = cli_options(
    RESCUE_CFG,
    {
        "max_secondary": 20,
        "min_read_mapq": 10,
        "min_read_aligned_bp": 2000,
        "read_anchor_bp": 1000,
        "max_read_event_discordance_bp": 30,
        "min_selected_reads": 8,
        "min_source_allele_reads": 3,
        "contig_anchor_bp": 10000,
        "min_contig_mapq": 20,
        "min_contig_identity": 0.99,
        "discordance_window_bp": 500,
        "max_event_discordance_bp": 20,
        "min_rescue_length_bp": 20000,
        "hifiasm_extra": "",
    },
)

QC_REPORT_OPTIONS = cli_options(
    QC_CFG,
    {
        "breakpoint_distance": PLAN_CFG.get("breakpoint_distance", 500),
        "length_similarity": PLAN_CFG.get("length_similarity", 0.70),
        "min_allele_reads_warn": 5,
        "min_allele_fraction_warn": 0.10,
        "max_selected_reads_warn": 100000,
        "max_local_contigs_warn": 10000,
        "min_target_coverage_warn": 0.80,
    },
)


def current_plan():
    return str(checkpoints.plan_local_v2_rescues.get().output.plan)


def plan_rows():
    return read_table(current_plan())


def plan_row(wildcards):
    matches = [
        row for row in plan_rows() if row["rescue_id"] == wildcards.rescue
    ]
    if len(matches) != 1 or matches[0]["plan_status"] != "READY":
        raise WorkflowError(f"No unique READY plan row for {wildcards.rescue}")
    return matches[0]


def ready_results(wildcards):
    rescue_ids = sorted(
        row["rescue_id"]
        for row in plan_rows()
        if row["plan_status"] == "READY"
    )
    return expand(RESCUE_RESULT, rescue=rescue_ids)


def rescue_reads(wildcards):
    sample = plan_row(wildcards)["sample_id"]
    return FASTQ.format(sample=sample)


def rescue_assembly(wildcards):
    return plan_row(wildcards)["assembly_path"]


localrules: local_sv_graph_v2, write_local_v2_hifi_manifest


rule local_sv_graph_v2:
    input:
        READ_MANIFEST,
        PLAN,
        RESCUE_SUMMARY,
        RESCUE_FASTA,
        GRAPH_INPUTS,
        GRAPH,
        GRAPH_SUMMARY,
        GRAPH_VERSIONS,
        GRAPH_PROVENANCE,
        GRAPH_BUBBLES,
        QC_INTEGRITY,
        QC_DELTA_SEGMENTS,
        QC_DELTA_LINKS,
        QC_RESCUES,
        QC_COHORT,
        QC_REPORT


rule write_local_v2_hifi_manifest:
    input:
        assembly_manifest=ASSEMBLY_MANIFEST,
        centenarians=CENTENARIAN_SAMPLES
    output:
        READ_MANIFEST
    params:
        script=ARCHIVE_SCRIPT,
        archive_root=HIFI["archive_root"],
        archive_pattern=lambda wildcards: HIFI["archive_bam_pattern"],
        fastq_pattern=lambda wildcards: FASTQ
    log:
        f"{OUT}/logs/write_hifi_archive_manifest.log"
    conda:
        RESCUE_ENV
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {params.script:q} manifest \
          --assembly-manifest {input.assembly_manifest:q} \
          --centenarian-samples {input.centenarians:q} \
          --archive-root {params.archive_root:q} \
          --archive-pattern {params.archive_pattern:q} \
          --fastq-pattern {params.fastq_pattern:q} \
          --output {output:q} > {log:q} 2>&1
        """


checkpoint plan_local_v2_rescues:
    input:
        catalog=CATALOG,
        graph_candidates=GRAPH_CANDIDATES,
        assembly_manifest=ASSEMBLY_MANIFEST,
        read_manifest=READ_MANIFEST,
        event_ids=PLAN_CFG.get("event_ids") or []
    output:
        plan=PLAN
    params:
        script=RESCUE_SCRIPT,
        options=PLAN_OPTIONS,
        event_ids=lambda wildcards, input: (
            f"--event-ids {shlex.quote(str(input.event_ids))}"
            if input.event_ids
            else ""
        )
    log:
        f"{OUT}/logs/plan_rescues.log"
    conda:
        RESCUE_ENV
    resources:
        mem_mb=int(PLAN_CFG.get("mem_mb", 16000)),
        runtime=int(PLAN_CFG.get("runtime", 360))
    shell:
        r"""
        mkdir -p "$(dirname {output.plan:q})" "$(dirname {log:q})"
        python {params.script:q} plan \
          --catalog {input.catalog:q} \
          --graph-candidates {input.graph_candidates:q} \
          --assembly-manifest {input.assembly_manifest:q} \
          --read-manifest {input.read_manifest:q} \
          --output {output.plan:q} \
          {params.options} {params.event_ids} > {log:q} 2>&1
        """


rule stage_local_v2_hifi_bams:
    input:
        plan=lambda wildcards: current_plan(),
        read_manifest=READ_MANIFEST
    output:
        temp(directory(STAGED))
    params:
        script=ARCHIVE_SCRIPT,
        daget=HIFI["daget"],
        dals=HIFI["dals"],
        poll_seconds=int(HIFI.get("poll_seconds", 30)),
        max_wait=int(HIFI.get("max_wait_min", 1440))
    log:
        f"{OUT}/logs/archive/daget.log"
    conda:
        RESCUE_ENV
    resources:
        partition=HIFI.get("archive_partition", "archive"),
        mem_mb=int(HIFI.get("stage_mem_mb", 100)),
        runtime=int(HIFI.get("stage_runtime", 1500))
    shell:
        r"""
        python {params.script:q} stage \
          --plan {input.plan:q} \
          --read-manifest {input.read_manifest:q} \
          --output-dir {output:q} \
          --daget {params.daget:q} \
          --dals {params.dals:q} \
          --poll-seconds {params.poll_seconds} \
          --max-wait-minutes {params.max_wait} \
          --log {log:q}
        """


rule copy_local_v2_hifi_bam:
    input:
        staged=STAGED
    output:
        bam=temp(STAGED_BAM)
    params:
        script=ARCHIVE_SCRIPT,
        source=lambda wildcards: archive_bam(wildcards.sample),
        marker=lambda wildcards: f"{STAGED}/{wildcards.sample}.staged"
    log:
        f"{OUT}/logs/archive/copy/{{sample}}.log"
    conda:
        RESCUE_ENV
    resources:
        partition=HIFI.get("archive_partition", "archive"),
        mem_mb=int(HIFI.get("copy_mem_mb", 100)),
        runtime=int(HIFI.get("copy_runtime", 720))
    shell:
        r"""
        python {params.script:q} copy \
          --source {params.source:q} \
          --marker {params.marker:q} \
          --output {output.bam:q} \
          --log {log:q}
        """


rule convert_local_v2_hifi_bam:
    input:
        staged_bam=STAGED_BAM,
        header_template=HEADER_TEMPLATE
    output:
        reads=protected(FASTQ),
        header_qc=HEADER_QC
    params:
        archive_bam=lambda wildcards: archive_bam(wildcards.sample),
        cohort=lambda wildcards: cohort(wildcards.sample)
    log:
        f"{OUT}/logs/hifi_conversion/{{sample}}.log"
    benchmark:
        f"{OUT}/benchmarks/hifi_conversion/{{sample}}.tsv"
    conda:
        RESCUE_ENV
    threads:
        int(HIFI.get("conversion_threads", 8))
    resources:
        mem_mb=int(HIFI.get("conversion_mem_mb", 16000)),
        runtime=int(HIFI.get("conversion_runtime", 720)),
        disk_mb=int(HIFI.get("conversion_disk_mb", 1000))
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.reads:q})" \
                 "$(dirname {output.header_qc:q})" \
                 "$(dirname {log:q})"
        test -s {input.staged_bam:q}
        test -r {input.staged_bam:q}

        header_work=$(mktemp -d -p "${{TMPDIR:-/tmp}}" "local-v2-hifi-header.XXXXXX")
        project_work=$(mktemp -d -p "$(dirname {output.reads:q})" \
          ".{wildcards.sample}.hifi-fastq.XXXXXX")
        trap 'rm -rf -- "$header_work" "$project_work"' EXIT
        samtools quickcheck -v {input.staged_bam:q}
        samtools view -H {input.staged_bam:q} > "$header_work/header.sam"
        awk '$1 == "@SQ"' {input.header_template:q} > "$header_work/template.sq"
        awk '$1 == "@SQ"' "$header_work/header.sam" > "$header_work/actual.sq"
        cmp "$header_work/template.sq" "$header_work/actual.sq"
        awk -F '\t' '
          $1 == "@RG" {{
            pacbio = ccs = 0
            for (i = 2; i <= NF; i++) {{
              if ($i == "PL:PACBIO") pacbio = 1
              if ($i ~ /^DS:/ && index($i, "READTYPE=CCS")) ccs = 1
            }}
            if (pacbio && ccs) found = 1
          }}
          END {{exit !found}}
        ' "$header_work/header.sam"

        samtools fastq -@ {threads} -F 0x900 -n {input.staged_bam:q} \
          > "$project_work/reads.fastq" 2> {log:q}
        test -s "$project_work/reads.fastq"
        fastq_uncompressed_size=$(stat -c '%s' "$project_work/reads.fastq")
        bgzip -@ {threads} -c "$project_work/reads.fastq" \
          > "$project_work/reads.fastq.gz" 2>> {log:q}
        gzip -t "$project_work/reads.fastq.gz" 2>> {log:q}
        mv "$project_work/reads.fastq.gz" {output.reads:q}

        {{
          printf 'metric\tvalue\n'
          printf 'sample_id\t%s\n' {wildcards.sample:q}
          printf 'cohort\t%s\n' {params.cohort:q}
          printf 'archive_bam\t%s\n' {params.archive_bam:q}
          printf 'staged_bam\t%s\n' {input.staged_bam:q}
          printf 'sq_count\t%s\n' "$(wc -l < "$header_work/actual.sq")"
          printf 'sq_sha256\t%s\n' "$(sha256sum "$header_work/actual.sq" | cut -d' ' -f1)"
          printf 'fastq_uncompressed_size_bytes\t%s\n' "$fastq_uncompressed_size"
          printf 'fastq_size_bytes\t%s\n' "$(stat -c '%s' {output.reads:q})"
        }} > {output.header_qc:q}
        """


rule run_local_v2_rescue:
    input:
        plan=lambda wildcards: current_plan(),
        reads=rescue_reads,
        assembly=rescue_assembly
    output:
        result=directory(RESCUE_RESULT)
    params:
        script=RESCUE_SCRIPT,
        results_dir=lambda wildcards, output: str(Path(output.result).parent),
        options=RESCUE_OPTIONS
    log:
        f"{OUT}/logs/rescue/{{rescue}}.log"
    benchmark:
        f"{OUT}/benchmarks/rescue/{{rescue}}.tsv"
    conda:
        RESCUE_ENV
    threads:
        int(RESCUE_CFG.get("threads", 16))
    resources:
        mem_mb=int(RESCUE_CFG.get("mem_mb", 64000)),
        runtime=int(RESCUE_CFG.get("runtime", 720)),
        disk_mb=int(RESCUE_CFG.get("disk_mb", 100000))
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.results_dir:q} "$(dirname {log:q})"
        python {params.script:q} run \
          --plan {input.plan:q} \
          --rescue-id {wildcards.rescue:q} \
          --output-dir {params.results_dir:q} \
          --threads {threads} \
          {params.options} > {log:q} 2>&1
        test -s {output.result:q}/qc.tsv
        test -s {output.result:q}/run_status.json
        """


rule summarize_local_v2_rescues:
    input:
        plan=lambda wildcards: current_plan(),
        results=ready_results
    output:
        summary=RESCUE_SUMMARY,
        fasta=RESCUE_FASTA,
        graph_inputs=GRAPH_INPUTS
    params:
        script=RESCUE_SCRIPT,
        results_dir=RESCUE_RESULTS
    log:
        f"{OUT}/logs/summarize_rescues.log"
    conda:
        RESCUE_ENV
    resources:
        mem_mb=4000,
        runtime=60
    shell:
        r"""
        python {params.script:q} summarize \
          --plan {input.plan:q} \
          --results-dir {params.results_dir:q} \
          --summary {output.summary:q} \
          --fasta {output.fasta:q} \
          --graph-inputs {output.graph_inputs:q} > {log:q} 2>&1
        """


rule build_local_sv_graph_v2:
    input:
        frozen_graph=FROZEN_GRAPH,
        frozen_summary=FROZEN_SUMMARY,
        frozen_inventory=FROZEN_INVENTORY,
        rescue_fasta=RESCUE_FASTA,
        graph_inputs=GRAPH_INPUTS,
        rescue_summary=RESCUE_SUMMARY,
        plan=PLAN,
        read_manifest=READ_MANIFEST,
        header_template=HEADER_TEMPLATE,
        graph_environment=GRAPH_ENV_INPUT
    output:
        graph=GRAPH,
        summary=GRAPH_SUMMARY,
        versions=GRAPH_VERSIONS,
        provenance=GRAPH_PROVENANCE
    params:
        script=GRAPH_SCRIPT,
        preset=GRAPH_CFG.get("minigraph_preset", "ggs")
    log:
        f"{OUT}/logs/build_local_sv_graph_v2.log"
    benchmark:
        f"{OUT}/benchmarks/build_local_sv_graph_v2.tsv"
    conda:
        GRAPH_ENV
    threads:
        int(GRAPH_CFG.get("threads", 16))
    resources:
        mem_mb=int(GRAPH_CFG.get("mem_mb", 240000)),
        runtime=int(GRAPH_CFG.get("runtime", 2880)),
        disk_mb=int(GRAPH_CFG.get("disk_mb", 20000))
    shell:
        r"""
        python {params.script:q} \
          --frozen-graph {input.frozen_graph:q} \
          --frozen-summary {input.frozen_summary:q} \
          --frozen-inventory {input.frozen_inventory:q} \
          --rescue-fasta {input.rescue_fasta:q} \
          --graph-inputs {input.graph_inputs:q} \
          --rescue-summary {input.rescue_summary:q} \
          --plan {input.plan:q} \
          --read-manifest {input.read_manifest:q} \
          --header-template {input.header_template:q} \
          --graph-environment {input.graph_environment:q} \
          --output-graph {output.graph:q} \
          --output-summary {output.summary:q} \
          --output-versions {output.versions:q} \
          --output-provenance {output.provenance:q} \
          --preset {params.preset:q} \
          --threads {threads} \
          --scratch-dir "${{TMPDIR:-/tmp}}" \
          --log {log:q}
        """


rule call_local_v2_graph_bubbles:
    input:
        graph=GRAPH
    output:
        bubbles=GRAPH_BUBBLES
    log:
        f"{OUT}/logs/call_local_v2_graph_bubbles.log"
    benchmark:
        f"{OUT}/benchmarks/call_local_v2_graph_bubbles.tsv"
    conda:
        GRAPH_ENV
    resources:
        mem_mb=int(GRAPH_CFG.get("mem_mb", 240000)),
        runtime=int(GRAPH_CFG.get("runtime", 2880)),
        disk_mb=int(GRAPH_CFG.get("disk_mb", 20000))
    shell:
        r"""
        mkdir -p "$(dirname {output.bubbles:q})" "$(dirname {log:q})"
        work=$(mktemp -d -p "${{TMPDIR:-/tmp}}" "local-v2-bubbles.XXXXXX")
        trap 'rm -rf -- "$work"' EXIT
        gfatools bubble {input.graph:q} > "$work/bubbles.tsv" 2> {log:q}
        mv "$work/bubbles.tsv" {output.bubbles:q}
        """


rule qc_local_v2_graph_structure:
    input:
        frozen_graph=FROZEN_GRAPH,
        graph=GRAPH,
        graph_summary=GRAPH_SUMMARY,
        frozen_inventory=FROZEN_INVENTORY,
        provenance=GRAPH_PROVENANCE
    output:
        integrity=QC_INTEGRITY,
        delta_segments=QC_DELTA_SEGMENTS,
        delta_links=QC_DELTA_LINKS,
        metrics=QC_GRAPH_METRICS
    params:
        script=QC_SCRIPT,
        build_log=f"{OUT}/logs/build_local_sv_graph_v2.log"
    log:
        f"{OUT}/logs/qc/graph_structure.log"
    conda:
        GRAPH_ENV
    threads: 1
    resources:
        mem_mb=int(QC_CFG.get("graph_mem_mb", 32000)),
        runtime=int(QC_CFG.get("graph_runtime", 360))
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.integrity:q})" "$(dirname {log:q})"
        python {params.script:q} graph \
          --frozen-graph {input.frozen_graph:q} \
          --rescued-graph {input.graph:q} \
          --graph-summary {input.graph_summary:q} \
          --frozen-inventory {input.frozen_inventory:q} \
          --provenance {input.provenance:q} \
          --build-log {params.build_log:q} \
          --integrity {output.integrity:q} \
          --delta-segments {output.delta_segments:q} \
          --delta-links {output.delta_links:q} \
          --metrics {output.metrics:q} > {log:q} 2>&1
        test -s {output.integrity:q}
        """


rule report_local_v2_qc:
    input:
        plan=PLAN,
        rescue_summary=RESCUE_SUMMARY,
        graph_inputs=GRAPH_INPUTS,
        catalog=CATALOG,
        catalog_evidence=CATALOG_EVIDENCE,
        assembly_manifest=ASSEMBLY_MANIFEST,
        integrity=QC_INTEGRITY,
        frozen_bubbles=FROZEN_BUBBLES,
        rescued_bubbles=GRAPH_BUBBLES
    output:
        bubble_delta=QC_BUBBLE_DELTA,
        rescues=QC_RESCUES,
        cohort=QC_COHORT,
        report=QC_REPORT
    params:
        script=QC_SCRIPT,
        min_sv=int(DISCOVERY_SCREEN.get("min_sv_size", 50)),
        options=QC_REPORT_OPTIONS
    log:
        f"{OUT}/logs/qc/report.log"
    conda:
        GRAPH_ENV
    threads: 1
    resources:
        mem_mb=int(QC_CFG.get("report_mem_mb", 16000)),
        runtime=int(QC_CFG.get("report_runtime", 360))
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.report:q})" "$(dirname {log:q})"
        python {params.script:q} report \
          --plan {input.plan:q} \
          --rescue-summary {input.rescue_summary:q} \
          --graph-inputs {input.graph_inputs:q} \
          --catalog {input.catalog:q} \
          --catalog-evidence {input.catalog_evidence:q} \
          --assembly-manifest {input.assembly_manifest:q} \
          --integrity {input.integrity:q} \
          --frozen-bubbles {input.frozen_bubbles:q} \
          --rescued-bubbles {input.rescued_bubbles:q} \
          --bubble-delta {output.bubble_delta:q} \
          --rescue-qc {output.rescues:q} \
          --cohort-qc {output.cohort:q} \
          --report {output.report:q} \
          --min-sv-size {params.min_sv} \
          {params.options} > {log:q} 2>&1
        test -s {output.report:q}
        """


wildcard_constraints:
    sample=r"[A-Za-z0-9_.-]+",
    rescue=r"[A-Za-z0-9_.-]+"
