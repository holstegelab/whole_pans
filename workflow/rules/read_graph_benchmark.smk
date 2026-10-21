"""Standalone linear/pangenome read-mapping benchmark."""

configfile: workflow.source_path("../../config/read_graph_benchmark.yaml")


import csv
import re
import shlex
from pathlib import Path


PROJECT_ROOT = Path(workflow.current_basedir).parents[1]
workdir: str(PROJECT_ROOT)
shell.executable("/bin/bash")

CFG_PATH = str(workflow.source_path("../../config/read_graph_benchmark.yaml"))
OUT = config["results"]
GRAPHS = {row["graph_id"]: row for row in config["graph_registry"]}
if len(GRAPHS) != len(config["graph_registry"]):
    raise WorkflowError("Duplicate graph IDs are present; run read_graph_preflight for details")

PY_ENV = str(workflow.source_path("../envs/read_graph_benchmark_python.yaml"))
LONG_ENV = str(workflow.source_path("../envs/read_graph_benchmark_long.yaml"))
VG_ENV = str(workflow.source_path("../envs/read_graph_benchmark_vg.yaml"))
REGISTRY_SCRIPT = str(workflow.source_path("../scripts/read_graph_registry.py"))
TMEM_SELECT_SCRIPT = str(workflow.source_path("../scripts/select_tmem_pilot.py"))
READ_SELECT_SCRIPT = str(workflow.source_path("../scripts/select_fastq_reads.py"))
METRICS_SCRIPT = str(workflow.source_path("../scripts/mapping_metrics.py"))
SIMULATE_SCRIPT = str(workflow.source_path("../scripts/simulate_hifi.py"))
RESCUE_SCRIPT = str(workflow.source_path("../scripts/rescue_capture.py"))
TMEM_SCRIPT = str(workflow.source_path("../scripts/tmem_evidence.py"))
SV_SCRIPT = str(workflow.source_path("../scripts/normalize_sv_catalog.py"))
REPORT_SCRIPT = str(workflow.source_path("../scripts/read_graph_benchmark_report.py"))
ARCHIVE_SCRIPT = str(workflow.source_path("../scripts/archive_hifi_reads.py"))
SV_PILOT_SCRIPT = str(workflow.source_path("../scripts/sv_short_read_pilot.py"))
SV_GAF_SCRIPT = str(workflow.source_path("../scripts/sv_gaf_metrics.py"))
SV_TARGET_SCRIPT = str(workflow.source_path("../scripts/tmem_whole_graph_target.py"))
SV_REPORT_SCRIPT = str(workflow.source_path("../scripts/sv_short_read_report.py"))

TABLES = f"{OUT}/tables"
LOGS = f"{OUT}/logs"
BENCHMARKS = f"{OUT}/benchmarks"
REPORTS = f"{OUT}/report"
GRAPH_REGISTRY = f"{TABLES}/graph_registry.tsv"
GRAPH_PREFLIGHT = f"{TABLES}/graph_preflight.tsv"
GRAPH_PREFLIGHT_REPORT = f"{REPORTS}/graph_preflight.md"
TMEM_MANIFEST = f"{OUT}/manifests/tmem_pilot_samples.tsv"
SAMPLES = f"{TABLES}/samples.tsv"
MEMBERSHIP = f"{TABLES}/sample_graph_membership.tsv"
LEAKAGE = f"{TABLES}/leakage_audit.tsv"
LONG_GATE = f"{OUT}/preflight/long_reads.ready"
SHORT_GATE = f"{OUT}/preflight/tmem_short_reads.ready"
RESCUE_GATE = f"{OUT}/preflight/v2_rescue.ready"
MAPPING_METRICS = f"{TABLES}/mapping_metrics_per_sample.tsv"
PAIRED_CONTRASTS = f"{TABLES}/paired_graph_contrasts.tsv"
SIM_ACCURACY = f"{TABLES}/simulated_mapping_accuracy.tsv"
RESCUE_CAPTURE = f"{TABLES}/v2_rescue_capture.tsv"
TMEM_EVIDENCE = f"{TABLES}/tmem_genotype_evidence.tsv"
NORMALIZED_SV = f"{TABLES}/normalized_sv_catalog.tsv"
SV_CLUSTERS = f"{TABLES}/normalized_sv_clusters.tsv"
TOOL_VERSIONS = f"{TABLES}/tool_versions.tsv"
PROVENANCE = f"{TABLES}/provenance.tsv"
FINAL_REPORT = f"{REPORTS}/read_graph_benchmark_report.md"

SV_CFG = config.get("sv_short_read_pilot", {})
SV_ENABLED = bool(SV_CFG.get("enabled", False))
SV_OUT = f"{OUT}/sv_short_reads"
SV_TABLES = f"{SV_OUT}/tables"
SV_LOGS = f"{SV_OUT}/logs"
SV_BENCHMARKS = f"{SV_OUT}/benchmarks"
SV_SELECTED = f"{SV_OUT}/manifests/selected_samples.tsv"
SV_MEMBERSHIP = f"{SV_TABLES}/sample_membership.tsv"
SV_INDEX_PREFLIGHT_JOB = f"{SV_OUT}/index_preflight/{{graph}}.tsv"
SV_INDEX_PREFLIGHT = f"{SV_TABLES}/index_preflight.tsv"
SV_INDEX_DIR = f"{SV_OUT}/indexes/vg-{SV_CFG.get('vg_version', 'unknown')}/{{graph}}"
SV_INDEX_MANIFEST = f"{SV_OUT}/index_manifests/{{graph}}.tsv"
SV_INDEX_VALIDATION = f"{SV_OUT}/index_manifests/{{graph}}.validated"
SV_FASTQ_QC = f"{SV_OUT}/fastq_qc/{{sample}}.tsv"
SV_SMOKE_R1 = f"{SV_OUT}/smoke_reads/{{sample}}.R1.fastq.gz"
SV_SMOKE_R2 = f"{SV_OUT}/smoke_reads/{{sample}}.R2.fastq.gz"
SV_SMOKE_STATS = f"{SV_OUT}/smoke_reads/{{sample}}.selection.tsv"
SV_QUERY_FASTA = f"{SV_OUT}/tmem/allele_queries.fa"
SV_QUERY_MANIFEST = f"{SV_OUT}/tmem/allele_queries.tsv"
SV_TARGET_GAF = f"{SV_OUT}/tmem/targets/{{graph}}.alleles.gaf"
SV_TARGET_JOB = f"{SV_OUT}/tmem/targets/{{graph}}.tsv"
SV_TARGET_PRESENCE = f"{SV_TABLES}/tmem_target_presence.tsv"
SV_CORE_METRIC = f"{SV_OUT}/mapping/{{mode}}/core_metrics/{{sample}}.{{graph}}.tsv"
SV_RAW_METRIC = f"{SV_OUT}/mapping/{{mode}}/raw_metrics/{{sample}}.{{graph}}.tsv"
SV_JOB_EVIDENCE = f"{SV_OUT}/tmem/evidence/{{mode}}/{{sample}}.{{graph}}.tsv"
SV_SMOKE_METRICS = f"{SV_TABLES}/smoke_mapping_metrics_per_sample.tsv"
SV_SMOKE_CONTRASTS = f"{SV_TABLES}/smoke_paired_graph_contrasts.tsv"
SV_SMOKE_EVIDENCE = f"{SV_TABLES}/smoke_tmem_evidence.tsv"
SV_FULL_METRICS = f"{SV_TABLES}/mapping_metrics_per_sample.tsv"
SV_FULL_CONTRASTS = f"{SV_TABLES}/paired_graph_contrasts.tsv"
SV_FULL_EVIDENCE = f"{SV_TABLES}/tmem_evidence.tsv"
SV_TOOL_VERSIONS = f"{SV_TABLES}/tool_versions.tsv"
SV_PROVENANCE = f"{SV_TABLES}/provenance.tsv"
SV_SMOKE_READY = f"{SV_OUT}/preflight/smoke.ready"
SV_SMOKE_REPORT = f"{SV_OUT}/report/sv_short_read_smoke.md"
SV_FULL_REPORT = f"{SV_OUT}/report/sv_short_read_pilot.md"

SAMPLE_CONFIG = {
    str(row["sample_id"]): row for row in config.get("samples", [])
    if row.get("include", True)
}
LONG_SAMPLES = sorted(
    sample for sample, row in SAMPLE_CONFIG.items()
    if row.get("data_type") in {"hifi", "ont"} and row.get("long_read")
)
ARCHIVE_MEMBER_SAMPLES = sorted(
    sample for sample, row in SAMPLE_CONFIG.items()
    if row.get("stage_from_archive", False)
)
HEADLINE_GRAPHS = [
    row["graph_id"] for row in config["graph_registry"]
    if row.get("comparison_scope") == "headline" and "long_reads" in row.get("mapping_capabilities", [])
]
LONG_JOBS = [(sample, "minigraph", graph) for sample in LONG_SAMPLES for graph in HEADLINE_GRAPHS]
LONG_JOBS += [(sample, "minimap2", "chm13_linear") for sample in LONG_SAMPLES]
SIM_SAMPLES = sorted(
    sample for sample, row in SAMPLE_CONFIG.items()
    if config.get("simulation", {}).get("enabled", False)
    and row.get("simulation", False)
    and row.get("truth_assembly_hap1") and row.get("truth_assembly_hap2")
)
SIM_JOBS = [(sample, graph) for sample in SIM_SAMPLES for graph in HEADLINE_GRAPHS]

LONG_ALIGNMENT = f"{OUT}/long_reads/alignments/{{sample}}/{{mapper}}.{{graph}}.gaf"
LONG_JOB_METRIC = f"{OUT}/long_reads/metrics/{{sample}}/{{mapper}}.{{graph}}.tsv"
SIM_FASTQ = f"{OUT}/simulation/reads/{{sample}}.simulated.hifi.fastq.gz"
SIM_ALIGNMENT = f"{OUT}/simulation/alignments/{{sample}}/minigraph.{{graph}}.gaf"
SIM_JOB_ACCURACY = f"{OUT}/simulation/accuracy/{{sample}}/minigraph.{{graph}}.tsv"
TMEM_GAF = f"{OUT}/tmem/alignments/{{sample}}/{{graph}}.gaf"
TMEM_JOB_EVIDENCE = f"{OUT}/tmem/evidence/{{sample}}/{{graph}}.tsv"

METADATA_RESOURCES = config.get("resources", {})
LONG_CFG = config["long_reads"]
ARCHIVE_CFG = config.get("long_read_archive_staging", {})
TMEM_CFG = config["tmem_pilot"]
SIM_CFG = config.get("simulation", {})
RESCUE_CFG = config["v2_rescue_validation"]
SV_CATALOG_CFG = config["sv_catalog"]
TMEM_GRAPH_IDS = tuple(str(graph) for graph in TMEM_CFG["graph_ids"])
TMEM_VG_VERSION = str(TMEM_CFG["vg_version"])
TMEM_INDEX_DIR = f"{OUT}/tmem/indexes/vg-{TMEM_VG_VERSION}"
TMEM_INDEX_PATTERNS = {
    "dist": f"{TMEM_INDEX_DIR}/{{graph}}.dist",
    "minimizer": f"{TMEM_INDEX_DIR}/{{graph}}.shortread.withzip.min",
    "zipcodes": f"{TMEM_INDEX_DIR}/{{graph}}.shortread.zipcodes",
}
TMEM_GRAPH_REGEX = "(?:" + "|".join(map(re.escape, TMEM_GRAPH_IDS)) + ")"
SV_GRAPH_IDS = tuple(str(graph) for graph in SV_CFG.get("graph_ids", [])) if SV_ENABLED else tuple()
SV_GRAPH_REGEX = "(?:" + "|".join(map(re.escape, SV_GRAPH_IDS)) + ")" if SV_GRAPH_IDS else "(?!)"
RESCUE_GRAPHS = (
    RESCUE_CFG.get("v1_graph_id", "local_v1"),
    RESCUE_CFG.get("v2_graph_id", "local_v2"),
)

if ARCHIVE_MEMBER_SAMPLES and not ARCHIVE_CFG:
    raise WorkflowError("stage_from_archive samples require long_read_archive_staging config")

for graph in TMEM_GRAPH_IDS:
    if graph not in GRAPHS:
        raise WorkflowError(f"Unknown TMEM graph: {graph}")
    configured_indexes = GRAPHS[graph].get("indexes", {})
    for index_name, pattern in TMEM_INDEX_PATTERNS.items():
        expected = pattern.format(graph=graph)
        if configured_indexes.get(index_name) != expected:
            raise WorkflowError(
                f"{graph}: indexes.{index_name} must be {expected}"
            )

if SV_ENABLED:
    if not SV_GRAPH_IDS:
        raise WorkflowError("sv_short_read_pilot.graph_ids is empty")
    if SV_CFG.get("baseline_graph") not in SV_GRAPH_IDS:
        raise WorkflowError("sv_short_read_pilot.baseline_graph must be in graph_ids")
    if SV_CFG.get("mapper") != "vg_giraffe":
        raise WorkflowError("sv_short_read_pilot currently supports mapper=vg_giraffe only")
    if SV_CFG.get("index_mode") != "rgfa_artificial_path_cover":
        raise WorkflowError("SV-rGFA Giraffe indexing mode must be rgfa_artificial_path_cover")
    if SV_CFG.get("mode", "full") not in {"smoke", "full"}:
        raise WorkflowError("sv_short_read_pilot.mode must be smoke or full")
    for graph in SV_GRAPH_IDS:
        if graph not in GRAPHS:
            raise WorkflowError(f"Unknown sv_short_read_pilot graph: {graph}")
        if GRAPHS[graph].get("graph_type") != "minigraph_rgfa":
            raise WorkflowError(f"{graph}: SV short-read targets must be minigraph_rgfa")

ARCHIVE_FASTQ = ARCHIVE_CFG.get("fastq_pattern", "")
ARCHIVE_READ_MANIFEST = f"{OUT}/staging/hifi/hifi_archive_reads.tsv"
ARCHIVE_STAGED = f"{OUT}/staging/hifi/archive/staged"
ARCHIVE_BAM = f"{OUT}/staging/hifi/compute_bams/{{sample}}.bam"
ARCHIVE_HEADER_QC = f"{OUT}/staging/hifi/header_qc/{{sample}}.header_qc.tsv"
ARCHIVE_SAMPLE_REGEX = "(?:" + "|".join(map(re.escape, ARCHIVE_MEMBER_SAMPLES)) + ")"
ARCHIVE_FASTQS = [SAMPLE_CONFIG[sample]["long_read"] for sample in ARCHIVE_MEMBER_SAMPLES]

for sample in ARCHIVE_MEMBER_SAMPLES:
    expected_fastq = ARCHIVE_FASTQ.format(sample=sample)
    if SAMPLE_CONFIG[sample]["long_read"] != expected_fastq:
        raise WorkflowError(
            f"{sample}: long_read must match long_read_archive_staging.fastq_pattern: {expected_fastq}"
        )
    if not SAMPLE_CONFIG[sample].get("archive_cohort"):
        raise WorkflowError(f"{sample}: stage_from_archive requires archive_cohort")


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


SV_SOURCE_MANIFEST_PATH = (PROJECT_ROOT / TMEM_MANIFEST).resolve()
SV_SELECTION_ROLES = tuple(str(role) for role in SV_CFG.get("sample_roles", []))
if SV_ENABLED and SV_CFG.get("include_matched_controls", False):
    SV_SELECTION_ROLES = tuple(dict.fromkeys(SV_SELECTION_ROLES + ("tmem_matched_control",)))
if SV_ENABLED:
    if not SV_SOURCE_MANIFEST_PATH.is_file():
        raise WorkflowError(
            "The frozen TMEM pilot manifest is required and will not be regenerated for the "
            f"SV short-read pilot: {SV_SOURCE_MANIFEST_PATH}"
        )
    SV_FROZEN_ROWS = [
        row for row in read_tsv(SV_SOURCE_MANIFEST_PATH)
        if row.get("role") in SV_SELECTION_ROLES
    ]
else:
    SV_FROZEN_ROWS = []
SV_SAMPLE_ROWS = {row["sample_id"]: row for row in SV_FROZEN_ROWS}
if len(SV_SAMPLE_ROWS) != len(SV_FROZEN_ROWS):
    raise WorkflowError("Duplicate selected sample IDs in the frozen TMEM pilot manifest")
SV_SAMPLES = tuple(sorted(SV_SAMPLE_ROWS))
SV_SAMPLE_REGEX = "(?:" + "|".join(map(re.escape, SV_SAMPLES)) + ")" if SV_SAMPLES else "(?!)"
SV_JOBS = tuple((sample, graph) for sample in SV_SAMPLES for graph in SV_GRAPH_IDS)


def sv_sample_row(wildcards):
    if wildcards.sample not in SV_SAMPLE_ROWS:
        raise WorkflowError(f"Unsupported SV short-read pilot sample: {wildcards.sample}")
    return SV_SAMPLE_ROWS[wildcards.sample]


def sv_paths(wildcards, mate):
    return split_paths(sv_sample_row(wildcards)[mate])


def sv_qc_inputs(wildcards):
    return sv_paths(wildcards, "r1") + sv_paths(wildcards, "r2")


def sv_quoted_paths(wildcards, mate):
    return " ".join(
        "--input " + shlex.quote(path) for path in sv_paths(wildcards, mate)
    )


def sv_fastq_stats(wildcards):
    if wildcards.mode == "smoke":
        return SV_SMOKE_STATS.format(sample=wildcards.sample)
    if wildcards.mode == "full":
        return SV_FASTQ_QC.format(sample=wildcards.sample)
    raise WorkflowError(f"Unsupported SV short-read run mode: {wildcards.mode}")


def sv_mapping_reads(wildcards):
    if wildcards.mode == "smoke":
        return [
            SV_SMOKE_R1.format(sample=wildcards.sample),
            SV_SMOKE_R2.format(sample=wildcards.sample),
        ]
    if wildcards.mode == "full":
        return sv_paths(wildcards, "r1") + sv_paths(wildcards, "r2")
    raise WorkflowError(f"Unsupported SV short-read run mode: {wildcards.mode}")


def sv_mapping_fastq_args(wildcards, mate):
    if wildcards.mode == "smoke":
        pattern = SV_SMOKE_R1 if mate == "r1" else SV_SMOKE_R2
        values = [pattern.format(sample=wildcards.sample)]
    else:
        values = sv_paths(wildcards, mate)
    return " ".join("--input " + shlex.quote(value) for value in values)


def sv_full_gate(wildcards):
    if wildcards.mode == "full" and SV_CFG.get("require_smoke_before_full", True):
        return [SV_SMOKE_READY]
    return []


def sv_metric_inputs(mode):
    return [
        SV_RAW_METRIC.format(mode=mode, sample=sample, graph=graph)
        for sample, graph in SV_JOBS
    ]


def sv_evidence_inputs(mode):
    return [
        SV_JOB_EVIDENCE.format(mode=mode, sample=sample, graph=graph)
        for sample, graph in SV_JOBS
    ]


def sv_index_file(wildcards, suffix):
    directory = SV_INDEX_DIR.format(graph=wildcards.graph)
    return f"{directory}/{wildcards.graph}{suffix}"


def include_inputs(flag, paths):
    return paths if flag else []


def sample_coverage(wildcards):
    return str(SAMPLE_CONFIG[wildcards.sample].get(
        "selected_coverage", LONG_CFG.get("default_selected_coverage", 20)
    ))


def sample_reads(wildcards):
    return SAMPLE_CONFIG[wildcards.sample]["long_read"]


def archive_bam(wildcards):
    sample = wildcards.sample
    return ARCHIVE_CFG["archive_bam_pattern"].format(
        archive_root=ARCHIVE_CFG["archive_root"],
        cohort=SAMPLE_CONFIG[sample]["archive_cohort"],
        sample=sample,
    )


def graph_path(wildcards):
    if wildcards.graph not in GRAPHS:
        raise WorkflowError(f"Unknown graph: {wildcards.graph}")
    return GRAPHS[wildcards.graph]["path"]


def validate_long_job(wildcards):
    job = (wildcards.sample, wildcards.mapper, wildcards.graph)
    if job not in LONG_JOBS:
        raise WorkflowError(f"Unsupported long-read mapping combination: {job}")
    return GRAPHS[wildcards.graph]["path"]


def long_mapping_gate(wildcards):
    # Rescue-only runs need only V1/V2 preflight. Other headline mappings keep
    # the stricter gate that validates every headline graph before they start.
    if wildcards.mapper == "minigraph" and wildcards.graph in RESCUE_GRAPHS:
        return RESCUE_GATE
    return LONG_GATE


def minimap2_preset(wildcards):
    return "map-ont" if SAMPLE_CONFIG[wildcards.sample]["data_type"] == "ont" else LONG_CFG.get("minimap2_preset", "map-hifi")


def tmem_rows(wildcards=None):
    path = str(checkpoints.select_tmem_pilot.get().output.manifest)
    return read_tsv(path)


def tmem_row(wildcards):
    matches = [row for row in tmem_rows(wildcards) if row["sample_id"] == wildcards.sample]
    if len(matches) != 1:
        raise WorkflowError(f"No unique TMEM manifest row for {wildcards.sample}")
    return matches[0]


def split_paths(value):
    return [path for path in str(value).split(",") if path]


def tmem_r1(wildcards):
    return split_paths(tmem_row(wildcards)["r1"])


def tmem_r2(wildcards):
    return split_paths(tmem_row(wildcards)["r2"])


def quoted_tmem_files(wildcards, mate):
    return " ".join(shlex.quote(path) for path in split_paths(tmem_row(wildcards)[mate]))


def tmem_targets(wildcards):
    samples = [row["sample_id"] for row in tmem_rows(wildcards)]
    return expand(TMEM_JOB_EVIDENCE, sample=samples, graph=TMEM_CFG["graph_ids"])


def long_metric_inputs():
    return [LONG_JOB_METRIC.format(sample=sample, mapper=mapper, graph=graph) for sample, mapper, graph in LONG_JOBS]


def sim_accuracy_inputs():
    return [SIM_JOB_ACCURACY.format(sample=sample, graph=graph) for sample, graph in SIM_JOBS]


def rescue_samples():
    samples = set()
    for event in RESCUE_CFG["events"]:
        samples.add(str(event["sample_id"]))
        samples.update(str(sample) for sample in event.get("independent_carriers", []))
        samples.update(str(sample) for sample in event.get("noncarrier_controls", []))
    missing = samples - set(LONG_SAMPLES)
    if missing:
        raise WorkflowError(
            "Rescue evidence samples must be enabled long-read sample entries: " + ",".join(sorted(missing))
        )
    return sorted(samples)


def rescue_alignment_inputs():
    return [
        LONG_ALIGNMENT.format(sample=sample, mapper="minigraph", graph=graph)
        for sample in rescue_samples()
        for graph in RESCUE_GRAPHS
    ]


def rescue_alignment_args():
    return " ".join(
        "--alignment " + shlex.quote(
            f"{graph}={sample}="
            f"{LONG_ALIGNMENT.format(sample=sample, mapper='minigraph', graph=graph)}"
        )
        for sample in rescue_samples()
        for graph in RESCUE_GRAPHS
    )


localrules: read_graph_benchmark_all, read_graph_preflight, prepare_long_read_member_controls, prepare_tmem_giraffe_indexes, prepare_sv_short_read_indexes, sv_short_read_smoke, sv_short_read_pilot, long_read_pilot, v2_rescue_validation, tmem_pilot, simulation_benchmark, unified_sv_catalog


rule read_graph_benchmark_all:
    input:
        GRAPH_REGISTRY,
        GRAPH_PREFLIGHT,
        SAMPLES,
        MEMBERSHIP,
        LEAKAGE,
        MAPPING_METRICS,
        PAIRED_CONTRASTS,
        SIM_ACCURACY,
        RESCUE_CAPTURE,
        NORMALIZED_SV,
        TMEM_EVIDENCE,
        TOOL_VERSIONS,
        PROVENANCE,
        FINAL_REPORT


rule read_graph_preflight:
    input:
        GRAPH_REGISTRY,
        GRAPH_PREFLIGHT,
        GRAPH_PREFLIGHT_REPORT,
        SAMPLES,
        MEMBERSHIP,
        LEAKAGE


rule prepare_long_read_member_controls:
    input:
        ARCHIVE_FASTQS,
        expand(ARCHIVE_HEADER_QC, sample=ARCHIVE_MEMBER_SAMPLES)


rule prepare_tmem_giraffe_indexes:
    input:
        expand(list(TMEM_INDEX_PATTERNS.values()), graph=TMEM_GRAPH_IDS)


rule prepare_sv_short_read_indexes:
    input:
        SV_INDEX_PREFLIGHT,
        SV_TOOL_VERSIONS,
        expand(SV_INDEX_DIR, graph=SV_GRAPH_IDS),
        expand(SV_INDEX_MANIFEST, graph=SV_GRAPH_IDS)


rule sv_short_read_smoke:
    input:
        SV_SELECTED,
        SV_MEMBERSHIP,
        SV_INDEX_PREFLIGHT,
        expand(SV_INDEX_MANIFEST, graph=SV_GRAPH_IDS),
        expand(SV_FASTQ_QC, sample=SV_SAMPLES),
        expand(SV_SMOKE_STATS, sample=SV_SAMPLES),
        SV_TARGET_PRESENCE,
        SV_SMOKE_METRICS,
        SV_SMOKE_CONTRASTS,
        SV_SMOKE_EVIDENCE,
        SV_SMOKE_READY,
        SV_TOOL_VERSIONS,
        SV_PROVENANCE,
        SV_SMOKE_REPORT


rule sv_short_read_pilot:
    input:
        SV_SELECTED,
        SV_MEMBERSHIP,
        SV_INDEX_PREFLIGHT,
        expand(SV_INDEX_MANIFEST, graph=SV_GRAPH_IDS),
        expand(SV_FASTQ_QC, sample=SV_SAMPLES),
        SV_TARGET_PRESENCE,
        SV_FULL_METRICS,
        SV_FULL_CONTRASTS,
        SV_FULL_EVIDENCE,
        SV_TOOL_VERSIONS,
        SV_PROVENANCE,
        SV_FULL_REPORT


rule long_read_pilot:
    input:
        MAPPING_METRICS,
        PAIRED_CONTRASTS


rule v2_rescue_validation:
    input:
        RESCUE_CAPTURE


rule tmem_pilot:
    input:
        tmem_targets,
        TMEM_EVIDENCE


rule simulation_benchmark:
    input:
        SIM_ACCURACY


rule unified_sv_catalog:
    input:
        NORMALIZED_SV,
        SV_CLUSTERS


rule build_tmem_giraffe_indexes:
    input:
        gbz=graph_path,
        vg_env=VG_ENV
    output:
        dist=TMEM_INDEX_PATTERNS["dist"],
        minimizer=TMEM_INDEX_PATTERNS["minimizer"],
        zipcodes=TMEM_INDEX_PATTERNS["zipcodes"]
    log:
        f"{LOGS}/tmem/{{graph}}.build_indexes.log"
    benchmark:
        f"{BENCHMARKS}/tmem/{{graph}}.build_indexes.tsv"
    conda:
        VG_ENV
    threads:
        int(TMEM_CFG.get("index_threads", TMEM_CFG.get("threads", 24)))
    resources:
        mem_mb=int(TMEM_CFG.get("index_mem_mb", TMEM_CFG.get("mem_mb", 32000))),
        runtime_min=int(TMEM_CFG.get("index_runtime_min", TMEM_CFG.get("runtime_min", 480))),
        disk_mb=int(TMEM_CFG.get("index_disk_mb", TMEM_CFG.get("disk_mb", 10000)))
    wildcard_constraints:
        graph=TMEM_GRAPH_REGEX
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.dist:q})" "$(dirname {log:q})"
        index_work=$(mktemp -d -p "$(dirname {output.dist:q})" ".{wildcards.graph}.indexes.XXXXXX")
        trap 'rm -rf -- "$index_work"' EXIT

        vg version > {log:q} 2>&1
        vg index -j "$index_work/graph.dist" {input.gbz:q} >> {log:q} 2>&1
        vg minimizer -p -t {threads} \
          -d "$index_work/graph.dist" \
          -o "$index_work/graph.shortread.withzip.min" \
          -z "$index_work/graph.shortread.zipcodes" \
          {input.gbz:q} >> {log:q} 2>&1

        test -s "$index_work/graph.dist"
        test -s "$index_work/graph.shortread.withzip.min"
        test -s "$index_work/graph.shortread.zipcodes"
        mv "$index_work/graph.dist" {output.dist:q}
        mv "$index_work/graph.shortread.withzip.min" {output.minimizer:q}
        mv "$index_work/graph.shortread.zipcodes" {output.zipcodes:q}
        trap - EXIT
        rm -rf -- "$index_work"
        """


rule select_sv_short_read_samples:
    input:
        config=CFG_PATH
    output:
        SV_SELECTED
    params:
        manifest=str(SV_SOURCE_MANIFEST_PATH)
    log:
        f"{SV_LOGS}/select_samples.log"
    benchmark:
        f"{SV_BENCHMARKS}/select_samples.tsv"
    conda:
        PY_ENV
    threads: 1
    resources:
        mem_mb=int(METADATA_RESOURCES.get("metadata_mem_mb", 4000)),
        runtime_min=int(METADATA_RESOURCES.get("metadata_runtime_min", 60))
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {SV_PILOT_SCRIPT:q} select \
          --config {input.config:q} --manifest {params.manifest:q} \
          --output {output:q} > {log:q} 2>&1
        """


rule sv_short_read_fastq_qc:
    input:
        selected=SV_SELECTED,
        fastqs=sv_qc_inputs
    output:
        SV_FASTQ_QC
    params:
        r1=lambda wildcards: " ".join("--r1 " + shlex.quote(value) for value in sv_paths(wildcards, "r1")),
        r2=lambda wildcards: " ".join("--r2 " + shlex.quote(value) for value in sv_paths(wildcards, "r2"))
    log:
        f"{SV_LOGS}/fastq_qc/{{sample}}.log"
    benchmark:
        f"{SV_BENCHMARKS}/fastq_qc/{{sample}}.tsv"
    conda:
        PY_ENV
    threads: 1
    resources:
        mem_mb=int(SV_CFG.get("metrics_mem_mb", 4000)),
        runtime_min=int(SV_CFG.get("runtime_min", 1440)),
        disk_mb=1000
    wildcard_constraints:
        sample=SV_SAMPLE_REGEX
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {SV_PILOT_SCRIPT:q} fastq-qc \
          --sample {wildcards.sample:q} {params.r1} {params.r2} \
          --output {output:q} > {log:q} 2>&1
        """


rule subsample_sv_short_reads:
    input:
        qc=SV_FASTQ_QC,
        fastqs=sv_qc_inputs
    output:
        r1=SV_SMOKE_R1,
        r2=SV_SMOKE_R2,
        stats=SV_SMOKE_STATS
    params:
        r1=lambda wildcards: " ".join("--r1 " + shlex.quote(value) for value in sv_paths(wildcards, "r1")),
        r2=lambda wildcards: " ".join("--r2 " + shlex.quote(value) for value in sv_paths(wildcards, "r2")),
        max_pairs=int(SV_CFG.get("smoke_max_pairs", 5000000)),
        seed=int(SV_CFG.get("subsampling_seed", 20260901))
    log:
        f"{SV_LOGS}/subsample/{{sample}}.log"
    benchmark:
        f"{SV_BENCHMARKS}/subsample/{{sample}}.tsv"
    conda:
        PY_ENV
    threads: 1
    resources:
        mem_mb=int(SV_CFG.get("metrics_mem_mb", 4000)),
        runtime_min=int(SV_CFG.get("runtime_min", 1440)),
        disk_mb=max(5000, int(SV_CFG.get("smoke_max_pairs", 5000000)) // 1000)
    wildcard_constraints:
        sample=SV_SAMPLE_REGEX
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.r1:q})" "$(dirname {log:q})"
        python {SV_PILOT_SCRIPT:q} subsample \
          --sample {wildcards.sample:q} {params.r1} {params.r2} \
          --qc {input.qc:q} --max-pairs {params.max_pairs} --seed {params.seed} \
          --output-r1 {output.r1:q} --output-r2 {output.r2:q} \
          --stats {output.stats:q} > {log:q} 2>&1
        """


rule sv_short_read_membership_audit:
    input:
        config=CFG_PATH,
        selected=SV_SELECTED
    output:
        SV_MEMBERSHIP
    log:
        f"{SV_LOGS}/sample_membership.log"
    benchmark:
        f"{SV_BENCHMARKS}/sample_membership.tsv"
    conda:
        PY_ENV
    threads: 1
    resources:
        mem_mb=int(METADATA_RESOURCES.get("metadata_mem_mb", 4000)),
        runtime_min=int(METADATA_RESOURCES.get("metadata_runtime_min", 60))
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {SV_PILOT_SCRIPT:q} membership \
          --config {input.config:q} --selected {input.selected:q} \
          --output {output:q} > {log:q} 2>&1
        """


rule sv_short_read_index_preflight_job:
    input:
        source=graph_path
    output:
        SV_INDEX_PREFLIGHT_JOB
    params:
        qc=lambda wildcards: GRAPHS[wildcards.graph].get("qc_status", ""),
        experimental=lambda wildcards: str(GRAPHS[wildcards.graph].get("experimental", False)).lower(),
        mode=SV_CFG.get("index_mode", "rgfa_artificial_path_cover")
    log:
        f"{SV_LOGS}/indexes/{{graph}}.preflight.log"
    benchmark:
        f"{SV_BENCHMARKS}/indexes/{{graph}}.preflight.tsv"
    conda:
        PY_ENV
    threads: 1
    resources:
        mem_mb=int(SV_CFG.get("metrics_mem_mb", 4000)),
        runtime_min=int(SV_CFG.get("index_runtime_min", 2880)),
        disk_mb=1000
    wildcard_constraints:
        graph=SV_GRAPH_REGEX
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {SV_PILOT_SCRIPT:q} index-preflight \
          --graph {wildcards.graph:q} --source {input.source:q} \
          --qc-status {params.qc:q} --experimental {params.experimental:q} \
          --indexing-mode {params.mode:q} --output {output:q} > {log:q} 2>&1
        """


rule merge_sv_short_read_index_preflight:
    input:
        expand(SV_INDEX_PREFLIGHT_JOB, graph=SV_GRAPH_IDS)
    output:
        SV_INDEX_PREFLIGHT
    params:
        inputs=lambda wildcards, input: " ".join("--input " + shlex.quote(str(value)) for value in input)
    log:
        f"{SV_LOGS}/indexes/merge_preflight.log"
    conda:
        PY_ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {SV_PILOT_SCRIPT:q} merge-index-preflight {params.inputs} \
          --output {output:q} > {log:q} 2>&1
        """


rule sv_short_read_tool_versions:
    input:
        environment=LONG_ENV
    output:
        SV_TOOL_VERSIONS
    log:
        f"{SV_LOGS}/tool_versions.log"
    benchmark:
        f"{SV_BENCHMARKS}/tool_versions.tsv"
    conda:
        LONG_ENV
    threads: 1
    resources:
        mem_mb=int(METADATA_RESOURCES.get("metadata_mem_mb", 4000)),
        runtime_min=int(METADATA_RESOURCES.get("metadata_runtime_min", 60))
    params:
        expected_vg_version=str(SV_CFG.get("vg_version", "1.76.1"))
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {SV_PILOT_SCRIPT:q} tools \
          --expected-vg-version {params.expected_vg_version:q} \
          --output {output:q} > {log:q} 2>&1
        """


rule build_sv_short_read_giraffe_index:
    input:
        source=graph_path,
        preflight=SV_INDEX_PREFLIGHT_JOB,
        tools=SV_TOOL_VERSIONS
    output:
        index_dir=directory(SV_INDEX_DIR),
        manifest=SV_INDEX_MANIFEST
    params:
        expected_version=str(SV_CFG.get("vg_version", "1.76.1")),
        indexing_mode=SV_CFG.get("index_mode", "rgfa_artificial_path_cover"),
        command=lambda wildcards: (
            f"vg autoindex --workflow sr-giraffe --gfa graph.gfa --prefix {wildcards.graph} "
            "--threads THREADS --tmp-dir TMPDIR"
        )
    log:
        f"{SV_LOGS}/indexes/{{graph}}.autoindex.log"
    benchmark:
        f"{SV_BENCHMARKS}/indexes/{{graph}}.autoindex.tsv"
    conda:
        VG_ENV
    threads:
        int(SV_CFG.get("index_threads", 32))
    resources:
        mem_mb=int(SV_CFG.get("index_mem_mb", 256000)),
        runtime_min=int(SV_CFG.get("index_runtime_min", 2880)),
        disk_mb=int(SV_CFG.get("index_disk_mb", 500000))
    wildcard_constraints:
        graph=SV_GRAPH_REGEX
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.index_dir:q})" "$(dirname {output.manifest:q})" "$(dirname {log:q})"
        work=$(mktemp -d -p "${{TMPDIR:-/tmp}}" "sv-giraffe-{wildcards.graph}.XXXXXX")
        stage=""
        trap 'rm -rf -- "$work"; if [ -n "$stage" ]; then rm -rf -- "$stage"; fi' EXIT
        mkdir -p "$work/products" "$work/autoindex-tmp"

        vg version > "$work/vg.version.txt" 2>&1
        cat "$work/vg.version.txt" > {log:q}
        grep -F {params.expected_version:q} "$work/vg.version.txt" >/dev/null
        vg autoindex --help >> {log:q} 2>&1 || true
        vg giraffe --help >> {log:q} 2>&1 || true

        case {input.source:q} in
          *.gz) gzip -cd -- {input.source:q} > "$work/graph.gfa" ;;
          *) cp -- {input.source:q} "$work/graph.gfa" ;;
        esac
        test -s "$work/graph.gfa"

        vg autoindex \
          --workflow sr-giraffe \
          --gfa "$work/graph.gfa" \
          --prefix "$work/products/{wildcards.graph}" \
          --threads {threads} \
          --tmp-dir "$work/autoindex-tmp" >> {log:q} 2>&1

        gbz="$work/products/{wildcards.graph}.giraffe.gbz"
        dist="$work/products/{wildcards.graph}.dist"
        minimizer="$work/products/{wildcards.graph}.shortread.withzip.min"
        zipcodes="$work/products/{wildcards.graph}.shortread.zipcodes"
        for product in "$gbz" "$dist" "$minimizer" "$zipcodes"; do
          test -s "$product"
          vg describe "$product" >> {log:q} 2>&1
        done
        vg validate "$gbz" >> {log:q} 2>&1

        final_dir={output.index_dir:q}
        if [ -e "$final_dir" ]; then
          echo "Refusing to replace an existing index directory: $final_dir" >> {log:q}
          exit 1
        fi
        stage=$(mktemp -d -p "$(dirname "$final_dir")" ".{wildcards.graph}.stage.XXXXXX")
        mv -- "$gbz" "$dist" "$minimizer" "$zipcodes" "$stage/"
        cp -- "$work/vg.version.txt" "$stage/vg.version.txt"
        python {SV_PILOT_SCRIPT:q} index-manifest \
          --graph {wildcards.graph:q} --source {input.source:q} \
          --preflight {input.preflight:q} --index-dir "$stage" \
          --reported-index-dir "$final_dir" --prefix-name {wildcards.graph:q} \
          --vg-version-file "$stage/vg.version.txt" --log {log:q} \
          --indexing-mode {params.indexing_mode:q} --command {params.command:q} \
          --output "$stage/completion_manifest.tsv" >> {log:q} 2>&1
        mv -- "$stage" "$final_dir"
        stage=""
        cp -- "$final_dir/completion_manifest.tsv" {output.manifest:q}.tmp
        mv -- {output.manifest:q}.tmp {output.manifest:q}
        trap - EXIT
        rm -rf -- "$work"
        """


rule validate_sv_short_read_giraffe_index:
    input:
        source=graph_path,
        manifest=SV_INDEX_MANIFEST
    output:
        temp(SV_INDEX_VALIDATION)
    log:
        f"{SV_LOGS}/indexes/{{graph}}.validate_reuse.log"
    benchmark:
        f"{SV_BENCHMARKS}/indexes/{{graph}}.validate_reuse.tsv"
    conda:
        PY_ENV
    threads: 1
    resources:
        mem_mb=int(SV_CFG.get("metrics_mem_mb", 4000)),
        runtime_min=int(SV_CFG.get("metrics_runtime_min", 120)),
        disk_mb=1000
    wildcard_constraints:
        graph=SV_GRAPH_REGEX
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {SV_PILOT_SCRIPT:q} validate-index \
          --graph {wildcards.graph:q} --source {input.source:q} \
          --manifest {input.manifest:q} --output {output:q} \
          > {log:q} 2>&1
        """


rule build_tmem_whole_graph_queries:
    input:
        vcf=SV_CFG.get("regional_graph_vcf", ""),
        targets=SV_CFG.get("target_insertions", ""),
        reference=SV_CFG.get("regional_reference_fasta", "")
    output:
        fasta=SV_QUERY_FASTA,
        manifest=SV_QUERY_MANIFEST
    params:
        vcf_position_offset=int(SV_CFG.get("regional_vcf_position_offset", 150000)),
        flank_bp=int(SV_CFG.get("target_flank_bp", 1000))
    log:
        f"{SV_LOGS}/tmem/build_allele_queries.log"
    benchmark:
        f"{SV_BENCHMARKS}/tmem/build_allele_queries.tsv"
    conda:
        PY_ENV
    threads: 1
    resources:
        mem_mb=int(SV_CFG.get("metrics_mem_mb", 4000)),
        runtime_min=int(SV_CFG.get("metrics_runtime_min", 120))
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.fasta:q})" "$(dirname {log:q})"
        python {SV_TARGET_SCRIPT:q} queries \
          --vcf {input.vcf:q} --targets {input.targets:q} --reference {input.reference:q} \
          --vcf-position-offset {params.vcf_position_offset} \
          --flank-bp {params.flank_bp} \
          --output-fasta {output.fasta:q} --output-manifest {output.manifest:q} \
          > {log:q} 2>&1
        """


rule map_tmem_alleles_to_whole_graph:
    input:
        graph=graph_path,
        queries=SV_QUERY_FASTA
    output:
        SV_TARGET_GAF
    log:
        f"{SV_LOGS}/tmem/{{graph}}.allele_mapping.log"
    benchmark:
        f"{SV_BENCHMARKS}/tmem/{{graph}}.allele_mapping.tsv"
    conda:
        LONG_ENV
    threads:
        int(SV_CFG.get("target_threads", 16))
    resources:
        mem_mb=int(SV_CFG.get("target_mem_mb", 64000)),
        runtime_min=int(SV_CFG.get("target_runtime_min", 240)),
        disk_mb=10000
    wildcard_constraints:
        graph=SV_GRAPH_REGEX
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        temporary={output:q}.tmp
        trap 'rm -f -- "$temporary"' EXIT
        minigraph -c -xasm --vc --secondary=yes -t {threads} {input.graph:q} {input.queries:q} \
          > "$temporary" 2> {log:q}
        test -s "$temporary"
        mv -- "$temporary" {output:q}
        trap - EXIT
        """


rule audit_tmem_whole_graph_target:
    input:
        gaf=SV_TARGET_GAF,
        queries=SV_QUERY_MANIFEST
    output:
        SV_TARGET_JOB
    params:
        min_coverage=float(SV_CFG.get("target_min_query_coverage", 0.95)),
        min_identity=float(SV_CFG.get("target_min_identity", 0.90)),
        min_mapq=int(SV_CFG.get("target_min_mapq", 20)),
        min_flank_aligned=int(SV_CFG.get("target_min_flank_aligned_bp", 900))
    log:
        f"{SV_LOGS}/tmem/{{graph}}.target_audit.log"
    benchmark:
        f"{SV_BENCHMARKS}/tmem/{{graph}}.target_audit.tsv"
    conda:
        PY_ENV
    threads: 1
    resources:
        mem_mb=int(SV_CFG.get("metrics_mem_mb", 4000)),
        runtime_min=int(SV_CFG.get("metrics_runtime_min", 120))
    wildcard_constraints:
        graph=SV_GRAPH_REGEX
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {SV_TARGET_SCRIPT:q} audit \
          --graph {wildcards.graph:q} --gaf {input.gaf:q} \
          --query-manifest {input.queries:q} \
          --min-coverage {params.min_coverage} \
          --min-identity {params.min_identity} \
          --min-mapq {params.min_mapq} \
          --min-flank-aligned {params.min_flank_aligned} \
          --output {output:q} > {log:q} 2>&1
        """


rule merge_tmem_whole_graph_targets:
    input:
        expand(SV_TARGET_JOB, graph=SV_GRAPH_IDS)
    output:
        SV_TARGET_PRESENCE
    params:
        inputs=lambda wildcards, input: " ".join("--input " + shlex.quote(str(value)) for value in input)
    log:
        f"{SV_LOGS}/tmem/merge_target_presence.log"
    conda:
        PY_ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {SV_TARGET_SCRIPT:q} merge {params.inputs} \
          --output {output:q} > {log:q} 2>&1
        """


rule map_sv_short_reads:
    input:
        gate=sv_full_gate,
        selected=SV_SELECTED,
        stats=sv_fastq_stats,
        reads=sv_mapping_reads,
        index_dir=SV_INDEX_DIR,
        index_validation=SV_INDEX_VALIDATION,
        target=SV_TARGET_JOB
    output:
        core=temp(SV_CORE_METRIC),
        evidence=SV_JOB_EVIDENCE
    params:
        role=lambda wildcards: sv_sample_row(wildcards).get("role", ""),
        batch=lambda wildcards: sv_sample_row(wildcards).get("batch", ""),
        r1=lambda wildcards: sv_mapping_fastq_args(wildcards, "r1"),
        r2=lambda wildcards: sv_mapping_fastq_args(wildcards, "r2"),
        gbz=lambda wildcards: sv_index_file(wildcards, ".giraffe.gbz"),
        dist=lambda wildcards: sv_index_file(wildcards, ".dist"),
        minimizer=lambda wildcards: sv_index_file(wildcards, ".shortread.withzip.min"),
        zipcodes=lambda wildcards: sv_index_file(wildcards, ".shortread.zipcodes"),
        baseline=SV_CFG.get("baseline_graph", "hprc_chm13_sv"),
        indexing_mode=SV_CFG.get("index_mode", "rgfa_artificial_path_cover"),
        min_mapq=int(SV_CFG.get("min_mapq", 5)),
        min_junction=int(SV_CFG.get("min_junction_reads_per_side", 3)),
        retain_debug=str(bool(SV_CFG.get("retain_debug_gaf", False))).lower(),
        debug_path=lambda wildcards: f"{SV_OUT}/mapping/debug_gaf/{wildcards.mode}/{wildcards.sample}.{wildcards.graph}.gaf.gz"
    log:
        f"{SV_LOGS}/mapping/{{mode}}/{{sample}}.{{graph}}.giraffe.log"
    benchmark:
        f"{SV_BENCHMARKS}/mapping/{{mode}}/{{sample}}.{{graph}}.giraffe.tsv"
    conda:
        VG_ENV
    threads:
        int(SV_CFG.get("threads", 24))
    resources:
        mem_mb=int(SV_CFG.get("mem_mb", 64000)),
        runtime_min=int(SV_CFG.get("runtime_min", 1440)),
        disk_mb=int(SV_CFG.get("disk_mb", 100000)),
        gaf_writers=1
    wildcard_constraints:
        sample=SV_SAMPLE_REGEX,
        graph=SV_GRAPH_REGEX,
        mode="smoke|full"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.core:q})" "$(dirname {output.evidence:q})" "$(dirname {log:q})"
        work=$(mktemp -d -p "${{TMPDIR:-/tmp}}" "sv-giraffe-map.XXXXXX")
        trap 'rm -rf -- "$work"' EXIT
        temporary_dist="$work/graph.dist"
        cp -- {params.dist:q} "$temporary_dist"
        touch --reference={params.minimizer:q} "$temporary_dist"
        test -s {params.gbz:q}
        test -s "$temporary_dist"
        test -s {params.minimizer:q}
        test -s {params.zipcodes:q}

        debug_args=()
        if [ {params.retain_debug:q} = true ]; then
          debug_args=(--debug-gaf "$work/debug.gaf.gz")
        fi

        vg giraffe -p -t {threads} -b default \
          -Z {params.gbz:q} -d "$temporary_dist" \
          -m {params.minimizer:q} -z {params.zipcodes:q} \
          -f <(python {SV_PILOT_SCRIPT:q} cat-fastq {params.r1} 2>> {log:q}) \
          -f <(python {SV_PILOT_SCRIPT:q} cat-fastq {params.r2} 2>> {log:q}) \
          -o GAF 2> {log:q} | \
        python {SV_GAF_SCRIPT:q} \
          --sample {wildcards.sample:q} --role {params.role:q} --batch {params.batch:q} \
          --graph {wildcards.graph:q} --baseline-graph {params.baseline:q} \
          --indexing-mode {params.indexing_mode:q} --run-mode {wildcards.mode:q} \
          --fastq-stats {input.stats:q} --target {input.target:q} \
          --min-mapq {params.min_mapq} --min-junction {params.min_junction} \
          --metrics-output "$work/core.tsv" --evidence-output "$work/evidence.tsv" \
          "${{debug_args[@]}}" 2>> {log:q}

        test -s "$work/core.tsv"
        test -s "$work/evidence.tsv"
        mv -- "$work/core.tsv" {output.core:q}
        mv -- "$work/evidence.tsv" {output.evidence:q}
        if [ {params.retain_debug:q} = true ]; then
          mkdir -p "$(dirname {params.debug_path:q})"
          test -s "$work/debug.gaf.gz"
          mv -- "$work/debug.gaf.gz" {params.debug_path:q}
        fi
        trap - EXIT
        rm -rf -- "$work"
        """


rule finalize_sv_short_read_mapping_metrics:
    input:
        core=SV_CORE_METRIC,
        stats=sv_fastq_stats,
        membership=SV_MEMBERSHIP,
        giraffe_log=f"{SV_LOGS}/mapping/{{mode}}/{{sample}}.{{graph}}.giraffe.log",
        mapping_benchmark=f"{SV_BENCHMARKS}/mapping/{{mode}}/{{sample}}.{{graph}}.giraffe.tsv"
    output:
        SV_RAW_METRIC
    params:
        qc=lambda wildcards: GRAPHS[wildcards.graph].get("qc_status", "")
    log:
        f"{SV_LOGS}/mapping/{{mode}}/{{sample}}.{{graph}}.finalize.log"
    benchmark:
        f"{SV_BENCHMARKS}/mapping/{{mode}}/{{sample}}.{{graph}}.finalize.tsv"
    conda:
        PY_ENV
    threads: 1
    resources:
        mem_mb=int(SV_CFG.get("metrics_mem_mb", 4000)),
        runtime_min=int(SV_CFG.get("metrics_runtime_min", 120))
    wildcard_constraints:
        sample=SV_SAMPLE_REGEX,
        graph=SV_GRAPH_REGEX,
        mode="smoke|full"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {SV_PILOT_SCRIPT:q} finalize-metrics \
          --core {input.core:q} --fastq-stats {input.stats:q} \
          --membership {input.membership:q} --benchmark {input.mapping_benchmark:q} \
          --giraffe-log {input.giraffe_log:q} --sample {wildcards.sample:q} \
          --graph {wildcards.graph:q} --graph-qc-status {params.qc:q} \
          --output {output:q} > {log:q} 2>&1
        """


rule merge_sv_smoke_mapping_metrics:
    input:
        sv_metric_inputs("smoke")
    output:
        SV_SMOKE_METRICS
    params:
        inputs=lambda wildcards, input: " ".join("--input " + shlex.quote(str(value)) for value in input)
    log:
        f"{SV_LOGS}/mapping/smoke/merge_metrics.log"
    conda:
        PY_ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {SV_PILOT_SCRIPT:q} merge-metrics {params.inputs} \
          --output {output:q} > {log:q} 2>&1
        """


rule merge_sv_full_mapping_metrics:
    input:
        sv_metric_inputs("full")
    output:
        SV_FULL_METRICS
    params:
        inputs=lambda wildcards, input: " ".join("--input " + shlex.quote(str(value)) for value in input)
    log:
        f"{SV_LOGS}/mapping/full/merge_metrics.log"
    conda:
        PY_ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {SV_PILOT_SCRIPT:q} merge-metrics {params.inputs} \
          --output {output:q} > {log:q} 2>&1
        """


rule sv_smoke_paired_graph_contrasts:
    input:
        metrics=SV_SMOKE_METRICS
    output:
        SV_SMOKE_CONTRASTS
    log:
        f"{SV_LOGS}/mapping/smoke/contrasts.log"
    conda:
        PY_ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {SV_PILOT_SCRIPT:q} contrasts --metrics {input.metrics:q} \
          --output {output:q} > {log:q} 2>&1
        """


rule sv_full_paired_graph_contrasts:
    input:
        metrics=SV_FULL_METRICS
    output:
        SV_FULL_CONTRASTS
    log:
        f"{SV_LOGS}/mapping/full/contrasts.log"
    conda:
        PY_ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {SV_PILOT_SCRIPT:q} contrasts --metrics {input.metrics:q} \
          --output {output:q} > {log:q} 2>&1
        """


rule merge_sv_smoke_tmem_evidence:
    input:
        sv_evidence_inputs("smoke")
    output:
        SV_SMOKE_EVIDENCE
    params:
        inputs=lambda wildcards, input: " ".join("--input " + shlex.quote(str(value)) for value in input)
    log:
        f"{SV_LOGS}/tmem/smoke.merge_evidence.log"
    conda:
        PY_ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {SV_PILOT_SCRIPT:q} merge-evidence {params.inputs} \
          --output {output:q} > {log:q} 2>&1
        """


rule merge_sv_full_tmem_evidence:
    input:
        sv_evidence_inputs("full")
    output:
        SV_FULL_EVIDENCE
    params:
        inputs=lambda wildcards, input: " ".join("--input " + shlex.quote(str(value)) for value in input)
    log:
        f"{SV_LOGS}/tmem/full.merge_evidence.log"
    conda:
        PY_ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {SV_PILOT_SCRIPT:q} merge-evidence {params.inputs} \
          --output {output:q} > {log:q} 2>&1
        """


rule gate_sv_short_read_smoke:
    input:
        metrics=SV_SMOKE_METRICS
    output:
        SV_SMOKE_READY
    params:
        expected=len(SV_JOBS),
        min_fraction=float(SV_CFG.get("smoke_min_mapped_fraction", 0.80))
    log:
        f"{SV_LOGS}/mapping/smoke/gate.log"
    conda:
        PY_ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {SV_PILOT_SCRIPT:q} smoke-gate --metrics {input.metrics:q} \
          --expected-rows {params.expected} --min-mapped-fraction {params.min_fraction} \
          --output {output:q} > {log:q} 2>&1
        """


rule sv_short_read_provenance:
    input:
        config=CFG_PATH,
        selected=SV_SELECTED,
        index_preflight=SV_INDEX_PREFLIGHT,
        tools=SV_TOOL_VERSIONS
    output:
        SV_PROVENANCE
    log:
        f"{SV_LOGS}/provenance.log"
    conda:
        PY_ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {SV_PILOT_SCRIPT:q} provenance --config {input.config:q} \
          --selected {input.selected:q} --index-preflight {input.index_preflight:q} \
          --tools {input.tools:q} --output {output:q} > {log:q} 2>&1
        """


rule render_sv_short_read_smoke_report:
    input:
        preflight=SV_INDEX_PREFLIGHT,
        selected=SV_SELECTED,
        membership=SV_MEMBERSHIP,
        manifests=expand(SV_INDEX_MANIFEST, graph=SV_GRAPH_IDS),
        index_benchmarks=expand(
            f"{SV_BENCHMARKS}/indexes/{{graph}}.autoindex.tsv", graph=SV_GRAPH_IDS
        ),
        metrics=SV_SMOKE_METRICS,
        contrasts=SV_SMOKE_CONTRASTS,
        targets=SV_TARGET_PRESENCE,
        evidence=SV_SMOKE_EVIDENCE,
        tools=SV_TOOL_VERSIONS,
        provenance=SV_PROVENANCE
    output:
        SV_SMOKE_REPORT
    params:
        manifests=lambda wildcards, input: " ".join("--index-manifest " + shlex.quote(str(value)) for value in input.manifests),
        index_benchmarks=lambda wildcards, input: " ".join(
            "--index-benchmark " + shlex.quote(f"{graph}={path}")
            for graph, path in zip(SV_GRAPH_IDS, input.index_benchmarks)
        )
    log:
        f"{SV_LOGS}/report/smoke.log"
    benchmark:
        f"{SV_BENCHMARKS}/report/smoke.tsv"
    conda:
        PY_ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {SV_REPORT_SCRIPT:q} --mode smoke \
          --index-preflight {input.preflight:q} --selected {input.selected:q} \
          --membership {input.membership:q} {params.manifests} {params.index_benchmarks} \
          --metrics {input.metrics:q} --contrasts {input.contrasts:q} \
          --target-presence {input.targets:q} --evidence {input.evidence:q} \
          --tools {input.tools:q} --provenance {input.provenance:q} \
          --output {output:q} > {log:q} 2>&1
        """


rule render_sv_short_read_full_report:
    input:
        preflight=SV_INDEX_PREFLIGHT,
        selected=SV_SELECTED,
        membership=SV_MEMBERSHIP,
        manifests=expand(SV_INDEX_MANIFEST, graph=SV_GRAPH_IDS),
        index_benchmarks=expand(
            f"{SV_BENCHMARKS}/indexes/{{graph}}.autoindex.tsv", graph=SV_GRAPH_IDS
        ),
        metrics=SV_FULL_METRICS,
        contrasts=SV_FULL_CONTRASTS,
        targets=SV_TARGET_PRESENCE,
        evidence=SV_FULL_EVIDENCE,
        tools=SV_TOOL_VERSIONS,
        provenance=SV_PROVENANCE
    output:
        SV_FULL_REPORT
    params:
        manifests=lambda wildcards, input: " ".join("--index-manifest " + shlex.quote(str(value)) for value in input.manifests),
        index_benchmarks=lambda wildcards, input: " ".join(
            "--index-benchmark " + shlex.quote(f"{graph}={path}")
            for graph, path in zip(SV_GRAPH_IDS, input.index_benchmarks)
        )
    log:
        f"{SV_LOGS}/report/full.log"
    benchmark:
        f"{SV_BENCHMARKS}/report/full.tsv"
    conda:
        PY_ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {SV_REPORT_SCRIPT:q} --mode full \
          --index-preflight {input.preflight:q} --selected {input.selected:q} \
          --membership {input.membership:q} {params.manifests} {params.index_benchmarks} \
          --metrics {input.metrics:q} --contrasts {input.contrasts:q} \
          --target-presence {input.targets:q} --evidence {input.evidence:q} \
          --tools {input.tools:q} --provenance {input.provenance:q} \
          --output {output:q} > {log:q} 2>&1
        """


rule preflight_graphs:
    input:
        config=CFG_PATH,
        tmem_indexes=expand(list(TMEM_INDEX_PATTERNS.values()), graph=TMEM_GRAPH_IDS)
    output:
        registry=GRAPH_REGISTRY,
        preflight=GRAPH_PREFLIGHT,
        report=GRAPH_PREFLIGHT_REPORT
    log:
        f"{LOGS}/preflight_graphs.log"
    benchmark:
        f"{BENCHMARKS}/preflight_graphs.tsv"
    conda:
        PY_ENV
    threads: 1
    resources:
        mem_mb=int(METADATA_RESOURCES.get("metadata_mem_mb", 4000)),
        runtime_min=int(METADATA_RESOURCES.get("metadata_runtime_min", 60))
    shell:
        r"""
        mkdir -p "$(dirname {output.registry:q})" "$(dirname {output.report:q})" "$(dirname {log:q})"
        python {REGISTRY_SCRIPT:q} preflight \
          --config {input.config:q} \
          --registry {output.registry:q} \
          --preflight {output.preflight:q} \
          --report {output.report:q} > {log:q} 2>&1
        """


checkpoint select_tmem_pilot:
    input:
        config=CFG_PATH
    output:
        manifest=TMEM_MANIFEST
    log:
        f"{LOGS}/select_tmem_pilot.log"
    benchmark:
        f"{BENCHMARKS}/select_tmem_pilot.tsv"
    conda:
        PY_ENV
    threads: 1
    resources:
        mem_mb=int(METADATA_RESOURCES.get("metadata_mem_mb", 4000)),
        runtime_min=int(METADATA_RESOURCES.get("metadata_runtime_min", 60))
    shell:
        r"""
        mkdir -p "$(dirname {output.manifest:q})" "$(dirname {log:q})"
        python {TMEM_SELECT_SCRIPT:q} \
          --config {input.config:q} \
          --output {output.manifest:q} > {log:q} 2>&1
        """


rule sample_registry_and_leakage:
    input:
        config=CFG_PATH,
        tmem=TMEM_MANIFEST
    output:
        samples=SAMPLES,
        membership=MEMBERSHIP,
        leakage=LEAKAGE
    log:
        f"{LOGS}/sample_registry_and_leakage.log"
    benchmark:
        f"{BENCHMARKS}/sample_registry_and_leakage.tsv"
    conda:
        PY_ENV
    threads: 1
    resources:
        mem_mb=int(METADATA_RESOURCES.get("metadata_mem_mb", 4000)),
        runtime_min=int(METADATA_RESOURCES.get("metadata_runtime_min", 60))
    shell:
        r"""
        mkdir -p "$(dirname {output.samples:q})" "$(dirname {log:q})"
        python {REGISTRY_SCRIPT:q} samples \
          --config {input.config:q} \
          --tmem-manifest {input.tmem:q} \
          --samples {output.samples:q} \
          --membership {output.membership:q} \
          --leakage {output.leakage:q} > {log:q} 2>&1
        """


rule gate_long_reads:
    input:
        preflight=GRAPH_PREFLIGHT,
        leakage=LEAKAGE
    output:
        temp(LONG_GATE)
    params:
        graphs=" ".join("--graph-id " + shlex.quote(graph) for graph in HEADLINE_GRAPHS)
    log:
        f"{LOGS}/gate_long_reads.log"
    conda:
        PY_ENV
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {REGISTRY_SCRIPT:q} gate \
          --preflight {input.preflight:q} \
          --leakage {input.leakage:q} \
          --mode long {params.graphs} \
          --output {output:q} > {log:q} 2>&1
        """


rule gate_tmem_short_reads:
    input:
        preflight=GRAPH_PREFLIGHT,
        leakage=LEAKAGE
    output:
        temp(SHORT_GATE)
    params:
        graphs=" ".join("--graph-id " + shlex.quote(graph) for graph in TMEM_CFG["graph_ids"])
    log:
        f"{LOGS}/gate_tmem_short_reads.log"
    conda:
        PY_ENV
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {REGISTRY_SCRIPT:q} gate \
          --preflight {input.preflight:q} \
          --leakage {input.leakage:q} \
          --mode short {params.graphs} \
          --output {output:q} > {log:q} 2>&1
        """


rule gate_v2_rescue:
    input:
        preflight=GRAPH_PREFLIGHT,
        leakage=LEAKAGE
    output:
        temp(RESCUE_GATE)
    params:
        graphs=" ".join(
            "--graph-id " + shlex.quote(graph)
            for graph in RESCUE_GRAPHS
        )
    log:
        f"{LOGS}/gate_v2_rescue.log"
    conda:
        PY_ENV
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {REGISTRY_SCRIPT:q} gate \
          --preflight {input.preflight:q} \
          --leakage {input.leakage:q} \
          --mode rescue {params.graphs} \
          --output {output:q} > {log:q} 2>&1
        """


rule write_benchmark_hifi_archive_manifest:
    input:
        assembly_manifest=ARCHIVE_CFG["assembly_manifest"],
        centenarians=ARCHIVE_CFG["centenarian_samples"]
    output:
        ARCHIVE_READ_MANIFEST
    params:
        archive_root=ARCHIVE_CFG["archive_root"],
        archive_pattern=lambda wildcards: ARCHIVE_CFG["archive_bam_pattern"],
        fastq_pattern=lambda wildcards: ARCHIVE_CFG["fastq_pattern"]
    log:
        f"{LOGS}/staging/write_hifi_archive_manifest.log"
    conda:
        PY_ENV
    threads: 1
    resources:
        mem_mb=int(METADATA_RESOURCES.get("metadata_mem_mb", 4000)),
        runtime_min=int(METADATA_RESOURCES.get("metadata_runtime_min", 60))
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {ARCHIVE_SCRIPT:q} manifest \
          --assembly-manifest {input.assembly_manifest:q} \
          --centenarian-samples {input.centenarians:q} \
          --archive-root {params.archive_root:q} \
          --archive-pattern {params.archive_pattern:q} \
          --fastq-pattern {params.fastq_pattern:q} \
          --output {output:q} > {log:q} 2>&1
        """


rule stage_benchmark_hifi_bams:
    input:
        read_manifest=ARCHIVE_READ_MANIFEST
    output:
        temp(directory(ARCHIVE_STAGED))
    params:
        samples=" ".join(
            "--sample " + shlex.quote(sample) for sample in ARCHIVE_MEMBER_SAMPLES
        ),
        daget=ARCHIVE_CFG["daget"],
        dals=ARCHIVE_CFG["dals"],
        poll_seconds=int(ARCHIVE_CFG.get("poll_seconds", 30)),
        max_wait=int(ARCHIVE_CFG.get("max_wait_min", 1440))
    log:
        f"{LOGS}/staging/archive/daget.log"
    conda:
        PY_ENV
    threads: 1
    resources:
        partition=ARCHIVE_CFG.get("archive_partition", "archive"),
        mem_mb=int(ARCHIVE_CFG.get("stage_mem_mb", 100)),
        runtime_min=int(ARCHIVE_CFG.get("stage_runtime_min", 1500))
    shell:
        r"""
        mkdir -p "$(dirname {log:q})"
        python {ARCHIVE_SCRIPT:q} stage \
          {params.samples} \
          --read-manifest {input.read_manifest:q} \
          --output-dir {output:q} \
          --daget {params.daget:q} \
          --dals {params.dals:q} \
          --poll-seconds {params.poll_seconds} \
          --max-wait-minutes {params.max_wait} \
          --log {log:q}
        """


rule copy_benchmark_hifi_bam:
    input:
        staged=ARCHIVE_STAGED
    output:
        bam=temp(ARCHIVE_BAM)
    params:
        source=archive_bam,
        marker=lambda wildcards: f"{ARCHIVE_STAGED}/{wildcards.sample}.staged"
    log:
        f"{LOGS}/staging/archive/copy/{{sample}}.log"
    conda:
        PY_ENV
    threads: 1
    resources:
        partition=ARCHIVE_CFG.get("archive_partition", "archive"),
        mem_mb=int(ARCHIVE_CFG.get("copy_mem_mb", 100)),
        runtime_min=int(ARCHIVE_CFG.get("copy_runtime_min", 720))
    wildcard_constraints:
        sample=ARCHIVE_SAMPLE_REGEX
    shell:
        r"""
        mkdir -p "$(dirname {output.bam:q})" "$(dirname {log:q})"
        python {ARCHIVE_SCRIPT:q} copy \
          --source {params.source:q} \
          --marker {params.marker:q} \
          --output {output.bam:q} \
          --log {log:q}
        """


rule convert_benchmark_hifi_bam:
    input:
        staged_bam=ARCHIVE_BAM,
        header_template=ARCHIVE_CFG["header_template"]
    output:
        reads=protected(ARCHIVE_FASTQ),
        header_qc=ARCHIVE_HEADER_QC
    params:
        archive_bam=archive_bam,
        cohort=lambda wildcards: SAMPLE_CONFIG[wildcards.sample]["archive_cohort"]
    log:
        f"{LOGS}/staging/hifi_conversion/{{sample}}.log"
    benchmark:
        f"{BENCHMARKS}/staging/hifi_conversion/{{sample}}.tsv"
    conda:
        LONG_ENV
    threads:
        int(ARCHIVE_CFG.get("conversion_threads", 8))
    resources:
        mem_mb=int(ARCHIVE_CFG.get("conversion_mem_mb", 16000)),
        runtime_min=int(ARCHIVE_CFG.get("conversion_runtime_min", 720)),
        disk_mb=int(ARCHIVE_CFG.get("conversion_disk_mb", 250000))
    wildcard_constraints:
        sample=ARCHIVE_SAMPLE_REGEX
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.reads:q})" \
                 "$(dirname {output.header_qc:q})" \
                 "$(dirname {log:q})"
        temporary_fastq={output.reads:q}.tmp
        temporary_qc={output.header_qc:q}.tmp
        header_work=$(mktemp -d -p "${{TMPDIR:-/tmp}}" "benchmark-hifi-header.XXXXXX")
        trap 'rm -f -- "$temporary_fastq" "$temporary_qc"; rm -rf -- "$header_work"' EXIT

        test -s {input.staged_bam:q}
        test -r {input.staged_bam:q}
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

        samtools fastq -@ {threads} -F 0x900 -n {input.staged_bam:q} 2> {log:q} \
          | bgzip -@ {threads} -c > "$temporary_fastq" 2>> {log:q}
        test -s "$temporary_fastq"
        gzip -t "$temporary_fastq" 2>> {log:q}

        {{
          printf 'metric\tvalue\n'
          printf 'sample_id\t%s\n' {wildcards.sample:q}
          printf 'cohort\t%s\n' {params.cohort:q}
          printf 'archive_bam\t%s\n' {params.archive_bam:q}
          printf 'sq_count\t%s\n' "$(wc -l < "$header_work/actual.sq")"
          printf 'sq_sha256\t%s\n' "$(sha256sum "$header_work/actual.sq" | cut -d' ' -f1)"
          printf 'fastq_size_bytes\t%s\n' "$(stat -c '%s' "$temporary_fastq")"
        }} > "$temporary_qc"

        mv "$temporary_fastq" {output.reads:q}
        mv "$temporary_qc" {output.header_qc:q}
        trap - EXIT
        rm -rf -- "$header_work"
        """


rule select_long_reads:
    input:
        reads=sample_reads
    output:
        names=temp(f"{OUT}/long_reads/selection/{{sample}}.read_names.txt"),
        stats=f"{OUT}/long_reads/selection/{{sample}}.selection.tsv"
    params:
        coverage=sample_coverage,
        genome_size=int(LONG_CFG["genome_size_bp"]),
        seed=int(LONG_CFG["selection_seed"])
    log:
        f"{LOGS}/long_reads/{{sample}}.select.log"
    benchmark:
        f"{BENCHMARKS}/long_reads/{{sample}}.select.tsv"
    conda:
        LONG_ENV
    threads: 1
    resources:
        mem_mb=int(LONG_CFG.get("metrics_mem_mb", 8000)),
        runtime_min=int(LONG_CFG.get("metrics_runtime_min", 120)),
        disk_mb=1000
    shell:
        r"""
        mkdir -p "$(dirname {output.names:q})" "$(dirname {log:q})"
        python {READ_SELECT_SCRIPT:q} \
          --fastq {input.reads:q} \
          --sample {wildcards.sample:q} \
          --coverage {params.coverage:q} \
          --genome-size {params.genome_size} \
          --seed {params.seed} \
          --read-names {output.names:q} \
          --stats {output.stats:q} > {log:q} 2>&1
        """


rule map_long_reads:
    input:
        gate=long_mapping_gate,
        reads=sample_reads,
        names=f"{OUT}/long_reads/selection/{{sample}}.read_names.txt"
    output:
        temp(LONG_ALIGNMENT)
    params:
        graph=validate_long_job,
        minigraph_preset=LONG_CFG.get("minigraph_preset", "lr"),
        minigraph_extra=LONG_CFG.get("minigraph_extra", ""),
        minimap2_preset=minimap2_preset,
        minimap2_extra=LONG_CFG.get("minimap2_extra", "--secondary=yes")
    log:
        f"{LOGS}/long_reads/{{sample}}.{{mapper}}.{{graph}}.log"
    benchmark:
        f"{BENCHMARKS}/long_reads/{{sample}}.{{mapper}}.{{graph}}.tsv"
    conda:
        LONG_ENV
    threads:
        int(LONG_CFG.get("threads", 24))
    resources:
        mem_mb=int(LONG_CFG.get("mem_mb", 64000)),
        runtime_min=int(LONG_CFG.get("runtime_min", 720)),
        disk_mb=int(LONG_CFG.get("disk_mb", 20000)),
        gaf_writers=1
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        temporary={output:q}.tmp
        trap 'rm -f -- "$temporary"' EXIT
        if [ {wildcards.mapper:q} = minigraph ]; then
            minigraph -cx{params.minigraph_preset} -t {threads} {params.minigraph_extra} \
              {params.graph:q} <(seqtk subseq {input.reads:q} {input.names:q}) \
              > "$temporary" 2> {log:q}
        elif [ {wildcards.mapper:q} = minimap2 ]; then
            minimap2 -cx {params.minimap2_preset:q} -t {threads} {params.minimap2_extra} \
              {params.graph:q} <(seqtk subseq {input.reads:q} {input.names:q}) \
              > "$temporary" 2> {log:q}
        else
            echo "Unsupported mapper: {wildcards.mapper}" >&2
            exit 1
        fi
        test -s "$temporary"
        mv "$temporary" {output:q}
        trap - EXIT
        """


rule summarize_long_mapping:
    input:
        alignment=LONG_ALIGNMENT,
        selection=f"{OUT}/long_reads/selection/{{sample}}.selection.tsv",
        samples=SAMPLES,
        membership=MEMBERSHIP,
        benchmark=f"{BENCHMARKS}/long_reads/{{sample}}.{{mapper}}.{{graph}}.tsv"
    output:
        temp(LONG_JOB_METRIC)
    log:
        f"{LOGS}/long_reads/{{sample}}.{{mapper}}.{{graph}}.metrics.log"
    benchmark:
        f"{BENCHMARKS}/long_reads/{{sample}}.{{mapper}}.{{graph}}.metrics.tsv"
    conda:
        PY_ENV
    threads: 1
    resources:
        mem_mb=int(LONG_CFG.get("metrics_mem_mb", 8000)),
        runtime_min=int(LONG_CFG.get("metrics_runtime_min", 120))
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {METRICS_SCRIPT:q} summarize \
          --alignment {input.alignment:q} \
          --selection-stats {input.selection:q} \
          --samples {input.samples:q} \
          --membership {input.membership:q} \
          --benchmark {input.benchmark:q} \
          --sample {wildcards.sample:q} \
          --graph {wildcards.graph:q} \
          --mapper {wildcards.mapper:q} \
          --output {output:q} > {log:q} 2>&1
        """


rule merge_long_mapping_metrics:
    input:
        long_metric_inputs()
    output:
        MAPPING_METRICS
    params:
        inputs=lambda wildcards, input: " ".join("--input " + shlex.quote(str(path)) for path in input)
    log:
        f"{LOGS}/merge_long_mapping_metrics.log"
    conda:
        PY_ENV
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {METRICS_SCRIPT:q} merge {params.inputs} --output {output:q} > {log:q} 2>&1
        """


rule paired_graph_contrasts:
    input:
        metrics=MAPPING_METRICS
    output:
        PAIRED_CONTRASTS
    log:
        f"{LOGS}/paired_graph_contrasts.log"
    conda:
        PY_ENV
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {METRICS_SCRIPT:q} contrasts \
          --metrics {input.metrics:q} \
          --output {output:q} > {log:q} 2>&1
        """


rule simulate_phased_hifi:
    input:
        hap1=lambda wildcards: SAMPLE_CONFIG[wildcards.sample]["truth_assembly_hap1"],
        hap2=lambda wildcards: SAMPLE_CONFIG[wildcards.sample]["truth_assembly_hap2"]
    output:
        fastq=temp(SIM_FASTQ),
        stats=f"{OUT}/simulation/reads/{{sample}}.simulation.tsv"
    params:
        profile_arg=(
            "--empirical-profile-json " + shlex.quote(SIM_CFG["empirical_profile_json"])
            if SIM_CFG.get("empirical_profile_json") else ""
        ),
        coverage=SIM_CFG.get("coverage", 5),
        seed=SIM_CFG.get("seed", 1),
        read_length_mean=SIM_CFG.get("read_length_mean", 15000),
        read_length_sd=SIM_CFG.get("read_length_sd", 4000),
        substitution_rate=SIM_CFG.get("substitution_rate", 0.0015),
        insertion_rate=SIM_CFG.get("insertion_rate", 0.0005),
        deletion_rate=SIM_CFG.get("deletion_rate", 0.0005)
    log:
        f"{LOGS}/simulation/{{sample}}.simulate.log"
    benchmark:
        f"{BENCHMARKS}/simulation/{{sample}}.simulate.tsv"
    conda:
        PY_ENV
    threads: 1
    resources:
        mem_mb=int(LONG_CFG.get("mem_mb", 64000)),
        runtime_min=int(LONG_CFG.get("runtime_min", 720)),
        disk_mb=int(LONG_CFG.get("disk_mb", 20000))
    shell:
        r"""
        mkdir -p "$(dirname {output.fastq:q})" "$(dirname {log:q})"
        python {SIMULATE_SCRIPT:q} \
          --sample {wildcards.sample:q} --hap1 {input.hap1:q} --hap2 {input.hap2:q} \
          --coverage {params.coverage} --seed {params.seed} \
          --read-length-mean {params.read_length_mean} \
          --read-length-sd {params.read_length_sd} \
          --substitution-rate {params.substitution_rate} \
          --insertion-rate {params.insertion_rate} \
          --deletion-rate {params.deletion_rate} \
          {params.profile_arg} \
          --fastq {output.fastq:q} --stats {output.stats:q} > {log:q} 2>&1
        """


rule map_simulated_hifi:
    input:
        gate=LONG_GATE,
        fastq=SIM_FASTQ
    output:
        temp(SIM_ALIGNMENT)
    params:
        graph=graph_path,
        preset=LONG_CFG.get("minigraph_preset", "lr"),
        extra=LONG_CFG.get("minigraph_extra", "")
    log:
        f"{LOGS}/simulation/{{sample}}.minigraph.{{graph}}.log"
    benchmark:
        f"{BENCHMARKS}/simulation/{{sample}}.minigraph.{{graph}}.tsv"
    conda:
        LONG_ENV
    threads:
        int(LONG_CFG.get("threads", 24))
    resources:
        mem_mb=int(LONG_CFG.get("mem_mb", 64000)),
        runtime_min=int(LONG_CFG.get("runtime_min", 720)),
        disk_mb=int(LONG_CFG.get("disk_mb", 20000)),
        gaf_writers=1
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        temporary={output:q}.tmp
        trap 'rm -f -- "$temporary"' EXIT
        minigraph -cx{params.preset} -t {threads} {params.extra} \
          {params.graph:q} {input.fastq:q} > "$temporary" 2> {log:q}
        test -s "$temporary"
        mv "$temporary" {output:q}
        trap - EXIT
        """


rule assess_simulated_mapping:
    input:
        alignment=SIM_ALIGNMENT
    output:
        temp(SIM_JOB_ACCURACY)
    params:
        tolerance=SIM_CFG.get("locus_tolerance_bp", 10000)
    log:
        f"{LOGS}/simulation/{{sample}}.minigraph.{{graph}}.accuracy.log"
    conda:
        PY_ENV
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {METRICS_SCRIPT:q} sim-accuracy \
          --alignment {input.alignment:q} --sample {wildcards.sample:q} \
          --graph {wildcards.graph:q} --tolerance {params.tolerance} \
          --output {output:q} > {log:q} 2>&1
        """


rule merge_simulated_mapping_accuracy:
    input:
        sim_accuracy_inputs()
    output:
        SIM_ACCURACY
    params:
        inputs=lambda wildcards, input: " ".join("--input " + shlex.quote(str(path)) for path in input)
    log:
        f"{LOGS}/merge_simulated_mapping_accuracy.log"
    conda:
        PY_ENV
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {METRICS_SCRIPT:q} merge-accuracy {params.inputs} --output {output:q} > {log:q} 2>&1
        """


rule assess_v2_rescue_capture:
    input:
        gate=RESCUE_GATE,
        alignments=rescue_alignment_inputs()
    output:
        RESCUE_CAPTURE
    params:
        v1=GRAPHS[RESCUE_CFG.get("v1_graph_id", "local_v1")]["path"],
        v2=GRAPHS[RESCUE_CFG.get("v2_graph_id", "local_v2")]["path"],
        alignments=rescue_alignment_args()
    log:
        f"{LOGS}/assess_v2_rescue_capture.log"
    benchmark:
        f"{BENCHMARKS}/assess_v2_rescue_capture.tsv"
    conda:
        PY_ENV
    threads: 1
    resources:
        mem_mb=int(LONG_CFG.get("mem_mb", 64000)),
        runtime_min=int(LONG_CFG.get("runtime_min", 720))
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {RESCUE_SCRIPT:q} \
          --config {CFG_PATH:q} --v1-graph {params.v1:q} --v2-graph {params.v2:q} \
          {params.alignments} --output {output:q} > {log:q} 2>&1
        """


rule map_tmem_short_reads:
    input:
        gate=SHORT_GATE,
        manifest=lambda wildcards: str(checkpoints.select_tmem_pilot.get().output.manifest),
        r1=tmem_r1,
        r2=tmem_r2,
        gbz=graph_path,
        dist=TMEM_INDEX_PATTERNS["dist"],
        minimizer=TMEM_INDEX_PATTERNS["minimizer"],
        zipcodes=TMEM_INDEX_PATTERNS["zipcodes"]
    output:
        temp(TMEM_GAF)
    params:
        r1_files=lambda wildcards: quoted_tmem_files(wildcards, "r1"),
        r2_files=lambda wildcards: quoted_tmem_files(wildcards, "r2")
    log:
        f"{LOGS}/tmem/{{sample}}.{{graph}}.giraffe.log"
    benchmark:
        f"{BENCHMARKS}/tmem/{{sample}}.{{graph}}.giraffe.tsv"
    conda:
        VG_ENV
    threads:
        int(TMEM_CFG.get("threads", 24))
    resources:
        mem_mb=int(TMEM_CFG.get("mem_mb", 32000)),
        runtime_min=int(TMEM_CFG.get("runtime_min", 480)),
        disk_mb=int(TMEM_CFG.get("disk_mb", 10000)),
        gaf_writers=1
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        temporary={output:q}.tmp
        temporary_dist=$(mktemp "${{TMPDIR:-/tmp}}/tmem.dist.XXXXXX")
        trap 'rm -f -- "$temporary" "$temporary_dist"' EXIT
        cp {input.dist:q} "$temporary_dist"
        touch --reference={input.minimizer:q} "$temporary_dist"
        vg giraffe -p -t {threads} -b default \
          -Z {input.gbz:q} -d "$temporary_dist" \
          -m {input.minimizer:q} -z {input.zipcodes:q} \
          -f <(gzip -cd -- {params.r1_files}) \
          -f <(gzip -cd -- {params.r2_files}) \
          -o GAF > "$temporary" 2> {log:q}
        test -s "$temporary"
        mv "$temporary" {output:q}
        trap - EXIT
        rm -f -- "$temporary_dist"
        """


rule summarize_tmem_evidence:
    input:
        gaf=TMEM_GAF,
        manifest=lambda wildcards: str(checkpoints.select_tmem_pilot.get().output.manifest),
        r1=tmem_r1,
        r2=tmem_r2
    output:
        temp(TMEM_JOB_EVIDENCE)
    params:
        targets=TMEM_CFG["target_insertions"],
        vcf=TMEM_CFG["graph_vcf"],
        kmer_size=TMEM_CFG.get("kmer_size", 31),
        edge_window=TMEM_CFG.get("edge_window_bp", 250),
        min_mapq=TMEM_CFG.get("min_mapq", 5),
        min_junction=TMEM_CFG.get("min_junction_reads_per_side", 3),
        min_reference=TMEM_CFG.get("min_reference_reads", 3)
    log:
        f"{LOGS}/tmem/{{sample}}.{{graph}}.evidence.log"
    benchmark:
        f"{BENCHMARKS}/tmem/{{sample}}.{{graph}}.evidence.tsv"
    conda:
        PY_ENV
    threads: 1
    resources:
        mem_mb=int(TMEM_CFG.get("mem_mb", 32000)),
        runtime_min=int(TMEM_CFG.get("runtime_min", 480))
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {TMEM_SCRIPT:q} summarize \
          --gaf {input.gaf:q} --targets {params.targets:q} --vcf {params.vcf:q} \
          --sample-manifest {input.manifest:q} --sample {wildcards.sample:q} --graph {wildcards.graph:q} \
          --kmer-size {params.kmer_size} --edge-window {params.edge_window} \
          --min-mapq {params.min_mapq} \
          --min-junction {params.min_junction} \
          --min-reference {params.min_reference} \
          --output {output:q} > {log:q} 2>&1
        """


rule merge_tmem_evidence:
    input:
        tmem_targets
    output:
        TMEM_EVIDENCE
    params:
        inputs=lambda wildcards, input: " ".join("--input " + shlex.quote(str(path)) for path in input)
    log:
        f"{LOGS}/merge_tmem_evidence.log"
    conda:
        PY_ENV
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {TMEM_SCRIPT:q} merge {params.inputs} --output {output:q} > {log:q} 2>&1
        """


rule normalize_sv_catalog:
    input:
        config=CFG_PATH
    output:
        NORMALIZED_SV
    log:
        f"{LOGS}/normalize_sv_catalog.log"
    benchmark:
        f"{BENCHMARKS}/normalize_sv_catalog.tsv"
    conda:
        PY_ENV
    threads: 1
    resources:
        mem_mb=int(METADATA_RESOURCES.get("report_mem_mb", 8000)),
        runtime_min=int(METADATA_RESOURCES.get("report_runtime_min", 120))
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {SV_SCRIPT:q} normalize --config {input.config:q} --output {output:q} > {log:q} 2>&1
        """


rule cluster_normalized_sv_catalog:
    input:
        catalog=NORMALIZED_SV
    output:
        SV_CLUSTERS
    params:
        breakpoint_distance=SV_CATALOG_CFG.get("breakpoint_distance", 500),
        length_similarity=SV_CATALOG_CFG.get("length_similarity", 0.70),
        sequence_similarity=SV_CATALOG_CFG.get("sequence_similarity", 0.80)
    log:
        f"{LOGS}/cluster_normalized_sv_catalog.log"
    benchmark:
        f"{BENCHMARKS}/cluster_normalized_sv_catalog.tsv"
    conda:
        PY_ENV
    threads: 1
    resources:
        mem_mb=int(METADATA_RESOURCES.get("report_mem_mb", 8000)),
        runtime_min=int(METADATA_RESOURCES.get("report_runtime_min", 120))
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {SV_SCRIPT:q} cluster --catalog {input.catalog:q} \
          --breakpoint-distance {params.breakpoint_distance} \
          --length-similarity {params.length_similarity} \
          --sequence-similarity {params.sequence_similarity} \
          --output {output:q} > {log:q} 2>&1
        """


rule benchmark_tool_versions:
    input:
        environment=LONG_ENV
    output:
        TOOL_VERSIONS
    log:
        f"{LOGS}/benchmark_tool_versions.log"
    conda:
        LONG_ENV
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {REGISTRY_SCRIPT:q} versions --output {output:q} > {log:q} 2>&1
        """


rule benchmark_provenance:
    input:
        config=CFG_PATH,
        preflight=GRAPH_PREFLIGHT
    output:
        PROVENANCE
    log:
        f"{LOGS}/benchmark_provenance.log"
    conda:
        PY_ENV
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {REGISTRY_SCRIPT:q} provenance \
          --config {input.config:q} --preflight {input.preflight:q} \
          --output {output:q} > {log:q} 2>&1
        """


rule render_read_graph_benchmark_report:
    input:
        preflight=GRAPH_PREFLIGHT,
        leakage=LEAKAGE,
        metrics=MAPPING_METRICS,
        contrasts=PAIRED_CONTRASTS,
        simulation=SIM_ACCURACY,
        rescue=RESCUE_CAPTURE,
        tmem=TMEM_EVIDENCE
    output:
        FINAL_REPORT
    log:
        f"{LOGS}/render_read_graph_benchmark_report.log"
    benchmark:
        f"{BENCHMARKS}/render_read_graph_benchmark_report.tsv"
    conda:
        PY_ENV
    threads: 1
    resources:
        mem_mb=int(METADATA_RESOURCES.get("report_mem_mb", 8000)),
        runtime_min=int(METADATA_RESOURCES.get("report_runtime_min", 120))
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {REPORT_SCRIPT:q} \
          --preflight {input.preflight:q} --leakage {input.leakage:q} \
          --metrics {input.metrics:q} --contrasts {input.contrasts:q} \
          --simulation {input.simulation:q} --rescue {input.rescue:q} \
          --tmem {input.tmem:q} --output {output:q} > {log:q} 2>&1
        """


wildcard_constraints:
    sample="[A-Za-z0-9_.-]+",
    graph="[A-Za-z0-9_.-]+",
    mapper="minigraph|minimap2"
