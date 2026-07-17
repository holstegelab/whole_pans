NOVEL_SV_CONFIG = workflow.source_path("../../config/config.yaml")
# This file is a standalone stage entry point, like the other files in
# workflow/rules/.  Do not include it from whole_pans/Snakefile: configfile and
# workdir are workflow-global directives.
configfile: NOVEL_SV_CONFIG


import csv
import gzip
import json
import os
import re
import shlex
from collections import defaultdict
from pathlib import Path
from pathlib import PurePosixPath


NOVEL = config.get("novel_sv_discovery", {})
NOVEL_SCREEN = NOVEL.get("graph_screen", {})
NOVEL_REFERENCE_CALLING = NOVEL.get("reference_calling", {})
NOVEL_CALLERS = NOVEL_REFERENCE_CALLING.get("callers", {})
NOVEL_CLUSTERING = NOVEL.get("clustering", {})
NOVEL_HIFI = NOVEL.get("hifi_recommendations", {})
NOVEL_PROVENANCE = NOVEL.get("graph_provenance", {})
NOVEL_PROJECT_ROOT = Path(workflow.current_basedir).resolve().parents[1]
NOVEL_ALLOWED_ROOT = NOVEL_PROJECT_ROOT.parent
workdir: str(NOVEL_PROJECT_ROOT)


NOVEL_PATHS = {
    "results": NOVEL.get("results", "results/novel_sv_discovery"),
    "frozen_graph": NOVEL.get("frozen_graph", "results/graphs/sv_pangenome.minigraph.gfa"),
    "graph_assemblies": NOVEL.get(
        "graph_assemblies", "results/metadata/sv_pangenome.ordered_assemblies.tsv"
    ),
    "cleaned_manifest": NOVEL.get(
        "cleaned_manifest", "results/post_decontamination_qc/resources/all_cleaned_assemblies.tsv"
    ),
    "original_qc": NOVEL.get("original_qc", "results/qc/summary/assembly_qc.tsv"),
    "post_qc": NOVEL.get("post_qc", "results/post_qc/summary/assembly_qc.tsv"),
    "contamination": NOVEL.get(
        "contamination", "results/decontamination/summary/contamination_summary.tsv"
    ),
    "previous_feasibility": NOVEL.get(
        "previous_feasibility", "results/qc_analysis/tables/excluded_haplotype_feasibility.csv"
    ),
    "graph_provenance.tool_versions": NOVEL_PROVENANCE.get(
        "tool_versions", "results/metadata/tool_versions.tsv"
    ),
    "graph_provenance.graph_summary": NOVEL_PROVENANCE.get(
        "graph_summary", "results/metadata/sv_pangenome.graph_summary.tsv"
    ),
    "graph_provenance.build_log": NOVEL_PROVENANCE.get(
        "build_log", "results/logs/build_sv_pangenome_graph.log"
    ),
    "graph_provenance.qc_report": NOVEL_PROVENANCE.get(
        "qc_report", "qc_analysis/report/pangenome_qc_report.md"
    ),
    "references.CHM13": NOVEL.get("references", {}).get("CHM13", "../refs/CHM13.fa.gz"),
    "references.hg38": NOVEL.get("references", {}).get(
        "hg38", "../refs/GRCh38_primary_chr1-22_XY.fa"
    ),
}
def novel_validate_config_path(key, value):
    if os.path.isabs(str(value)):
        raise WorkflowError(
            f"novel_sv_discovery.{key} must be relative to whole_pans/Snakefile, got: {value}"
        )
    # Check the configured path lexically.  resolve() follows existing symlinks
    # and incorrectly rejects project-local links to server-side reference/data
    # storage outside the checkout.
    normalized = Path(os.path.abspath(NOVEL_PROJECT_ROOT / str(value)))
    if not normalized.is_relative_to(NOVEL_ALLOWED_ROOT):
        raise WorkflowError(
            f"novel_sv_discovery.{key} escapes the pangenome project tree: {value}"
        )


for key, value in NOVEL_PATHS.items():
    novel_validate_config_path(key, value)
for key in ("sample_metadata", "read_manifest"):
    value = NOVEL.get(key, "")
    if value:
        novel_validate_config_path(key, value)

NOVEL_OUTDIR = str(PurePosixPath(str(NOVEL_PATHS["results"])))
NOVEL_GRAPH = str(PurePosixPath(str(NOVEL_PATHS["frozen_graph"])))
NOVEL_GRAPH_ASSEMBLIES = str(PurePosixPath(str(NOVEL_PATHS["graph_assemblies"])))
NOVEL_CLEANED_MANIFEST = str(PurePosixPath(str(NOVEL_PATHS["cleaned_manifest"])))
NOVEL_ORIGINAL_QC = str(PurePosixPath(str(NOVEL_PATHS["original_qc"])))
NOVEL_POST_QC = str(PurePosixPath(str(NOVEL_PATHS["post_qc"])))
NOVEL_CONTAMINATION = str(PurePosixPath(str(NOVEL_PATHS["contamination"])))
NOVEL_PREVIOUS_FEASIBILITY = str(PurePosixPath(str(NOVEL_PATHS["previous_feasibility"])))
NOVEL_BUILD_TOOL_VERSIONS = str(
    PurePosixPath(str(NOVEL_PATHS["graph_provenance.tool_versions"]))
)
NOVEL_BUILD_GRAPH_SUMMARY = str(
    PurePosixPath(str(NOVEL_PATHS["graph_provenance.graph_summary"]))
)
NOVEL_BUILD_LOG = str(PurePosixPath(str(NOVEL_PATHS["graph_provenance.build_log"])))
NOVEL_QC_REPORT = str(PurePosixPath(str(NOVEL_PATHS["graph_provenance.qc_report"])))
NOVEL_REFERENCES = {
    reference: str(PurePosixPath(str(NOVEL_PATHS[f"references.{reference}"])))
    for reference in ("CHM13", "hg38")
}
NOVEL_GRAPH_ENV = str(workflow.source_path("../envs/sv_graph.yaml"))
NOVEL_REFERENCE_ENV = str(workflow.source_path("../envs/sv_reference.yaml"))
NOVEL_CATALOG_ENV = str(workflow.source_path("../envs/sv_catalog.yaml"))
# PAV is a source workflow with a broad native dependency set.  Keep it in a
# dedicated environment rather than coupling its solve to graph screening.
NOVEL_PAV_ENV = str(workflow.source_path("../envs/sv_discovery.yaml"))
NOVEL_SCREEN_SCRIPT = str(workflow.source_path("../scripts/screen_novel_graph_svs.py"))
NOVEL_MANIFEST_SCRIPT = str(workflow.source_path("../scripts/build_novel_sv_manifest.py"))
NOVEL_PAV_SCRIPT = str(workflow.source_path("../scripts/prepare_pav_inputs.py"))
NOVEL_MERGE_SCRIPT = str(workflow.source_path("../scripts/merge_novel_sv_catalog.py"))
NOVEL_HIFI_SCRIPT = str(workflow.source_path("../scripts/recommend_hifi_samples.py"))

NOVEL_DISCOVERY_MANIFEST = f"{NOVEL_OUTDIR}/manifest/assembly_discovery_manifest.tsv"
NOVEL_REFRESHED_FEASIBILITY = f"{NOVEL_OUTDIR}/manifest/refreshed_excluded_haplotype_feasibility.tsv"
NOVEL_PROVISIONAL_HIFI = f"{NOVEL_OUTDIR}/hifi/provisional_hifi_samples.tsv"
NOVEL_GRAPH_SEGMENTS = f"{NOVEL_OUTDIR}/frozen_graph/graph_segments.tsv.gz"
NOVEL_GRAPH_TASKS = f"{NOVEL_OUTDIR}/graph_screen/tasks.tsv"
NOVEL_GRAPH_SUMMARY = f"{NOVEL_OUTDIR}/graph_screen/summary/all_assembly_novel_sv_summary.tsv"
NOVEL_GRAPH_CANDIDATES = f"{NOVEL_OUTDIR}/graph_screen/summary/all_residual_sv_candidates.tsv.gz"
NOVEL_GRAPH_COMPLEX = f"{NOVEL_OUTDIR}/graph_screen/summary/all_complex_alignments.tsv.gz"
NOVEL_GRAPH_ALIGNMENTS = f"{NOVEL_OUTDIR}/graph_screen/summary/all_contig_alignments.tsv.gz"
NOVEL_CALL_MANIFEST = f"{NOVEL_OUTDIR}/reference_calls/reference_call_manifest.tsv"
NOVEL_CATALOG = (
    f"{NOVEL_OUTDIR}/catalog/provisional_per_coordinate_frame_assembly_sv_catalog.tsv"
)
NOVEL_EVIDENCE = f"{NOVEL_OUTDIR}/catalog/all_assembly_evidence.tsv.gz"
NOVEL_RECOMMENDED_HIFI = f"{NOVEL_OUTDIR}/hifi/recommended_hifi_samples.tsv"
NOVEL_TOOL_VERSIONS = f"{NOVEL_OUTDIR}/provenance/discovery_tool_versions.tsv"


def novel_local_path(path):
    candidate = Path(path)
    return candidate if candidate.is_absolute() else NOVEL_PROJECT_ROOT / candidate


def novel_open_table(path):
    local = novel_local_path(path)
    opener = gzip.open if str(local).endswith(".gz") else open
    delimiter = "," if str(local).lower().endswith(".csv") else "\t"
    with opener(local, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def novel_normalize_assembly_id(value):
    name = os.path.basename(str(value))
    for suffix in (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name[: -len(".clean")] if name.endswith(".clean") else name


NOVEL_ASSEMBLIES = {}
NOVEL_ASSEMBLY_METADATA = {}
NOVEL_ASSEMBLY_PATTERN = re.compile(
    r"^(?P<sample>.+?)\.hifi\.hifiasm\.bp\.(?P<haplotype>hap[12])\.p_ctg(?:\..*)?$"
)
for row in novel_open_table(NOVEL_CLEANED_MANIFEST):
    assembly_id = novel_normalize_assembly_id(row["assembly_id"])
    if assembly_id in NOVEL_ASSEMBLIES:
        raise WorkflowError(f"Duplicate assembly ID in {NOVEL_CLEANED_MANIFEST}: {assembly_id}")
    path = row.get("cleaned_fasta_path") or row.get("path") or row.get("cleaned_path")
    if not path:
        raise WorkflowError(f"No FASTA path for {assembly_id} in {NOVEL_CLEANED_MANIFEST}")
    match = NOVEL_ASSEMBLY_PATTERN.match(assembly_id)
    if not match:
        raise WorkflowError(f"Unexpected assembly ID in discovery manifest: {assembly_id}")
    NOVEL_ASSEMBLIES[assembly_id] = path
    NOVEL_ASSEMBLY_METADATA[assembly_id] = {
        "sample_id": match.group("sample"),
        "haplotype": match.group("haplotype"),
    }

NOVEL_ASSEMBLY_IDS = sorted(NOVEL_ASSEMBLIES)
NOVEL_SAMPLE_HAPS = defaultdict(dict)
for assembly_id, values in NOVEL_ASSEMBLY_METADATA.items():
    sample = values["sample_id"]
    haplotype = values["haplotype"]
    if haplotype in NOVEL_SAMPLE_HAPS[sample]:
        raise WorkflowError(f"Duplicate {sample} {haplotype} assembly")
    NOVEL_SAMPLE_HAPS[sample][haplotype] = assembly_id
NOVEL_PAIRED_SAMPLES = sorted(
    sample for sample, haplotypes in NOVEL_SAMPLE_HAPS.items() if set(haplotypes) == {"hap1", "hap2"}
)
NOVEL_UNPAIRED_SAMPLES = sorted(set(NOVEL_SAMPLE_HAPS) - set(NOVEL_PAIRED_SAMPLES))

NOVEL_BATCH_SIZE = int(NOVEL_SCREEN.get("batch_size", 1))
if NOVEL_BATCH_SIZE < 1:
    raise WorkflowError("novel_sv_discovery.graph_screen.batch_size must be at least 1")

NOVEL_SVIM_ENABLED = bool(NOVEL_CALLERS.get("svim_asm", True))
NOVEL_DIPCALL_ENABLED = bool(NOVEL_CALLERS.get("dipcall", True))
NOVEL_PAV_ENABLED = bool(NOVEL_CALLERS.get("pav", False))
if not any((NOVEL_SVIM_ENABLED, NOVEL_DIPCALL_ENABLED, NOVEL_PAV_ENABLED)):
    raise WorkflowError(
        "Enable at least one novel_sv_discovery.reference_calling caller; "
        "a graph-only catalog must use a separate, explicitly named target."
    )
NOVEL_PAV_SNAKEFILE = str(
    PurePosixPath(str(NOVEL_REFERENCE_CALLING.get("pav_snakefile", "../tools/pav/Snakefile")))
)
if os.path.isabs(NOVEL_PAV_SNAKEFILE):
    raise WorkflowError("novel_sv_discovery.reference_calling.pav_snakefile must be relative")
novel_validate_config_path("reference_calling.pav_snakefile", NOVEL_PAV_SNAKEFILE)

NOVEL_SEX_CONFIG = NOVEL_REFERENCE_CALLING.get("sex_chromosomes", {})
NOVEL_SEX_MODE = str(NOVEL_SEX_CONFIG.get("mode", "autosomes_only"))
if NOVEL_SEX_MODE not in {"autosomes_only", "sex_aware"}:
    raise WorkflowError(
        "novel_sv_discovery.reference_calling.sex_chromosomes.mode must be "
        "autosomes_only or sex_aware"
    )
NOVEL_SEX_CONTIG_REGEX = str(
    NOVEL_SEX_CONFIG.get("exclude_contig_regex", r"(?:chr)?[XY]")
)
if not NOVEL_SEX_CONTIG_REGEX:
    raise WorkflowError(
        "novel_sv_discovery.reference_calling.sex_chromosomes.exclude_contig_regex "
        "must not be empty"
    )
try:
    re.compile(NOVEL_SEX_CONTIG_REGEX)
except re.error as error:
    raise WorkflowError(
        "Invalid sex_chromosomes.exclude_contig_regex: "
        f"{NOVEL_SEX_CONTIG_REGEX!r}: {error}"
    )
NOVEL_EXCLUDE_CONTIG_REGEX = (
    NOVEL_SEX_CONTIG_REGEX
    if NOVEL_SEX_MODE == "autosomes_only"
    else ""
)
NOVEL_PAR_BEDS = {
    reference: str(NOVEL_SEX_CONFIG.get("par_bed", {}).get(reference, ""))
    for reference in NOVEL_REFERENCES
}
NOVEL_CHRY_PAR_HARD_MASKED = {
    reference: (
        NOVEL_SEX_CONFIG.get("chrY_par_hard_masked", {}).get(reference, False)
        is True
    )
    for reference in NOVEL_REFERENCES
}
for reference, path in NOVEL_PAR_BEDS.items():
    if path:
        novel_validate_config_path(f"reference_calling.sex_chromosomes.par_bed.{reference}", path)

NOVEL_SAMPLE_SEX = {}
NOVEL_PATERNAL_HAPLOTYPE = {}
if NOVEL_SEX_MODE == "sex_aware":
    metadata_path = NOVEL.get("sample_metadata", "")
    if not metadata_path:
        raise WorkflowError(
            "sex_aware calling requires novel_sv_discovery.sample_metadata with "
            "sample_id, sex, and paternal_haplotype for male samples"
        )
    for row in novel_open_table(metadata_path):
        sample = row.get("sample_id") or row.get("sample")
        if not sample:
            raise WorkflowError(f"Missing sample_id in {metadata_path}")
        sex = str(row.get("sex", row.get("reported_sex", ""))).strip().lower()
        if sex in {"m", "male", "1", "xy"}:
            NOVEL_SAMPLE_SEX[sample] = "male"
        elif sex in {"f", "female", "2", "xx"}:
            NOVEL_SAMPLE_SEX[sample] = "female"
        else:
            raise WorkflowError(f"Missing or unrecognized sex for {sample}: {sex!r}")
        paternal = str(
            row.get("paternal_haplotype", row.get("paternal_hap", ""))
        ).strip().lower()
        if NOVEL_SAMPLE_SEX[sample] == "male":
            if paternal not in {"hap1", "hap2"}:
                raise WorkflowError(
                    f"Male sample {sample} needs paternal_haplotype=hap1 or hap2 "
                    "for dipcall's paternal-first convention"
                )
            NOVEL_PATERNAL_HAPLOTYPE[sample] = paternal
    missing_sex = sorted(set(NOVEL_PAIRED_SAMPLES) - set(NOVEL_SAMPLE_SEX))
    if missing_sex:
        raise WorkflowError(
            "Missing sex metadata for paired samples: " + ", ".join(missing_sex[:20])
        )
    if NOVEL_DIPCALL_ENABLED and any(
        sex == "male" for sex in NOVEL_SAMPLE_SEX.values()
    ):
        missing_par = [reference for reference, path in NOVEL_PAR_BEDS.items() if not path]
        if missing_par:
            raise WorkflowError(
                "Male dipcall requires reference-specific PAR BED files for: "
                + ", ".join(missing_par)
            )
        unmasked_references = [
            reference
            for reference, is_masked in NOVEL_CHRY_PAR_HARD_MASKED.items()
            if not is_masked
        ]
        if unmasked_references:
            raise WorkflowError(
                "Male dipcall requires PAR sequence to be hard-masked on reference "
                "chrY. Point novel_sv_discovery.references at verified masked FASTAs "
                "and set reference_calling.sex_chromosomes.chrY_par_hard_masked=true "
                "for: "
                + ", ".join(unmasked_references)
            )


def novel_assembly_input(wildcards):
    return NOVEL_ASSEMBLIES[wildcards.assembly]


def novel_reference_input(wildcards):
    return NOVEL_REFERENCES[wildcards.reference]


def novel_prepared_reference(reference):
    return f"{NOVEL_OUTDIR}/reference_calls/{reference}/reference.fa"


def novel_prepared_reference_input(wildcards):
    return novel_prepared_reference(wildcards.reference)


def novel_haplotype_input(wildcards, haplotype):
    assembly_id = NOVEL_SAMPLE_HAPS[wildcards.sample][haplotype]
    return NOVEL_ASSEMBLIES[assembly_id]


def novel_haplotype_alignment(wildcards, haplotype):
    assembly_id = NOVEL_SAMPLE_HAPS[wildcards.sample][haplotype]
    return f"{NOVEL_OUTDIR}/reference_calls/{wildcards.reference}/alignments/{assembly_id}.bam"


def novel_dipcall_haplotype_input(wildcards, order):
    if NOVEL_SEX_MODE == "sex_aware" and NOVEL_SAMPLE_SEX[wildcards.sample] == "male":
        paternal = NOVEL_PATERNAL_HAPLOTYPE[wildcards.sample]
        haplotype = paternal if order == "first" else ("hap2" if paternal == "hap1" else "hap1")
    else:
        haplotype = "hap1" if order == "first" else "hap2"
    return NOVEL_ASSEMBLIES[NOVEL_SAMPLE_HAPS[wildcards.sample][haplotype]]


def novel_dipcall_par_input(wildcards):
    if (
        NOVEL_SEX_MODE == "sex_aware"
        and NOVEL_SAMPLE_SEX[wildcards.sample] == "male"
    ):
        return NOVEL_PAR_BEDS[wildcards.reference]
    return []


def novel_dipcall_sex_args(wildcards):
    par = novel_dipcall_par_input(wildcards)
    return f"-x {shlex.quote(str(par))}" if par else ""


def novel_svim_vcf(reference, sample):
    return f"{NOVEL_OUTDIR}/reference_calls/{reference}/svim_asm/{sample}.vcf.gz"


def novel_dipcall_vcf(reference, sample):
    return f"{NOVEL_OUTDIR}/reference_calls/{reference}/dipcall/{sample}.vcf.gz"


def novel_dipcall_bed(reference, sample):
    return f"{NOVEL_OUTDIR}/reference_calls/{reference}/dipcall/{sample}.callable.bed.gz"


def novel_pav_vcf(reference, sample):
    return f"{NOVEL_OUTDIR}/reference_calls/{reference}/pav/{sample}/pav_{sample}.vcf.gz"


def novel_caller_exclude_contig_regex(caller):
    """Keep sex-aware X/Y calls only from the PAR-aware caller."""
    if NOVEL_SEX_MODE == "autosomes_only" or caller != "dipcall":
        return NOVEL_SEX_CONTIG_REGEX
    return ""


NOVEL_REFERENCE_CALL_RECORDS = []
if NOVEL_SVIM_ENABLED:
    for reference in NOVEL_REFERENCES:
        for sample in NOVEL_PAIRED_SAMPLES:
            NOVEL_REFERENCE_CALL_RECORDS.append(
                {
                    "caller": "svim_asm",
                    "coordinate_system": reference,
                    "sample_id": sample,
                    "assembly_id": "",
                    "haplotype": "diploid",
                    "path": novel_svim_vcf(reference, sample),
                    "callable_bed": "",
                    "exclude_contig_regex": novel_caller_exclude_contig_regex("svim_asm"),
                }
            )
if NOVEL_DIPCALL_ENABLED:
    for reference in NOVEL_REFERENCES:
        for sample in NOVEL_PAIRED_SAMPLES:
            NOVEL_REFERENCE_CALL_RECORDS.append(
                {
                    "caller": "dipcall",
                    "coordinate_system": reference,
                    "sample_id": sample,
                    "assembly_id": "",
                    "haplotype": "diploid",
                    "path": novel_dipcall_vcf(reference, sample),
                    "callable_bed": novel_dipcall_bed(reference, sample),
                    "exclude_contig_regex": novel_caller_exclude_contig_regex("dipcall"),
                }
            )
if NOVEL_PAV_ENABLED:
    for reference in NOVEL_REFERENCES:
        for sample in NOVEL_PAIRED_SAMPLES:
            NOVEL_REFERENCE_CALL_RECORDS.append(
                {
                    "caller": "pav",
                    "coordinate_system": reference,
                    "sample_id": sample,
                    "assembly_id": "",
                    "haplotype": "diploid",
                    "path": novel_pav_vcf(reference, sample),
                    "callable_bed": "",
                    "exclude_contig_regex": novel_caller_exclude_contig_regex("pav"),
                }
            )
NOVEL_REFERENCE_OUTPUTS = [row["path"] for row in NOVEL_REFERENCE_CALL_RECORDS]
NOVEL_REFERENCE_OUTPUTS.extend(
    row["callable_bed"] for row in NOVEL_REFERENCE_CALL_RECORDS if row["callable_bed"]
)


localrules: novel_sv_discovery, write_novel_sv_reference_call_manifest

wildcard_constraints:
    reference="CHM13|hg38",
    task="[0-9]+"


rule novel_sv_discovery:
    input:
        graph_inventory=f"{NOVEL_OUTDIR}/frozen_graph/graph_inventory.tsv",
        graph_versions=f"{NOVEL_OUTDIR}/frozen_graph/tool_versions.tsv",
        graph_bubbles=f"{NOVEL_OUTDIR}/frozen_graph/existing_bubbles.tsv",
        graph_inputs=f"{NOVEL_OUTDIR}/frozen_graph/ordered_assemblies.tsv",
        build_versions=f"{NOVEL_OUTDIR}/frozen_graph/build_tool_versions.tsv",
        build_summary=f"{NOVEL_OUTDIR}/frozen_graph/build_graph_summary.tsv",
        build_log=f"{NOVEL_OUTDIR}/frozen_graph/build.log",
        qc_report=f"{NOVEL_OUTDIR}/frozen_graph/pangenome_qc_report.md",
        manifest=NOVEL_DISCOVERY_MANIFEST,
        refreshed_feasibility=NOVEL_REFRESHED_FEASIBILITY,
        graph_summary=NOVEL_GRAPH_SUMMARY,
        graph_candidates=NOVEL_GRAPH_CANDIDATES,
        graph_complex=NOVEL_GRAPH_COMPLEX,
        graph_alignments=NOVEL_GRAPH_ALIGNMENTS,
        call_manifest=NOVEL_CALL_MANIFEST,
        catalog=NOVEL_CATALOG,
        evidence=NOVEL_EVIDENCE,
        provisional_hifi=NOVEL_PROVISIONAL_HIFI,
        recommended_hifi=NOVEL_RECOMMENDED_HIFI,
        discovery_versions=NOVEL_TOOL_VERSIONS,


rule index_frozen_graph_for_novel_sv:
    input:
        gfa=NOVEL_GRAPH
    output:
        NOVEL_GRAPH_SEGMENTS
    params:
        script=NOVEL_SCREEN_SCRIPT
    log:
        f"{NOVEL_OUTDIR}/logs/index_frozen_graph.log"
    conda:
        NOVEL_GRAPH_ENV
    threads: 1
    resources:
        mem_mb=8000,
        runtime=120
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {params.script:q} index-graph \
          --gfa {input.gfa:q} \
          --output {output:q} > {log:q} 2>&1
        test -s {output:q}
        """


rule inventory_frozen_graph_for_novel_sv:
    input:
        gfa=NOVEL_GRAPH,
        ordered=NOVEL_GRAPH_ASSEMBLIES,
        build_versions=NOVEL_BUILD_TOOL_VERSIONS,
        build_summary=NOVEL_BUILD_GRAPH_SUMMARY,
        build_log=NOVEL_BUILD_LOG,
        qc_report=NOVEL_QC_REPORT
    output:
        inventory=f"{NOVEL_OUTDIR}/frozen_graph/graph_inventory.tsv",
        versions=f"{NOVEL_OUTDIR}/frozen_graph/tool_versions.tsv",
        bubbles=f"{NOVEL_OUTDIR}/frozen_graph/existing_bubbles.tsv",
        ordered=f"{NOVEL_OUTDIR}/frozen_graph/ordered_assemblies.tsv",
        build_versions=f"{NOVEL_OUTDIR}/frozen_graph/build_tool_versions.tsv",
        build_summary=f"{NOVEL_OUTDIR}/frozen_graph/build_graph_summary.tsv",
        build_log=f"{NOVEL_OUTDIR}/frozen_graph/build.log",
        qc_report=f"{NOVEL_OUTDIR}/frozen_graph/pangenome_qc_report.md"
    log:
        f"{NOVEL_OUTDIR}/logs/inventory_frozen_graph.log"
    conda:
        NOVEL_GRAPH_ENV
    threads: 1
    resources:
        mem_mb=16000,
        runtime=240
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.inventory:q})" "$(dirname {log:q})"
        {{
          printf 'metric\tvalue\n'
          printf 'path\t%s\n' {input.gfa:q}
          printf 'size_bytes\t%s\n' "$(stat -c %s {input.gfa:q})"
          printf 'sha256\t%s\n' "$(sha256sum {input.gfa:q} | awk '{{print $1}}')"
          printf 'ordered_assemblies_sha256\t%s\n' "$(sha256sum {input.ordered:q} | awk '{{print $1}}')"
          printf 'build_tool_versions_sha256\t%s\n' "$(sha256sum {input.build_versions:q} | awk '{{print $1}}')"
          printf 'build_graph_summary_sha256\t%s\n' "$(sha256sum {input.build_summary:q} | awk '{{print $1}}')"
          printf 'build_log_sha256\t%s\n' "$(sha256sum {input.build_log:q} | awk '{{print $1}}')"
          printf 'qc_report_sha256\t%s\n' "$(sha256sum {input.qc_report:q} | awk '{{print $1}}')"
          printf 'segments\t%s\n' "$(awk -F '\t' '$1 == "S" {{n++}} END {{print n+0}}' {input.gfa:q})"
          printf 'links\t%s\n' "$(awk -F '\t' '$1 == "L" {{n++}} END {{print n+0}}' {input.gfa:q})"
        }} > {output.inventory:q}
        minigraph_version=$(minigraph --version 2>&1 || true)
        gfatools_version=$({{ gfatools version 2>&1 || gfatools 2>&1 || true; }})
        {{
          printf 'tool\tversion\n'
          printf 'minigraph\t'; awk 'NF {{print; found=1; exit}} END {{if (!found) print "unavailable"}}' <<< "$minigraph_version"
          printf 'gfatools\t'; awk 'NF {{print; found=1; exit}} END {{if (!found) print "unavailable"}}' <<< "$gfatools_version"
        }} > {output.versions:q}
        gfatools bubble {input.gfa:q} > {output.bubbles:q} 2> {log:q}
        cp {input.ordered:q} {output.ordered:q}
        cp {input.build_versions:q} {output.build_versions:q}
        cp {input.build_summary:q} {output.build_summary:q}
        cp {input.build_log:q} {output.build_log:q}
        cp {input.qc_report:q} {output.qc_report:q}
        test -s {output.inventory:q}
        """


rule record_novel_sv_tool_versions:
    input:
        graph_env=NOVEL_GRAPH_ENV,
        reference_env=NOVEL_REFERENCE_ENV,
        catalog_env=NOVEL_CATALOG_ENV,
        pav_env=NOVEL_PAV_ENV,
        pav=lambda wildcards: [NOVEL_PAV_SNAKEFILE] if NOVEL_PAV_ENABLED else []
    output:
        NOVEL_TOOL_VERSIONS
    params:
        pav_enabled=str(NOVEL_PAV_ENABLED).lower(),
        pav_snakefile=NOVEL_PAV_SNAKEFILE
    log:
        f"{NOVEL_OUTDIR}/logs/record_discovery_tool_versions.log"
    conda:
        NOVEL_REFERENCE_ENV
    threads: 1
    resources:
        mem_mb=2000,
        runtime=30
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        {{
          printf 'tool\tversion\texecutable\texecutable_sha256\n'
          for tool in minigraph gfatools minimap2 samtools bcftools tabix svim-asm run-dipcall snakemake; do
            executable=$(command -v "$tool" || true)
            if [ -z "$executable" ]; then
              printf '%s\tunavailable\t\t\n' "$tool"
              continue
            fi
            case "$tool" in
              gfatools) version=$(gfatools version 2>&1 | awk 'NF {{print; exit}}' || true) ;;
              run-dipcall) version='0.3 (environment pin)' ;;
              *) version=$("$tool" --version 2>&1 | awk 'NF {{print; exit}}' || true) ;;
            esac
            [ -n "$version" ] || version='version flag unavailable'
            checksum=$(sha256sum "$executable" | awk '{{print $1}}')
            version=${{version//$'\t'/ }}
            printf '%s\t%s\t%s\t%s\n' "$tool" "$version" "$executable" "$checksum"
          done
          for environment in {input.graph_env:q} {input.reference_env:q} {input.catalog_env:q} {input.pav_env:q}; do
            printf 'conda_environment\tsha256:%s\t%s\t\n' \
              "$(sha256sum "$environment" | awk '{{print $1}}')" "$environment"
          done
          if [ {params.pav_enabled:q} = true ]; then
            printf 'pav_snakefile\tsha256:%s\t%s\t\n' \
              "$(sha256sum {params.pav_snakefile:q} | awk '{{print $1}}')" {params.pav_snakefile:q}
          fi
        }} > {output:q} 2> {log:q}
        test -s {output:q}
        """


rule build_novel_sv_discovery_manifest:
    input:
        cleaned=NOVEL_CLEANED_MANIFEST,
        graph_assemblies=NOVEL_GRAPH_ASSEMBLIES,
        original_qc=NOVEL_ORIGINAL_QC,
        post_qc=NOVEL_POST_QC,
        contamination=NOVEL_CONTAMINATION,
        previous_feasibility=NOVEL_PREVIOUS_FEASIBILITY,
        optional=lambda wildcards: [
            path for path in (NOVEL.get("sample_metadata", ""), NOVEL.get("read_manifest", ""))
            if path
        ]
    output:
        manifest=NOVEL_DISCOVERY_MANIFEST,
        refreshed=NOVEL_REFRESHED_FEASIBILITY,
        provisional_hifi=NOVEL_PROVISIONAL_HIFI
    params:
        script=NOVEL_MANIFEST_SCRIPT,
        optional_args=lambda wildcards: " ".join(
            argument
            for argument in (
                "--sample-metadata " + shlex.quote(str(NOVEL.get("sample_metadata")))
                if NOVEL.get("sample_metadata") else "",
                "--read-manifest " + shlex.quote(str(NOVEL.get("read_manifest")))
                if NOVEL.get("read_manifest") else "",
            )
            if argument
        )
    log:
        f"{NOVEL_OUTDIR}/logs/build_discovery_manifest.log"
    conda:
        NOVEL_CATALOG_ENV
    threads: 1
    resources:
        mem_mb=8000,
        runtime=120
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.manifest:q})" "$(dirname {log:q})"
        python {params.script:q} \
          --cleaned-manifest {input.cleaned:q} \
          --graph-assemblies {input.graph_assemblies:q} \
          --original-qc {input.original_qc:q} \
          --post-qc {input.post_qc:q} \
          --contamination {input.contamination:q} \
          --previous-feasibility {input.previous_feasibility:q} \
          --output {output.manifest:q} \
          --refreshed-feasibility {output.refreshed:q} \
          --provisional-hifi {output.provisional_hifi:q} \
          {params.optional_args} > {log:q} 2>&1
        test -s {output.manifest:q}
        """


checkpoint create_novel_sv_graph_tasks:
    input:
        manifest=NOVEL_DISCOVERY_MANIFEST,
        graph_assemblies=NOVEL_GRAPH_ASSEMBLIES
    output:
        tasks=NOVEL_GRAPH_TASKS
    params:
        script=NOVEL_SCREEN_SCRIPT,
        batch_size=NOVEL_BATCH_SIZE
    log:
        f"{NOVEL_OUTDIR}/logs/create_graph_tasks.log"
    conda:
        NOVEL_GRAPH_ENV
    threads: 1
    resources:
        mem_mb=4000,
        runtime=60
    shell:
        r"""
        set -euo pipefail
        python {params.script:q} tasks \
          --manifest {input.manifest:q} \
          --graph-assemblies {input.graph_assemblies:q} \
          --batch-size {params.batch_size} \
          --output {output.tasks:q} > {log:q} 2>&1
        test -s {output.tasks:q}
        """


def novel_graph_task_directories(wildcards):
    checkpoint_output = checkpoints.create_novel_sv_graph_tasks.get(**wildcards).output.tasks
    rows = novel_open_table(checkpoint_output)
    task_ids = sorted({row["task_id"] for row in rows})
    if not task_ids:
        raise WorkflowError(f"No graph-screen tasks found in {checkpoint_output}")
    return expand(
        f"{NOVEL_OUTDIR}/graph_screen/tasks/{{task}}",
        task=task_ids,
    )


rule screen_assembly_against_frozen_graph:
    input:
        graph=NOVEL_GRAPH,
        segment_index=NOVEL_GRAPH_SEGMENTS,
        manifest=NOVEL_DISCOVERY_MANIFEST,
        graph_assemblies=NOVEL_GRAPH_ASSEMBLIES,
        tasks=NOVEL_GRAPH_TASKS
    output:
        directory(f"{NOVEL_OUTDIR}/graph_screen/tasks/{{task}}")
    params:
        script=NOVEL_SCREEN_SCRIPT,
        primary_chain=int(NOVEL_SCREEN.get("primary_min_chain_score", 5000)),
        primary_secondary=int(NOVEL_SCREEN.get("primary_secondary", 5)),
        primary_extra=shlex.quote(str(NOVEL_SCREEN.get("minigraph_extra", ""))),
        rescue_fraction=float(NOVEL_SCREEN.get("rescue_callable_fraction", 0.85)),
        sensitivity_chain=int(NOVEL_SCREEN.get("sensitivity_min_chain_score", 1000)),
        sensitivity_secondary=int(NOVEL_SCREEN.get("sensitivity_secondary", 20)),
        sensitivity_extra=shlex.quote(
            str(NOVEL_SCREEN.get("sensitivity_minigraph_extra", ""))
        ),
        min_sv=int(NOVEL_SCREEN.get("min_sv_size", 50)),
        raw_min=int(NOVEL_SCREEN.get("raw_min_size", 30)),
        indel_gap=int(NOVEL_SCREEN.get("adjacent_indel_gap", 50)),
        min_alignment=int(NOVEL_SCREEN.get("min_alignment", 5000)),
        min_anchor=int(NOVEL_SCREEN.get("min_anchor", 2000)),
        min_mapq=int(NOVEL_SCREEN.get("min_mapq", 5)),
        min_identity=float(NOVEL_SCREEN.get("min_identity", 0.90))
    log:
        f"{NOVEL_OUTDIR}/graph_screen/task_logs/{{task}}.log"
    benchmark:
        f"{NOVEL_OUTDIR}/benchmarks/graph_screen/{{task}}.tsv"
    conda:
        NOVEL_GRAPH_ENV
    threads: int(NOVEL_SCREEN.get("threads", 16))
    resources:
        mem_mb=int(NOVEL_SCREEN.get("mem_mb", 96000)),
        runtime=int(NOVEL_SCREEN.get("runtime", NOVEL_SCREEN.get("runtime_min", 360)))
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {params.script:q} run \
          --graph {input.graph:q} \
          --segment-index {input.segment_index:q} \
          --manifest {input.manifest:q} \
          --graph-assemblies {input.graph_assemblies:q} \
          --tasks {input.tasks:q} \
          --output-dir {output:q} \
          --task-id {wildcards.task:q} \
          --threads {threads} \
          --primary-min-chain-score {params.primary_chain} \
          --primary-secondary {params.primary_secondary} \
          --minigraph-extra {params.primary_extra} \
          --rescue-callable-fraction {params.rescue_fraction} \
          --sensitivity-min-chain-score {params.sensitivity_chain} \
          --sensitivity-secondary {params.sensitivity_secondary} \
          --sensitivity-minigraph-extra {params.sensitivity_extra} \
          --min-sv-size {params.min_sv} \
          --raw-min-size {params.raw_min} \
          --adjacent-indel-gap {params.indel_gap} \
          --min-alignment {params.min_alignment} \
          --min-anchor {params.min_anchor} \
          --min-mapq {params.min_mapq} \
          --min-identity {params.min_identity} > {log:q} 2>&1
        test -s {output:q}/.complete
        test -s {output:q}/task_outputs.tsv
        """


rule summarize_novel_sv_graph_screen:
    input:
        manifest=NOVEL_DISCOVERY_MANIFEST,
        tasks=NOVEL_GRAPH_TASKS,
        task_dirs=novel_graph_task_directories
    output:
        summary=NOVEL_GRAPH_SUMMARY,
        candidates=NOVEL_GRAPH_CANDIDATES,
        complex=NOVEL_GRAPH_COMPLEX,
        alignments=NOVEL_GRAPH_ALIGNMENTS
    params:
        script=NOVEL_SCREEN_SCRIPT,
        output_dir=f"{NOVEL_OUTDIR}/graph_screen",
        task_dir_args=lambda wildcards, input: " ".join(
            f"--task-dir {shlex.quote(str(path))}" for path in input.task_dirs
        )
    log:
        f"{NOVEL_OUTDIR}/logs/summarize_graph_screen.log"
    conda:
        NOVEL_GRAPH_ENV
    threads: 1
    resources:
        mem_mb=16000,
        runtime=240
    shell:
        r"""
        set -euo pipefail
        python {params.script:q} summarize \
          --manifest {input.manifest:q} \
          --tasks {input.tasks:q} \
          {params.task_dir_args} \
          --output-dir {params.output_dir:q} > {log:q} 2>&1
        test -s {output.summary:q}
        """


rule prepare_novel_sv_linear_reference:
    input:
        novel_reference_input
    output:
        fasta=f"{NOVEL_OUTDIR}/reference_calls/{{reference}}/reference.fa",
        fai=f"{NOVEL_OUTDIR}/reference_calls/{{reference}}/reference.fa.fai",
        dict=f"{NOVEL_OUTDIR}/reference_calls/{{reference}}/reference.dict"
    log:
        f"{NOVEL_OUTDIR}/logs/reference_calls/{{reference}}.prepare_reference.log"
    benchmark:
        f"{NOVEL_OUTDIR}/benchmarks/reference_calls/{{reference}}.prepare_reference.tsv"
    conda:
        NOVEL_REFERENCE_ENV
    threads: 1
    resources:
        mem_mb=4000,
        runtime=120
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.fasta:q})" "$(dirname {log:q})"
        tmp={output.fasta:q}.tmp.$$
        trap 'rm -f "$tmp" "$tmp.fai" "$tmp.dict"' EXIT
        case {input:q} in
          *.gz) gzip -dc {input:q} > "$tmp" ;;
          *) cp {input:q} "$tmp" ;;
        esac
        test -s "$tmp"
        samtools faidx "$tmp" > {log:q} 2>&1
        samtools dict -o "$tmp.dict" "$tmp" >> {log:q} 2>&1
        mv "$tmp" {output.fasta:q}
        mv "$tmp.fai" {output.fai:q}
        mv "$tmp.dict" {output.dict:q}
        trap - EXIT
        test -s {output.fasta:q}
        test -s {output.fai:q}
        test -s {output.dict:q}
        """


rule index_novel_sv_linear_reference:
    input:
        fasta=novel_prepared_reference_input,
        fai=f"{NOVEL_OUTDIR}/reference_calls/{{reference}}/reference.fa.fai"
    output:
        f"{NOVEL_OUTDIR}/reference_calls/{{reference}}/reference.mmi"
    params:
        preset=str(NOVEL_REFERENCE_CALLING.get("alignment_preset", "asm5"))
    log:
        f"{NOVEL_OUTDIR}/logs/reference_calls/{{reference}}.index.log"
    benchmark:
        f"{NOVEL_OUTDIR}/benchmarks/reference_calls/{{reference}}.index.tsv"
    conda:
        NOVEL_REFERENCE_ENV
    threads: int(NOVEL_REFERENCE_CALLING.get("threads", 16))
    resources:
        mem_mb=int(NOVEL_REFERENCE_CALLING.get("alignment_mem_mb", 64000)),
        runtime=int(
            NOVEL_REFERENCE_CALLING.get(
                "alignment_runtime", NOVEL_REFERENCE_CALLING.get("alignment_runtime_min", 360)
            )
        )
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        minimap2 -x {params.preset:q} -t {threads} -d {output:q} {input.fasta:q} 2> {log:q}
        test -s {output:q}
        """


rule align_assembly_for_svim_asm:
    input:
        assembly=novel_assembly_input,
        reference=novel_prepared_reference_input,
        fai=f"{NOVEL_OUTDIR}/reference_calls/{{reference}}/reference.fa.fai",
        index=f"{NOVEL_OUTDIR}/reference_calls/{{reference}}/reference.mmi"
    output:
        bam=temp(f"{NOVEL_OUTDIR}/reference_calls/{{reference}}/alignments/{{assembly}}.bam"),
        bai=temp(f"{NOVEL_OUTDIR}/reference_calls/{{reference}}/alignments/{{assembly}}.bam.bai")
    log:
        f"{NOVEL_OUTDIR}/logs/reference_calls/{{reference}}.{{assembly}}.align.log"
    benchmark:
        f"{NOVEL_OUTDIR}/benchmarks/reference_calls/{{reference}}.{{assembly}}.align.tsv"
    params:
        preset=str(NOVEL_REFERENCE_CALLING.get("alignment_preset", "asm5")),
        sort_mem_mb=int(NOVEL_REFERENCE_CALLING.get("sort_mem_mb_per_thread", 2048))
    conda:
        NOVEL_REFERENCE_ENV
    threads: int(NOVEL_REFERENCE_CALLING.get("threads", 16))
    resources:
        mem_mb=int(NOVEL_REFERENCE_CALLING.get("alignment_mem_mb", 64000)),
        runtime=int(
            NOVEL_REFERENCE_CALLING.get(
                "alignment_runtime", NOVEL_REFERENCE_CALLING.get("alignment_runtime_min", 360)
            )
        )
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.bam:q})" "$(dirname {log:q})"
        minimap2 -a -x {params.preset:q} --cs -r2k -t {threads} {input.index:q} {input.assembly:q} 2> {log:q} \
          | samtools sort -m {params.sort_mem_mb}M -@ {threads} -o {output.bam:q} - 2>> {log:q}
        samtools index -@ {threads} {output.bam:q} {output.bai:q} 2>> {log:q}
        """


rule call_svs_with_svim_asm:
    input:
        hap1_bam=lambda wildcards: novel_haplotype_alignment(wildcards, "hap1"),
        hap1_bai=lambda wildcards: novel_haplotype_alignment(wildcards, "hap1") + ".bai",
        hap2_bam=lambda wildcards: novel_haplotype_alignment(wildcards, "hap2"),
        hap2_bai=lambda wildcards: novel_haplotype_alignment(wildcards, "hap2") + ".bai",
        reference=novel_prepared_reference_input,
        fai=f"{NOVEL_OUTDIR}/reference_calls/{{reference}}/reference.fa.fai"
    output:
        vcf=f"{NOVEL_OUTDIR}/reference_calls/{{reference}}/svim_asm/{{sample}}.vcf.gz",
        tbi=f"{NOVEL_OUTDIR}/reference_calls/{{reference}}/svim_asm/{{sample}}.vcf.gz.tbi"
    log:
        f"{NOVEL_OUTDIR}/logs/reference_calls/{{reference}}.{{sample}}.svim_asm.log"
    benchmark:
        f"{NOVEL_OUTDIR}/benchmarks/reference_calls/{{reference}}.{{sample}}.svim_asm.tsv"
    conda:
        NOVEL_REFERENCE_ENV
    threads: 1
    resources:
        mem_mb=int(NOVEL_REFERENCE_CALLING.get("caller_mem_mb", 64000)),
        runtime=int(
            NOVEL_REFERENCE_CALLING.get(
                "caller_runtime", NOVEL_REFERENCE_CALLING.get("caller_runtime_min", 720)
            )
        )
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.vcf:q})" "$(dirname {log:q})"
        scratch_root="${{TMPDIR:-$(dirname {output.vcf:q})}}"
        workdir=$(mktemp -d -p "$scratch_root" "svim-asm.{wildcards.reference}.{wildcards.sample}.XXXXXX")
        trap 'rm -rf "$workdir"' EXIT
        svim-asm diploid "$workdir" {input.hap1_bam:q} {input.hap2_bam:q} {input.reference:q} > {log:q} 2>&1
        bcftools sort "$workdir/variants.vcf" -Oz -o {output.vcf:q} >> {log:q} 2>&1
        tabix -f -p vcf {output.vcf:q}
        """


rule call_svs_with_dipcall:
    input:
        reference=novel_prepared_reference_input,
        fai=f"{NOVEL_OUTDIR}/reference_calls/{{reference}}/reference.fa.fai",
        index=f"{NOVEL_OUTDIR}/reference_calls/{{reference}}/reference.mmi",
        first=lambda wildcards: novel_dipcall_haplotype_input(wildcards, "first"),
        second=lambda wildcards: novel_dipcall_haplotype_input(wildcards, "second"),
        par=novel_dipcall_par_input
    output:
        vcf=f"{NOVEL_OUTDIR}/reference_calls/{{reference}}/dipcall/{{sample}}.vcf.gz",
        bed=f"{NOVEL_OUTDIR}/reference_calls/{{reference}}/dipcall/{{sample}}.callable.bed.gz"
    log:
        f"{NOVEL_OUTDIR}/logs/reference_calls/{{reference}}.{{sample}}.dipcall.log"
    benchmark:
        f"{NOVEL_OUTDIR}/benchmarks/reference_calls/{{reference}}.{{sample}}.dipcall.tsv"
    params:
        sex_args=novel_dipcall_sex_args
    conda:
        NOVEL_REFERENCE_ENV
    threads: int(NOVEL_REFERENCE_CALLING.get("threads", 16))
    resources:
        mem_mb=int(NOVEL_REFERENCE_CALLING.get("caller_mem_mb", 64000)),
        runtime=int(
            NOVEL_REFERENCE_CALLING.get(
                "caller_runtime", NOVEL_REFERENCE_CALLING.get("caller_runtime_min", 720)
            )
        )
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.vcf:q})" "$(dirname {log:q})"
        scratch_root="${{TMPDIR:-$(dirname {output.vcf:q})}}"
        workdir=$(mktemp -d -p "$scratch_root" "dipcall.{wildcards.reference}.{wildcards.sample}.XXXXXX")
        trap 'rm -rf "$workdir"' EXIT
        prefix="$workdir/call"
        run-dipcall \
          -t {threads} \
          -d {input.index:q} \
          {params.sex_args} \
          "$prefix" \
          {input.reference:q} \
          {input.first:q} \
          {input.second:q} > "$workdir/call.mak" 2> {log:q}
        # Each independent Dipcall Make target can launch a Minimap2 process
        # using all declared threads. Run targets serially to avoid multiplying
        # both CPU and reference-index memory beyond this rule's allocation.
        make -j1 -f "$workdir/call.mak" >> {log:q} 2>&1
        cp "$prefix.dip.vcf.gz" {output.vcf:q}
        if [ -s "$prefix.dip.bed.gz" ]; then
            cp "$prefix.dip.bed.gz" {output.bed:q}
        else
            bgzip -c "$prefix.dip.bed" > {output.bed:q}
        fi
        test -s {output.vcf:q}
        test -s {output.bed:q}
        """


rule prepare_pav_sample_inputs:
    input:
        reference=novel_prepared_reference_input,
        fai=f"{NOVEL_OUTDIR}/reference_calls/{{reference}}/reference.fa.fai",
        hap1=lambda wildcards: novel_haplotype_input(wildcards, "hap1"),
        hap2=lambda wildcards: novel_haplotype_input(wildcards, "hap2")
    output:
        config=f"{NOVEL_OUTDIR}/reference_calls/{{reference}}/pav/{{sample}}/config.json",
        assemblies=f"{NOVEL_OUTDIR}/reference_calls/{{reference}}/pav/{{sample}}/assemblies.tsv"
    params:
        script=NOVEL_PAV_SCRIPT,
        extra=json.dumps(NOVEL_REFERENCE_CALLING.get("pav_extra_config", {}), sort_keys=True)
    log:
        f"{NOVEL_OUTDIR}/logs/reference_calls/{{reference}}.{{sample}}.prepare_pav.log"
    conda:
        NOVEL_PAV_ENV
    threads: 1
    resources:
        mem_mb=4000,
        runtime=60
    shell:
        r"""
        set -euo pipefail
        python {params.script:q} \
          --sample {wildcards.sample:q} \
          --reference {input.reference:q} \
          --hap1 {input.hap1:q} \
          --hap2 {input.hap2:q} \
          --extra-config-json {params.extra:q} \
          --config-output {output.config:q} \
          --assemblies-output {output.assemblies:q} > {log:q} 2>&1
        """


rule call_svs_with_pav:
    input:
        config=f"{NOVEL_OUTDIR}/reference_calls/{{reference}}/pav/{{sample}}/config.json",
        assemblies=f"{NOVEL_OUTDIR}/reference_calls/{{reference}}/pav/{{sample}}/assemblies.tsv",
        snakefile=NOVEL_PAV_SNAKEFILE
    output:
        vcf=f"{NOVEL_OUTDIR}/reference_calls/{{reference}}/pav/{{sample}}/pav_{{sample}}.vcf.gz"
    log:
        f"{NOVEL_OUTDIR}/logs/reference_calls/{{reference}}.{{sample}}.pav.log"
    benchmark:
        f"{NOVEL_OUTDIR}/benchmarks/reference_calls/{{reference}}.{{sample}}.pav.tsv"
    conda:
        NOVEL_PAV_ENV
    threads: int(NOVEL_REFERENCE_CALLING.get("threads", 16))
    resources:
        mem_mb=int(NOVEL_REFERENCE_CALLING.get("caller_mem_mb", 64000)),
        runtime=int(
            NOVEL_REFERENCE_CALLING.get(
                "caller_runtime", NOVEL_REFERENCE_CALLING.get("caller_runtime_min", 720)
            )
        )
    shell:
        r"""
        set -euo pipefail
        analysis_dir=$(dirname {input.config:q})
        pav_snakefile=$(readlink -f {input.snakefile:q})
        cd "$analysis_dir"
        OPENBLAS_NUM_THREADS=1 snakemake \
          --snakefile "$pav_snakefile" \
          --cores {threads} \
          --rerun-incomplete \
          "pav_{wildcards.sample}.vcf.gz" > {log:q} 2>&1
        test -s "pav_{wildcards.sample}.vcf.gz"
        """


rule write_novel_sv_reference_call_manifest:
    input:
        NOVEL_REFERENCE_OUTPUTS
    output:
        NOVEL_CALL_MANIFEST
    run:
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        fields = [
            "caller",
            "coordinate_system",
            "sample_id",
            "assembly_id",
            "haplotype",
            "path",
            "callable_bed",
            "exclude_contig_regex",
        ]
        with open(output[0], "w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fields,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(NOVEL_REFERENCE_CALL_RECORDS)


rule merge_assembly_only_novel_sv_catalog:
    input:
        graph_candidates=NOVEL_GRAPH_CANDIDATES,
        call_manifest=NOVEL_CALL_MANIFEST,
        assembly_manifest=NOVEL_DISCOVERY_MANIFEST
    output:
        catalog=NOVEL_CATALOG,
        evidence=NOVEL_EVIDENCE
    params:
        script=NOVEL_MERGE_SCRIPT,
        min_sv=int(NOVEL_SCREEN.get("min_sv_size", 50)),
        distance=int(NOVEL_CLUSTERING.get("breakpoint_distance", 500)),
        similarity=float(NOVEL_CLUSTERING.get("length_similarity", 0.70)),
        exclude_contig_regex=NOVEL_EXCLUDE_CONTIG_REGEX
    log:
        f"{NOVEL_OUTDIR}/logs/merge_assembly_only_catalog.log"
    conda:
        NOVEL_CATALOG_ENV
    threads: 1
    resources:
        mem_mb=32000,
        runtime=360
    shell:
        r"""
        set -euo pipefail
        python {params.script:q} \
          --graph-candidates {input.graph_candidates:q} \
          --call-manifest {input.call_manifest:q} \
          --assembly-manifest {input.assembly_manifest:q} \
          --catalog-output {output.catalog:q} \
          --evidence-output {output.evidence:q} \
          --min-sv-size {params.min_sv} \
          --breakpoint-distance {params.distance} \
          --length-similarity {params.similarity} \
          --exclude-contig-regex {params.exclude_contig_regex:q} > {log:q} 2>&1
        test -s {output.catalog:q}
        """


rule rank_recommended_hifi_samples:
    input:
        manifest=NOVEL_DISCOVERY_MANIFEST,
        screen=NOVEL_GRAPH_SUMMARY,
        catalog=NOVEL_CATALOG,
        evidence=NOVEL_EVIDENCE
    output:
        NOVEL_RECOMMENDED_HIFI
    params:
        script=NOVEL_HIFI_SCRIPT,
        callable_threshold=float(NOVEL_HIFI.get("callable_threshold", 0.85)),
        controls=int(NOVEL_HIFI.get("control_count", 20))
    log:
        f"{NOVEL_OUTDIR}/logs/rank_recommended_hifi.log"
    conda:
        NOVEL_CATALOG_ENV
    threads: 1
    resources:
        mem_mb=8000,
        runtime=120
    shell:
        r"""
        set -euo pipefail
        python {params.script:q} \
          --assembly-manifest {input.manifest:q} \
          --screen-summary {input.screen:q} \
          --catalog {input.catalog:q} \
          --evidence {input.evidence:q} \
          --callable-threshold {params.callable_threshold} \
          --control-count {params.controls} \
          --output {output:q} > {log:q} 2>&1
        test -s {output:q}
        """
