PANGENOME_QC_CONFIG = workflow.source_path("../../config/config.yaml")
configfile: PANGENOME_QC_CONFIG


import os
from pathlib import Path
from pathlib import PurePosixPath


PANGENOME_QC = config.get("pangenome", {})
PANGENOME_QC_PROJECT_ROOT = Path(workflow.current_basedir).resolve().parents[1]
workdir: str(PANGENOME_QC_PROJECT_ROOT)

PANGENOME_QC_PATHS = {
    "results": PANGENOME_QC.get("results", "results/sv_pangenome"),
    "qc_results": PANGENOME_QC.get(
        "qc_results",
        f"{PANGENOME_QC.get('results', 'results/sv_pangenome')}/qc_analysis",
    ),
    "vg": PANGENOME_QC.get("vg", "../vg"),
}
for key, value in PANGENOME_QC_PATHS.items():
    if os.path.isabs(str(value)):
        raise WorkflowError(
            f"pangenome.{key} must be relative to the Snakefile, got: {value}"
        )

PANGENOME_QC_OUTDIR = str(PurePosixPath(str(PANGENOME_QC_PATHS["results"])))
PANGENOME_QC_RESULT_DIR = str(PurePosixPath(str(PANGENOME_QC_PATHS["qc_results"])))
PANGENOME_QC_VG = str(PurePosixPath(str(PANGENOME_QC_PATHS["vg"])))
PANGENOME_QC_GFA = f"{PANGENOME_QC_OUTDIR}/graphs/sv_pangenome.minigraph.gfa"
PANGENOME_QC_GRAPH_SUMMARY = f"{PANGENOME_QC_OUTDIR}/metadata/sv_pangenome.graph_summary.tsv"
PANGENOME_QC_ORDERED = f"{PANGENOME_QC_OUTDIR}/metadata/sv_pangenome.ordered_assemblies.tsv"
PANGENOME_QC_MASH = f"{PANGENOME_QC_OUTDIR}/metadata/chm13.mash_distances.tsv"
PANGENOME_QC_TOOL_VERSIONS = f"{PANGENOME_QC_OUTDIR}/metadata/tool_versions.tsv"
PANGENOME_QC_BUILD_LOG = f"{PANGENOME_QC_OUTDIR}/logs/build_sv_pangenome_graph.log"
PANGENOME_QC_BUILD_BENCHMARK = (
    f"{PANGENOME_QC_OUTDIR}/benchmarks/build_sv_pangenome_graph.tsv"
)
PANGENOME_QC_POST_SUMMARY = (
    "../whole_pangenome/assembly_qc_decontaminated/results/summary"
)
PANGENOME_QC_DECONTAM_SUMMARY = (
    "../whole_pangenome/assembly_decontamination/results/summary"
)
PANGENOME_QC_SCRIPT = str(workflow.source_path("../scripts/pangenome_qc_analysis.py"))
PANGENOME_QC_ENV = str(workflow.source_path("../envs/pangenome.yaml"))


rule pangenome_qc:
    input:
        f"{PANGENOME_QC_RESULT_DIR}/data/seg_records.parquet",
        f"{PANGENOME_QC_RESULT_DIR}/data/rank_tally.tsv",
        f"{PANGENOME_QC_RESULT_DIR}/tables/integrity_checks.csv",
        f"{PANGENOME_QC_RESULT_DIR}/tables/graph_overview_stats.csv",
        f"{PANGENOME_QC_RESULT_DIR}/tables/per_source_contribution.csv",
        f"{PANGENOME_QC_RESULT_DIR}/tables/per_chromosome_nonref.csv",
        f"{PANGENOME_QC_RESULT_DIR}/tables/mash_outliers.csv",
        f"{PANGENOME_QC_RESULT_DIR}/figures/segment_length_distribution.png",
        f"{PANGENOME_QC_RESULT_DIR}/figures/pangenome_growth_curve.png",
        f"{PANGENOME_QC_RESULT_DIR}/figures/per_chromosome_nonref.png",
        f"{PANGENOME_QC_RESULT_DIR}/figures/mash_distance_diversity.png",
        f"{PANGENOME_QC_RESULT_DIR}/report/pangenome_qc_report.md",
        f"{PANGENOME_QC_RESULT_DIR}/metadata/qc_tool_versions.tsv",
        f"{PANGENOME_QC_RESULT_DIR}/metadata/run_manifest.tsv",


rule analyze_sv_pangenome_graph:
    input:
        gfa=PANGENOME_QC_GFA,
        graph_summary=PANGENOME_QC_GRAPH_SUMMARY,
        ordered=PANGENOME_QC_ORDERED,
        mash=PANGENOME_QC_MASH,
        tool_versions=PANGENOME_QC_TOOL_VERSIONS,
        build_log=PANGENOME_QC_BUILD_LOG,
        build_benchmark=PANGENOME_QC_BUILD_BENCHMARK,
        post_qc_assembly=f"{PANGENOME_QC_POST_SUMMARY}/assembly_qc.tsv",
        post_qc_sample=f"{PANGENOME_QC_POST_SUMMARY}/sample_qc.tsv",
        post_qc_included=f"{PANGENOME_QC_POST_SUMMARY}/graph_included_assemblies.txt",
        contamination=f"{PANGENOME_QC_DECONTAM_SUMMARY}/contamination_summary.tsv"
    output:
        seg_records=f"{PANGENOME_QC_RESULT_DIR}/data/seg_records.parquet",
        rank_tally=f"{PANGENOME_QC_RESULT_DIR}/data/rank_tally.tsv",
        integrity=f"{PANGENOME_QC_RESULT_DIR}/tables/integrity_checks.csv",
        overview=f"{PANGENOME_QC_RESULT_DIR}/tables/graph_overview_stats.csv",
        per_source=f"{PANGENOME_QC_RESULT_DIR}/tables/per_source_contribution.csv",
        per_chrom=f"{PANGENOME_QC_RESULT_DIR}/tables/per_chromosome_nonref.csv",
        mash_outliers=f"{PANGENOME_QC_RESULT_DIR}/tables/mash_outliers.csv",
        length_fig=f"{PANGENOME_QC_RESULT_DIR}/figures/segment_length_distribution.png",
        growth_fig=f"{PANGENOME_QC_RESULT_DIR}/figures/pangenome_growth_curve.png",
        chrom_fig=f"{PANGENOME_QC_RESULT_DIR}/figures/per_chromosome_nonref.png",
        mash_fig=f"{PANGENOME_QC_RESULT_DIR}/figures/mash_distance_diversity.png",
        report=f"{PANGENOME_QC_RESULT_DIR}/report/pangenome_qc_report.md",
        versions=f"{PANGENOME_QC_RESULT_DIR}/metadata/qc_tool_versions.tsv",
        manifest=f"{PANGENOME_QC_RESULT_DIR}/metadata/run_manifest.tsv"
    params:
        script=PANGENOME_QC_SCRIPT,
        output_dir=PANGENOME_QC_RESULT_DIR,
        vg=PANGENOME_QC_VG,
        ln_sample_limit=int(PANGENOME_QC.get("ln_sample_limit", 100000)),
        mad_k=float(PANGENOME_QC.get("mash_mad_k", 3.5)),
        low_matching_fraction=float(
            PANGENOME_QC.get("low_matching_hash_fraction", 0.90)
        )
    log:
        f"{PANGENOME_QC_RESULT_DIR}/logs/analyze_sv_pangenome_graph.log"
    benchmark:
        f"{PANGENOME_QC_RESULT_DIR}/metadata/analyze_sv_pangenome_graph.benchmark.tsv"
    conda:
        PANGENOME_QC_ENV
    threads: 1
    resources:
        mem_mb=int(PANGENOME_QC.get("analysis_mem_mb", 16000)),
        runtime_min=int(PANGENOME_QC.get("analysis_runtime_min", 240))
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})"
        python {params.script:q} \
          --gfa {input.gfa:q} \
          --graph-summary {input.graph_summary:q} \
          --ordered-assemblies {input.ordered:q} \
          --mash-distances {input.mash:q} \
          --tool-versions {input.tool_versions:q} \
          --build-log {input.build_log:q} \
          --build-benchmark {input.build_benchmark:q} \
          --post-qc-assembly {input.post_qc_assembly:q} \
          --post-qc-sample {input.post_qc_sample:q} \
          --post-qc-included {input.post_qc_included:q} \
          --contamination-summary {input.contamination:q} \
          --vg {params.vg:q} \
          --output-dir {params.output_dir:q} \
          --ln-sample-limit {params.ln_sample_limit} \
          --mad-k {params.mad_k} \
          --low-matching-fraction {params.low_matching_fraction} > {log:q} 2>&1
        test -s {output.report:q}
        """
