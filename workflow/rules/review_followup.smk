"""Immediate follow-up to the September 15 review; also usable standalone."""

configfile: workflow.source_path("../../config/config.yaml")
configfile: workflow.source_path("../../config/review_followup.yaml")

import json
import shlex
from pathlib import Path

RF_ROOT = Path(workflow.current_basedir).parents[1]
workdir: str(RF_ROOT)
shell.executable("/bin/bash")

RF = config["review_followup"]
RF_OUT = RF["results"]
RF_PANEL = RF["panel"]
RF_RES = RF["resources"]
RF_V2 = config["local_sv_graph_v2"]
RF_V2_OUT = RF_V2["results"]
RF_SCRIPT = str(workflow.source_path("../scripts/review_followup.py"))
RF_QC_SCRIPT = str(workflow.source_path("../scripts/qc_local_sv_graph_v2.py"))
RF_ENV = str(workflow.source_path("../envs/review_followup.yaml"))
RF_PYTHON_ENV = str(workflow.source_path("../envs/sv_catalog.yaml"))
RF_EXAMPLE_NAMES = ("ranked_candidates.tsv", "ranked_genic_candidates.tsv",
                    "selected_candidates.tsv", "selected_genic_candidates.tsv")
RF_REFERENCES = config["novel_sv_discovery"]["references"]
RF_QC = f"{RF_OUT}/v2_qc"
RF_DONORS = f"{RF_OUT}/hprc_provenance"
RF_ALLELES = f"{RF_OUT}/panel"
RF_FROZEN = f"{RF_OUT}/tmem_frozen"


rule review_followup:
    input:
        f"{RF_OUT}/report.md",
        f"{RF_OUT}/checksums.sha256"


# Supplies the existing corrected analysis_v2 targets and completes missing
# sample/graph evidence before freezing. It does not change sample selection.
include: "personalized_sv_comparison.smk"


rule review_copy_v2_graph_qc:
    input:
        integrity=f"{RF_V2_OUT}/qc/tables/graph_integrity.tsv",
        segments=f"{RF_V2_OUT}/qc/tables/delta_segments.tsv",
        links=f"{RF_V2_OUT}/qc/tables/delta_links.tsv",
        metrics=f"{RF_V2_OUT}/qc/tables/graph_qc_metrics.tsv",
        provenance=f"{RF_V2_OUT}/metadata/sv_pangenome.local.v2.provenance.tsv",
        inventory=RF_V2["frozen_graph_inventory"]
    output:
        integrity=f"{RF_QC}/tables/graph_integrity.tsv",
        segments=f"{RF_QC}/tables/delta_segments.tsv",
        links=f"{RF_QC}/tables/delta_links.tsv",
        metrics=f"{RF_QC}/tables/graph_qc_metrics.tsv",
        provenance=f"{RF_QC}/metadata/v2.provenance.tsv",
        inventory=f"{RF_QC}/metadata/v1.graph_inventory.tsv"
    log: f"{RF_OUT}/logs/copy_v2_graph_qc.log"
    conda: RF_PYTHON_ENV
    threads: 1
    resources:
        mem_mb=1000,
        runtime_min=30
    shell:
        r"""
        mkdir -p "$(dirname {output.integrity:q})" "$(dirname {output.provenance:q})" "$(dirname {log:q})"
        {{
          cp -- {input.integrity:q} {output.integrity:q}
          cp -- {input.segments:q} {output.segments:q}
          cp -- {input.links:q} {output.links:q}
          cp -- {input.metrics:q} {output.metrics:q}
          cp -- {input.provenance:q} {output.provenance:q}
          cp -- {input.inventory:q} {output.inventory:q}
        }} > {log:q} 2>&1
        """


rule review_corrected_v2_qc:
    input:
        script=RF_QC_SCRIPT,
        plan=f"{RF_V2_OUT}/rescue/plan.tsv",
        summary=f"{RF_V2_OUT}/rescue/rescue_summary.tsv",
        graph_inputs=f"{RF_V2_OUT}/rescue/clean_rescue_graph_inputs.tsv",
        catalog=RF_V2["catalog"],
        evidence=RF_V2["catalog_evidence"],
        assemblies=RF_V2["assembly_manifest"],
        integrity=f"{RF_QC}/tables/graph_integrity.tsv",
        frozen_bubbles=RF_V2["frozen_graph_bubbles"],
        rescued_bubbles=f"{RF_V2_OUT}/metadata/sv_pangenome.local.v2.bubbles.tsv"
    output:
        bubbles=f"{RF_QC}/tables/bubble_delta.tsv",
        rescues=f"{RF_QC}/tables/rescued_sv_graph_qc.tsv",
        cohort=f"{RF_QC}/tables/cohort_completeness.tsv",
        report=f"{RF_QC}/report/local_sv_graph_v2_qc_report.md"
    params:
        options=" ".join(
            "--" + name.replace("_", "-") + " " + shlex.quote(str(value))
            for name, value in {
                "min_sv_size": config["novel_sv_discovery"]["graph_screen"].get("min_sv_size", 50),
                "breakpoint_distance": RF_V2.get("qc", {}).get("breakpoint_distance", 500),
                "length_similarity": RF_V2.get("qc", {}).get("length_similarity", 0.70),
                "min_allele_reads_warn": RF_V2.get("qc", {}).get("min_allele_reads_warn", 5),
                "min_allele_fraction_warn": RF_V2.get("qc", {}).get("min_allele_fraction_warn", 0.10),
                "max_selected_reads_warn": RF_V2.get("qc", {}).get("max_selected_reads_warn", 100000),
                "max_local_contigs_warn": RF_V2.get("qc", {}).get("max_local_contigs_warn", 10000),
                "min_target_coverage_warn": RF_V2.get("qc", {}).get("min_target_coverage_warn", 0.80),
            }.items()
        )
    log: f"{RF_OUT}/logs/v2_qc.log"
    conda: RF_PYTHON_ENV
    threads: 1
    resources:
        mem_mb=RF_RES["qc_mem_mb"],
        runtime_min=RF_RES["qc_runtime_min"]
    shell:
        r"""
        mkdir -p "$(dirname {output.report:q})" "$(dirname {log:q})"
        python {input.script:q} report --plan {input.plan:q} \
          --rescue-summary {input.summary:q} --graph-inputs {input.graph_inputs:q} \
          --catalog {input.catalog:q} --catalog-evidence {input.evidence:q} \
          --assembly-manifest {input.assemblies:q} --integrity {input.integrity:q} \
          --frozen-bubbles {input.frozen_bubbles:q} --rescued-bubbles {input.rescued_bubbles:q} \
          --bubble-delta {output.bubbles:q} --rescue-qc {output.rescues:q} \
          --cohort-qc {output.cohort:q} --report {output.report:q} \
          {params.options} > {log:q} 2>&1
        """


rule review_hprc_donor_provenance:
    input:
        script=RF_SCRIPT,
        manifest=RF["hprc_build_manifest"],
        graph=PSC["source_graphs"]["hprc_plus_local"],
        bubbles=f"{PSC_OUT}/sv_catalogs/hprc_plus_local.bubbles.tsv",
        tables=[f"{RF['hprc_examples']}/{name}" for name in RF_EXAMPLE_NAMES]
    output:
        candidates=f"{RF_DONORS}/all_candidates.tsv",
        alleles=f"{RF_DONORS}/insertion_alleles.fa",
        order=f"{RF_DONORS}/build_order.tsv",
        tables=[f"{RF_DONORS}/{name}" for name in RF_EXAMPLE_NAMES]
    params:
        output=RF_DONORS,
        first_rank=int(RF["first_local_source_rank"])
    log: f"{RF_OUT}/logs/hprc_provenance.log"
    conda: RF_PYTHON_ENV
    threads: 1
    resources:
        mem_mb=RF_RES["preparation_mem_mb"],
        runtime_min=RF_RES["preparation_runtime_min"]
    shell:
        r"""
        mkdir -p {params.output:q} "$(dirname {log:q})"
        python {input.script:q} provenance --manifest {input.manifest:q} \
          --graph {input.graph:q} --bubbles {input.bubbles:q} --tables {input.tables:q} \
          --first-rank {params.first_rank} --output {params.output:q} > {log:q} 2>&1
        """


rule review_freeze_tmem:
    input:
        script=RF_SCRIPT,
        evidence=PSC_SV_EVIDENCE,
        inventory=PSC_INVENTORY,
        targets=PSC_GRAPH_TARGETS,
        report=PSC_REPORT,
        samples=PSC_SAMPLES_TSV,
        membership=PSC_MEMBERSHIP,
        metrics=PSC_MAPPING_METRICS,
        comparison=PSC_COMPARISON,
        catalog=PSC_SV_CATALOG,
        validation=PSC_TMEM_VALIDATION,
        queries=PSC_TMEM_FASTA,
        query_manifest=PSC_TMEM_TARGETS,
        stats_script=PSC_STATS,
        target_script=PSC_TMEM,
        rules=str(workflow.source_path("personalized_sv_comparison.smk")),
        followup_rules=str(workflow.source_path("review_followup.smk")),
        environment=PSC_ENV,
        config=str(workflow.source_path("../../config/config.yaml")),
        followup_config=str(workflow.source_path("../../config/review_followup.yaml"))
    output:
        snapshot=protected(directory(RF_FROZEN))
    params:
        samples=list(PSC_SAMPLES),
        graphs=list(PSC_GRAPHS),
        effective_config=lambda wc: json.dumps(config, sort_keys=True),
        artifacts=lambda wc, input: " ".join(
            "--artifact " + shlex.quote(f"{name}={input[field]}") for name, field in (
                ("samples.tsv", "samples"), ("sample_graph_membership.tsv", "membership"),
                ("mapping_metrics.tsv", "metrics"), ("paired_comparison.tsv", "comparison"),
                ("sv_catalogs.tsv", "catalog"), ("tmem_pilot_validation.tsv", "validation"),
                ("allele_queries.fa", "queries"), ("tmem_targets.tsv", "query_manifest"),
                ("personalized_sv_statistics.py", "stats_script"), ("tmem_whole_graph_target.py", "target_script"),
                ("personalized_sv_comparison.smk", "rules"), ("review_followup.smk", "followup_rules"),
                ("personalized_sv_comparison.yaml", "environment"), ("config.yaml", "config"),
                ("review_followup.yaml", "followup_config"), ("review_followup.py", "script"),
            )
        )
    log: f"{RF_OUT}/logs/freeze_tmem.log"
    conda: RF_PYTHON_ENV
    threads: 1
    resources:
        mem_mb=2000,
        runtime_min=30
    shell:
        r"""
        mkdir -p "$(dirname {log:q})"
        python {input.script:q} freeze --evidence {input.evidence:q} \
          --inventory {input.inventory:q} --targets {input.targets:q} --report {input.report:q} \
          --samples {params.samples:q} --graphs {params.graphs:q} \
          --effective-config {params.effective_config:q} {params.artifacts} \
          --output {output.snapshot:q} > {log:q} 2>&1
        """


rule review_sequence_resolved_panel:
    input:
        script=RF_SCRIPT,
        candidates=f"{RF_DONORS}/all_candidates.tsv",
        alleles=f"{RF_DONORS}/insertion_alleles.fa",
        reference=RF_REFERENCES["CHM13"],
        tmem_fasta=PSC_TMEM_FASTA,
        tmem_targets=PSC_TMEM_TARGETS,
        selection=[RF_PANEL["selection"]] if RF_PANEL.get("selection") else [],
        repeats=RF_PANEL.get("repeat_beds", []),
        segdups=RF_PANEL.get("segdup_beds", [])
    output:
        panel=f"{RF_ALLELES}/panel.tsv",
        queries=f"{RF_ALLELES}/query_manifest.tsv",
        chm13=f"{RF_ALLELES}/CHM13.queries.fa",
        hg38=f"{RF_ALLELES}/hg38.queries.fa",
        worksheet=f"{RF_ALLELES}/validation_worksheet.tsv"
    params:
        output=RF_ALLELES,
        max_loci=int(RF_PANEL["max_loci"]),
        genic=int(RF_PANEL["genic_loci"]),
        repeat_loci=int(RF_PANEL["repeat_review_loci"]),
        flank=int(RF_PANEL["flank_bp"]),
        optional=lambda wc, input: " ".join(
            flag + " " + shlex.quote(str(path))
            for flag, paths in (("--selection", input.selection), ("--repeat-bed", input.repeats),
                                ("--segdup-bed", input.segdups)) for path in paths
        )
    log: f"{RF_OUT}/logs/prepare_panel.log"
    conda: RF_PYTHON_ENV
    threads: 1
    resources:
        mem_mb=RF_RES["preparation_mem_mb"],
        runtime_min=RF_RES["preparation_runtime_min"]
    shell:
        r"""
        mkdir -p {params.output:q} "$(dirname {log:q})"
        python {input.script:q} panel --candidates {input.candidates:q} --alleles {input.alleles:q} \
          --reference {input.reference:q} --tmem-fasta {input.tmem_fasta:q} \
          --tmem-targets {input.tmem_targets:q} --max-loci {params.max_loci} \
          --genic-loci {params.genic} --repeat-loci {params.repeat_loci} --flank {params.flank} \
          {params.optional} --output {params.output:q} > {log:q} 2>&1
        """


rule review_allele_genome_placement:
    input:
        queries=f"{RF_ALLELES}/{{reference}}.queries.fa",
        reference=lambda wc: RF_REFERENCES[wc.reference]
    output:
        paf=f"{RF_ALLELES}/placement/{{reference}}.paf",
        version=f"{RF_ALLELES}/placement/{{reference}}.minimap2.version.txt"
    log: f"{RF_OUT}/logs/placement/{{reference}}.log"
    benchmark: f"{RF_OUT}/benchmarks/placement/{{reference}}.tsv"
    conda: RF_ENV
    threads: RF_RES["placement_threads"]
    resources:
        mem_mb=RF_RES["placement_mem_mb"],
        runtime_min=RF_RES["placement_runtime_min"]
    wildcard_constraints:
        reference="CHM13|hg38"
    shell:
        r"""
        mkdir -p "$(dirname {output.paf:q})" "$(dirname {log:q})"
        minimap2 --version > {output.version:q}
        minimap2 -x asm5 -c --cs=long --secondary=yes -N 20 -I 8G -t {threads} \
          {input.reference:q} {input.queries:q} > {output.paf:q} 2> {log:q}
        """


rule review_followup_report:
    input:
        script=RF_SCRIPT,
        panel=f"{RF_ALLELES}/panel.tsv",
        queries=f"{RF_ALLELES}/query_manifest.tsv",
        qc=f"{RF_QC}/report/local_sv_graph_v2_qc_report.md",
        frozen=RF_FROZEN,
        pafs=expand(f"{RF_ALLELES}/placement/{{reference}}.paf", reference=["CHM13", "hg38"])
    output:
        placement=f"{RF_ALLELES}/placement_review.tsv",
        report=f"{RF_OUT}/report.md"
    params:
        flank=int(RF_PANEL["flank_bp"]),
        pafs=" ".join("--paf " + shlex.quote(f"{ref}={RF_ALLELES}/placement/{ref}.paf")
                      for ref in ("CHM13", "hg38"))
    log: f"{RF_OUT}/logs/report.log"
    conda: RF_PYTHON_ENV
    threads: 1
    resources:
        mem_mb=2000,
        runtime_min=30
    shell:
        r"""
        mkdir -p "$(dirname {log:q})"
        python {input.script:q} report --panel {input.panel:q} --queries {input.queries:q} \
          {params.pafs} --flank {params.flank} --placement {output.placement:q} \
          --output {output.report:q} > {log:q} 2>&1
        """


rule review_followup_checksums:
    input:
        rules.review_followup_report.output,
        rules.review_copy_v2_graph_qc.output,
        rules.review_corrected_v2_qc.output,
        rules.review_hprc_donor_provenance.output,
        rules.review_sequence_resolved_panel.output,
        expand(f"{RF_ALLELES}/placement/{{reference}}.paf", reference=["CHM13", "hg38"]),
        expand(f"{RF_ALLELES}/placement/{{reference}}.minimap2.version.txt", reference=["CHM13", "hg38"]),
        RF["review"], RF_SCRIPT, RF_QC_SCRIPT,
        RF["hprc_build_manifest"],
        [f"{RF['hprc_examples']}/{name}" for name in RF_EXAMPLE_NAMES],
        str(workflow.source_path("review_followup.smk")),
        str(workflow.source_path("../../config/review_followup.yaml")),
        RF_ENV
    output: f"{RF_OUT}/checksums.sha256"
    log: f"{RF_OUT}/logs/checksums.log"
    conda: RF_PYTHON_ENV
    threads: 1
    resources:
        mem_mb=1000,
        runtime_min=30
    shell:
        r"""
        mkdir -p "$(dirname {log:q})"
        sha256sum {input:q} > {output:q} 2> {log:q}
        """
