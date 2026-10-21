"""Local SV alleles -> resource comparison -> measurement -> gene/phenotype evidence.

Run this standalone from whole_pans. The default performs inventory only.
"""

configfile: workflow.source_path("../../config/sv_gene_phenotype.yaml")

import csv
import shlex
from pathlib import Path

SGP_ROOT = Path(workflow.current_basedir).parents[1]
workdir: str(SGP_ROOT)
shell.executable("/bin/bash")

SGP = config["sv_gene_phenotype"]
SGP_OUT = SGP["results"].rstrip("/")
SGP_TABLE_ENV = str(SGP_ROOT / "workflow/envs/svgp_tables.yaml")
SGP_PLACE_ENV = str(SGP_ROOT / "workflow/envs/svgp_placement.yaml")
SGP_CATALOGUE = str(SGP_ROOT / "workflow/scripts/svgp_catalogue.py")
SGP_PROVENANCE = str(SGP_ROOT / "workflow/scripts/review_followup.py")
SGP_PREP = SGP["catalogue"]
SGP_COMPARE = SGP["representation"]["resources"]
SGP_COMPUTE = SGP["compute"]
SGP_TABLES = ("ranked_candidates.tsv", "ranked_genic_candidates.tsv",
              "selected_candidates.tsv", "selected_genic_candidates.tsv")
SGP_PHENO = SGP.get("phenotypes", {})
SGP_SAMPLE_INPUT = (SGP_OUT + "/metadata/samples_with_phenotypes.tsv"
                    if SGP_PHENO.get("enabled", False) else SGP["inventory"]["samples"])


rule inventory:
    input:
        expand(SGP_OUT + "/metadata/{name}.tsv", name=[
            "samples", "reads", "resources", "sample_flow", "phenotype_by_batch_depth",
            "donor_overlap", "ont_pairing", "missing_inputs", "exclusions"]),
        SGP_OUT + "/metadata/reported_ont_runs.tsv",
        SGP_OUT + "/existing_calls/source_files.tsv",
        SGP_OUT + "/existing_calls/source_summary.tsv"
    default_target: True


rule phenotypes:
    input:
        SGP_OUT + "/metadata/phenotypes.tsv",
        SGP_OUT + "/metadata/phenotype_join_report.tsv",
        SGP_OUT + "/metadata/phenotype_import_summary.tsv"


rule svgp_import_phenotypes:
    input:
        script=str(SGP_ROOT / "workflow/scripts/svgp_phenotypes.py"),
        common=str(SGP_ROOT / "workflow/scripts/svgp_inventory.py"),
        workbook=SGP_PHENO.get("workbook", ""),
        samples=SGP["inventory"]["samples"],
        crosswalk=SGP_PHENO.get("crosswalk", "config/svgp_phenotype_crosswalk.example.tsv")
    output:
        phenotypes=SGP_OUT + "/metadata/phenotypes.tsv",
        joins=SGP_OUT + "/metadata/phenotype_join_report.tsv",
        summary=SGP_OUT + "/metadata/phenotype_import_summary.tsv",
        samples=SGP_OUT + "/metadata/samples_with_phenotypes.tsv"
    params:
        sheet=SGP_PHENO.get("sheet", "Blad1"),
        id_column=SGP_PHENO.get("id_column", "Sample_ID_WGS_given"),
        phenotype_column=SGP_PHENO.get("phenotype_column", "Diagnosis_svPPA"),
        out=SGP_OUT + "/metadata"
    threads: 1
    resources:
        mem_mb=2000,
        runtime_min=15
    log: SGP_OUT + "/logs/phenotypes.log"
    conda: str(SGP_ROOT / "workflow/envs/svgp_phenotypes.yaml")
    shell:
        "python {input.script:q} --workbook {input.workbook:q} --sheet {params.sheet:q} "
        "--samples {input.samples:q} --crosswalk {input.crosswalk:q} "
        "--id-column {params.id_column:q} --phenotype-column {params.phenotype_column:q} "
        "--output {params.out:q} > {log:q} 2>&1"


rule catalogue:
    input:
        SGP_OUT + "/catalogue/alleles.tsv",
        SGP_OUT + "/catalogue/coordinate_crosswalk.tsv",
        SGP_OUT + "/catalogue/rejections.tsv"


rule representation:
    input:
        SGP_OUT + "/representation/resource_matrix.tsv",
        SGP_OUT + "/representation/evidence.tsv",
        SGP_OUT + "/representation/flow.tsv"


rule panel:
    input:
        SGP_OUT + "/panel/panel.tsv"


checkpoint svgp_inventory:
    input:
        script=str(SGP_ROOT / "workflow/scripts/svgp_inventory.py"),
        samples=SGP_SAMPLE_INPUT,
        reads=SGP["inventory"]["reads"],
        resources=SGP["inventory"]["resources"],
        crosswalk=[SGP["inventory"]["crosswalk"]] if SGP["inventory"].get("crosswalk") else []
    output:
        expand(SGP_OUT + "/metadata/{name}.tsv", name=[
            "samples", "reads", "resources", "sample_flow", "phenotype_by_batch_depth",
            "donor_overlap", "ont_pairing", "missing_inputs", "exclusions"])
    params:
        out=SGP_OUT,
        crosswalk=lambda wc, input: "--crosswalk " + shlex.quote(input.crosswalk[0]) if input.crosswalk else ""
    threads: 1
    resources:
        mem_mb=2000,
        runtime_min=30
    log: SGP_OUT + "/logs/inventory.log"
    conda: SGP_TABLE_ENV
    shell:
        "python {input.script:q} --samples {input.samples:q} --reads {input.reads:q} "
        "--resources {input.resources:q} {params.crosswalk} --output {params.out:q} > {log:q} 2>&1"


# Reuse only the existing provenance helper, without including the old broad DAG.
rule svgp_correct_provenance:
    input:
        script=SGP_PROVENANCE,
        manifest=SGP_PREP["build_manifest"],
        graph=SGP_PREP["graph"],
        bubbles=SGP_PREP["bubbles"],
        tables=[SGP_PREP["examples"] + "/" + name for name in SGP_TABLES]
    output:
        candidates=SGP_OUT + "/provenance/all_candidates.tsv",
        alleles=SGP_OUT + "/provenance/insertion_alleles.fa",
        order=SGP_OUT + "/provenance/build_order.tsv",
        tables=[SGP_OUT + "/provenance/" + name for name in SGP_TABLES]
    params:
        out=SGP_OUT + "/provenance",
        first_rank=SGP_PREP["first_local_source_rank"]
    threads: 1
    resources:
        mem_mb=SGP_COMPUTE["preparation_mem_mb"],
        runtime_min=SGP_COMPUTE["preparation_runtime_min"]
    log: SGP_OUT + "/logs/provenance.log"
    benchmark: SGP_OUT + "/benchmarks/provenance.tsv"
    conda: SGP_TABLE_ENV
    shell:
        "python {input.script:q} provenance --manifest {input.manifest:q} "
        "--graph {input.graph:q} --bubbles {input.bubbles:q} --tables {input.tables:q} "
        "--first-rank {params.first_rank} --output {params.out:q} > {log:q} 2>&1"


rule svgp_catalogue:
    input:
        script=SGP_CATALOGUE,
        candidates=SGP_OUT + "/provenance/all_candidates.tsv",
        alleles=SGP_OUT + "/provenance/insertion_alleles.fa",
        reference=SGP_PREP["reference"],
        optional=[SGP_PREP[key] for key in ("source_review", "tmem_fasta", "tmem_targets", "controls", "control_fasta") if SGP_PREP.get(key)]
    output:
        alleles=SGP_OUT + "/catalogue/alleles.tsv",
        fasta=SGP_OUT + "/catalogue/alleles.fa",
        queries=SGP_OUT + "/catalogue/locus_queries.fa",
        manifest=SGP_OUT + "/catalogue/query_manifest.tsv",
        crosswalk=SGP_OUT + "/catalogue/coordinate_crosswalk.tsv",
        rejections=SGP_OUT + "/catalogue/rejections.tsv",
        flow=SGP_OUT + "/catalogue/candidate_flow.tsv",
        totals=SGP_OUT + "/catalogue/descriptive_totals.tsv"
    params:
        graph_version=SGP_PREP["graph_version"],
        graph_sha256=SGP_PREP["graph_sha256"],
        flank=SGP_PREP["flank_bp"],
        out=SGP_OUT + "/catalogue",
        optional=" ".join("--" + key.replace("_", "-") + " " + shlex.quote(SGP_PREP[key])
                          for key in ("source_review", "tmem_fasta", "tmem_targets", "controls", "control_fasta") if SGP_PREP.get(key))
    threads: 1
    resources:
        mem_mb=SGP_COMPUTE["preparation_mem_mb"],
        runtime_min=SGP_COMPUTE["preparation_runtime_min"]
    log: SGP_OUT + "/logs/catalogue.log"
    benchmark: SGP_OUT + "/benchmarks/catalogue.tsv"
    conda: SGP_TABLE_ENV
    shell:
        "python {input.script:q} catalogue --candidates {input.candidates:q} --alleles {input.alleles:q} "
        "--reference {input.reference:q} --graph-version {params.graph_version:q} "
        "--graph-sha256 {params.graph_sha256:q} --flank {params.flank} {params.optional} "
        "--output {params.out:q} > {log:q} 2>&1"


rule svgp_linear_placement:
    input:
        reference=lambda wc: SGP_COMPARE[wc.resource]["path"],
        queries=SGP_OUT + "/catalogue/locus_queries.fa"
    output:
        paf=SGP_OUT + "/representation/alignments/{resource}.paf",
        version=SGP_OUT + "/representation/alignments/{resource}.version.txt"
    wildcard_constraints:
        resource="GRCh38|CHM13"
    threads: SGP_COMPUTE["placement_threads"]
    resources:
        mem_mb=SGP_COMPUTE["placement_mem_mb"],
        runtime_min=SGP_COMPUTE["placement_runtime_min"]
    log: SGP_OUT + "/logs/placement_{resource}.log"
    benchmark: SGP_OUT + "/benchmarks/placement_{resource}.tsv"
    conda: SGP_PLACE_ENV
    shell:
        "minimap2 --version > {output.version:q}; "
        "minimap2 -x asm5 -c --cs=long --secondary=yes -t {threads} "
        "{input.reference:q} {input.queries:q} > {output.paf:q} 2> {log:q}"


rule svgp_graph_placement:
    input:
        graph=SGP_COMPARE["HPRC"]["path"],
        queries=SGP_OUT + "/catalogue/locus_queries.fa"
    output:
        gaf=SGP_OUT + "/representation/alignments/HPRC.gaf",
        version=SGP_OUT + "/representation/alignments/HPRC.version.txt"
    threads: SGP_COMPUTE["placement_threads"]
    resources:
        mem_mb=SGP_COMPUTE["placement_mem_mb"],
        runtime_min=SGP_COMPUTE["placement_runtime_min"]
    log: SGP_OUT + "/logs/placement_HPRC.log"
    benchmark: SGP_OUT + "/benchmarks/placement_HPRC.tsv"
    conda: SGP_PLACE_ENV
    shell:
        "minigraph --version > {output.version:q}; "
        "minigraph -c -xasm --vc --secondary=yes -t {threads} "
        "{input.graph:q} {input.queries:q} > {output.gaf:q} 2> {log:q}"


rule svgp_search_manifest:
    input:
        SGP_OUT + "/representation/alignments/HPRC.gaf",
        SGP_OUT + "/representation/alignments/GRCh38.paf",
        SGP_OUT + "/representation/alignments/CHM13.paf"
    output: SGP_OUT + "/representation/search_manifest.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime_min=5
    log: SGP_OUT + "/logs/search_manifest.log"
    conda: SGP_TABLE_ENV
    run:
        with open(output[0], "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["resource", "alignment_path", "resource_version", "resource_sha256", "search_complete"])
            for resource, path in zip(("HPRC", "GRCh38", "CHM13"), input):
                writer.writerow([resource, path, SGP_COMPARE[resource]["version"], SGP_COMPARE[resource]["sha256"], "true"])
        Path(log[0]).write_text("Whole-resource mapper jobs completed; adequacy still requires anchored review.\n")


rule svgp_representation:
    input:
        script=SGP_CATALOGUE,
        catalogue=SGP_OUT + "/catalogue/alleles.tsv",
        queries=SGP_OUT + "/catalogue/query_manifest.tsv",
        searches=SGP_OUT + "/representation/search_manifest.tsv",
        review=SGP["representation"]["review"],
        candidate_flow=SGP_OUT + "/catalogue/candidate_flow.tsv"
    output:
        matrix=SGP_OUT + "/representation/resource_matrix.tsv",
        evidence=SGP_OUT + "/representation/evidence.tsv",
        counts=SGP_OUT + "/representation/resource_counts.tsv",
        flow=SGP_OUT + "/representation/flow.tsv"
    params: out=SGP_OUT + "/representation"
    threads: 1
    resources:
        mem_mb=8000,
        runtime_min=60
    log: SGP_OUT + "/logs/representation.log"
    conda: SGP_TABLE_ENV
    shell:
        "python {input.script:q} representation --catalogue {input.catalogue:q} --queries {input.queries:q} "
        "--search-manifest {input.searches:q} --review {input.review:q} --output {params.out:q} > {log:q} 2>&1"


rule svgp_panel:
    input:
        script=SGP_CATALOGUE,
        catalogue=SGP_OUT + "/catalogue/alleles.tsv",
        queries=SGP_OUT + "/catalogue/query_manifest.tsv",
        fasta=SGP_OUT + "/catalogue/locus_queries.fa",
        matrix=SGP_OUT + "/representation/resource_matrix.tsv",
        selection=SGP["panel"]["selection"]
    output:
        panel=SGP_OUT + "/panel/panel.tsv",
        queries=SGP_OUT + "/panel/locus_queries.fa",
        manifest=SGP_OUT + "/panel/query_manifest.tsv",
        counts=SGP_OUT + "/panel/counts.tsv"
    params: out=SGP_OUT + "/panel"
    threads: 1
    resources:
        mem_mb=4000,
        runtime_min=30
    log: SGP_OUT + "/logs/panel.log"
    conda: SGP_TABLE_ENV
    shell:
        "python {input.script:q} panel --catalogue {input.catalogue:q} --queries {input.queries:q} "
        "--fasta {input.fasta:q} --resource-matrix {input.matrix:q} --selection {input.selection:q} "
        "--output {params.out:q} > {log:q} 2>&1"


include: "svgp_existing_calls.smk"
include: "svgp_measurement.smk"
include: "svgp_annotation.smk"
include: "svgp_association.smk"


rule pilot:
    input: SGP_PILOT_EVIDENCE, SGP_PILOT_SUMMARY


rule ont_validation:
    input: SGP_ONT_EVIDENCE, SGP_ONT_SUMMARY


rule cohort:
    input:
        SGP_COHORT_EVIDENCE,
        SGP_COHORT_QC,
        SGP_OUT + "/genotypes/carrier_matrix.tsv.gz",
        SGP_OUT + "/genotypes/dosage_matrix.tsv.gz",
        SGP_OUT + "/genotypes/callability.tsv.gz"


rule association:
    input:
        SGP_OUT + "/association/effects.tsv",
        SGP_OUT + "/association/counts.tsv",
        SGP_OUT + "/association/status.tsv"


include: "svgp_reported_ont.smk"
include: "svgp_prediction.smk"
include: "svgp_report.smk"
