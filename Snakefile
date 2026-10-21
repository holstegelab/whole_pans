MAIN_CONFIG = workflow.source_path("config/config.yaml")
configfile: MAIN_CONFIG


rule all:
    input:
        f"{config['results']['decontamination']}/summary/contamination_summary.tsv",
        f"{config['results']['decontamination']}/summary/contig_actions.tsv",
        f"{config['results']['decontamination']}/summary/review_candidates.tsv",
        f"{config['results']['decontamination']}/summary/graph_cleaned_assemblies.txt",
        f"{config['results']['post_decontamination_qc']}/summary/assembly_qc.tsv",
        f"{config['results']['post_decontamination_qc']}/summary/sample_qc.tsv",
        f"{config['results']['post_decontamination_qc']}/summary/graph_included_assemblies.txt",
        f"{config['results']['post_decontamination_qc']}/summary/graph_excluded_assemblies.tsv",


include: "workflow/rules/post_decontamination_QC.smk"
include: "workflow/rules/pangenome.smk"
include: "workflow/rules/hprc_plus_cohort.smk"
include: "workflow/rules/pangenome_qc.smk"
include: "workflow/rules/review_followup.smk"
include: "workflow/rules/minigraph_cactus.smk"
include: "workflow/rules/minigraph_cactus_reuse.smk"
