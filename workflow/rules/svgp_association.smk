"""One frozen binary contrast; a small, explicit family of carrier tests."""

SGP_ASSOC = SGP.get("association", {})


def svgp_association_reviews(wildcards):
    paths = []
    for table, fields in ((SGP_ASSOC.get("design", "config/svgp_association_design.example.tsv"), ("review_evidence",)),
                          (SGP_ASSOC.get("tests", "config/svgp_association_tests.example.tsv"), ("selection_evidence", "snp_evidence"))):
        with open(table, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                paths.extend(row[field] for field in fields
                             if row.get(field) not in (None, "", ".", "NA", "unknown"))
    return sorted(set(paths))


rule svgp_association:
    input:
        samples=SGP_OUT + "/metadata/samples.tsv",
        evidence=SGP_OUT + "/genotypes/evidence.tsv.gz",
        panel=SGP_OUT + "/panel/panel.tsv",
        lock=SGP["measurement"]["measurement_lock"],
        release=SGP_OUT + "/gates/cohort_release.tsv",
        design=SGP_ASSOC.get("design", "config/svgp_association_design.example.tsv"),
        tests=SGP_ASSOC.get("tests", "config/svgp_association_tests.example.tsv"),
        reviews=svgp_association_reviews
    output:
        design=SGP_OUT + "/association/design.md",
        tests=SGP_OUT + "/association/test_manifest.tsv",
        counts=SGP_OUT + "/association/counts.tsv",
        effects=SGP_OUT + "/association/effects.tsv",
        missingness=SGP_OUT + "/association/missingness.tsv",
        power=SGP_OUT + "/association/power_grid.tsv",
        plot=SGP_OUT + "/association/effects.svg",
        status=SGP_OUT + "/association/status.tsv",
        versions=SGP_OUT + "/association/software_versions.txt"
    params:
        graph_id=SGP["measurement"].get("main_graph", "hprc_plus_local")
    threads: 1
    resources:
        mem_mb=SGP_ASSOC.get("mem_mb", 8000),
        runtime_min=SGP_ASSOC.get("runtime_min", 120)
    log: SGP_OUT + "/logs/association.log"
    conda: str(SGP_ROOT / "workflow/envs/svgp_association.yaml")
    script: "../scripts/svgp_association.R"
