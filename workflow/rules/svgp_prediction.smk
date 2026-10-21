"""Optional, prespecified nested prediction comparison after the evidence chain."""

SGP_PRED = SGP.get("prediction", {})


def svgp_prediction_reviews(wildcards):
    paths = []
    for table in (SGP_PRED["design"], SGP_PRED["features"]):
        with open(table, newline="") as handle:
            paths.extend(row["review_evidence"] for row in csv.DictReader(handle, delimiter="\t")
                         if row.get("review_evidence") not in (None, "", ".", "NA", "unknown"))
    return sorted(set(paths))


rule prediction:
    input: SGP_OUT + "/prediction/status.tsv"


if SGP_PRED.get("enabled", False):
    rule svgp_prediction:
        input:
            samples=SGP_OUT + "/metadata/samples.tsv",
            evidence=SGP_OUT + "/genotypes/evidence.tsv.gz",
            panel=SGP_OUT + "/panel/panel.tsv",
            lock=SGP["measurement"]["measurement_lock"],
            release=SGP_OUT + "/gates/cohort_release.tsv",
            pilot=SGP_OUT + "/pilot/summary.tsv",
            ont=SGP_OUT + "/validation/ont_by_sample_locus.tsv",
            annotation=SGP_OUT + "/annotation/allele_gene_context.tsv",
            association=SGP_OUT + "/association/status.tsv",
            association_design=SGP.get("association", {}).get("design", "config/svgp_association_design.example.tsv"),
            design=SGP_PRED["design"],
            features=SGP_PRED["features"],
            reviews=svgp_prediction_reviews
        output:
            status=SGP_OUT + "/prediction/status.tsv",
            splits=SGP_OUT + "/prediction/split_manifest.tsv",
            oof=SGP_OUT + "/prediction/out_of_fold.tsv",
            repeats=SGP_OUT + "/prediction/stability_predictions.tsv",
            metrics=SGP_OUT + "/prediction/metrics.tsv",
            uncertainty=SGP_OUT + "/prediction/conditional_uncertainty.tsv",
            features=SGP_OUT + "/prediction/feature_manifest.tsv",
            report=SGP_OUT + "/prediction/report.md",
            versions=SGP_OUT + "/prediction/software_versions.txt"
        params:
            settings=SGP_PRED,
            graph_id=SGP["measurement"]["main_graph"]
        threads: SGP_PRED.get("threads", 4)
        resources:
            mem_mb=SGP_PRED.get("mem_mb", 8000),
            runtime_min=SGP_PRED.get("runtime_min", 120)
        log: SGP_OUT + "/logs/prediction.log"
        benchmark: SGP_OUT + "/benchmarks/prediction.tsv"
        conda: str(SGP_ROOT / "workflow/envs/svgp_prediction.yaml")
        script: "../scripts/svgp_prediction.py"
else:
    rule svgp_prediction_skipped:
        output: SGP_OUT + "/prediction/status.tsv"
        threads: 1
        resources:
            mem_mb=100,
            runtime_min=1
        log: SGP_OUT + "/logs/prediction.log"
        run:
            Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
            Path(log[0]).parent.mkdir(parents=True, exist_ok=True)
            message = "stage\tstatus\treason\nprediction\tSKIPPED\tDisabled by default; requires completed earlier stages and a frozen prediction design\n"
            Path(output[0]).write_text(message)
            Path(log[0]).write_text(message)
