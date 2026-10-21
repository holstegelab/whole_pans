"""Small fixed-feature prediction comparison; Snakemake owns scheduling and gates.

Only sklearn fitting requires Python. Features are fixed before phenotype release;
imputation, encoding, scaling and ridge tuning are refitted within training folds.
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from joblib import parallel_backend
from sklearn.compose import ColumnTransformer
from sklearn.exceptions import ConvergenceWarning
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, log_loss, roc_auc_score
from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from threadpoolctl import threadpool_limits


def read(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False,
                       na_values=["", ".", "NA", "unknown", "not_assessed"])


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def write(rows, path, columns=None):
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False, na_rep=".")


def existing(path):
    return isinstance(path, str) and Path(path).is_file() and Path(path).stat().st_size > 0


def columns(value):
    return [] if pd.isna(value) else [x.strip() for x in value.split(";") if x.strip()]


def skip(reason):
    for name, path in sm.output.items():
        if name not in {"report", "versions"}:
            write([], path, ["stage", "status", "reason"])
    write([dict(stage="prediction", status="SKIPPED", reason=reason)], sm.output.status)
    Path(sm.output.report).write_text("Prediction was not evaluated: " + reason + "\n")
    sys.exit(0)


def splits(y, groups, folds, seed):
    cv = StratifiedGroupKFold(n_splits=folds, shuffle=True, random_state=seed)
    result = list(cv.split(np.zeros(len(y)), y, groups))
    if any(len(np.unique(y[a])) != 2 or len(np.unique(y[b])) != 2
           or set(groups[a]) & set(groups[b]) for a, b in result):
        skip("Insufficient class/group structure for the prespecified nested splits")
    return result


def model(numeric, categorical):
    transforms = []
    if numeric:
        transforms.append(("numeric", make_pipeline(
            SimpleImputer(strategy="median", keep_empty_features=True), StandardScaler()), numeric))
    if categorical:
        transforms.append(("categorical", make_pipeline(
            SimpleImputer(strategy="constant", fill_value="MISSING", keep_empty_features=True),
            OneHotEncoder(handle_unknown="ignore", sparse_output=False)), categorical))
    return make_pipeline(ColumnTransformer(transforms), LogisticRegression(
        l1_ratio=0.0, solver="lbfgs", max_iter=2000, tol=1e-7))


def scores(y, probability):
    p = np.clip(probability, np.finfo(float).eps, 1 - np.finfo(float).eps)
    intercept = slope = np.nan
    calibration_status = "UNESTIMABLE"
    separated = max(p[y == 0]) <= min(p[y == 1]) or max(p[y == 1]) <= min(p[y == 0])
    if separated:
        calibration_status = "SEPARATED_OR_CONSTANT"
    if np.ptp(p) > 0 and not separated:
        try:
            fit = LogisticRegression(C=np.inf, solver="lbfgs", max_iter=2000, tol=1e-8)
            fit.fit(np.log(p / (1 - p)).reshape(-1, 1), y)
            intercept, slope = fit.intercept_[0], fit.coef_[0, 0]
            calibration_status = "DESCRIPTIVE_IN_SAMPLED_CONTRAST"
        except (ValueError, ConvergenceWarning):
            pass
    return dict(log_loss=log_loss(y, p, labels=[0, 1]), auc=roc_auc_score(y, p),
                brier=brier_score_loss(y, p), calibration_intercept=intercept,
                calibration_slope=slope, calibration_status=calibration_status)


sm = snakemake
Path(sm.log[0]).parent.mkdir(parents=True, exist_ok=True)
sys.stdout = sys.stderr = open(sm.log[0], "w")
for path in sm.output:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
Path(sm.output.versions).write_text("\n".join([
    "Python " + sys.version.split()[0], "scikit-learn " + sklearn.__version__,
    "numpy " + np.__version__, "pandas " + pd.__version__,
    "Stable scikit-learn 1.9.1 checked 2026-10-20",
    "https://scikit-learn.org/stable/whats_new/v1.9.html",
    "https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html"]) + "\n")
warnings.filterwarnings("error", category=ConvergenceWarning)
settings = sm.params.settings
design = read(sm.input.design)
association_design = read(sm.input.association_design)
features = read(sm.input.features)
require(len(design) == len(association_design) == 1 and len(features) > 0,
        "Supply one frozen prediction design, one association design and explicit features")
d, a = design.iloc[0], association_design.iloc[0]
require(design.notna().all().all() and d.status == "APPROVED" and existing(d.review_evidence),
        "Prediction design must be complete, approved and documented")
require(a.status == "APPROVED" and d.contrast_id == a.contrast_id and d.phenotype_blind_features == "yes",
        "Reuse the approved contrast and freeze phenotype-blind features")
frozen = pd.to_datetime(d.design_frozen_utc, utc=True)
release_time = pd.to_datetime(d.phenotype_release_utc, utc=True)
require(frozen < release_time and d.phenotype_release_utc == a.phenotype_release_utc,
        "Prediction design must predate the same documented phenotype release")
require(features.feature_id.notna().all() and not features.feature_id.duplicated().any()
        and features.status.eq("APPROVED").all() and features.review_evidence.map(existing).all()
        and features.selection_basis.isin(["external_prespecified", "phenotype_blind_panel"]).all()
        and (pd.to_datetime(features.frozen_utc, utc=True) <= frozen).all(),
        "Every feature requires a documented, prespecified definition; association-hit selection is forbidden")
require(features.feature_type.isin(["SNP", "SV_CARRIER"]).all(), "Only supplied SNP information and calibrated carriers are supported")
snps = features[features.feature_type == "SNP"]
svs = features[features.feature_type == "SV_CARRIER"]
require(len(snps) > 0 and snps.selection_basis.eq("external_prespecified").all()
        and svs.selection_basis.eq("phenotype_blind_panel").all()
        and 5 <= svs.locus_id.nunique() <= 10 and len(svs) <= 10,
        "Use externally specified SNP columns and five to ten phenotype-blind SV features")
require(snps.metadata_column.notna().all() and not snps.metadata_column.duplicated().any()
        and svs[["locus_id", "allele_id"]].notna().all().all()
        and not svs.duplicated(["locus_id", "allele_id"]).any(), "Incomplete or duplicated feature definitions")
samples = read(sm.input.samples)
require(not set(svs.feature_id) & set(samples.columns), "SV feature identifiers must differ from all metadata column names")
require(samples.person_id.notna().all() and not samples.person_id.duplicated().any(), "Expected exactly one row per person")
eligible = (samples.analysis_eligible.eq("yes") & samples.discovery_donor.eq("no")
            & samples.related_to_discovery_donor.eq("no") & samples.unrelated_primary.eq("yes")
            & samples[a.phenotype_column].isin([a.case_label, a.control_label]))
samples = samples[eligible].sort_values("person_id").reset_index(drop=True)
if not len(samples):
    skip("No eligible independent unrelated people in the frozen contrast")
y = samples[a.phenotype_column].eq(a.case_label).astype(int).to_numpy()
require(samples[d.group_column].notna().all(), "Family/person grouping must be complete")
groups = samples[d.group_column].to_numpy()
if len(np.unique(y)) != 2 or min(np.bincount(y)) < 5 or len(np.unique(groups)) < 5:
    skip("Insufficient people/classes/groups for outer five-fold evaluation")
require(a.age_overlap_review == "PASS" and samples.age_definition.eq(a.age_definition).all(), "Unreviewed age definitions/overlap")
age = pd.to_numeric(samples[a.age_column], errors="raise")
if age.isna().any() or max(age[y == 0].min(), age[y == 1].min()) > min(age[y == 0].max(), age[y == 1].max()):
    skip("Age missingness or nonoverlap makes the sampled contrast unsuitable for prediction")
numeric, categorical = columns(a.numeric_covariates), columns(a.categorical_covariates)
require(len(numeric + categorical) > 0 and len(set(numeric + categorical)) == len(numeric + categorical), "Declare distinct baseline covariates")
require(not set(numeric + categorical + snps.metadata_column.tolist()) & {a.phenotype_column, "person_id", d.group_column}, "Outcome/identity/group cannot be predictors")
require(not set(numeric + categorical) & set(snps.metadata_column), "SNP columns must add to the covariate baseline")
for name in set(["batch"] + categorical):
    if samples[name].isna().any() or samples.groupby(name)[a.phenotype_column].nunique().max() < 2:
        skip("Outcome is confounded with, or missing, batch/categorical covariates")
graph = sm.params.graph_id
release = read(sm.input.release)
require(((release.graph_id == graph) & (release.status == "APPROVED")).sum() == 1, "No approved cohort gate")
if not read(sm.input.association).status.eq("ESTIMATED").any():
    skip("Earlier association design produced no estimable tests; retain the main evidence deliverables")
locks = read(sm.input.lock)
panel, annotation, ont = read(sm.input.panel), read(sm.input.annotation), read(sm.input.ont)
evidence = read(sm.input.evidence)
evidence = evidence[evidence.graph_id.eq(graph) & evidence.person_id.isin(samples.person_id)]
require(not evidence.duplicated(["person_id", "locus_id", "allele_id"]).any(), "Duplicate assay observations")
for row in svs.itertuples():
    same = lambda table: table.locus_id.eq(row.locus_id) & table.allele_id.eq(row.allele_id)
    lock = locks[same(locks) & locks.graph_id.eq(graph)]
    require(len(lock) == 1 and lock.iloc[0].status == "LOCKED" and lock.iloc[0].technical_status == "PASS"
            and lock.iloc[0].locked_before_phenotypes == "yes", "SV feature lacks a passing phenotype-blind lock")
    require((same(panel) & panel.panel_role.eq("candidate_addition")).any()
            and (same(annotation) & annotation.status.eq("ANNOTATED_CONTEXT")).any(), "SV feature lacks passing panel/annotation evidence")
    require(panel.loc[same(panel), "chromosome"].str.removeprefix("chr").isin(
        [str(chromosome) for chromosome in range(1, 23)]).all(), "Primary SV prediction supports autosomal loci only")
    if not (same(ont) & ont.status.isin(["RESOLVED_CARRIER", "RESOLVED_GENOTYPE"]) & ont.carrier_status.eq("ALT_SUPPORTED")).any():
        skip("A prespecified SV feature has no resolved ONT ALT corroboration: " + row.feature_id)
    ev = evidence[same(evidence)].set_index("person_id").reindex(samples.person_id)
    require(ev.graph_sha256.dropna().eq(lock.iloc[0].graph_sha256).all(), "Prediction evidence/locked graph mismatch")
    call = ev.callable.eq("yes") & ev.measurement_type.eq("carrier")
    samples[row.feature_id] = ev.carrier_status.map({"ALT_SUPPORTED": 1.0, "REF_SUPPORTED": 0.0}).where(call).to_numpy()
    if not samples[row.feature_id].notna().any():
        skip("A prespecified SV feature has no callable eligible people: " + row.feature_id)
    features.loc[features.feature_id.eq(row.feature_id), "callable_people"] = samples[row.feature_id].notna().sum()
sv_columns, snp_columns = svs.feature_id.tolist(), snps.metadata_column.tolist()
require(not set(sv_columns) & set(numeric + categorical + snp_columns), "SV feature names collide with metadata predictors")
for name in numeric + snp_columns:
    samples[name] = pd.to_numeric(samples[name], errors="raise")
    require(np.isfinite(samples[name].dropna()).all(), "Nonfinite predictor values")
    if not samples[name].notna().any():
        skip("A frozen numeric predictor has no observed values: " + name)
numeric_models = {"M0": numeric, "M1": numeric + snp_columns, "M2": numeric + snp_columns + sv_columns}
seeds = [int(seed) for seed in settings.get("seeds", [20261020, 20261021, 20261022])]
grid = [float(c) for c in settings.get("c_grid", [0.1, 1.0, 10.0])]
require(len(seeds) >= 2 and len(set(seeds)) == len(seeds) and 0 < len(grid) <= 5
        and all(np.isfinite(c) and c > 0 for c in grid), "Prespecify distinct stability seeds and a small positive C grid")
predictions, manifests, metrics = [], [], []
for repeat, seed in enumerate(seeds):
    outer = splits(y, groups, 5, seed)
    probabilities = {name: np.full(len(y), np.nan) for name in numeric_models}
    for fold, (train, test) in enumerate(outer, 1):
        inner = splits(y[train], groups[train], 3, seed + fold)
        for level, partitions, indices in (("outer", [(train, test)], np.arange(len(y))), ("inner", inner, train)):
            for inner_fold, (_, validation) in enumerate(partitions, 1):
                for idx in indices[validation]:
                    manifests.append(dict(seed=seed, outer_fold=fold, level=level, inner_fold=inner_fold if level == "inner" else ".",
                                          person_id=samples.person_id[idx], group_id=groups[idx], role="validation"))
        for name, nums in numeric_models.items():
            search = GridSearchCV(model(nums, categorical), {"logisticregression__C": grid},
                                  scoring="neg_log_loss", cv=inner, n_jobs=sm.threads, error_score="raise")
            # Threaded fits use the job's allocated cores; BLAS gets one thread.
            with parallel_backend("threading"), threadpool_limits(limits=1):
                search.fit(samples.iloc[train], y[train])
                probabilities[name][test] = search.predict_proba(samples.iloc[test])[:, 1]
            for idx in test:
                predictions.append(dict(seed=seed, person_id=samples.person_id[idx], group_id=groups[idx],
                                        outer_fold=fold, outcome=y[idx], model=name,
                                        probability=probabilities[name][idx], selected_C=search.best_params_["logisticregression__C"]))
    for name, p in probabilities.items():
        require(np.isfinite(p).all(), "Every model needs one held-out prediction per person")
        metrics.append(dict(seed=seed, model=name, n=len(y), **scores(y, p)))
    metrics.append(dict(seed=seed, model="M2_minus_M1", n=len(y),
                        log_loss=log_loss(y, probabilities["M2"]) - log_loss(y, probabilities["M1"]),
                        auc=roc_auc_score(y, probabilities["M2"]) - roc_auc_score(y, probabilities["M1"]),
                        brier=brier_score_loss(y, probabilities["M2"]) - brier_score_loss(y, probabilities["M1"])))
    if repeat == 0:
        primary = probabilities
predictions = pd.DataFrame(predictions)
predictions[predictions.seed.eq(seeds[0])].to_csv(sm.output.oof, sep="\t", index=False)
predictions.to_csv(sm.output.repeats, sep="\t", index=False)
write(manifests, sm.output.splits)
write(metrics, sm.output.metrics)
features.to_csv(sm.output.features, sep="\t", index=False, na_rep=".")
# Conditional scoring interval only: fixed predictions, group-level resampling.
losses = {}
for name in ["M1", "M2"]:
    p = np.clip(primary[name], np.finfo(float).eps, 1 - np.finfo(float).eps)
    losses[name] = -(y * np.log(p) + (1 - y) * np.log(1 - p))
delta = losses["M2"] - losses["M1"]
members = [np.flatnonzero(groups == group) for group in np.unique(groups)]
rng = np.random.default_rng(seeds[0])
replicates = int(settings.get("bootstrap_replicates", 1000))
require(replicates >= 1000, "Use at least 1000 conditional bootstrap replicates")
bootstrap = [delta[np.concatenate([members[g] for g in rng.integers(len(members), size=len(members))])].mean() for _ in range(replicates)]
low, high = np.quantile(bootstrap, [0.025, 0.975])
write([dict(metric="M2_minus_M1_log_loss", estimate=delta.mean(), conditional_ci95_lower=low,
            conditional_ci95_upper=high, unit="person/family group", replicates=replicates,
            status="CONDITIONAL_ONLY", reason="Fixed OOF predictions; excludes training and split uncertainty")], sm.output.uncertainty)
Path(sm.output.report).write_text("\n".join([
    "# Preliminary internal prediction comparison", "Baseline: " + d.baseline_description,
    "SNP definition: " + d.snp_definition, "APOE: " + d.apoe_description,
    "Existing linear SV comparison: " + d.linear_sv_comparison,
    "Seeds: " + str(seeds) + "; first seed is primary. Ridge C grid: " + str(grid),
    "Same eligible independent unrelated people and identical outer 5/inner 3 grouped stratified splits for all models.",
    "M0: frozen covariates; M1: M0 plus supplied SNP information; M2: M1 plus frozen calibrated SV carrier features.",
    "Negative M2-minus-M1 log loss/Brier is improvement; positive AUC difference is improvement.",
    "Missing predictors remain missing until training-fold imputation. No global feature screening or association-hit selection.",
    "Repeated nested partitions describe instability, with the same people; folds/repeats are not independent observations.",
    "The bootstrap interval conditions on primary fixed OOF predictions. It excludes training/split uncertainty and is not a full model-performance interval.",
    "Full retraining uncertainty and batch-held-out sensitivity are not implemented; do not claim a definitive performance benefit from this module.",
    "Calibration/Brier/probabilities describe the sampled contrast, not population disease risk. No clinical utility or external validation is established.",
    "Discovery ascertainment remains despite excluding discovery donors/relatives. ONT source support is not independent donor validation.",
    "https://scikit-learn.org/stable/auto_examples/model_selection/plot_nested_cross_validation_iris.html"]) + "\n")
write([dict(stage="prediction", status="EVALUATED_CONDITIONAL_UNCERTAINTY", reason="Preliminary internal comparison; see report limitations")], sm.output.status)
