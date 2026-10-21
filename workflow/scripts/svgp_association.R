# Firth inference needs an established statistical library; Snakemake schedules
# the job and tracks its inputs. No phenotype-driven model/locus selection occurs.
# Official manual: https://georgheinze.r-universe.dev/logistf/doc/manual.html
sink(snakemake@log[[1]], split = TRUE)
suppressPackageStartupMessages(library(logistf))
missing_values <- c("", ".", "NA", "unknown", "not_assessed")
read_tsv <- function(path) {
    con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
    on.exit(close(con))
    read.delim(con, colClasses = "character", check.names = FALSE,
               quote = "", comment.char = "", na.strings = missing_values)
}
require_columns <- function(x, columns) {
    if (anyDuplicated(names(x)) || !all(columns %in% names(x)))
        stop("Missing/duplicate required TSV columns: ", paste(setdiff(columns, names(x)), collapse = ", "))
}
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE,
                                          row.names = FALSE, na = ".")
split_columns <- function(x) if (is.na(x)) character() else trimws(strsplit(x, ";", fixed = TRUE)[[1]])
has_file <- function(x) !is.na(x) && file.exists(x) && !is.na(file.info(x)$size) && !file.info(x)$isdir && file.info(x)$size > 0
numeric_column <- function(x, name) {
    v <- suppressWarnings(as.numeric(x))
    if (any(!is.na(x) & (!is.finite(v)))) stop("Non-numeric values in ", name)
    v
}
samples <- read_tsv(snakemake@input[["samples"]])
evidence <- read_tsv(snakemake@input[["evidence"]])
panel <- read_tsv(snakemake@input[["panel"]])
locks <- read_tsv(snakemake@input[["lock"]])
release <- read_tsv(snakemake@input[["release"]])
design <- read_tsv(snakemake@input[["design"]])
tests <- read_tsv(snakemake@input[["tests"]])
design_fields <- c("contrast_id", "phenotype_column", "case_label", "control_label", "contrast_description",
                   "age_column", "age_definition", "numeric_covariates", "categorical_covariates",
                   "sensitivity_numeric_covariates", "sensitivity_categorical_covariates",
                   "min_pooled_carriers", "age_overlap_review", "ascertainment_description",
                   "design_frozen_utc", "technical_lock_utc", "phenotype_release_utc", "reviewer", "status", "review_evidence")
require_columns(design, design_fields)
require_columns(tests, c("test_id", "locus_id", "allele_id", "family", "panel_role", "analysis",
                         "conditional_snp", "selection_evidence", "snp_evidence"))
if (nrow(design) != 1 || !nrow(tests)) stop("Supply one frozen design and a nonempty test manifest")
d <- as.list(design[1, ])
required_design <- setdiff(design_fields, c("numeric_covariates", "categorical_covariates",
                                            "sensitivity_numeric_covariates", "sensitivity_categorical_covariates"))
if (anyNA(unlist(d[required_design])) || d$status != "APPROVED" || !has_file(d$review_evidence))
    stop("Design requires an approved, documented review and all required fields")
times <- as.POSIXct(unlist(d[c("design_frozen_utc", "technical_lock_utc", "phenotype_release_utc")]),
                    format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
if (anyNA(times) || any(times[1:2] >= times[3])) stop("Freeze design and technical rules before phenotype release")
if (d$case_label == d$control_label) stop("Case and comparator labels must differ")
minimum <- numeric_column(d$min_pooled_carriers, "min_pooled_carriers")
if (minimum < 10 || minimum != floor(minimum)) stop("Prespecify at least ten pooled informative carriers")
require_columns(samples, c("person_id", "analysis_eligible", "discovery_donor", "related_to_discovery_donor",
                           "unrelated_primary", "batch", "depth", "age_definition", d$phenotype_column, d$age_column))
require_columns(evidence, c("person_id", "locus_id", "allele_id", "graph_id", "graph_sha256", "reference_build",
                            "carrier_status", "callable", "measurement_type"))
require_columns(panel, c("locus_id", "allele_id", "panel_role", "chromosome"))
require_columns(locks, c("graph_id", "graph_sha256", "locus_id", "allele_id", "status", "technical_status", "locked_before_phenotypes"))
require_columns(release, c("graph_id", "status"))
graph <- snakemake@params[["graph_id"]]
if (sum(release$graph_id == graph & release$status == "APPROVED", na.rm = TRUE) != 1)
    stop("No approved cohort release for the assay graph")
if (anyNA(samples$person_id) || anyDuplicated(samples$person_id) || anyNA(tests$test_id) || anyDuplicated(tests$test_id))
    stop("Missing or duplicate person/test identifiers")
key <- function(x) paste(x$locus_id, x$allele_id, sep = "::")
if (anyNA(tests[c("locus_id", "allele_id", "family", "panel_role", "analysis", "selection_evidence")]))
    stop("Incomplete test manifest")
p <- match(key(tests), key(panel))
if (anyNA(p) || any(tests$panel_role != panel$panel_role[p]) || !all(vapply(tests$selection_evidence, has_file, logical(1))))
    stop("Every test must match the frozen panel and documented phenotype-blind selection")
if (any(!sub("^chr", "", panel$chromosome[p]) %in% as.character(1:22)))
    stop("This primary analysis supports autosomal loci only; sex chromosomes require a separate design")
candidate <- tests$panel_role == "candidate_addition"
if (!length(unique(tests$locus_id[candidate])) %in% 5:10) stop("Prespecify five to ten candidate loci")
if (length(unique(tests$family[candidate])) != 1) stop("Keep the candidate tests in one declared multiplicity family")
if (length(intersect(tests$family[candidate], tests$family[!candidate]))) stop("Keep controls in separate test families")
if (!all(unique(key(tests)[candidate]) %in% key(tests[tests$analysis == "carrier", ])))
    stop("Every selected candidate allele needs its primary carrier test")
if (anyDuplicated(paste(key(tests), tests$analysis, tests$conditional_snp))) stop("Duplicate prespecified analyses")
analyses <- c("carrier", "technical_sensitivity", "sv_matched", "snp_only", "joint")
if (any(!tests$analysis %in% analyses)) stop("Unknown prespecified analysis")
snp_rows <- tests$analysis %in% c("sv_matched", "snp_only", "joint")
if (anyNA(tests$conditional_snp[snp_rows]) || !all(vapply(tests$snp_evidence[snp_rows], has_file, logical(1))))
    stop("Matched SNP comparisons require a prespecified metadata SNP column and evidence")
for (group in unique(paste(key(tests)[snp_rows], tests$conditional_snp[snp_rows]))) {
    rows <- snp_rows & paste(key(tests), tests$conditional_snp) == group
    if (!all(c("sv_matched", "snp_only", "joint") %in% tests$analysis[rows]))
        stop("Declare SV-only, SNP-only and joint comparisons on identical complete cases")
}
keep <- with(samples, analysis_eligible == "yes" & discovery_donor == "no" &
             related_to_discovery_donor == "no" & unrelated_primary == "yes")
keep <- !is.na(keep) & keep & samples[[d$phenotype_column]] %in% c(d$case_label, d$control_label)
samples <- samples[keep, , drop = FALSE]
if (!nrow(samples)) stop("No independent, unrelated eligible people in this contrast")
samples$y <- as.integer(samples[[d$phenotype_column]] == d$case_label)
age <- numeric_column(samples[[d$age_column]], d$age_column)
samples$depth <- numeric_column(samples$depth, "depth")
base_numeric <- split_columns(d$numeric_covariates)
base_factor <- split_columns(d$categorical_covariates)
extra_numeric <- split_columns(d$sensitivity_numeric_covariates)
extra_factor <- split_columns(d$sensitivity_categorical_covariates)
all_covariates <- c(base_numeric, base_factor, extra_numeric, extra_factor)
require_columns(samples, unique(c(all_covariates, tests$conditional_snp[snp_rows])))
if (!length(c(base_numeric, base_factor)) || any(c("y", d$phenotype_column) %in% all_covariates))
    stop("Declare justified baseline covariates; the outcome cannot be a covariate")
if (length(intersect(c(base_numeric, extra_numeric), c(base_factor, extra_factor))))
    stop("A covariate cannot be both numeric and categorical")
for (name in unique(c(base_numeric, extra_numeric, tests$conditional_snp[snp_rows])))
    samples[[name]] <- numeric_column(samples[[name]], name)
design_reason <- ""
if (length(unique(samples$y)) != 2) design_reason <- "MISSING_OUTCOME_GROUP"
if (anyNA(age) || anyNA(samples$age_definition) || any(samples$age_definition != d$age_definition))
    design_reason <- "AGE_MISSING_OR_DEFINITION_MISMATCH"
if (d$age_overlap_review != "PASS") design_reason <- "SEVERE_OR_UNREVIEWED_AGE_NONOVERLAP"
if (!nzchar(design_reason) && max(min(age[samples$y == 0]), min(age[samples$y == 1])) >
    min(max(age[samples$y == 0]), max(age[samples$y == 1]))) design_reason <- "AGE_RANGES_DO_NOT_OVERLAP"
for (name in unique(c("batch", base_factor))) {
    cells <- split(samples$y, samples[[name]])
    if (!length(cells) || anyNA(samples[[name]]) || all(vapply(cells, function(z) length(unique(z)) == 1, logical(1))))
        design_reason <- paste0("PHENOTYPE_CONFOUNDED_OR_MISSING_", name)
}
evidence <- evidence[evidence$graph_id == graph & evidence$person_id %in% samples$person_id, , drop = FALSE]
if (anyNA(evidence[c("person_id", "locus_id", "allele_id", "graph_id", "graph_sha256")])) stop("Evidence has missing identifiers")
if (anyDuplicated(paste(evidence$person_id, key(evidence)))) stop("Duplicate person/allele evidence")
effects <- counts <- missingness <- power <- list()
for (i in seq_len(nrow(tests))) {
    t <- tests[i, ]; label <- paste(t$locus_id, t$allele_id, sep = "::")
    ev <- evidence[key(evidence) == label, , drop = FALSE]
    ev <- ev[match(samples$person_id, ev$person_id), , drop = FALSE]
    call <- !is.na(ev$callable) & ev$callable == "yes" & ev$measurement_type == "carrier" &
            ev$carrier_status %in% c("ALT_SUPPORTED", "REF_SUPPORTED")
    call[is.na(call)] <- FALSE
    carrier <- ifelse(call, as.integer(ev$carrier_status == "ALT_SUPPORTED"), NA_integer_)
    ncase <- sum(call & samples$y == 1); ncontrol <- sum(call & samples$y == 0)
    case_carriers <- sum(carrier[samples$y == 1] == 1, na.rm = TRUE)
    control_carriers <- sum(carrier[samples$y == 0] == 1, na.rm = TRUE)
    case_ci <- if (ncase) binom.test(case_carriers, ncase)$conf.int else c(NA_real_, NA_real_)
    control_ci <- if (ncontrol) binom.test(control_carriers, ncontrol)$conf.int else c(NA_real_, NA_real_)
    counts[[i]] <- data.frame(test_id = t$test_id, locus_id = t$locus_id, allele_id = t$allele_id,
        graph_id = graph, reference_build = paste(unique(na.omit(ev$reference_build)), collapse = ";"),
        attempted_cases = sum(samples$y == 1), attempted_controls = sum(samples$y == 0),
        callable_cases = ncase, callable_controls = ncontrol,
        carrier_cases = case_carriers, carrier_controls = control_carriers,
        case_carrier_freq_ci95_lower = case_ci[1], case_carrier_freq_ci95_upper = case_ci[2],
        control_carrier_freq_ci95_lower = control_ci[1], control_carrier_freq_ci95_upper = control_ci[2],
        uncertain_cases = sum(!call & samples$y == 1), uncertain_controls = sum(!call & samples$y == 0),
        measurement_type = "carrier", status = "DESCRIPTIVE", reason = "Genotype and allele count not inferred")
    for (group in c("phenotype", "batch", "depth_quartile")) {
        strata <- switch(group, phenotype = as.character(samples$y), batch = samples$batch,
                         depth_quartile = ifelse(is.na(samples$depth), "unknown", paste0("Q", pmin(4, ceiling(rank(samples$depth, na.last = "keep", ties.method = "min") / sum(!is.na(samples$depth)) * 4)))))
        for (level in unique(strata[!is.na(strata)])) {
            select <- !is.na(strata) & strata == level
            missingness[[length(missingness) + 1]] <- data.frame(test_id = t$test_id, group = group,
                stratum = level, attempted = sum(select), callable = sum(call[select]),
                missing_fraction = mean(!call[select]), status = "DESCRIPTIVE", reason = "Post-lock bias diagnostic")
        }
    }
    reason <- design_reason
    lock <- locks[locks$graph_id == graph & key(locks) == label, , drop = FALSE]
    if (nrow(lock) != 1 || is.na(lock$status[1]) || lock$status[1] != "LOCKED" ||
        is.na(lock$technical_status[1]) || lock$technical_status[1] != "PASS" ||
        is.na(lock$locked_before_phenotypes[1]) || lock$locked_before_phenotypes[1] != "yes") reason <- "NO_PHENOTYPE_BLIND_LOCK"
    if (nrow(lock) == 1 && (is.na(lock$graph_sha256[1]) || any(ev$graph_sha256[call] != lock$graph_sha256[1])))
        reason <- "MEASUREMENT_GRAPH_MISMATCH"
    use_numeric <- unique(c(base_numeric, if (t$analysis == "technical_sensitivity") extra_numeric))
    use_factor <- unique(c(base_factor, if (t$analysis == "technical_sensitivity") extra_factor))
    model <- data.frame(y = samples$y, carrier = carrier)
    for (j in seq_along(c(use_numeric, use_factor))) {
        name <- c(use_numeric, use_factor)[j]
        model[[paste0("cov", j)]] <- if (name %in% use_numeric) samples[[name]] else factor(samples[[name]])
    }
    if (snp_rows[i]) model$snp <- samples[[t$conditional_snp]]
    model <- droplevels(model[complete.cases(model), , drop = FALSE])
    pooled_carriers <- sum(model$carrier == 1)
    if (pooled_carriers < minimum) reason <- "SPARSE_POOLED_INFORMATIVE_CARRIER_COUNT"
    if (!sum(model$carrier == 0)) reason <- "NO_INFORMATIVE_NONCARRIER"
    if (length(unique(model$y)) < 2) reason <- "NO_CALLABLE_OUTCOME_GROUP"
    for (name in grep("^cov", names(model), value = TRUE)) if (is.factor(model[[name]])) {
        cells <- split(model$y, model[[name]])
        if (length(cells) && all(vapply(cells, function(z) length(unique(z)) == 1, logical(1)))) reason <- "PHENOTYPE_CONFOUNDED_WITH_COVARIATE"
    }
    predictors <- c(if (t$analysis != "snp_only") "carrier", grep("^cov", names(model), value = TRUE),
                    if (t$analysis %in% c("snp_only", "joint")) "snp")
    formula <- reformulate(predictors, response = "y")
    x <- tryCatch(model.matrix(formula, model), error = function(e) NULL)
    if (is.null(x) || nrow(x) <= ncol(x) || qr(x)$rank < ncol(x)) reason <- "CONFOUNDED_OR_RANK_DEFICIENT_MODEL"
    estimate <- lower <- upper <- pvalue <- NA_real_
    term <- if (t$analysis == "snp_only") "snp" else "carrier"
    if (!nzchar(reason)) {
        fit <- tryCatch(withCallingHandlers(logistf(formula, data = model, firth = TRUE, pl = TRUE,
                        plconf = match(term, colnames(x)), control = logistf.control(maxit = 100)),
                        warning = function(w) stop(conditionMessage(w))),
                        error = function(e) e)
        if (inherits(fit, "error")) reason <- paste("FIT_FAILED", conditionMessage(fit)) else {
            estimate <- exp(fit$coefficients[term]); lower <- exp(fit$ci.lower[term]); upper <- exp(fit$ci.upper[term]); pvalue <- fit$prob[term]
            if (any(!is.finite(c(estimate, lower, upper, pvalue))) || lower <= 0 || upper < lower || pvalue < 0 || pvalue > 1) reason <- "INVALID_PROFILE_RESULT"
        }
    }
    if (nzchar(reason)) estimate <- lower <- upper <- pvalue <- NA_real_
    correlation <- if (snp_rows[i] && nrow(model) > 2 && sd(model$carrier) > 0 && sd(model$snp) > 0) cor(model$carrier, model$snp) else NA_real_
    missing_table <- table(factor(samples$y, levels = 0:1), factor(call, levels = c(FALSE, TRUE)))
    missing_p <- if (all(rowSums(missing_table) > 0) && all(colSums(missing_table) > 0)) fisher.test(missing_table)$p.value else NA_real_
    effects[[i]] <- data.frame(t, contrast_id = d$contrast_id, graph_id = graph, measurement_type = "carrier",
        effect_unit = if (t$analysis == "snp_only") "per supplied SNP dosage unit" else "ALT-supported carrier versus REF-supported",
        n = nrow(model), model_cases = sum(model$y == 1), model_controls = sum(model$y == 0),
        model_carrier_cases = sum(model$carrier == 1 & model$y == 1),
        model_carrier_controls = sum(model$carrier == 1 & model$y == 0),
        pooled_carriers = pooled_carriers, pooled_noncarriers = sum(model$carrier == 0),
        odds_ratio = estimate, ci95_lower = lower, ci95_upper = upper, p_raw = pvalue,
        missingness_fisher_p = missing_p, carrier_snp_correlation = correlation,
        independent_carrier_snp_states = if (snp_rows[i]) nrow(unique(model[c("carrier", "snp")])) else NA_integer_,
        status = if (nzchar(reason)) "SKIPPED" else "ESTIMATED", reason = if (nzchar(reason)) reason else "PASS")
    # An approximate, unadjusted 2x2 power grid, using pooled carrier frequency
    # and the actual complete-case split. This is not adjusted Firth model power.
    n1 <- sum(model$y == 1); n0 <- sum(model$y == 0); freq <- if (nrow(model)) mean(model$carrier) else NA_real_
    burden <- sum(tests$family == t$family); alpha <- 0.05 / burden
    ors <- c(seq(1, 5, by = 0.25), 6, 8, 10, 15, 20)
    values <- rep(NA_real_, length(ors))
    if (n0 > 0 && n1 > 0 && is.finite(freq) && freq > 0 && freq < 1) for (z in seq_along(ors)) {
        odds <- ors[z]
        p0 <- uniroot(function(q) (n0 * q + n1 * odds * q / (1 - q + odds * q)) / (n0 + n1) - freq, c(0, 1))$root
        p1 <- odds * p0 / (1 - p0 + odds * p0)
        ncp <- (p1 - p0)^2 / (freq * (1 - freq) * (1 / n0 + 1 / n1))
        values[z] <- pchisq(qchisq(1 - alpha, 1), 1, ncp = ncp, lower.tail = FALSE)
    }
    detectable <- ors[which(values >= 0.8)[1]]
    power[[i]] <- data.frame(test_id = t$test_id, n_cases = n1, n_controls = n0,
        pooled_carrier_frequency = freq, family_tests = burden, alpha_bonferroni = alpha,
        assumed_or = ors, approximate_power = values, minimum_grid_or_80pct_power = detectable,
        status = "APPROXIMATE", reason = "Unadjusted carrier 2x2 noncentral chi-square; not covariate-adjusted Firth power; larger-OR direction only")
}
effects <- do.call(rbind, effects)
# Failed/sparse tests remain in the declared family with p=1 for correction.
effects$p_holm <- ave(ifelse(is.na(effects$p_raw), 1, effects$p_raw), effects$family, FUN = function(x) p.adjust(x, "holm"))
write_tsv(effects, snakemake@output[["effects"]])
write_tsv(do.call(rbind, counts), snakemake@output[["counts"]])
write_tsv(do.call(rbind, missingness), snakemake@output[["missingness"]])
write_tsv(do.call(rbind, power), snakemake@output[["power"]])
write_tsv(tests, snakemake@output[["tests"]])
write_tsv(effects[c("test_id", "locus_id", "allele_id", "graph_id", "status", "reason")], snakemake@output[["status"]])
writeLines(c("# Prespecified SV carrier association", paste(names(d), unlist(d), sep = ": "),
    "", paste("Eligible independent unrelated people:", nrow(samples)), paste("Design diagnostic:", if (nzchar(design_reason)) design_reason else "PASS"),
    "Primary analysis excludes graph donors and their close relatives. Technical sensitivities are explicit manifest rows.",
    "No genotype dosage or allele frequency is inferred from carrier status. No carrier-in-both-groups filter is used.",
    "Holm correction includes all declared tests within each family, including skipped tests as p=1.",
    "Missingness P values are descriptive diagnostics, not additional association claims.",
    "Discovery ascertainment persists despite phenotype-blind selection. Conditional associations do not prove causality.",
    "Power is an approximate unadjusted carrier-table grid, ignoring covariates, Firth penalization and ascertainment.",
    "No local SNP scan, prediction or phenotype-driven covariate/model selection is performed."), snakemake@output[["design"]])
writeLines(c(paste("R", getRversion()), paste("logistf", packageVersion("logistf")),
    "Official stable release checked 2026-10-20: logistf 1.26.1",
    "https://georgheinze.r-universe.dev/logistf/doc/manual.html"), snakemake@output[["versions"]])
svg(snakemake@output[["plot"]], width = 8, height = max(3, 0.35 * nrow(effects) + 1.5))
shown <- effects[effects$status == "ESTIMATED", , drop = FALSE]
if (!nrow(shown)) { plot.new(); text(0.5, 0.5, "No estimable prespecified association tests") } else {
    par(mar = c(5, 12, 2, 1))
    plot(shown$odds_ratio, seq_len(nrow(shown)), log = "x", yaxt = "n", ylab = "", xlab = "Odds ratio (profile likelihood 95% interval)",
         xlim = range(c(shown$ci95_lower, shown$ci95_upper)), pch = 19)
    segments(shown$ci95_lower, seq_len(nrow(shown)), shown$ci95_upper, seq_len(nrow(shown)))
    axis(2, at = seq_len(nrow(shown)), labels = shown$test_id, las = 2); abline(v = 1, lty = 2)
}
dev.off()
sink()
