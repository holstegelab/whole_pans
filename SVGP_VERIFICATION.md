# Server verification

No tests, syntax checks, environment solves or dry runs were executed locally.
Repository `AGENTS.md` requires all validation on the server. Static review cannot
establish that the installed server tools, conda package builds, real input
schemas and graph indexes are compatible.

Run from the deployed repository on the server, using its existing Snakemake
environment and submission profile. Preserve existing inputs and historical
results. Do not run the broad historical `review_followup` aggregate.

```bash
cd /gpfs/work3/0/qtholstg/pangenome/whole_pans
python -m py_compile workflow/scripts/svgp_*.py tests/test_sv*.py
Rscript -e 'invisible(parse(file="workflow/scripts/svgp_association.R"))'
snakemake -s workflow/rules/sv_gene_phenotype.smk -n -p inventory
```

For both fixture suites, activate the existing server environment used by the
phenotype import, providing `openpyxl` (specified in
`workflow/envs/svgp_phenotypes.yaml`), then run:

```bash
python -m unittest discover -s tests -p 'test_sv*.py' -v
```

This pattern includes both `test_sv_gene_phenotype.py` and
`test_svgp_supplied_inputs.py`. The original fixtures use Python's standard
library and the existing `review_followup.py` helper; the supplied-input fixtures
also use `openpyxl` for tiny temporary workbooks. Both create tiny temporary
files. The graph reader
is mocked, so they do not invoke `vg` or read research datasets. They cover
explicit/ambiguous sample aliases and unknown eligibility;
fingerprint confirmation with unavailable evidence; release review binding to
the actual pilot/ONT inputs and equivalent relative path spellings; both source
orientations; zero-length reference insertion boundaries and VCF anchors;
unmapped-versus-unresolved resource decisions; controls excluded from local
addition denominators; duplicate molecules and overlapping target intervals;
separate multiallelic carrier no-calls and missing diploid dosage; insertion
feature boundaries and exact TSS distances. These are regression checks, not
biological assay validation.

Populate the supplied manifests and complete stage-specific reviews before
requesting later targets. Header-only templates record missing inputs; they do
not establish a usable cohort, panel, measurement lock or analysis design.

The 21 October update adds the supplied clinical workbook and existing-call
sources. Keep the workbook unchanged and inspect both clinical identifier joins
and raw caller sample labels. The default inventory now includes workbook import
and VCF header inventory, with no caller reruns or raw-read mapping. On the server:

```bash
snakemake -s workflow/rules/sv_gene_phenotype.smk -n -p phenotypes
snakemake -s workflow/rules/sv_gene_phenotype.smk -n -p reported_ont
snakemake -s workflow/rules/sv_gene_phenotype.smk -n -p existing_calls_inventory
```

After reviewing explicit crosswalks, selecting one gCNV callset and freezing a
small ExpansionHunter repeat list, dry-run `gcnv_import` and `repeat_import`
separately. Check that only existing single-sample VCF imports appear. Inspect one
file from each source before the batch import. Confirm raw gCNV CN/GT/FILTER and
EH repeat-unit/CI/missing states survive; caller evidence is not a validated graph
allele. GATK segment intervals and ExpansionHunter repeat intervals use different
VCF anchoring rules; see their input guides. The user confirmed svPPA versus
source-labelled controls and age at sampling (current age in this new cohort).
Confirm the imported `source_age`, numeric `age_at_sampling` and
`age_definition` fields preserve missing values and only enrich explicitly
matched people. Check that `other_dementia` remains in metadata and is excluded
only by the primary contrast. The configured `svgp_association_design.tsv` remains
`DRAFT`: review age overlap, ascertainment, covariates and design-freeze evidence
before approval. Missing ages continue to block inference under the current
conservative age gate; the workbook does not verify healthy or centenarian
control status. Prediction remains disabled.

```bash
snakemake -s workflow/rules/sv_gene_phenotype.smk -n -p catalogue
snakemake -s workflow/rules/sv_gene_phenotype.smk -n -p representation
snakemake -s workflow/rules/sv_gene_phenotype.smk -n -p panel
snakemake -s workflow/rules/sv_gene_phenotype.smk -n -p svgp_pilot_mapping
snakemake -s workflow/rules/sv_gene_phenotype.smk -n -p pilot
snakemake -s workflow/rules/sv_gene_phenotype.smk -n -p ont_validation
snakemake -s workflow/rules/sv_gene_phenotype.smk -n -p annotation
snakemake -s workflow/rules/sv_gene_phenotype.smk -n -p association
snakemake -s workflow/rules/sv_gene_phenotype.smk -n -p prediction
```

Check printed dependencies before execution: no Cactus build, V2 rescans,
unrequested historical mapping, local mount prefixes or cohort-wide two-graph
mapping should appear. Only the bounded pilot maps each configured person to
two graphs; new cohort people use one main graph. The default `inventory`
target must schedule no mapping.

For the first native mapping check, keep the full frozen pilot selection and
request the two mapping directories for one explicitly selected pilot person.
Replace `PILOT_PERSON_ID` below, and adapt the results prefix if configured
differently. Use the existing server profile for execution after inspecting
this dry run:

```bash
snakemake -s workflow/rules/sv_gene_phenotype.smk -n -p \
  ../whole_pangenome/sv_gene_phenotype_2026-10-20/mapping/hprc_only/PILOT_PERSON_ID \
  ../whole_pangenome/sv_gene_phenotype_2026-10-20/mapping/hprc_plus_local/PILOT_PERSON_ID
```

Inspect
mapping counts, assigned molecule completeness, junction matches, depth,
reference/ALT support, no-call reasons and benchmark files before requesting
`pilot`. A complete calibrated panel needs the planned full pilot and ONT;
one person only checks mechanics. Do not approve cohort release from that
single-person check.

After full pilot/ONT review, inspect `cohort` with `-n -p` before releasing
cohort jobs. Check that configured simultaneous CPU/memory/disk demands fit the
existing profile; benchmark extrapolation is an estimate. Temporary-disk use
marked `NOT_MEASURED` needs separate server measurement before budget approval.
Keep no-call cells and callable carrier denominators visible. This carrier
assay does not estimate diploid dosage or allele frequency.

Review generated association status, carrier counts, profile-likelihood
intervals and the complete adjusted test family. Prediction is optional and
disabled by default. Its uncertainty interval is conditional on fixed held-out
predictions; repeated split results describe additional stability, not new
independent people. No numerical result or positive biological claim is implied
by having prepared this workflow.

Refresh the read-only evidence snapshot after running each stage:

```bash
snakemake -s workflow/rules/sv_gene_phenotype.smk --cores 1 --forcerun svgp_study_status study_status
```

An available output is not proof of passed scientific gates. Read the row-level
statuses, rejection ledger, no-call reasons and frozen review evidence.
