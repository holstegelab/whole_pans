# SV alleles, gene context and phenotype workflow

This standalone workflow implements the 20 October 2026 coding plan. Start in
`/gpfs/work3/0/qtholstg/pangenome/whole_pans` on the server. It reuses the existing
SV graph and provenance helper; it does not depend on Cactus, V2 rescue jobs,
EdgeDepth or a four-reference cohort remapping project.

The default target is **inventory only**. No genomic analysis has been run as
part of writing this implementation. Header-only example manifests report
missing inputs; they are not a cohort of zero observed people.

## Inputs supplied on 21 October

The configuration now points to the clinical workbook and existing gCNV and
ExpansionHunter results supplied by the user. Their source paths are configurable
and are translated to repository-relative server locations. The original files
and existing analyses remain read-only inputs.

- [Phenotype input details](SVGP_PHENOTYPE_INPUTS.md): `Blad1` has 397 sequencing
  records. `Sample_ID_WGS_given` is the unique workbook key; the other identifier
  columns have conflicting duplicate groups. Fill the explicit phenotype-to-person
  crosswalk before joining clinical rows to curated people. Both diagnosis columns
  and `DNA_from` are retained. `Diagnosis_svPPA` keeps its three original labels.
  The confirmed primary contrast is `svPPA` versus `control`; `other_dementia`
  remains in metadata and is excluded only from that contrast. The user confirmed
  age at sampling (current age in this new cohort); the importer preserves raw
  `source_age` and records numeric `age_at_sampling` without imputing missing ages.
- [gCNV input details](SVGP_GCNV_INPUTS.md): standard and WGS-autosome analyses
  overlap. Select **one** callset after reviewing its QC. Copy-number estimates,
  filters and quality fields remain raw caller evidence; they are not exact
  insertion-allele genotypes. Both inspected runs use a custom masked/decoy GRCh38
  reference and an older recorded GATK snapshot.
- [ExpansionHunter input details](SVGP_EXPANSIONHUNTER_INPUTS.md): choose a small
  phenotype-blind repeat subset. Repeat-unit estimates, confidence bounds, support
  and low-depth/missing states are retained. Exact VCF sample names, including
  `.markdup`, require explicit person crosswalks. These are repeat counts, not
  structural-allele copy dosage.
- [Reported ONT runs](SVGP_ONT_INPUTS.md): the supplied assignments retain seven
  exact sample labels across five runs. `2310032p01` has two runs. The small
  manifest records reported availability; person identities, read paths and SNP
  fingerprints remain to be checked before ONT measurement.

The default inventory additionally reads VCF headers and imports workbook metadata.
It does not rerun gCNV/ExpansionHunter, scan BAM/CRAM data or start graph mapping.
`phenotypes` imports the workbook and enriches only existing curated people.
`existing_calls_inventory` lists both caller sources without extracting calls.
`reported_ont` records the supplied ONT assignments without reading sequence data.

After filling the shared call crosswalk, choosing the gCNV callset and providing
the prespecified repeat list, inspect these separate import targets on the server:

```bash
snakemake -s workflow/rules/sv_gene_phenotype.smk -n -p phenotypes
snakemake -s workflow/rules/sv_gene_phenotype.smk -n -p reported_ont
snakemake -s workflow/rules/sv_gene_phenotype.smk -n -p existing_calls_inventory
snakemake -s workflow/rules/sv_gene_phenotype.smk -n -p gcnv_import
snakemake -s workflow/rules/sv_gene_phenotype.smk -n -p repeat_import
```

Imports use one existing VCF per job; a selected repeat catalog is prepared once.
For a bounded import, the existing profile can limit `existing_call_jobs`, for
example `--resources existing_call_jobs=16 --cores 16`. Tables are streamed into
separate compressed caller evidence files. Missing calls never become normal
alleles, and technical replicates never become independent people. No imported
CNV or repeat feature enters association or prediction automatically. Same-person,
same-locus/build and measurement/QC review is required first; these calls from the
same Illumina reads are not an independent ONT truth set. The caller reference
FASTA checksum is still unknown and must be frozen before stronger comparisons.

XLSX parsing uses a dedicated openpyxl 3.1.5 environment. The existing-call imports
use Python's standard library; the original caller environments are unchanged.
These additions were inspected statically, not executed as local workflows.

## Start here

Edit `config/sv_gene_phenotype.yaml`. Copy each needed `config/svgp_*.example.tsv`
to a study input filename before filling it, then point the configuration at the
copy. All examples contain headers only. The configured association design,
`config/svgp_association_design.tsv`, contains the confirmed contrast and sampling
age as a proposed baseline adjustment; it is a draft requiring the remaining
study review. Keep historical tables and inputs intact.

```bash
cd /gpfs/work3/0/qtholstg/pangenome/whole_pans
snakemake -s workflow/rules/sv_gene_phenotype.smk -n -p inventory
snakemake -s workflow/rules/sv_gene_phenotype.smk --use-conda --cores 1 inventory
```

Use the existing submission profile for substantial jobs. Paths in the config
are relative to this checkout, not the local mounted path. Missing local tools
are not evidence of a server problem. No local environment installation is needed.

| Target | Work and decision |
|---|---|
| `inventory` | Reconcile explicit people, lanes, resources, donor overlap and ONT links. |
| `catalogue` | Correct donor provenance, verify source insertions, reconstruct REF/ALT and independent flanks. |
| `representation` | Search the whole frozen HPRC graph and both linear FASTAs; retain reviewed, unresolved and unavailable states. |
| `panel` | Freeze exactly 20 candidate, 5 represented and 5 technical-control loci, including TMEM106B as represented. |
| `svgp_pilot_mapping` | Map the selected people to HPRC-only and HPRC+local, without requiring calibrated thresholds. |
| `pilot` | Extract molecule evidence and apply available carrier locks; uncalibrated cells stay no-calls. |
| `pilot_compute` | Record actual mapping benchmarks/output sizes; temporary-disk peak still needs server monitoring. |
| `ont_validation` | Whole-genome ONT force calls plus separately adjudicated locus evidence and paired concordance. |
| `cohort` | Reviewed release; one main graph per eligible person, reusing main-graph pilot mappings. |
| `annotation` | Build-aware genes, transcripts, TSS windows, reference repeats and regulatory context. |
| `association` | One prespecified binary contrast and a declared family of Firth carrier tests. |
| `prediction` | Optional nested held-out baseline/plus-SV comparison; disabled by default. |
| `study_status` | Refresh the stage/output inventory without starting missing analyses. |
| `dossiers` | Assemble actual evidence for 2–3 explicitly selected candidate alleles. |

Request one stage at a time and inspect its proposed dependencies first:

```bash
snakemake -s workflow/rules/sv_gene_phenotype.smk -n -p catalogue
snakemake -s workflow/rules/sv_gene_phenotype.smk -n -p representation
snakemake -s workflow/rules/sv_gene_phenotype.smk -n -p svgp_pilot_mapping
snakemake -s workflow/rules/sv_gene_phenotype.smk -n -p cohort
```

Do not use the old broad `review_followup` aggregate for this study. Only its
provenance Python helper is reused; the old mapping/V2 dependency graph is not
included. The corrected tables are written under the new results root.

## Inputs that need actual study decisions

`svgp_samples.example.tsv` is one row per **person**. `svgp_reads.example.tsv`
has explicit library/lane/mate records. The optional crosswalk uses
`namespace,alias,person_id`; it never strips sample suffixes. Additional numeric
covariate/SNP columns in the sample table are preserved. Inventory checks path
availability, not every FASTQ record. Declare manifest completeness only after
checking all expected lanes and native pairing. The old five batch manifests
are starting points, not proof that all approximately 500 people are covered.

Unknown donor/relative/unrelatedness, phenotype or covariate fields exclude a
person from primary association. They may still enter a separately labelled
technical pilot. A declared SNP fingerprint and its evidence are needed before
scoring ONT–Illumina concordance. Neither a shared label nor a force call verifies
identity. Ascertainment and discovery-donor overlap remain study limitations.

Freeze resource release, build, scope and existing verified checksum in the
resource manifest. Copy the verified discovery GFA checksum into the catalogue
config. The mapping GBZ has its **own** checksum; target/lock crosswalks link it
to the catalogue GFA. GRCh38 here is primary chromosomes only. HPRC claims refer
to the configured frozen SV-only graph, not every HPRC assembly or later release.

The catalogue is insertion-first. Unsupported or unrecoverable candidates remain
in the rejection/attempted denominator. Source offsets and orientations are
checked against actual source FASTAs. Real source flanks are used in ALT queries;
nonmatching source/reference flanks need a recorded placement review. Exact
flank sequence alone does not establish unique whole-genome placement.

`representation` saves PAF/GAF with alignment differences and independent-flank
hits. Fill the representation review from these alignments and local graph or
assembly inspection. A definitive decision requires resolved two-sided synteny,
coverage, adequate search and full structural review. A whole-resource repeat
hit is separate from representation at this locus. Unmapped always remains
unresolved until adequate locus evidence exists. Strict additional counts require
all three resources to be `NOT_REPRESENTED_AT_RESOLVED_LOCUS`.

The known-control table and FASTA can supply represented controls outside the
local candidate catalogue. FASTA names are `control_id|REF` and `control_id|ALT`;
both include matching flanks. TMEM queries can instead reuse the existing TMEM
FASTA/target manifest. Resource review must support represented-control labels.
The panel does not estimate genome-wide novelty or enrichment.

See [SVGP_MEASUREMENT_INPUTS.md](SVGP_MEASUREMENT_INPUTS.md) for graph-bundle,
oriented-target, molecule, threshold-lock and ONT schemas. These are substantive
inputs: this workflow does not invent graph-aware PCR duplicate marking or ONT
truth. Genomic competitors, conversion/chopping and node translations must be
documented in target evidence. A per-person map preserves library boundaries;
mates, supplementary observations and duplicate fragments do not become extra
molecules.

For example, the initial main-graph bundle can point at existing index files
after their provenance and compatibility have been checked:

```yaml
# Under sv_gene_phenotype.measurement; add the HPRC-only bundle similarly.
graphs:
  hprc_plus_local:
    gbz: ../whole_pangenome/direct_sv_graph_comparison/indexes/hprc_plus_local/source.giraffe.gbz
    dist: ../whole_pangenome/direct_sv_graph_comparison/indexes/hprc_plus_local/source.dist
    minimizer: ../whole_pangenome/direct_sv_graph_comparison/indexes/hprc_plus_local/source.shortread.withzip.min
    zipcodes: ../whole_pangenome/direct_sv_graph_comparison/indexes/hprc_plus_local/source.shortread.zipcodes
    manifest: config/local_graph_indexes.tsv
    targets: config/local_assay_targets.tsv
    graph_sha256: "FILL_WITH_VERIFIED_GBZ_SHA256"
    node_space_id: "CHOOSE_A_FROZEN_BUNDLE_ID"
    reference_build: CHM13
    reference_paths: []
```

These historical filenames do not certify the bundle. The manifest records the
four exact component checksums and node space. The shared distance index must
already be read-only; mapping does not change its permissions. Header-only
examples must be replaced with actual verified records before mapping.

The implemented calibrated measurement is **carrier status**. Dosage, diploid
genotypes and allele frequency remain missing. Optional audited `vg pack/call`
provides a supporting VCF; it is not promoted to validated diploid calls. The
carrier matrix uses `1`, `0`, `.` for supported ALT, supported REF, no-call.
Allele frequencies or HWE cannot be inferred from those carrier values.

## Annotation and phenotype interpretation

The normalized feature TSV retains source build/release, gene/transcript IDs,
strand, cell type and evidence. Prepare it from the frozen gene, repeat/segdup
and SCREEN resources. Supply explicit splice-region intervals and record their
definition. `assessed_features` lists only feature types whose relevant coverage
is complete. Empty resources mean `NOT_ASSESSED`, never repeat-free/intergenic.
Regulatory target links need functional/interaction evidence; nearest-gene
assignment is descriptive. Justify neuronal/glial cell choices for dementia.

Coordinates are **0-based, half-open**. An insertion has `start=end`; novel
inserted bases do not become a reference interval. BED search proxies are used
only to find neighboring features; exact boundaries are restored in the output.
Cross-build annotation requires both unique flanks, orientation and synteny.
The initial catalogue crosswalk leaves unverified counterpart coordinates
unresolved. Optional RepeatMasker needs a supplied library/release and existing
server environment. Inserted-sequence repeats stay separate from flank repeats.

The configured draft design uses `phenotype=svPPA` as case and `phenotype=control`
as comparator, with `age_at_sampling` as the age field and proposed numeric
baseline adjustment. Source `control` labels have not established cognitively
healthy or centenarian status. `other_dementia` rows remain in the cohort metadata
but do not enter this primary comparison. No other contrast or phenotype recoding
is selected.

Before association can run, the design TSV must have exactly one approved row;
timestamps use
`YYYY-MM-DDTHH:MM:SSZ`. Freeze design and technical rules before phenotype
release. Covariate lists are semicolon separated; numeric values must be actual
numeric columns. The review must discuss ascertainment, age definitions and
overlap, batch/recruitment confounding and covariate choice. The supplied design
retains `DRAFT`, unreviewed age overlap and empty reviewer, evidence and freeze
dates until those decisions are made. Missing ages remain missing; the current
association gate skips inference when any eligible person in the contrast lacks
age or has a mismatched age definition. Adjustment does not fix nonoverlapping
ages or phenotype/batch confounding; the code skips nonestimable designs.

Declare 5–10 phenotype-blind candidate loci in the test manifest. Every allele
needs a primary `carrier` row; represented controls use separate families.
This first association and SV-prediction design accepts autosomal loci only.
Sparse tests retain counts and intervals and stay in the family-wide Holm
denominator. Ten pooled carriers is a stability screen, not guaranteed power;
carriers need not occur in both outcomes. The power grid is an approximate
unadjusted carrier comparison, not covariate-adjusted Firth-model power.
Optional SNP-only/SV-only/joint comparisons require an externally justified SNP
column, all three prespecified rows and identical complete-case people.

Prediction never uses full-cohort association hits as a selection step. Its
SNP baseline must be supplied and described exactly; covariates plus APOE are
not all SNP-based prediction. All learned transforms/tuning occur inside
training folds. Scores describe the sampled contrast, not absolute clinical
disease risk. A null association, no predictive gain or wide uncertainty is an
acceptable outcome.

## Parallel work and resource budget

Independent people/graphs and linear-resource searches run in parallel through
Snakemake. Libraries within one person run sequentially; read lanes stream in
matching mate order without duplicating the cohort. Each expensive job has a
benchmark. This avoids every-person-by-every-graph cohort expansion.

The initial mapping request is 16 threads/80 GB per job, a conservative allocation,
not a measured new-assay cost. On a 336-GB node, start with at most three such
jobs plus overhead; a 1,440-GB node may allow more after the pilot. For example,
after the dry run and profile review:

```bash
snakemake -s workflow/rules/sv_gene_phenotype.smk --use-conda \
  --cores 48 --resources mem_mb=280000 svgp_pilot_mapping
```

Use the existing server submission profile when allocating cluster jobs; this
example alone is not a scheduler submission command. Benchmark actual CPU/wall
time, peak memory, temporary disk and output size before approving cohort
release. No automatic deletion of existing alignments or source data is added.

## Verification and current documentation

Repository instructions forbid local tests, syntax checks and dry runs. All
implementation reviews here were static. Run the supplied regression fixtures
and targeted dry runs **on the server**, then one complete pilot person before
releasing a batch. Check the DAG for unexpected legacy or whole-cohort jobs.
Environment solving and actual tool-flag/index compatibility remain server checks.

Current official releases checked on 20 October 2026: Snakemake 9.27.0,
vg 1.77.0, minimap2 2.31, minigraph 0.21, Sniffles 2.8.1,
SAMtools/BCFtools 1.24, bedtools 2.31.1, RepeatMasker 4.2.4, logistf 1.26.1
and scikit-learn 1.9.1.
Dedicated environments pin the analysis tools; the existing Snakemake installation
is reused without an implicit upgrade. Gfatools has no confirmed current release
pin in this implementation; retain its actual server version with the frozen
provenance outputs. Do not claim that an unexecuted conda solution is validated.

Official methods and documentation:

- [Snakemake rules](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html), [deployment](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html), [release](https://github.com/snakemake/snakemake/releases/tag/v9.27.0).
- [Minigraph bubble/source format](https://github.com/lh3/minigraph#calling-structural-variations), [Minigraph paper](https://doi.org/10.1186/s13059-020-02168-z).
- [Minimap2 assembly alignment](https://github.com/lh3/minimap2#full-genomeassembly-alignment), [paper](https://doi.org/10.1093/bioinformatics/bty191).
- [Giraffe mapping](https://github.com/vgteam/vg/wiki/Mapping-short-reads-with-Giraffe), [graph SV genotyping](https://github.com/vgteam/vg/wiki/SV-Genotyping-and-variant-calling), [genotyping paper](https://doi.org/10.1186/s13059-020-1941-7).
- [Sniffles](https://github.com/fritzsedlazeck/Sniffles), [Sniffles2 paper](https://doi.org/10.1038/s41587-023-02024-y).
- [bedtools closest](https://bedtools.readthedocs.io/en/latest/content/tools/closest.html), [RepeatMasker](https://www.repeatmasker.org/RepeatMasker/), [SCREEN](https://screen.wenglab.org/).
- [logistf documentation](https://search.r-project.org/CRAN/refmans/logistf/html/logistf.html), [Firth inference basis](https://doi.org/10.1002/sim.1047).
- [scikit-learn release](https://scikit-learn.org/stable/whats_new/v1.9.html), [nested evaluation](https://scikit-learn.org/stable/auto_examples/model_selection/plot_nested_cross_validation_iris.html).

Biological claims must follow the evidence: resource-additional, independently
corroborated, callable carrier, gene-context candidate, association and predictive
gain are separate conclusions. No functional or clinical result is implied by
running this workflow.
