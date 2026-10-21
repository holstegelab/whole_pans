# Clinical workbook input

The supplied `Sequencing_analysisfile_version20260619.xlsx` contains worksheet
`Blad1` with 397 data rows. Read-only inspection found these source fields:

| Field | Use |
| --- | --- |
| `Sample_ID_WGS_given` | Unique workbook-row key in all 397 rows |
| `SequencingID_final` | Retained source alias; repeated values prevent a unique join |
| `Sample_ID_Biobank_raw` | Retained source alias; repeated values prevent a unique join |
| `Diagnosis_svPPA` | Source phenotype labels: `control` 175, `other_dementia` 164, `svPPA` 58 |
| `Diagnosis_complete` | Retained detailed diagnosis with 34 source categories |
| `sex` | Source values `F`, `M`, `NA` |
| `age` | Raw value retained as `source_age`; numeric `age_at_sampling` for 380 rows, 16 `NA` and one blank remain missing |
| `DNA_from` | Retained tissue source, including blood, brain and buccal material |

The workbook is not a 397-person cohort manifest. `Sample_ID_Biobank_raw` has
386 distinct values, with conflicting clinical fields in its duplicate groups.
Do not collapse those rows or infer people by removing sample suffixes.
`SequencingID_final` has 360 distinct values and also contains conflicting
duplicate groups. The unique WGS source key identifies a workbook row, not a
verified person or FASTQ identity.

Copy `config/svgp_phenotype_crosswalk.example.tsv` to a curated TSV and fill an
explicit `Sample_ID_WGS_given` to `person_id` mapping after checking the source
records. Each ID and person must occur once; unresolved rows may remain unmapped.
The workflow does not create people from workbook labels. If multiple samples
belong to one person, choose one reviewed clinical row for that person in this
crosswalk. Keep the other source rows visible in the join report.

The importer records both diagnosis columns, tissue source, row number and the
small workbook's SHA-256. It exports every workbook row in `phenotypes.tsv` and
reports unmapped rows in `phenotype_join_report.tsv`. With curated samples
provided, `samples_with_phenotypes.tsv` retains every curated person. The importer
fills unknown sex, phenotype, phenotype ID and age at sampling only for explicit
crosswalk matches.
Conflicting known sex or phenotype values stop the import. An existing phenotype
ID can remain in its own namespace because the crosswalk records the linkage.
Sex labels are compared literally: a curated `female` and source `F` require
explicit harmonization of the curated metadata before import. This avoids
silently mixing equivalent but differently spelled model-factor levels.

`Diagnosis_svPPA` is imported as its three original categories for availability
summaries. The user selected **svPPA versus control** as the primary contrast.
`config/svgp_association_design.tsv` therefore specifies `svPPA` as case and
`control` as comparator using the imported `phenotype` column. `other_dementia`
remains in the metadata and is excluded only from this primary contrast. Source
`control` labels do not establish cognitively healthy or centenarian status.
`Diagnosis_complete` includes heterogeneous diagnoses and clinical notes; no
diagnostic recoding is automatic.

The user confirmed that `age` means **age at sampling**, also described as current
age in this new cohort. The importer preserves `source_age`, creates the numeric
`age_at_sampling` field and records `age_definition=age_at_sampling` only for
explicitly matched people. Unknown ages remain missing; conflicting known ages
or definitions require review. The draft design proposes age at sampling as a
baseline numeric adjustment. It remains `DRAFT`: age overlap, ascertainment,
other covariates, frozen design dates and reviewer evidence still require study
review. The current conservative association gate skips inference if an eligible
person in the contrast has missing age or a mismatched age definition. Age
adjustment does not resolve nonoverlapping age distributions.

The workbook does not establish ancestry, sequencing batch, relatedness,
discovery-donor status or analysis eligibility. Those fields still require the
curated cohort metadata specified in the plan. Prediction remains disabled.

This workbook is only read. Neither the author workbook nor existing manifests
are overwritten. The XLSX reader requires `openpyxl` in the Snakemake rule's conda
environment. No workflow, importer, tests or validation commands were run locally,
following `AGENTS.md`; use the server verification instructions.
