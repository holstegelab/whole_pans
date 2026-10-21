# Supplied-input review

Static review on 21 October 2026 covered the clinical workbook import, existing
gCNV segment import, selected ExpansionHunter repeat import and their Snakemake
dependencies. No local tests, syntax checks, dry runs, environment solves or
workflow executions were performed, as required by `AGENTS.md`.

Concrete corrections made:

- All three helpers now reject output paths that alias an input or another
  output. The small shared guard compares normalized path spellings and existing
  filesystem identities without `realpath`. Imported shared-helper files are
  explicit Snakemake inputs, so changing that code invalidates downstream jobs.
- Existing-call source folders are tracked checkpoint parameters. Changing roots
  remains visible even when both the old and new folders contain no VCF files.
  The checkpoint inventories headers; aggregation then discovers independent
  per-file imports through its native checkpoint dependency. No upstream cycle
  or caller/mapping dependency was introduced.
- VCF imports/header inventory reject malformed sample-column headers. EH strips
  CRLF line endings before exact identity comparisons. Aggregation rejects bad
  table headers, ragged rows, schema changes and missing/extra per-file status rows.
- Missing new optional config sections no longer cause existing-call rule
  definitions to fail while being read. Default shipped config still enables
  workbook availability import; disabling `phenotypes.enabled` removes that
  workbook branch from `inventory` dependencies.
- Removed an inaccurate EH source-VCF checksum claim. The evidence retains its
  source path and selected catalog checksum. A single streaming VCF pass is kept.
- Retained literal clinical sex comparisons and documented that `female` versus
  `F` requires explicit curated harmonization. No silent sex/diagnosis/age
  recoding or new model field was added.

Plan fidelity and scientific limits were retained: workbook rows do not create
people; joins require explicit source-to-person mappings; age retains an
unspecified time point; diagnosis retains three source labels. gCNV overlapping
runs require one explicit selection and retain CN even when GT is missing.
FILTER `.` remains unassessed. Only reported reference segments are labelled
reference. EH requires a small prespecified repeat list and preserves count/CI,
LowDepth, missing and mismatched-coordinate states. Repeat counts do not become
copy number or graph dosage. Neither caller becomes an independent truth set or
an automatic association feature.

The official EH v5 catalog specification was checked for default component IDs
and 0-based half-open catalog regions. The upstream GATK segment VCF composer
was checked for direct segment-start/end output, consistent with the import's
`POS-1, END` BED conversion. Existing source calls retain their actual old caller
provenance; no caller update or reinterpretation was performed.

The new `tests/test_svgp_supplied_inputs.py` contains small **unexecuted** fixtures
for explicit workbook joins, duplicate people, unknown age meaning, sex conflicts,
input preservation, gCNV gain/reference/missing-CN states, selected catalog IDs,
EH anchor coordinates, partial uncertainty, LowDepth and exact `.markdup` sample
labels with CRLF input. These fixtures need `openpyxl` and belong in the server
environment used for phenotype import. They are mechanical regression checks,
not biological or caller calibration.

## Fresh final simplicity and optimization review

A separate reviewer reread only this new input integration after the corrections.
No further code change was necessary: VCF imports remain independent native
Snakemake jobs, the EH catalog subset is prepared once, and aggregation streams
tables without loading a cohort matrix. The extra gCNV checksum pass was retained
for clear provenance rather than complicating its parser. No new tool, phenotype
definition, identity alias, graph-equivalence rule or association was introduced.

The server verification command now discovers both fixture suites in the
phenotype environment with `openpyxl`. The EH manual examples now match shipped
manifest/output paths and distinguish optional manual log parsing from the
configured rules, which leave unestablished caller version unknown. These were
documentation corrections. No local execution or validation was performed.

## Confirmed phenotype and reported ONT follow-up — 21 October 2026

A separate reviewer inspected the later clarification changes. The user has now
confirmed **age at sampling** (current age in this new cohort) and **svPPA versus
control**. This supersedes the earlier statements above that the age time point
and primary contrast were unknown. The supplied-input fixtures now cover sampling
age rather than unknown age meaning.

The workbook importer retains raw `source_age` and exports numeric
`age_at_sampling` without filling missing values. Only explicitly matched curated
people receive the age field/definition; existing equal ages are retained, and
conflicting ages or definitions stop the import. Source diagnosis labels are
unchanged. The draft association design matches the existing R interface:
`svPPA` is the case, `control` is the comparator and age at sampling is the
proposed numeric adjustment. `other_dementia` remains in metadata and is filtered
out only for that comparison. `DRAFT`, unreviewed age overlap and absent study
review/freeze evidence continue to prevent inference; the user clarification
does not establish control health status or solve age/batch confounding.

The ONT manifest preserves all seven exact sample labels, eight sample/run
assignments and five run labels, including both runs for `2310032p01`. Its small
native Snakemake rule records reported availability and source provenance only.
Person IDs, read paths, identity and fingerprint status remain unresolved; the
default inventory does not start ONT measurement, alter the pilot selection or
mark reported samples as independent validation. Manifest parsing and duplicate
checks are direct and require no new helper or framework.

No additional implementation change was needed for simplicity or plan fidelity.
This was static source/interface review only: no tests, syntax checks, Snakemake
commands, dry runs, environment changes or workflow execution were performed.
The updated age fixtures and reported-ONT dry run remain server checks described
in `SVGP_VERIFICATION.md`.

## Fresh final optimization review of the clarification follow-up

A new reviewer reread the clarified phenotype import, its fixtures, the draft
association design, the native reported-ONT rule and their workflow interfaces.
No further code change was warranted. Sampling age is retained through the
inventory into the existing R design interface; missing ages and the draft study
status retain their existing gates. The seven exact ONT labels and eight run
assignments are handled by one small manifest rule with no additional parser or
mapping job. Existing per-file caller imports remain independently schedulable,
and this small metadata branch adds no useful parallelism to extract. Simplicity
was preferred over speculative optimization. This review introduced no new tool,
identity inference, sample selection, age imputation or study approval, and no
local tests, syntax checks, dry runs or workflow execution were performed.
