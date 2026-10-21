# Independent simplification and error review

Static review against `SV_gene_phenotype_coding_plan_2026-10-20.md`, 20 October
2026. Repository instructions prohibit local tests, syntax checks, dry runs,
environment solving and workflow execution; none were performed.

The workflow keeps Snakemake scheduling, per-person/per-graph parallelism,
native mapping/interval commands and the existing provenance helper. Separate
Python/R scripts handle only data joins, graph-edit evidence or statistical
models. No new biological analysis was added.

Corrections made:

- Replaced unchecked decompression process substitutions with checked streams.
  Shared distance indices must already be read-only; inputs are not changed.
- Made raw FASTQs tracked dependencies after the native inventory checkpoint
  resolves explicit identities. Empty pilot/cohort lists fail when requested,
  without blocking inventory; full pilot targets enforce 12–20 people/two graphs.
- Separated technical measurement from independent-primary phenotype
  eligibility. Evidence retains analysis groups and frequencies have separate
  denominators for those groups.
- Merged overlapping assay intervals consistently in target lengths and depth.
- Required confirmed ONT fingerprints and positive, finite resolved evidence.
  Cohort release now connects ONT evidence and callable locked alleles on the
  same main graph.
- Required passing technical locks for association and one declared candidate
  multiplicity family. Declared generated count tables as Snakemake outputs.
- Removed the duplicate measurement-default configuration and unused gfatools
  environment dependency. Kept one main configuration and the existing helper.

Integration inspected: annotation consumes passing source-GFA-linked locks;
prediction consumes main-graph carrier evidence and keeps fitting within
training splits; reports retain candidate-record/allele/locus denominators and
separate donor overlap. The current official Snakemake R script preamble uses
the bundled `methods` package, not `jsonlite`, so no unused dependency was added.

Remaining execution gates are intentional: populate and freeze identities,
resources/checksums, exact control/panel sequences, source/resource reviews,
graph/node-space targets, duplicate-resolved molecules, measurement calibration,
ONT adjudication and phenotype design. No absence, genotype, association or
prediction result is asserted by preparation code. Server verification is
still required, including a single-person run before the full bounded pilot.
