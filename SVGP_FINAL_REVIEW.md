# Final independent review

Static review against `SV_gene_phenotype_coding_plan_2026-10-20.md`, after
implementation, independent simplification and the subsequent error pass.
Prepared 20 October 2026. No tests, syntax checks, workflow commands, dry runs,
environment solves or genomic analyses were executed locally, as required by
`AGENTS.md`. Fixtures are supplied for server execution and remain unexecuted.

The implemented sequence follows the plan: explicit people/resources; corrected
source alleles; whole-resource searches with unresolved states; the 20/5/5-locus
panel; paired bounded Illumina pilot; independently adjudicated ONT; a reviewed
single-graph cohort; build-aware gene context; one prespecified carrier contrast;
optional nested prediction. Inventory remains the default and requires no
mapping. No additional biological method or study hypothesis was introduced.

Concrete corrections in this pass:

- Bound the cohort-release review to the actual pilot and ONT evidence inputs.
  Relative path spellings are normalized lexically, without `realpath`, symlink
  resolution or local/server path translation.
- Required an existing nonempty SNP-fingerprint evidence file before confirming
  ONT–Illumina identity. Source-anchor, definitive resource and cross-build
  coordinate reviews also require existing evidence files.
- Collected both FASTQ-stream exit statuses before accepting successful mapping,
  retaining a failed decompression as a failed job.
- Enforced the plan's autosome-only scope for primary association and SV
  prediction. A separate sex-chromosome design is not implemented implicitly.
- Added small server-only fixtures for missing fingerprint/resource evidence
  and release-input binding, including equivalent relative path spellings.

The changed branches, argument names, output fields and fixtures were reread
after editing. Carrier observations retain missing diploid dosage and allele
frequency; unresolved resource/coordinate states remain explicit. Source-GFA
and mapped-GBZ checksums retain their distinct roles. Sparse association tests
stay in their declared Holm denominator. Prediction preprocessing and tuning
remain inside training folds.

Optimization choices favor readability. Independent people/graphs and resource
searches remain Snakemake jobs; libraries within a person remain sequential.
The main cohort reuses pilot mappings and never expands to all graphs. FASTQs
stream without cohort copies; molecule/depth unions avoid counting duplicates
twice. Index checksums are audited once per graph bundle. Additional caching,
read recruitment, assignment-table sharding and multi-step job grouping were
not added: they require measured server bottlenecks and complicate provenance.
Use the pilot's actual benchmarks before changing concurrency or memory.

Official tool pins and references from the earlier documentation review were
retained. The [current Snakemake deployment documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html)
was checked again: conda settings on `run` rules apply to their `shell()` calls.
No global Snakemake upgrade is implied; actual server flag/index compatibility
and availability of the pinned conda builds remain unverified.

Scientific inputs remain substantive prerequisites: source/resource inspection,
graph/node translations, duplicate-resolved molecules, calibrated locks,
fingerprints, ONT adjudication and a usable phenotype design. Existing-file
checks do not judge scientific adequacy or replace that review. No diploid
genotype accuracy, independent ONT ALT validation for every candidate, causal
regulatory effect or predictive benefit is asserted. Candidate dossiers and
the flow graphic summarize available evidence; publication-ready locus diagrams
still require the verified study coordinates/features and scientific review.
Prediction reports explicitly omit full retraining uncertainty and a separate
batch-held-out sensitivity. `SVGP_VERIFICATION.md` supplies server checks before
any cohort release.
