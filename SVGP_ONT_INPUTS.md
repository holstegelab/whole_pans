# Reported ONT samples and runs

The user supplied these assignments on 21 October 2026. The source is
`config/svgp_ont_runs.tsv`; one row means one reported sample/run assignment.

| Sample label | Reported ONT runs |
| --- | --- |
| 2310032p01 | ONT00096, ONT00098 |
| 2410018p03 | ONT00145 |
| 2410020p01 | ONT00144 |
| 2410024p01 | ONT00144 |
| 2410024p02 | ONT00145 |
| 2410026p01 | ONT00136 |
| 2410037p01 | ONT00096 |

There are seven sample labels, eight assignments and five distinct run labels.
Runs can contain several samples, and one sample can appear in several runs.
These counts do not establish the number of people or independent validations.
Keep the sample labels exactly as supplied, including the `p01` and `p02` suffixes.

Configure the manifest using `sv_gene_phenotype.inventory.ont_runs`. The
`reported_ont` target writes `metadata/reported_ont_runs.tsv` below the configured
results directory. The default inventory also includes this output. The rule
only reads the small assignment manifest and records its path, checksum and
reported provenance; it does not search or scan sequencing files. It rejects
duplicate sample/run rows and malformed or incomplete manifest rows.

The output marks availability `USER_REPORTED` and retains `person_id` and
`read_path` as `unknown`, with identity
`UNRESOLVED`, reads `NOT_PROVIDED`, fingerprints `NOT_ASSESSED` and measurement
`NOT_READY`. It is separate from the curated `ont_pairing.tsv` inventory. The
reported assignments do not release samples for ONT genotyping or concordance.

To use these runs for the planned ONT corroboration, first supply the explicit
sample read or alignment paths and their reference/build or read preset, connect
them to curated people using an explicit crosswalk, and review ONT–Illumina SNP
fingerprints and discovery-donor overlap. Demultiplexed files and libraries must
be reviewed before combining the two runs reported for `2310032p01`. Do not infer
per-sample file paths from a shared run label or mark absent files as reference
calls. Existing `measurement.ont` settings and adjudication gates remain in force.

From the server `whole_pans` checkout, inspect the proposed inventory job with:

```bash
snakemake -s workflow/rules/sv_gene_phenotype.smk reported_ont --cores 1 -n
```

After reviewing the dry run, the same command without `-n` writes only this
reported assignment inventory. No command has been run locally, in accordance
with the repository instructions.
