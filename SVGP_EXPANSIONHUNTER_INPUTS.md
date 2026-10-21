# Existing ExpansionHunter evidence

Use the existing calls as a separate, explicitly linked evidence source for plan
steps 6, 8 and 9. They do not establish that a graph allele has been measured,
and they are not an independent truth set for calls from the same Illumina reads.
No ExpansionHunter rerun, BAM processing or genome-wide expansion search is added.

## Inspected inputs

The user supplied `/home/georgii/snellius/data/dementree_runs/expansionhunter`.
From the server checkout the corresponding starting location is configurable as
`../../dementree_runs/expansionhunter`; confirm it on the server.

- `results/` contains single-sample `.vcf`, `.json` and `_realigned.bam` files.
- `ExpansionHunter.smk` uses `variant_catalog.masked_reference.json` and the
  masked GRCh38 reference
  `/gpfs/work3/0/qtholstg/hg38_res_v2/hg38_phix/GRCh38_masked_v2_decoy_excludes_GPRIN2_DUSP22_FANCD2_decoy_HLA_PhiX.fa`.
  This is a distinct reference resource from the graph's CHM13 coordinate system.
- A representative log records **ExpansionHunter v5.0.0**. Its VCF sample label
  is `NL_DemenTree_batch2_2410006p04.markdup`, while its file stem lacks `.markdup`.
  Neither stem nor label proves a person identity.
- Historical `analysis/qc/run_completeness_summary.tsv` reports 448 expected,
  433 complete JSON/VCF pairs and 15 missing pairs. These are historical August
  counts, not a newly computed complete inventory. A realigned BAM does not prove
  that VCF/JSON generation completed.
- Historical `analysis/run_metadata.tsv` records catalog SHA256
  `2c9dcd7c5a95fd20cdfeb3c6ac46576b9c019bc5b41874b09b1c33717a22aec1`
  and no reference SHA256. The prepared subset records its actual catalog checksum
  when run; a filename or an old checksum alone does not freeze the reference.
- Existing named-locus QC covers 31 loci/38 repeat components. Sample approval
  fields inspected were `not_reviewed`, and historical documentation identifies
  unresolved IDs, a sex discrepancy and discordant putative controls. Existing
  `named_calls.parquet`, association and known-expansion screening outputs remain
  historical results, with their original limitations.

Only directory names, small reports, a representative VCF header/leading records,
the workflow and a log were inspected. Large JSON files, read files and whole
cohort genotypes were not rescanned.

## Small explicit inputs

Use a reviewed shared identity crosswalk, with these tab-separated columns:

```text
source_id	caller_sample_id	person_id	review_status
```

For this source use `source_id=expansionhunter`, the **exact VCF sample column**,
a person ID already reconciled with the phenotype/read inventory, and
`review_status=REVIEWED`. Other states remain `identity_status=UNRESOLVED`.
The importer never strips a batch prefix, `.markdup`, `p01` or any other suffix.
Each external label must occur once for the source; multiple reviewed technical
assays may map to one person and must not become independent observations.

Choose a small repeat subset before examining phenotype results. The selection
TSV has these columns, using **ExpansionHunter catalog** IDs:

```text
locus_id	variant_id
```

These are not the graph's `locus_id` or `allele_id`. A multi-component locus can
have several variant IDs. An empty selection deliberately prevents extraction;
the workflow must not silently select all genome-wide repeats or choose loci by
their existing association P values. An absent catalog component remains
`NOT_IN_CATALOG`.

## Import and parallel execution

Snakemake prepares the selected catalog once and schedules one existing VCF per
job. The small Python helper is needed for catalog/VCF field parsing and exact
identity joins; Snakemake handles files, dependencies and sample parallelism.
The import does not require ExpansionHunter itself in the environment.

The helper's server-side command contract is:

```bash
python workflow/scripts/svgp_expansionhunter.py catalog \
  --catalog ../../dementree_runs/expansionhunter/variant_catalog.masked_reference.json \
  --loci config/svgp_existing_repeat_loci.example.tsv \
  --output ../whole_pangenome/sv_gene_phenotype_2026-10-20/existing_calls/expansionhunter/selected_catalog.tsv

python workflow/scripts/svgp_expansionhunter.py extract \
  --vcf "$EH_VCF" \
  --catalog ../whole_pangenome/sv_gene_phenotype_2026-10-20/existing_calls/expansionhunter/selected_catalog.tsv \
  --crosswalk config/svgp_existing_call_crosswalk.example.tsv \
  --run-log "$EH_LOG" \
  --reference-build GRCh38_custom_masked_decoy \
  --reference-fasta "$EH_REFERENCE" \
  --out-evidence "$EH_EVIDENCE_TSV" \
  --out-status "$EH_STATUS_TSV"
```

Copy and curate the two header-only example manifests, then update their config
paths before import; the commands show the shipped paths for reference. The
variables denote explicitly selected server paths for manual use. The configured
Snakemake rules provide the equivalent VCF, reference and output arguments but
do not supply `--run-log`, so caller version remains `unknown`. Optional manual
`--run-log` use reads at most its first 20 lines for the caller's actual version.
Missing version/build/checksum information stays
`unknown`; `--reference-sha256` accepts a separately established checksum and
does not hash or validate a large reference. Check its provenance before use.
The source VCF is streamed once, keeping only selected records in memory. The
35.7 MB catalog is parsed once, not once per person. JSON sidecars, BAMs and
full genome-wide Parquet tables are unnecessary for this import.

## Output meaning and limits

Each selected repeat has one evidence row per source VCF, retaining exact caller
sample name/person linkage, source file path, catalog checksum, declared reference,
version evidence, motif, reference repeat count, raw VCF coordinates, `GT`,
`REPCN`, `REPCI`, `SO`, `ADSP`, `ADFL`, `ADIR` and locus coverage `LC`.

`allele1_repeat_units` and `allele2_repeat_units` are counts of motif copies on
each reported allele. They are **not** genomic copy number, expansion dosage,
graph SV dosage or phased haplotypes. Preserve caller order: allele 1 need not
be the shorter allele because reference alleles can be listed first. A haploid
call has one allele; expected ploidy and reported sex still need review.

Confidence interval bounds are retained as the caller's integer repeat-unit
endpoints, including both endpoints, with no BED-style conversion. Genomic
intervals are different: catalog `ReferenceRegion` is 0-based, half-open. For
EH repeat records, VCF `POS` is the 1-based preceding anchor, so the repeat
interval is **`[POS, END)`**, with `start0 = POS`, not `POS - 1`. The importer
requires exact chromosome/start/end agreement with the selected catalog.
It does not rename chromosomes or lift coordinates between assemblies.

`CALLER_PASS` means the source caller emitted PASS and complete repeat count/CI
fields. Every row remains `review_status=NOT_VALIDATED`. Other explicit states
include `LOW_DEPTH`, `PARTIAL_GENOTYPE`, `NO_CALL`, `FILTERED_OR_UNASSESSED`,
`CATALOG_COORDINATE_MISMATCH`, `RECORD_NOT_REPORTED`, `NOT_IN_CATALOG`,
`VCF_UNAVAILABLE` and header-only/missing-header states. Raw partial estimates
and LowDepth counts are retained for review. Missing, filtered or absent records
are never assigned a normal/reference allele.

The source header defines LowDepth as overall depth below 10x **or** fewer than
five reads spanning one/both breakends. LC alone cannot remove that filter.
PASS does not establish repeat structure, precision for long alleles, pathogenic
classification, resolved interruptions or calibrated carrier status. Ordinary
`Repeat` calls near fragment length can underestimate longer alleles; current
minimal import does not read JSON fragment/read-length metrics and therefore
does not clear this limitation. Use existing locus QC and selected read review
before interpreting candidates.

A VCF reaching EOF does not prove a completed caller run; the status table marks
this explicitly. Missing expected selected records remain in denominators, and
technical repeats must be reconciled before a per-person association matrix is
made. Incomplete identity/reference/measurement QC must be resolved before
adding a prespecified repeat feature to a conditional model or prediction
baseline. No feature is added automatically.

Repeat-locus overlap with a graph allele supports contextual review only.
Exact allele equivalence requires the same person, resolved build/locus, motif,
repeat configuration and compatible measurement evidence. Neither nearby repeat
length nor missing EH output proves whether a local graph insertion is present.

## Documentation checked 21 October 2026

- [Official releases](https://github.com/Illumina/ExpansionHunter/releases):
  standalone latest listed release v5.0.0, matching the inspected log.
- [Official catalog format](https://github.com/Illumina/ExpansionHunter/blob/master/docs/04_VariantCatalogFiles.md):
  component IDs, 0-based intervals, Repeat/RareRepeat and fragment-length limits.
- [Official VCF format](https://github.com/Illumina/ExpansionHunter/blob/master/docs/06_OutputVcfFiles.md):
  repeat-unit and read-support meaning. Its example uses older CN/CI names;
  the actual v5.0.0 headers are authoritative for REPCN/REPCI/ADSP/ADFL/ADIR/LC.
- [Versioned v5.0.0 VCF writer](https://github.com/Illumina/ExpansionHunter/blob/v5.0.0/ehunter/io/VcfWriter.cpp)
  and [allele encoder](https://github.com/Illumina/ExpansionHunter/blob/v5.0.0/ehunter/io/VcfWriterHelpers.cpp):
  preceding-anchor coordinates and reference-first allele ordering.
- [Official JSON format](https://github.com/Illumina/ExpansionHunter/blob/master/docs/05_OutputJsonFiles.md):
  available fragment length, read length, expected allele count and locus depth.
- [Dolzhenko et al., 2019](https://doi.org/10.1093/bioinformatics/btz431):
  the sequence-graph repeat-genotyping method.

No local test, syntax check, validation command or workflow run was performed,
following this repository's remote-only execution instructions. The commands
above are for the server and preserve existing analysis outputs.
