# Existing GATK-gCNV inputs

This optional branch imports existing postprocessed **segment VCFs** as linear
copy-number evidence. It does not rerun gCNV or add a discovery analysis. The
files below were inspected read-only on 21 October 2026; source data remain
untouched. Paths used by the workflow must be server paths supplied in config.
The supplied local mount root is
`/home/georgii/snellius/data/dementree_runs/CNV_1_4`; existing run summaries name
the corresponding server root `/gpfs/work3/0/qtholstg/dementree_runs/CNV_1_4`.
Confirm that location in the execution environment before importing.

The supplied run contains two overlapping analyses under `GATK_gCNV/`:

| Run | Segment directories | Existing overview sample count |
| --- | --- | ---: |
| Standard | `GENOTYPED_CALLS_segments_0`, `GENOTYPED_CALLS_segments_1` | 324 |
| WGS autosomes | `wgs_autosomes_all_samples/GENOTYPED_CALLS_segments_0` | 358 |

The counts come from the pre-existing `CNV_analysis/combined_overview.txt`,
not a fresh full scan. It reports 324 shared samples and 34 additional samples
in WGS autosomes. Select one run deliberately; they are not independent people
or independent validation callsets. The standard run's existing QC table flags
cohort 1 as a six-sample cohort. Existing QC tables also flag burden outliers.
Neither these summaries nor the caller's output automatically qualify a sample
or segment for the planned association analysis.

Both inspected segment headers identify GATK
`4.4.0.0-27-gabe8148-SNAPSHOT`, reference
`GRCh38_masked_v2_decoy_excludes_GPRIN2_DUSP22_FANCD2_decoy_HLA_PhiX.fa`, and
`--autosomal-ref-copy-number 2`. The command-line reference lives under
`/gpfs/work3/0/qtholstg/hg38_res_v2/hg38_phix/`.
Chromosome names use `chr`; the inspected chr1 length is 248956422. The custom
masked/decoy reference is not interchangeable with a CHM13 graph coordinate
system. A header path is provenance, not verification of sequence identity.

Sample names are taken from `#CHROM`, for example
`NL_DemenTree_batch1_10027`; filenames also encode cohort and sample index.
Use an explicit reviewed sample-to-person crosswalk. Do not strip batch prefixes
or assume that the suffix identifies the same biological person in the
phenotype workbook. Distinct sequencing identifiers may represent the same
person, including repeat libraries.

The segment VCF schema includes `END`, `GT`, `CN`, `NP`, `QA`, `QS`, `QSS` and
`QSE`. Preserve them and the original `QUAL` and `FILTER`. Coordinates imported
as `start0=POS-1`, `end0=END` describe the reported segment, including gaps between
the retained bins; they are not sequence-resolved breakpoints. `NP` counts bins
or targets and helps describe that limitation.
The [upstream segment VCF composer](https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/copynumber/gcnv/GermlineCNVSegmentVariantComposer.java)
writes segment start and end directly into VCF POS/END. This differs from the
preceding-anchor convention used by ExpansionHunter repeat records.

An observed duplication can have `CN=3` and `GT=./.`: gCNV does not necessarily
encode a conventional diploid genotype for a gain. Preserve the reported copy
number without manufacturing a genotype. The importer also retains explicit
reference-copy segments, clearly labelled; their presence does not establish
that every base between retained bins was callable. Missing VCFs, missing
records and missing CN values never become reference calls. `FILTER=.` means
unfiltered, not `PASS`; raw evidence is not a high-confidence callset.

The old downstream summaries used `FILTER in ./PASS`, `QS>=50`, `NP>=2`, length
at least 1000, excluding chrX/chrY/chrM. These were project-specific analysis
filters. This import does not silently adopt them or relabel them as calibrated
association filters. Gene-level effects and an equivalence to a graph allele
require reviewed coordinates, type, allelic structure and identity; simple
interval overlap does not establish equivalence. No automatic prediction
features or association tests are generated from this evidence.

The import is scheduled as independent files through Snakemake. It reads VCF
text sequentially and uses only Python's standard library for parsing and
provenance bookkeeping, which Snakemake itself cannot supply. Discovery,
parallelism, resources and dependencies remain in the Snakemake rules.

The helper's interface is `svgp_gcnv.py --vcf INPUT --source-id SOURCE
--crosswalk TSV --out-evidence TSV --out-status TSV`, with optional
`--expected-sample` and `--reference-build`. The shared crosswalk columns are
`source_id`, `caller_sample_id`, `person_id`, `review_status`; only the exact
status `REVIEWED` authorizes that explicit identity mapping. Use
`gcnv_standard` or `gcnv_wgs_autosomes` as source IDs. Unresolved identities
remain visible and do not become curated people. Evidence retains all reported
segments with `reported_state` describing CNV, explicit reference, missing CN,
or unsupported allele representation. The per-file status includes compressed
input SHA256, header digest, header version/reference, and record-state counts.

Current upstream checks: [GATK PostprocessGermlineCNVCalls documentation](https://gatk.broadinstitute.org/hc/en-us/articles/54698739237019-PostprocessGermlineCNVCalls)
(updated 20 August 2026) describes complementary interval and segment outputs,
ALT type relative to contig ploidy, and the special handling of declared sex
chromosomes. [Official GATK releases](https://github.com/broadinstitute/gatk/releases/tag/4.7.0.0)
lists 4.7.0.0 as the current release checked on 21 October 2026. The importer
records the version that actually generated each existing file; it does not
upgrade the caller or reinterpret these files as outputs of the current release.
The [GATK-gCNV paper](https://doi.org/10.1038/s41588-023-01449-0) benchmarks rare
CNV discovery from exome data; those performance results do not establish
accuracy for this WGS cohort, custom reference or a particular graph allele.

No local tests, Snakemake runs, dry runs, or syntax checks were performed,
following the repository's `AGENTS.md`. Perform the documented workflow checks
on the server before analysis.
