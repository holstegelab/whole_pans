# Personalized short-read comparison of four SV pangenome graphs

The implemented direct-mapping comparison now writes corrected downstream results
to `direct_sv_graph_comparison/analysis_v2/` while reusing the existing mapping
GAFs and indexes. See [the corrected comparison rerun guide](PERSONALIZED_SV_COMPARISON_RERUN.md)
for the target-coordinate fixes, interpretation changes and server commands.

## 1. Objective

Test whether personalized versions of population pangenome graphs improve
short-read SV representation and mapping for the DemenTree population.

This is a deliberately small first experiment. One Snakemake invocation should
perform the complete analysis:

```text
prepare source graphs and metadata
        -> personalize each graph for each sample
        -> build short-read indexes
        -> map real paired-end reads
        -> calculate coverage and SV evidence
        -> compare the four graph targets
```

Do not add a separate discovery workflow, synthetic read benchmark,
universal-versus-personalized comparison, or elaborate orchestration.

The experiment is SV-only. The graph targets are the following four source
graphs, each converted into a personalized per-sample reference before read
mapping:

| Target | Source graph | Purpose |
|---|---|---|
| `hprc` | HPRC CHM13 SV graph only | HPRC baseline |
| `local_v1` | Local SV graph V1 | Same-population graph |
| `local_v2` | Local rescued SV graph V2 | Updated local graph |
| `hprc_plus_local` | HPRC graph plus local assemblies | Combined graph |

The comparison answers: which source graph produces the best personalized
reference for these samples?

## 2. Inputs

### Real short reads

Read all five manifests:

```text
../DemenTree/NL_DemenTree_batch1.tsv
../DemenTree/NL_DemenTree_batch2.tsv
../DemenTree/NL_DemenTree_batch3.tsv
../DemenTree/NL_DemenTree_batch4.tsv
../DemenTree/NL_DemenTree_batch5.tsv
```

The currently observed files do not contain a header. They contain one sample
per row with comma-separated R1 and R2 FASTQ paths. Define the column names in
one place and validate the format before mapping.

Create one manifest with at least:

```text
sample_id  batch  read_type  sex  r1_paths  r2_paths
```

The five files contain approximately 90 rows each. Derive the final sample
count from the manifests, deduplicate only exact duplicate sample IDs, and
write the resulting manifest as an output.

### Initial 20-sample pilot

Run the first test on the following 20 short-read samples. These are test
samples from the same population; they must be checked against every
source-graph membership manifest and confirmed to be absent from all four
graphs before mapping.

#### High probability of the large insertion

These are the five highest-ranked `very_high` candidates in
`TMEM_UPD/shortread_19kb/results/summary/likely_19kb_insertion_candidates.tsv`.
They have strong support at both breakpoints and are included in the
edge-k-mer validation as `kmer_supported_alt_edges`.

```text
NL_DemenTree_batch4_2420004p01
NL_DemenTree_batch3_2310038p03
NL_DemenTree_batch4_2410005p01
NL_DemenTree_batch1_2310028p02
NL_DemenTree_batch1_2410062p02
```

#### Lower prior probability of the large insertion

These are the five lowest-ranked `moderate_candidate` samples in the TMEM
candidate table. They have weaker, but still non-zero, two-breakpoint support
and are lower-probability insertion candidates. They are not confirmed
genotypes.

```text
NL_DemenTree_batch3_2410039p04
NL_DemenTree_batch3_2410028p05
NL_DemenTree_batch4_2410079p01
NL_DemenTree_batch4_2410050p01
NL_DemenTree_batch4_2310010p03
```

#### High probability of insertion absence

These five samples are controls in the independent TMEM edge-k-mer validation
table. All have `no_alt_edge_kmer_signal`, the strongest available short-read
evidence for absence in this dataset. Do not use the broad insertion-node
genotype table for this category because it calls all 358 tested samples
insertion-present.

```text
NL_DemenTree_batch2_2410080p01
NL_DemenTree_batch2_2310037p02
NL_DemenTree_batch1_2011-054
NL_DemenTree_batch3_2410027p01
NL_DemenTree_batch3_2310024p02
```

#### Not tested on the TMEM local graph

The existing TMEM pilot manifest covers batches 1–4. These five batch-5
samples were therefore selected as samples outside that pilot. Confirm this
from the runtime TMEM manifest/result table; do not infer insertion absence
from this group.

```text
NL_DemenTree_batch5_1997-089A
NL_DemenTree_batch5_1998-037A
NL_DemenTree_batch5_2000-063A
NL_DemenTree_batch5_2002-075
NL_DemenTree_batch5_2003-027A
```

Set these explicit IDs in `sample_groups` for the first run. Do not select this
pilot randomly. The workflow must verify that every ID exists in the five
manifests, has paired reads, and is absent from all four graph membership
manifests. If a check fails, stop and replace the sample explicitly; do not
silently run an imbalanced pilot.

The TMEM results contain 358 tested samples from batches 1–4. The broad
`19kb_insertion_genotypes.tsv` result labels all of them insertion-present, so
it is not suitable for defining absence. The groups above use the stricter
breakpoint candidate table and independent edge-k-mer validation.

Do not assume that `/home/gozhegov/...` and `/gpfs/work3/...` paths are both
visible on the execution server. Keep optional path-prefix translation in the
config, or provide an execution-side manifest with valid paths.

### Source graph paths

Use the configured paths from `config/config.yaml`. The expected graph roles
are:

```yaml
personalized_sv_comparison:
  source_graphs:
    hprc: /gpfs/work3/0/qtholstg/pangenome/HPRC/graphs/hprc-v2.0-mc-chm13.sv.gfa.gz
    local_v1: ../whole_pangenome/sv_pangenome/results/graphs/sv_pangenome.minigraph.gfa
    local_v2: ../whole_pangenome/sv_pangenome/local_sv_graph_v2/graphs/sv_pangenome.local.v2.rescued.minigraph.gfa
    hprc_plus_local: ../whole_pangenome/hprc_plus_cohort/graphs/hprc-v2.0-mc-chm13.cohort.sv.gfa
```

The exact paths may be changed on the execution server. Fail clearly if a
source graph is missing; do not silently substitute a different graph.

The source graphs are made by:

```text
whole_pans/workflow/rules/pangenome.smk
whole_pans/workflow/rules/local_sv_graph_v2.smk
whole_pans/workflow/rules/hprc_plus_cohort.smk
```

For this experiment, these workflows provide the four source graph files. The
comparison workflow consumes them, records provenance, and creates personalized
derivatives. Rebuilding the population graphs inside the short-read comparison
is out of scope.

### TMEM insertion inputs

Use the existing TMEM outputs as a positive-control source:

```text
TMEM_UPD/tmem_out/pangenome_cactus_pipeline.vcf.gz
TMEM_UPD/shortread_19kb/results/graph/target_insertions.tsv
TMEM_UPD/shortread_19kb/results/summary/likely_19kb_insertion_candidates.tsv
TMEM_UPD/shortread_19kb/results/kmer_19kb_validation/allele_specific_kmer_analysis.tsv
```

Reuse existing sample-level TMEM evidence where available, but perform the new
comparison using the same four personalized whole-genome targets and the same
read-evidence thresholds.

## 3. Implemented direct-mapping mode

The local Minigraph rGFAs contain graph topology but no biological `P` or `W`
haplotype paths, and no path-preserving GBZs are available. The implemented
workflow therefore maps each sample directly to each complete source graph.
`vg autoindex --workflow sr-giraffe` builds one reusable Giraffe index per
source graph, using a greedy path cover where biological paths are absent. New
outputs are written to `whole_pangenome/direct_sv_graph_comparison/`; the old
incomplete personalization outputs are left untouched.

Mapping writes compressed GAF directly, without a whole-genome GAM
intermediate. The GAF is then surjected to the CHM13 backbone for BAM-based
coverage statistics. Shared distance indexes are kept read-only while Giraffe
jobs run concurrently.

This is a direct four-graph mapping comparison, not sample-specific haplotype
sampling. Path-cover threads are indexing aids and must not be interpreted as
observed haplotypes. The personalized design below can be restored later if
path-preserving source GBZs become available.

### Original personalized design

For every sample and each of the four source graphs, make a small subgraph
containing haplotypes selected using that sample’s own real-read k-mer counts.
Map the sample only to this personalized subgraph.

Use the method described by Sirén et al., *Personalized pangenome references*,
Nat Methods 2024, DOI `10.1038/s41592-024-02407-2`:

1. preprocess the source graph into blocks and haplotype/minimizer metadata;
2. count graph-informative k-mers in the sample’s reads;
3. classify their presence using estimated read k-mer coverage;
4. select a small number of compatible haplotypes per block;
5. construct a subgraph from selected haplotypes and required backbone;
6. create a Giraffe short-read index for that subgraph.

Preprocessing is done once per source graph. Personalization is done once per
sample/source-graph pair. The universal source graphs themselves are not
mapping targets in this experiment.

The personalization process cannot recover an SV absent from its source graph.
Therefore compare the four source graph inventories, but obtain all read mapping
and coverage results from the four personalized targets.

## 4. Leakage rule

Record graph membership for every source graph. If a sample contributed one or
both haplotypes to a source graph, label it as a graph member.

Use the 20 pilot samples above only if they are held out from all four source
graphs. If a selected sample is already present in a graph, stop the pilot and
replace it with an unrepresented sample. Do not treat a graph-member result as
independent evidence that the graph represents that sample.

The combined `hprc_plus_local` graph necessarily contains local assemblies.
This remains a valid requested target, but its membership and possible leakage
must be reported explicitly.

## 5. One-pass workflow design

Create one Snakefile, for example:

```text
whole_pans/workflow/rules/personalized_sv_comparison.smk
```

Include it from the existing `whole_pans/Snakefile`, or make it runnable
directly with a dedicated config. One normal Snakemake invocation should build
the final report and trigger all required upstream jobs.

The workflow should contain only the following practical stages.

### Stage A: prepare manifests and graph metadata

Rules:

```text
prepare_dementree_manifest
inventory_source_graphs
prepare_tmem_targets
```

`prepare_dementree_manifest`:

- read the five TSVs;
- add `batch` from the filename;
- validate R1/R2 pairing and sample IDs;
- check that input paths are non-empty and readable at runtime;
- write one deterministic sample manifest.

`inventory_source_graphs`:

- check all four source graph files;
- calculate checksums;
- record graph reference/build, source path, file size, and graph type;
- record source graph membership manifests where available;
- count segments, links, paths, and SV bubbles if supported.

`prepare_tmem_targets`:

- normalize the TMEM VCF and `target_insertions.tsv`;
- retain approximately 19-kb insertion candidates and flanks;
- record regional-to-genome coordinate conversion and reference build;
- create one stable target table used by all four comparisons.

### Stage B: personalize each graph

Rules:

```text
preprocess_personalization_graph  # 4 jobs, one per source graph
personalize_graph                 # 4 x number_of_samples jobs
```

`preprocess_personalization_graph` creates reusable metadata for one source
graph: graph blocks, haplotype paths, graph-informative k-mers/minimizers, and
the source graph inventory.

`personalize_graph` takes one sample’s R1/R2 reads and one preprocessed source
graph. It writes:

```text
personalized/{target}/{sample}.gfa
personalized/{target}/{sample}.selection.tsv
personalized/{target}/{sample}.kmer_stats.tsv
```

The selection table records selected haplotype/path IDs, block IDs, estimated
coverage, k-mer classification status, and tool parameters.

Suggested initial parameters should be configurable and identical for all four
targets:

```yaml
personalized_sv_comparison:
  kmer_size: 29
  minimizer_window: 11
  block_size_bp: 10000
  haplotypes_per_block: 8
  diploid_sampling: true
  min_kmer_count: 2
```

Use a fixed number of selected haplotypes so personalized graph sizes are
comparable. Flag samples with unusable k-mer coverage rather than silently
using a different personalization mode.

Pseudocode:

```text
for target in [hprc, local_v1, local_v2, hprc_plus_local]:
    preprocess(source_graph[target])

    for sample in samples:
        counts = count_relevant_kmers(sample.R1, sample.R2)
        coverage, status = estimate_coverage(counts)
        assert status in ["usable", "flagged"]

        selected = select_haplotypes(
            preprocessed_graph[target], counts, coverage,
            haplotypes_per_block=8,
            diploid=True,
        )
        write_subgraph(selected, include_backbone=True,
                       output=personalized[target][sample])
```

### Stage C: index and map

Rules:

```text
build_personalized_giraffe_index  # 4 x number_of_samples jobs
map_sample_to_personalized_graph  # 4 x number_of_samples jobs
```

For each personalized GFA, build a GBZ, distance index, minimizer index, and
zipcode index using one pinned `vg` version. Validate the index immediately in
the same rule or in a directly dependent rule.

Map paired-end reads with `vg giraffe`. Use identical mapper settings, thread
count, MAPQ handling, and output format for all four targets. The exact command
must match the installed `vg` version; the interface is:

```text
vg giraffe -p -t THREADS -b default -o GAF \
  -Z TARGET.giraffe.gbz \
  -m TARGET.shortread.withzip.min \
  -z TARGET.shortread.zipcodes \
  -f SAMPLE_R1.fastq.gz -f SAMPLE_R2.fastq.gz \
  > SAMPLE.TARGET.gaf
```

Each mapping job writes its own log and benchmark file. Do not concatenate all
samples into one mapping job: one sample/target is the comparison unit and
independent jobs allow the four targets to run in parallel.

### Stage D: analyze and compare

Rules:

```text
calculate_mapping_metrics
calculate_sv_evidence
summarize_19kb_insertion
compare_personalized_graphs
```

`calculate_mapping_metrics` summarizes each sample/target mapping. Use a
common linear backbone or common callable regions for the main coverage
comparison. Do not sum graph node lengths to create a cross-graph coverage
value.

`calculate_sv_evidence` evaluates support for SV alleles represented in each
personalized graph. Use graph paths/bubbles and GAF alignments to count:

- reference-path support;
- alternate-path support;
- split/junction support;
- supporting read pairs;
- MAPQ-qualified support;
- ambiguous or conflicting evidence.

`compare_personalized_graphs` joins results by sample ID and SV ID. It emits
paired differences with `hprc` as the baseline.

The final report rule depends on all four graph targets, all selected samples,
and all analysis tables. Running the final report target once runs the complete
workflow.

## 6. Metrics

### Mapping and coverage

For every sample and target report:

- input read and read-pair count;
- mapped-read fraction;
- MAPQ 0, MAPQ >=5, MAPQ >=20, and MAPQ >=40 fractions;
- mean/median MAPQ;
- mean alignment identity;
- fraction of read bases aligned;
- clipped bases per read;
- secondary/multimapping fraction;
- proper-pair or fragment-learning status;
- common-backbone breadth and depth summary;
- wall time, CPU time, peak memory, and output size.

The primary coverage endpoints are common-backbone mapped fraction and
MAPQ-qualified mapped fraction. Report graph-specific node/path coverage only
as a secondary diagnostic.

### SV representation

Before mapping, create one normalized SV catalog per source graph. For each
graph report:

- number of SV alleles by type and size bin;
- number of insertion, deletion, inversion, duplication, and complex alleles;
- number of alleles shared between the four source graphs;
- local-only and HPRC-only alleles;
- number of graph paths/haplotypes supporting each allele.

This describes what each source graph can represent. Read-level results still use
only personalized graph targets.

### Sample SV evidence

For each sample/target/SV record raw evidence counts and a simple call class:

```text
ALT_SUPPORTED
REF_SUPPORTED
HET_LIKE
NO_DECISION
CONFLICTING
```

Use identical thresholds for all four targets. `NO_DECISION` is not an SV
absence call; it can reflect low depth, repeats, missing unique anchors, or an
underspecified graph path.

## 7. Approximately 19-kb TMEM insertion test

This is the main targeted positive-control analysis.

1. Read the existing TMEM VCF and target table.
2. Identify all approximately 19-kb insertion candidates.
3. Normalize each candidate to contig, breakpoint, insertion length, alternate
   sequence, and left/right flanks.
4. Match each candidate against the four source graph SV catalogs using
   breakpoint and sequence similarity.
5. Map candidate carriers, known non-carriers, and available blinded samples to
   each sample-specific personalized target.
6. Count reads supporting both insertion junctions and reads covering the
   inserted sequence.
7. Compare carrier/non-carrier discrimination across the four targets.

Junction support is primary. Coverage over the inserted sequence alone is not
enough because reads may map to homologous sequence elsewhere.

| Result | Interpretation |
|---|---|
| All four support the insertion | HPRC is sufficient for this target; compare ambiguity and coverage |
| Local v1/v2 support it, HPRC does not | Evidence that local assemblies add a population-relevant allele |
| HPRC+local supports it, HPRC does not | Local alleles rescue the combined graph |
| No personalized target supports it | Investigate coordinate conversion, graph projection, anchors, and read depth |
| One target has excess non-carrier support | Possible graph ambiguity or false-positive mapping; penalize that target |

Use the existing TMEM coordinate offset only after checking it against the VCF
header/reference. Store converted whole-genome coordinates and source regional
coordinates in the target table.

## 8. Comparison and decision criteria

Use the sample as the paired unit. For every endpoint compare:

```text
local_v1        - hprc
local_v2        - hprc
hprc_plus_local - hprc
local_v2        - local_v1
```

Primary endpoints:

1. common-backbone MAPQ-qualified mapped fraction;
2. fraction of callable SVs with a confident allele decision;
3. number of supported non-reference SV alleles per sample;
4. correct support for the approximately 19-kb insertion.

Secondary endpoints:

- ambiguous/secondary mappings;
- false ALT support in known non-carriers;
- personalized graph size;
- index size and construction time;
- mapping time and peak memory;
- fraction of samples with unusable personalization.

Summarize by batch and read-depth bin. Do not pool reads across the cohort. A
simple paired bootstrap or paired Wilcoxon/permutation test is sufficient for
the first report.

The preferred graph should have better or equal SV evidence and coverage,
without a material increase in ambiguous mapping, and with acceptable total
runtime. A graph with more catalogued SVs but no additional confident read
support should be reported as having greater representation, not as a proven
mapping improvement.

## 9. Minimal output layout

```text
results/personalized_sv_comparison/
  samples.tsv
  graph_inventory.tsv
  sample_graph_membership.tsv
  tmem_targets.tsv
  preprocessed/{hprc,local_v1,local_v2,hprc_plus_local}/...
  personalized/{target}/{sample}.gfa
  personalized/{target}/{sample}.selection.tsv
  indexes/{target}/{sample}/...
  mapping/{target}/{sample}.gaf.gz
  mapping_metrics.tsv
  sv_catalogs.tsv
  sv_evidence.tsv
  tmem_19kb_evidence.tsv
  paired_comparison.tsv
  report.md
```

Every table should include, where relevant:

```text
sample_id  batch  target  source_graph  graph_checksum  mapper_version
```

The final report should directly answer:

- which personalized target has the best overall coverage;
- which target gives the most confident SV evidence;
- which target detects the approximately 19-kb insertion most reliably;
- whether local v1, local v2, or HPRC+local adds value over HPRC;
- the runtime and storage cost of personalization and mapping.

## 10. Implementation checklist

- Add one config section for the four source graph paths and parameters.
- Add one end-to-end comparison Snakefile/rule file with a final report target.
- Reuse existing graph indexing and SV evidence helpers where they support
  personalized targets.
- Do not add universal graph mapping targets.
- Do not add synthetic reads.
- Do not add long-read validation or a second discovery workflow.
- Keep graph preparation, mapping, and analysis as ordinary Snakemake rules
  with explicit inputs, outputs, logs, threads, memory, and runtime.
- Run the 20-sample pilot first using the same final workflow, then change only
  the sample-selection config to run all five batches.
- Verify the pilot with a Snakemake dry run before launching all samples.
