# Corrected four-graph comparison

Run these commands on Snellius, from `/gpfs/work3/0/qtholstg/pangenome/whole_pans`.
Use an allocated compute environment with the indicated CPU/memory resources.
No local tests, workflow runs or dry runs were performed while editing, in accordance
with this repository's `AGENTS.md`.

## Existing data and new outputs

The workflow continues to use `personalized_sv_comparison.results` for existing
GBZ/index files, sample metadata, bubble inventories, compressed read GAFs and
per-sample mapping metrics. Corrected derived results go into its `analysis_v2/`
subdirectory. An optional `personalized_sv_comparison.analysis_results` config key
can select a different directory, for example for a second MAPQ threshold.

Existing mapping results are sufficient for this rerun. The workflow regenerates
the small REF/ALT query alignments, exports segment-to-node translations from the
existing GBZs, recounts targeted evidence from the saved read GAFs, normalizes the
catalogue and produces the paired comparison and report. It does not need to remap
the sample FASTQs or rebuild indexes when those existing products are complete.

The evidence stage streams each callable sample/graph GAF once. Across 80 large
GAFs this is substantial I/O, even though only locus-associated read names are kept
in memory. Limit concurrent jobs according to the storage and compute allocation.
Graph targets that cannot be assessed produce blank evidence counts without a
whole-genome scan.

## Server checks and rerun

Run the focused synthetic checks with the server Python environment:

```bash
python -m unittest discover -s tests -p 'test_personalized_sv_comparison.py'
```

Inspect the complete proposed rerun:

```bash
snakemake -s workflow/rules/personalized_sv_comparison.smk \
  personalized_sv_comparison --use-conda --cores 8 \
  --resources mem_mb=32000 --rerun-triggers mtime --dry-run
```

`--rerun-triggers mtime` is useful for this migration: it prevents unrelated
recorded code/parameter metadata changes from forcing old mapping jobs. New
analysis outputs have new names and depend explicitly on their analysis scripts.
Check that the proposed jobs do not include `build_source_giraffe_index` or
`map_sample_to_source_graph`. If they do, resolve the missing/stale upstream
products before starting; do not force a new index over alignments from an old one.

First generate just the corrected graph targets:

```bash
snakemake -s workflow/rules/personalized_sv_comparison.smk \
  summarize_tmem_graph_targets --use-conda --cores 8 \
  --resources mem_mb=32000 --rerun-triggers mtime
```

Inspect `analysis_v2/tmem_graph_targets.tsv`. A usable graph should have
`callable_complete_graph_specific_target`, `edge_namespace=vg_node_ids`, and
nonempty left/right ALT and REF-defining edges. A failed audit is retained in the
report; it is not interpreted as allele absence. If translation fails, check the
source GFA/GBZ provenance before proceeding. Do not strip segment prefixes or
guess numeric IDs.

Then finish the evidence comparison and report:

```bash
snakemake -s workflow/rules/personalized_sv_comparison.smk \
  personalized_sv_comparison --use-conda --cores 8 \
  --resources mem_mb=32000 --rerun-triggers mtime
```

Primary outputs under `analysis_v2/`:

- `report.md`: mapping effects, graph callability, prior-group evidence, paired
  TMEM agreement, and per-sample junction counts.
- `tmem_graph_targets.tsv`: REF/ALT alignment qualification and translated edges.
- `tmem_19kb_evidence.tsv` and `sv_evidence.tsv`: the same targeted TMEM evidence;
  these are not a genome-wide genotyped SV callset.
- `mapping_metrics.tsv` and `paired_comparison.tsv`: reused mapping measurements,
  paired effects in original units and percentage points, direction counts,
  read-count changes, and BH-adjusted p-values for overall comparisons.
- `sv_catalogs.tsv`: normalized chromosome keys, graph membership for exact
  coordinate/type matches, and explicit `allele_equivalence=not_assessed`.

## What was corrected

Minigraph's default GAF can encode a reference traversal as `CHM13.chr7` and ALT
traversals as stable named intervals. Such coordinates cannot be compared directly
with numeric nodes in Giraffe GAF. The new target mapping requests `minigraph --vc`
and exports the actual mapping GBZ translation with
`vg gbwt -Z --translation FILE GRAPH.gbz`. Junctions are then translated using the
first/last indexed nodes of each oriented source segment. Node chopping and
reverse orientations are handled without guessing node IDs. Traversal strings in
the target table retain source segment coordinates; edge/node fields use vg IDs.

These options are documented in the
[Minigraph manual](https://lh3.github.io/minigraph/minigraph.html) and
[vg 1.76.1 GBWT implementation](https://github.com/vgteam/vg/blob/v1.76.1/src/subcommand/gbwt_main.cpp).

Flank qualification now counts bases aligned to graph sequence from the CIGAR,
excluding inserted/clipped bases. Missing MAPQ (255) is excluded. The prior-group
check also ignores samples outside the selected pilot; its old AWK lookup created
extra rows with empty groups.

Unavailable allele evidence is blank, not zero. Allele counts collapse `/1` and
`/2` mate names into fragments. Historical `*_reads` column names are retained,
with `counting_unit=read_name_collapsed_fragments` making their meaning explicit.
PCR duplicates are not removed. Fragments supporting incompatible REF and ALT
edges are excluded from allele counts; sufficiently many flag `CONFLICTING`.
Only ALT-specific internal nodes count as insertion-sequence evidence.

`CONFLICTING` and `NO_DECISION` are excluded from usable-decision denominators.
An entirely unresolved comparison reports N/A and has no TMEM discrimination
ranking. Positive prior groups and prior absence groups are not validated truth
genotypes. ALT support and REF support are evidence classes, not formal diploid
genotypes. The report does not estimate sensitivity, specificity or disease risk.

Catalogue keys normalize the configured reference's `CHM13.chr7` and
`CHM13#0#chr7` aliases while preserving the original names. Matching coordinates
and SV type do not establish identical allele sequences; nonmatching coordinates
can still describe equivalent variants. Novelty requires additional allele-level
comparison. Genome-wide mapping gains still need validation against independent
truth data before they can be called accuracy gains.

The existing MAPQ/junction thresholds are retained and included in the evidence
table and report. To assess threshold sensitivity, use a separate
`analysis_results` directory and repeat with `evidence_thresholds.min_mapq: 20`;
do not mix tables generated with different thresholds in one comparison.

The mapping-statistics implementation now uses a small MAPQ histogram instead
of keeping one MAPQ value per mapped read. Existing mapping metric files are
reused; the new missing-MAPQ handling applies when that rule is run again.
