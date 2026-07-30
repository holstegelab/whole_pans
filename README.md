# whole-pans

Snakemake workflow for human haplotype assembly QC, non-human sequence
decontamination, and QC of the cleaned assemblies. Original assemblies are
read-only. The default target runs decontamination plus the lightweight
post-decontamination sequence QC and produces the final assembly list for graph
construction. The full original-assembly QC remains available as target `qc`.

## Layout

```text
whole_pans/
├── Snakefile                         # pipeline entry point
├── config/config.yaml                # paths, thresholds, and resources
├── environment.yaml                 # Snakemake driver environment
└── workflow/
    ├── rules/
    │   ├── QC.smk
    │   ├── decontamination.smk
    │   └── post_decontamination_QC.smk
    ├── envs/
    │   ├── tools.yaml                # shared QC/decontamination tools
    │   └── compleasm.yaml
    └── scripts/
```

Runtime state, results, logs, and Python caches are excluded by `.gitignore`.

## Configuration

Edit `config/config.yaml` before running on another server. In particular,
check:

- `assemblies.directory` and `assemblies.patterns`;
- both reference FASTA paths under `references`;
- the Kraken2 database path and target taxids under `kraken`;
- all three result directories and per-rule resource requests.

Decontamination processes every assembly discovered from `assemblies.directory`
and `assemblies.patterns`, including assemblies that failed the original QC.
The second QC pass evaluates all cleaned assemblies and writes the final graph
inclusion list. It reuses the cleaned SeqKit statistics already produced by the
decontamination stage and applies only the thresholds that can materially change
after removing contigs: total length, contig N50, and N content. It does not run
Compleasm or new CHM13/hg38 alignments for the cleaned FASTAs.

In the current server configuration, workflow code remains in `whole_pans/`
while all result trees are under `whole_pangenome/`.

Keep the decontamination and post-decontamination QC directories separate so
the second QC pass cannot overwrite the original QC results.

The QC summarizer currently expects the reference keys `CHM13` and `hg38`.
Assembly names must follow the existing hifiasm convention, for example:

```text
SAMPLE.hifi.hifiasm.bp.hap1.p_ctg.fa.gz
SAMPLE.hifi.hifiasm.bp.hap2.p_ctg.fa.gz
```

When exactly one haplotype fails QC, the sample is marked `PARTIAL`: the
failing haplotype is excluded while its PASS/WARN mate remains eligible. Both
are excluded when both fail or when the pair is missing or duplicates a
haplotype. Change `workflow/scripts/summarize_qc.py` if the naming scheme
differs.

Decontamination uses Kraken2 with the same `pluspf_20230605` database as
`/home/georgii/AlzHub/short_read_analyzing_pipeline_Snakemake`. The database is
passed as a runtime parameter because it is only visible on the server. Only
contigs classified under the configured bacterial and viral taxids are removal
candidates; target calls with substantial human-reference support are retained
for review instead of being removed automatically.

The Kraken report parser accepts both the standard six-column report and the
eight-column report produced by `--report-minimizer-data`. Cleaned FASTA
identifiers are written as `<assembly_id>.<original_contig>` so every rGFA
source name is unique across haplotypes. The original FASTAs, removed FASTAs,
and review FASTAs retain their original identifiers; `split_maps/*.tsv` records
the original-to-cleaned name mapping.

## Run

Create the small driver environment once:

```bash
conda env create -f environment.yaml
conda activate whole-pans
```

Inspect the QC jobs and create the isolated rule environments:

```bash
snakemake -n --cores 1 qc
snakemake --cores 1 --use-conda --conda-create-envs-only
```

When running the optional original `qc` target, Compleasm needs its lineage
downloaded once from a node with internet access:

```bash
snakemake --cores 1 --use-conda compleasm_download
```

Run the complete pipeline:

```bash
snakemake --cores 32 \
  --resources mem_mb=120000 \
  --use-conda \
  --rerun-incomplete \
  --printshellcmds
```

The workflow contains no scheduler submission settings. Add the Snakemake
profile or executor used by the cluster. Compleasm uses `${TMPDIR:-/tmp}` for
its large temporary file tree when it runs on original assemblies; set `TMPDIR`
through the execution environment or cluster profile when node-local scratch
is available. Compleasm is part of the original QC target only.

Useful targets are:

```bash
snakemake --cores 32 --use-conda qc
snakemake --cores 32 --use-conda decontamination
snakemake --cores 32 --use-conda post_decontamination_qc
snakemake --cores 1 --use-conda pangenome_qc
snakemake --cores 32 --use-conda all
```

#### Future extension: map complete TMEM haplotype paths

A complementary, more comprehensive screen is to map every regional TMEM
haplotype sequence to the whole-genome graph. This tests complete observed
haplotypes and can identify complex or balanced graph paths that are difficult
to express as one REF/ALT length change in the regional VCF. Use this path-level
screen to discover candidate differences, then use the allele-centered workflow
above to confirm and deduplicate them.

Current-input observations that the future script must preserve:

- `TMEM_UPD/tmem_out/seqfile.tsv` has 984 entries in the master input list;
- the final TMEM subproblem seqfile has 380 entries: the hg38 reference plus
  379 non-reference haplotypes from 288 samples;
- the regional GFA contains 596 `W` records and no `P` records because some
  haplotypes are represented by multiple walk fragments.

Prefer the 379 non-reference FASTAs selected in this subproblem seqfile as the
mapping queries:

`TMEM_UPD/tmem_out/chrom-subproblems/seqfiles/chr7_11946976-12488798_sub_150000_391823.seqfile`

A preparation script should read this seqfile, skip the hg38 reference when
appropriate, discard empty inputs, and write one multi-FASTA with headers that
retain the sample, haplotype, original contig, and source interval. Using the
original subproblem FASTAs avoids turning one haplotype into several independent
queries at GFA walk breaks.

If the embedded graph walks themselves are required, they can instead be
extracted with ODGI. Run these commands from `whole_pans/`:

```bash
gzip -dc ../TMEM_UPD/tmem_out/pangenome_cactus_pipeline.gfa.gz \
  > /tmp/tmem_pangenome.gfa

odgi build -g /tmp/tmem_pangenome.gfa \
  -o /tmp/tmem_pangenome.og -t 16

odgi paths -i /tmp/tmem_pangenome.og -L \
  > tmem_graph_path_names.txt

odgi paths -i /tmp/tmem_pangenome.og -f \
  > tmem_graph_paths.fa
```

Map all regional sequences in one invocation so that the large whole-genome
graph is loaded only once:

```bash
minigraph -cxasm -l5000 -t 16 \
  ../whole_pangenome/sv_pangenome/results/graphs/sv_pangenome.minigraph.gfa \
  tmem_regional_haplotypes.fa \
  > tmem_haplotypes.whole_graph.gaf
```

The future GAF summarization script should:

1. select the best primary alignment while retaining secondary alignments for
   ambiguity checks;
2. require high query coverage and two-sided anchors around each internal
   difference;
3. parse `cg:Z` and flag internal insertions or deletions of at least 50 bp;
4. flag split alignments, inversions, unexpected orientation, and inconsistent
   graph traversal separately;
5. ignore terminal clipping caused by regional extraction boundaries unless it
   has independent two-sided support;
6. project flanking `SR:i:0` graph nodes through their `SN` and `SO` tags to
   CHM13 coordinates, while retaining node-plus-offset graph coordinates;
7. cluster repeated events across carrier haplotypes by CHM13 anchor pair, SV
   type and length, and inserted-sequence similarity;
8. match event clusters back to the regional VCF where possible and submit each
   candidate to the existing allele-centered confirmation workflow.

An end-to-end path alignment with no residual insertion/deletion of at least
50 bp supports representability of that regional haplotype in the frozen whole
graph. A residual event is a candidate missing allele, not an automatic absence
call: repetitive, low-MAPQ, split, and boundary alignments must remain
`uncertain` until local allele remapping confirms them.

Suggested outputs are one row per path alignment, one row per raw residual
event, one row per clustered candidate SV, and a summary report linking GRCh38,
CHM13, and graph node/offset coordinates.

Each stage file is also a standalone entry point. From the project root
(`/gpfs/work3/0/qtholstg/pangenome`), run:

```bash
# Original-assembly QC
snakemake --profile ~/.config/snakemake/zslurm/ \
  --snakefile whole_pans/workflow/rules/QC.smk \
  --use-conda --rerun-incomplete

# Decontamination of every discovered assembly
snakemake --profile ~/.config/snakemake/zslurm/ \
  --snakefile whole_pans/workflow/rules/decontamination.smk \
  --use-conda --rerun-incomplete

# Decontamination followed by QC of cleaned assemblies
snakemake --profile ~/.config/snakemake/zslurm/ \
  --snakefile whole_pans/workflow/rules/post_decontamination_QC.smk \
  --use-conda --rerun-incomplete
```

No target name is required: each file defines its own stage as the default
target. Add `--dry-run` to inspect the DAG before submission.

When expanding an existing QC-filtered decontamination run to all assemblies,
run the normal decontamination command above without `--forcerun`. Existing
Kraken, decision, cleaned FASTA, and stats outputs are reused; missing assembly
outputs are added to the new all-assembly manifest and generated automatically.
The first expanded run is still large because every newly added assembly needs
Kraken2 and CHM13/hg38 alignments.

After changing the Kraken report parser or cleaned FASTA naming, rebuild the
affected results in this order. The first command reuses the existing Kraken2
classification and report files; it does not force `kraken_contigs`.

```bash
# 1. Reclassify from existing Kraken reports and rewrite the cleaned FASTAs.
snakemake --profile ~/.config/snakemake/zslurm/ \
  --snakefile whole_pans/workflow/rules/decontamination.smk \
  --use-conda --rerun-incomplete \
  --forcerun classify_contigs clean_assembly cleaned_stats summarize_decontamination

# 2. Recalculate QC from the corrected cleaned FASTAs.
snakemake --profile ~/.config/snakemake/zslurm/ \
  --snakefile whole_pans/workflow/rules/post_decontamination_QC.smk \
  --use-conda --rerun-incomplete

# 3. Rebuild the graph with unique donor and reference sequence names.
snakemake --profile ~/.config/snakemake/zslurm/ \
  --snakefile whole_pans/workflow/rules/pangenome.smk \
  --use-conda --rerun-incomplete --force

# 4. Recreate the graph QC report.
snakemake --profile ~/.config/snakemake/zslurm/ \
  --snakefile whole_pans/workflow/rules/pangenome_qc.smk \
  --use-conda --rerun-incomplete
```

Run each command only after the previous command has completed successfully.
The graph builder now stops instead of publishing a GFA if Minigraph reports
inconsistent rGFA names or names associated with multiple source ranks.

`decontamination` does not filter its inputs using the first QC result.
`post_decontamination_qc` automatically runs missing decontamination work and
then applies the QC thresholds to all cleaned assemblies. This allows an
assembly rejected for excess non-human sequence or low query-aligned percentage
to become eligible after cleaning.

## Main outputs

QC outputs under the configured `results.qc/summary` directory:

- `assembly_qc.tsv`: metrics and PASS/WARN/FAIL reasons per assembly;
- `sample_qc.tsv`: paired-haplotype decisions, including `PARTIAL` samples;
- `graph_included_assemblies.txt`: assemblies accepted by the original QC;
- `graph_excluded_assemblies.tsv`: excluded assemblies and reasons.

Decontamination outputs under the configured `results.decontamination`
directory:

- `summary/contamination_summary.tsv`: contamination counts and bases;
- `summary/contig_actions.tsv`: REMOVE and REVIEW decisions; SPLIT fields are
  retained for compatibility and remain empty for Kraken-based decontamination;
- `summary/review_candidates.tsv`: ambiguous retained contigs;
- `summary/graph_cleaned_assemblies.txt`: cleaned FASTAs for graph construction;
- `cleaned/*.clean.fa.gz`: cleaned assemblies;
- `removed/*.nonhuman.fa.gz`: removed sequence;
- `review/*.review.fa.gz`: retained ambiguous contigs;
- `split_maps/*.split_map.tsv`: original-to-cleaned coordinates.

Post-decontamination QC outputs under the configured
`results.post_decontamination_qc/summary` directory:

- `assembly_qc.tsv` and `sample_qc.tsv`: cleaned-assembly sequence metrics and
  PASS/WARN/FAIL decisions based on total length, contig N50, and N content;
- `graph_included_assemblies.txt`: final cleaned FASTAs accepted for graph
  construction;
- `graph_excluded_assemblies.tsv`: cleaned FASTAs that fail the second QC pass.

Use the post-decontamination `graph_included_assemblies.txt` for graph
construction. The similarly named decontamination list contains every cleaned
assembly before the second QC pass.

SV pangenome outputs under the configured `pangenome.results` directory:

- `graphs/sv_pangenome.minigraph.gfa`: SV-level Minigraph graph;
- `metadata/sv_pangenome.ordered_assemblies.tsv`: CHM13, hg38, then cleaned
  assemblies ordered by Mash distance to CHM13;
- `metadata/chm13.mash_distances.tsv`: Mash distances used for ordering;
- `metadata/sv_pangenome.graph_summary.tsv`: minimal construction-time GFA
  summary.

### Screen cleaned assemblies for graph-missing SV candidates

`workflow/scripts/screen_novel_graph_svs.py` maps complete cleaned assemblies
back to the frozen Minigraph rGFA. It parses internal insertion/deletion
operations, retains 30--49 bp sensitivity evidence and nearby indel clusters,
and requires mapping quality, alignment length, and two-sided aligned anchors
before calling an event a high-confidence graph-missing SV. Residual events are candidates, not
confirmed variants; repeats, fragmented contigs, complex rearrangements, and
assembly errors still require local validation. A contig with multiple primary
alignments is reported as `REVIEW_SPLIT_OR_COMPLEX_ALIGNMENT` instead of being
declared free of candidates; balanced inversions and translocations require
follow-up breakpoint analysis beyond the insertion/deletion CIGAR screen.

The complete assembly-only first pass is implemented in
`workflow/rules/novel_sv_discovery.smk`. It freezes and inventories the current
graph, refreshes the excluded-assembly feasibility data, screens every cleaned
assembly, calls paired samples against CHM13 and hg38, builds a provisional
assembly-SV catalog, and emits both provisional and evidence-ranked
sample-level HiFi request lists. Unpaired assemblies remain in the graph
screen but are skipped by the diploid reference callers. From the directory
containing `whole_pans/`, run the environment solve and DAG checks before
submitting production jobs:

```bash
snakemake --snakefile whole_pans/workflow/rules/novel_sv_discovery.smk \
  --use-conda --conda-create-envs-only --cores 1 novel_sv_discovery

snakemake --snakefile whole_pans/workflow/rules/novel_sv_discovery.smk \
  --use-conda --dry-run novel_sv_discovery

snakemake --profile ~/.config/snakemake/zslurm/ \
  --snakefile whole_pans/workflow/rules/novel_sv_discovery.smk \
  --use-conda --rerun-incomplete novel_sv_discovery
```

All runtime paths are under `novel_sv_discovery` in `config/config.yaml` and
are relative to `whole_pans/Snakefile`. Do not replace them with `realpath`
values from a mounted checkout; the compute-server mount can differ.

SVIM-asm and dipcall are enabled by default. SVIM-asm runs once per phased pair
in diploid mode. Each input reference is copied to a writable uncompressed
`reference.fa`, indexed once as `.fai`, `.dict`, and `.mmi`, and the `.mmi` is
shared by minimap2 and dipcall. By default, known X/Y calls are excluded because
sample sex, paternal-haplotype origin, and reference-specific PAR intervals are
not yet available; raw X/Y caller VCF records are retained outside the merged
catalog. Unprojected graph rows are retained but marked
`UNRESOLVED_GRAPH_ORIGIN`, so they must not be assumed to be autosomal. Enable
`sex_aware` mode only after those metadata and PAR BED files
are populated and PAR sequence on reference chrY has been verified as
hard-masked, as described in `NOVEL_SV_DISCOVERY_WORKFLOW_GUIDE.md`. In that
mode, X/Y records from non-PAR-aware SVIM-asm and PAV remain excluded.

Runtime dependencies are separated by purpose:

- `sv_graph.yaml`: Minigraph/gfatools graph screening;
- `sv_reference.yaml`: reference preparation, SVIM-asm, and dipcall;
- `sv_catalog.yaml`: manifest, catalog, and ranking scripts;
- `sv_discovery.yaml`: isolated native dependencies for the optional PAV
  source workflow.

PAV is an upstream source Snakemake workflow rather than a Bioconda
executable. For the full three-caller production union, install a tagged,
recursive PAV checkout on the server, point
`reference_calling.pav_snakefile` to its `Snakefile`, and set
`reference_calling.callers.pav: true`.

The screening utility remains usable without Snakemake. First create the
reusable graph-coordinate index and deterministic task table:

```bash
python workflow/scripts/screen_novel_graph_svs.py index-graph \
  --gfa ../whole_pangenome/sv_pangenome/results/graphs/sv_pangenome.minigraph.gfa \
  --output ../whole_pangenome/sv_pangenome/novel_sv_discovery/frozen_graph/graph_segments.tsv.gz

python workflow/scripts/screen_novel_graph_svs.py tasks \
  --manifest ../whole_pangenome/sv_pangenome/novel_sv_discovery/manifest/assembly_discovery_manifest.tsv \
  --graph-assemblies ../whole_pangenome/sv_pangenome/results/metadata/sv_pangenome.ordered_assemblies.tsv \
  --batch-size 2 \
  --output ../whole_pangenome/sv_pangenome/novel_sv_discovery/graph_screen/tasks.tsv
```

Task IDs come directly from `tasks.tsv`; do not recalculate a task count from
the manifest. Each job atomically publishes one declared directory containing
all outputs for its batch. For a manual pilot, choose a task ID and write to a
new pilot directory:

```bash
TASK_ID=0001

python workflow/scripts/screen_novel_graph_svs.py run \
  --graph ../whole_pangenome/sv_pangenome/results/graphs/sv_pangenome.minigraph.gfa \
  --segment-index ../whole_pangenome/sv_pangenome/novel_sv_discovery/frozen_graph/graph_segments.tsv.gz \
  --manifest ../whole_pangenome/sv_pangenome/novel_sv_discovery/manifest/assembly_discovery_manifest.tsv \
  --graph-assemblies ../whole_pangenome/sv_pangenome/results/metadata/sv_pangenome.ordered_assemblies.tsv \
  --tasks ../whole_pangenome/sv_pangenome/novel_sv_discovery/graph_screen/tasks.tsv \
  --output-dir ../whole_pangenome/sv_pangenome/novel_sv_discovery/graph_screen/pilot_tasks/"$TASK_ID" \
  --task-id "$TASK_ID" \
  --threads 16
```

The output directory must not already exist. A successful task contains a
non-empty `.complete`, `task_outputs.tsv`, and all per-assembly outputs; failed
temporary directories are retained with a `.failed` suffix for diagnosis.

Before launching all tasks, pilot 5--10 graph members, 10--20 missing-mate or
`best_rescue` assemblies, and several fragmented assemblies by requesting their
task IDs from `tasks.tsv`. Add leave-one-out or synthetic positive controls for
SV-class/size recall. Use measured memory and runtime to update the configured
resources and batch size. The 3.5 GB graph is loaded independently by each
Minigraph process, so avoid excessive concurrency on a shared filesystem.

After all tasks finish, let the checkpoint pass the exact task-directory set to
the streamed aggregator. It rejects missing tasks and missing or empty
per-assembly outputs:

```bash
snakemake --snakefile workflow/rules/novel_sv_discovery.smk \
  --use-conda --rerun-incomplete --cores 32 \
  summarize_novel_sv_graph_screen
```

The main first-pass outputs are:

- `manifest/assembly_discovery_manifest.tsv`: all 982 haplotypes with graph
  membership, refreshed QC/rescue status, contamination, metadata, and read
  access fields;
- `graph_screen/summary/all_assembly_novel_sv_summary.tsv`: callable fractions
  and per-assembly decisions;
- `graph_screen/summary/all_residual_sv_candidates.tsv.gz`: raw, review, and
  high-confidence graph-residual evidence;
- `graph_screen/coordinate_qc/all_residual_sv_candidates.annotated.tsv.gz`:
  stable-coordinate repair of saved Minigraph residuals, accompanied by a
  fail-fast `coordinate_qc.tsv` report;
- `reference_calls/reference_call_manifest.tsv`: traceable PAV, SVIM-asm, and
  dipcall outputs for paired samples against CHM13 and hg38, according to
  enabled callers;
- `catalog/provisional_graph_residual_sv_catalog.tsv`: nonmember graph-residual
  clusters only; compatible CHM13 calls are validation evidence, while
  reference-only and graph-member-only clusters do not create catalog rows;
- `hifi/recommended_hifi_samples.tsv`: one ranked row per recommended sample,
  including poorly callable samples with zero assembly candidates.

The catalog merger streams VCFs into a temporary SQLite database instead of
loading the cohort union into RAM. Put `TMPDIR` on sufficiently large
node-local scratch for this job. The read-based discovery, read genotyping,
validated graph-v2 rebuild, and capture re-screen remain deferred second-pass
work: they require a populated read manifest or a reviewed graph-v2 input list and
must not be inferred while those inputs are unavailable. `minigraph --call`
only genotypes bubbles already represented in the graph; it is not used as the
missing-allele discovery method.

The separate `pangenome_qc` target analyzes an existing or newly built GFA:

```bash
snakemake --cores 1 --use-conda pangenome_qc
```

Graph QC follows `../pangenome_qc_plan.md` and writes a separate results tree
at `pangenome.qc_results`:

- `data/seg_records.parquet`: streamed rGFA `S`-record checkpoint with segment
  ID, `SR`, `LN`, `SN`, normalized chromosome, `SO`, and `is_ref`;
- `data/rank_tally.tsv`: per-rank segment count, segment bp, and link count;
- `tables/integrity_checks.csv`: parser-vs-summary, rank, link-integrity,
  tag-completeness, input-QC, and Mash-order checks;
- `tables/graph_overview_stats.csv`: topology, content, segment-length, and
  degree metrics;
- `tables/per_source_contribution.csv`: build-order non-reference sequence
  contribution and cumulative growth;
- `tables/per_chromosome_nonref.csv`: non-reference bp and bp/Mb by chromosome;
- `tables/mash_outliers.csv`: Mash-distance and low matching-hash outliers with
  assembly-QC/decontamination context;
- `figures/*.png`: segment-length, growth-curve, chromosome, and Mash-distance
  plots;
- `report/pangenome_qc_report.md`: written synthesis with links to tables and
  figures;
- `metadata/qc_tool_versions.tsv` and `metadata/run_manifest.tsv`: QC tool
  versions and output manifest.

Reference alignment is intended for gross outlier and contamination checks.
Real human structural variation and reference-specific sequence should not be
interpreted as assembly error without supporting read evidence.
