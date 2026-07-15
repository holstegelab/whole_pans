# whole-pans

Snakemake workflow for human haplotype assembly QC, non-human sequence
decontamination, and QC of the cleaned assemblies. Original assemblies are
read-only. The default target runs all three stages and produces a final,
post-decontamination assembly list for graph construction.

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
inclusion list. Missing CHM13/hg38 alignments needed for decontamination are
generated under `results.qc`; the original QC summaries remain unchanged.

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

Compleasm needs its lineage downloaded once from a node with internet access:

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
its large temporary file tree; set `TMPDIR` through the execution environment
or cluster profile when node-local scratch is available.

Useful targets are:

```bash
snakemake --cores 32 --use-conda qc
snakemake --cores 32 --use-conda decontamination
snakemake --cores 32 --use-conda post_decontamination_qc
snakemake --cores 1 --use-conda pangenome_qc
snakemake --cores 32 --use-conda all
```

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

- `assembly_qc.tsv` and `sample_qc.tsv`: QC metrics recalculated from cleaned
  FASTAs;
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
