# Read-to-graph benchmark

This standalone workflow benchmarks long- and short-read mapping without changing the graph-construction workflows. It treats all graph and read inputs as read-only and writes under the configurable `results` directory in `config/read_graph_benchmark.yaml`.

The five whole-genome long-read targets are:

1. plain T2T-CHM13 v2.0;
2. the HPRC v2.0 CHM13 SV graph;
3. local cohort SV graph V1;
4. local cohort SV graph V2 with rescued SVs;
5. the HPRC CHM13 SV graph augmented with the 196 local V1 haplotypes.

The regional TMEM GBZ is registered separately as `comparison_scope=tmem_pilot`. It is not a sixth whole-genome target and cannot enter headline contrasts.

## Scientific guardrails

- Only `held_out_primary` samples can contribute to the headline generalization comparison.
- Graph-member controls, V2 rescue-source controls and held-out samples remain separate sample roles and separate result rows.
- Real-read outputs are mapping quality/control metrics, not mapping accuracy. Correct-locus accuracy is available only in the optional truth-labelled simulation mode.
- MAPQ is summarized within a mapper. Minigraph and minimap2 MAPQ values must not be directly compared.
- Minigraph SV-level rGFA and Minigraph-Cactus GBZ are distinct graph representations. The workflow does not silently convert one into the other.
- V2 is experimental. Its structural integrity passed, but its overall cohort-accounting QC is currently **FAIL**. This is recorded as provenance and is not interpreted as a mapping result.
- Reads are observations, not independent biological replicates. The workflow reports per-sample effects and no read-level inferential p-values.

## Files

- `workflow/rules/read_graph_benchmark.smk`: standalone entry point and named targets.
- `config/read_graph_benchmark.yaml`: graph registry, samples, rescue events, TMEM pairs, seeds, mapper parameters and resources.
- `workflow/scripts/read_graph_registry.py`: graph preflight, sample membership, leakage gates, versions and provenance.
- `workflow/scripts/select_fastq_reads.py`: deterministic coverage selection by read name. It writes a name list and statistics, not another multi-gigabyte FASTQ.
- `workflow/scripts/mapping_metrics.py`: per-sample PAF/GAF metrics, explicit paired contrasts and simulated-read accuracy.
- `workflow/scripts/rescue_capture.py`: targeted V1/V2 rescue traversal and anchor evidence.
- `workflow/scripts/select_tmem_pilot.py`: fixed pairs plus phenotype-blind stratified random samples.
- `workflow/scripts/tmem_evidence.py`: strict junction, allele-edge k-mer and reference-path evidence.
- `workflow/scripts/simulate_hifi.py`: optional deterministic phased-assembly simulation with truth in read names.
- `workflow/scripts/normalize_sv_catalog.py`: normalized ingestion plus an isolated clustering stage.
- `workflow/scripts/read_graph_benchmark_report.py`: compact Markdown synthesis.
- `workflow/scripts/clean_read_graph_benchmark.py`: dry-run-first cleanup of completed, regenerable intermediates.
- `workflow/envs/read_graph_benchmark_*.yaml`: per-rule Conda environments.

## Expected inputs

Paths are interpreted from the repository root on Snellius. The default server checkout is:

```text
/gpfs/work3/0/qtholstg/pangenome/whole_pans
```

The defaults expect the existing layout:

```text
../refs/CHM13.fa.gz
../whole_pangenome/sv_pangenome/results/
../whole_pangenome/sv_pangenome/local_sv_graph_v2/
../whole_pangenome/hprc_plus_cohort/
/gpfs/work3/0/qtholstg/pangenome/HPRC/
/gpfs/work3/0/qtholstg/pangenome/TMEM_UPD/
```

Do not replace these with the local SSHFS mount path. If the server layout changes, edit only `config/read_graph_benchmark.yaml`.

### Graph registry schema

Each `graph_registry` list entry records:

- `graph_id`, label and comparison scope;
- `graph_type`: `linear_fasta`, `minigraph_rgfa` or `gbz`;
- backbone and exact reference/release description;
- graph path and explicit companion indexes;
- membership manifest and its parser format;
- checksum, checksum file and whether to recompute it;
- declared mapping capabilities;
- experimental and QC status;
- provenance source.

Keeping the registry as a YAML list allows the preflight helper to detect duplicate graph IDs. Large graph hashes are recorded but not recomputed by default. Set `verify_checksum: true` for a deliberate full-file verification run.

Preflight checks file existence and size, `.tmp` final targets, companion indexes, headline CHM13 compatibility, representation/mapper compatibility, membership auditability and configured checksum provenance. It always writes `graph_preflight.tsv` and `graph_preflight.md`. Separate gate rules then fail before mapping if a requested target is unavailable.

The HPRC-plus-local graph is marked `preflight_passed_membership_auditable` after the completed server preflight. A future empty or `.tmp` artifact is still reported as unavailable and must not be renamed or treated as a finished graph.

### Sample schema and roles

Long-read samples are configured in the `samples` list. The generated `samples.tsv` contains:

```text
sample_id, data_type, long_read, r1, r2, nominal_coverage,
selected_coverage, cohort, population, sex, batch,
archive_cohort, stage_from_archive, truth_assembly_hap1,
truth_assembly_hap2, truth_callset, callable_regions, role,
include, inclusion_reason, exclusion_reason
```

Allowed roles are:

- `local_v1_member_control`;
- `hprc_member_control`;
- `v2_rescue_source_control`;
- `held_out_primary`;
- `tmem_candidate`;
- `tmem_matched_control`;
- `tmem_random_blinded`.

Local membership is derived from `sv_pangenome.ordered_assemblies.tsv`. V2 membership is the V1 donor set plus donors with `PASS` rescue rows. HPRC membership is derived from its assemblies table. Two haplotypes from one donor count once because membership is stored as a set of donor/sample IDs.

For a held-out sample, every headline graph membership must be known and false. A match or an unauditable membership manifest produces a `FAIL` row in `leakage_audit.tsv`, and the mapping gate stops. The default config contains two local-V1 member controls and two V2 rescue-source controls, but no independent held-out samples, so it deliberately cannot produce a headline generalization estimate. Add genuinely absent samples as `held_out_primary` only after confirming their identifiers against all membership manifests.

The two local-V1 member controls are `19R0534` and `2348252`. They were selected reproducibly from 16 eligible AD-cohort donors with both haplotypes in local V1 and post-QC `PASS` status by sorting SHA-256 hashes of the exact UTF-8 string `20260806:<sample_id>` (without a trailing newline). This avoids choosing by graph order or distance to CHM13. Their raw HiFi FASTQs are recalled from the existing tape archive and converted by the benchmark preparation rules; the source BAMs and copied BAMs remain temporary, while final FASTQs are retained in the shared HiFi staging directory.

With two member controls and two rescue-source controls, this is a balanced descriptive pilot, not an adequately powered population-level comparison.

## Long-read benchmark

Coverage normalization is deterministic. The helper first records the input compressed-file SHA-256, total reads and bases, then selects reads using a seeded BLAKE2b threshold on read names. Mapping streams those selected names through `seqtk subseq`; it does not write a second selected FASTQ. Set `selected_coverage` per sample to `15`, `20` or `full`. Full input is therefore preserved unchanged for sensitivity runs.

All SV-level targets use:

```text
minigraph -cx lr
```

with the same configured extra parameters. The plain CHM13 FASTA is accepted by Minigraph as a single-path graph control. It is also mapped with the conventional `minimap2 -x map-hifi` baseline. The exact Minigraph sequence-to-graph and FASTA-control behavior follows the official Minigraph user guide; current `vg giraffe` GBZ/index arguments follow the official vg mapping guide. Tool versions are recorded at runtime.

Per-sample metrics include read/base mapping fractions, weighted primary identity, clipping, split/supplementary evidence, secondary/multimapping rate, alignment-length and MAPQ quantiles, runtime and peak RSS. Explicit Minigraph-only deltas are:

- HPRC minus linear;
- local V1 minus HPRC;
- local V2 minus local V1;
- HPRC-plus-local minus HPRC.

Each metric row records graph-specific membership. Each delta keeps the left/right membership states, a transition such as `nonmember_to_member`, `sample_role` and `headline_eligible`; control roles are never pooled with held-out samples. In particular:

- local-V1 member controls are members of local V1, local V2 and HPRC-plus-local, but not plain CHM13 or HPRC;
- the two rescue-source controls are nonmembers of local V1 and HPRC-plus-local but members of local V2;
- `HPRC_plus_local_minus_HPRC` therefore compares an out-to-in transition in the member controls against an out-to-out transition in the rescue controls;
- `local_V2_minus_local_V1` provides the complementary out-to-in transition for the rescue controls and an in-to-in control for local-V1 members.

Raw GAFs, per-job metric tables and read-name selection lists are workflow
intermediates and are marked `temp()`. Snakemake removes each one after every
consumer in the active DAG has completed. V1/V2 rescue validation reuses the
same canonical long-read GAFs as the mapping benchmark; it no longer writes a
second copy under `v2_rescue/alignments/`.

## V2 rescue validation

The default targeted events are the three requested `PASS` rescues:

- `PSV_49dd19af3e38911c`: sample 8062838, insertion 8,707 bp;
- `PSV_7f13e3c00c873932`: sample 8062838, insertion 2,895 bp;
- `PSV_9ba713bbb71cd157`: sample 7326950, deletion 6,555 bp.

The helper refuses any configured rescue whose summary status, type or approximate length does not match. It reports both-anchor records, V2-specific traversal records, residual maximum indel and clipping summaries for V1 and V2, ambiguity from secondary mappings, and one of:

- `captured_in_V2`;
- `not_captured`;
- `ambiguous_repetitive`.

Evidence is further labelled `source_only_evidence` or `independently_replicated`. Independent carriers and non-carrier controls must be added explicitly to each event and also enabled as long-read samples. The default events have no known independent carriers, so the workflow cannot claim independent replication.

## TMEM pilot selection and evidence

The six configured candidate/control pairs are fixed exactly as requested. The first sample in each pair is labelled `insertion_present` and the second `insertion_absent`, based on the prior strict junction screen. These expectations are not independent truth because the same read data supplied that screen. The evidence table therefore reports both the expectation basis and concordance rather than presenting this pilot as an unbiased genotype-accuracy benchmark.

The locked pilot threshold requires at least three mapped reads on each insertion junction. In the prior screen, candidate samples had at least 21 reads on their limiting junction and matched controls had at most two. Four additional samples are selected with a fixed seed from remaining paired-end samples and remain blinded. Selection uses batch and the configured mapping-depth proxy when available; it never opens or uses phenotype data. When the fixed-pair batch distribution permits, the default four-sample allocation selects across batches 1–4. The resulting manifest records the seed and reason.

Only the regional GBZ supports this pilot. The workflow builds its registered `.dist`, `.shortread.withzip.min` and `.shortread.zipcodes` products together with the pinned vg mapping version, leaving the original Minigraph-Cactus indexes unchanged. Giraffe receives all paths explicitly and a per-job temporary copy of the generated distance index. The scratch copy is timestamped to its associated minimizer, avoiding both shared read/write access and vg's dependency-mtime rejection. Whole-genome rGFA targets remain `unsupported` for short reads unless separately constructed, path-preserving GBZ targets are later registered. Such construction must be an explicit provenance-tracked stage; it is not hidden in this workflow.

TMEM GAFs and per-sample evidence tables are also temporary. The merged
`tmem_genotype_evidence.tsv`, logs and benchmarks remain persistent.

TMEM calls require strict support for both ALT junctions. The output also reports ALT left/right edge k-mers, reference-path support, allele balance, no-call/ambiguity, expected-allele concordance and agreement across supported registered mapping targets. Internal insertion sequence is not used as carrier evidence because it is repetitive/background-heavy. AluYb8 is also excluded because all 358 screened samples were called Alu-present. Short reads are never described as spanning the complete ~18.7 kb allele.

## Optional simulated mapping truth

Set `simulation.enabled: true` and `simulation: true` on selected samples with two phased truth FASTAs. The lightweight simulator uses the fixed seed and configured length/error profile. An `empirical_profile_json` may override the length and error-rate fields. Read names use:

```text
sim|sample|haplotype|contig:start-end|strand|ordinal
```

Simulated FASTQs, GAFs and per-job accuracy tables are temporary; the merged
accuracy table and simulation statistics are retained.

Results report correct chromosome/locus, wrong locus, aligned fraction, clipping and identity. Correct haplotype/path is reported only where the graph alignment path makes it assessable; Minigraph rGFA generally lacks embedded donor haplotype paths, so an explicit not-assessable value is retained.

## Unified SV catalog

`sv_catalog.inputs` accepts existing normalized provisional tables and VCF inputs. Output preserves CHM13 anchors, optional GRCh38 coordinates, graph node/offset, type/length, literal inserted sequence and SHA-256, source graph/sample/method, confidence, membership and projection status. Unprojected or unresolved rows are retained.

Clustering is a separate rule and output. It uses configurable CHM13 breakpoint distance, length similarity and, when both inserted sequences are available, sequence similarity. This keeps ingestion independent from heuristic equivalence decisions. Raw VCF record counts are never presented as biological SV counts.

## Outputs

The requested compact tables are under `../whole_pangenome/read_graph_benchmark/tables/` by default:

```text
graph_registry.tsv
graph_preflight.tsv
samples.tsv
sample_graph_membership.tsv
leakage_audit.tsv
mapping_metrics_per_sample.tsv
paired_graph_contrasts.tsv
simulated_mapping_accuracy.tsv
v2_rescue_capture.tsv
normalized_sv_catalog.tsv
tmem_genotype_evidence.tsv
tool_versions.tsv
provenance.tsv
```

The isolated clustering output is `normalized_sv_clusters.tsv`. Reports are under `report/`.

## Snellius commands

Run every command from the server checkout. No scheduler directives are embedded; resource requests are passed to the existing profile.

```bash
cd /gpfs/work3/0/qtholstg/pangenome/whole_pans
```

### Quota-conscious full run

Conda environments contain many small files. Keep those disposable files out
of the project checkout by giving Snakemake a scratch Conda prefix. The custom
`gaf_writers` resource limits the workflow to one multi-gigabyte GAF-producing
job at a time; without the command-line limit, the resource annotation has no
effect.

```bash
srun -p genoa -c 192 -t 5-0:0:0 \
  snakemake -c 192 \
  --snakefile workflow/rules/read_graph_benchmark.smk \
  --use-conda \
  --conda-prefix "$TMPDIR/whole-pans-read-graph-conda" \
  --resources gaf_writers=1 \
  --rerun-incomplete \
  --retries 1 \
  --rerun-triggers mtime \
  read_graph_benchmark_all
```

Run this command from `/gpfs/work3/0/qtholstg/pangenome/whole_pans`. While
resolving the quota problem, omit `--keep-going`: if a downstream summarization
fails, it is safer to stop than to let unrelated mapping jobs continue filling storage.
The retained scientific results still go to the configured project results
directory; only the reproducible Conda environments use `$TMPDIR`.

Check both space and inode usage on the server:

```bash
prjspc-quota /gpfs/work3/0/qtholstg
scratch-quota
du --inodes -d 2 .snakemake \
  ../whole_pangenome/read_graph_benchmark | sort -n
```

If `.snakemake/conda` is large, ask the installed Snakemake version to list
and remove only environments no longer referenced by this workflow:

```bash
snakemake -s workflow/rules/read_graph_benchmark.smk \
  --use-conda --list-conda-envs
snakemake -s workflow/rules/read_graph_benchmark.smk \
  --use-conda --conda-cleanup-envs
```

Do not remove a Conda prefix while a job is using it. The cleanup script below
deliberately does not touch Conda environments or Snakemake metadata.

### Cleanup of an existing run

First run the cleanup script without `--delete`. It inventories the results
and lists only intermediates whose non-empty downstream result exists and is
at least as new:

```bash
python3 workflow/scripts/clean_read_graph_benchmark.py --verbose
```

After confirming that no Snakemake process or cluster job is using this result
directory, delete those candidates:

```bash
python3 workflow/scripts/clean_read_graph_benchmark.py --delete
```

The script refuses deletion while entries exist in
`.snakemake/locks`. It never selects source reads or graphs, final
tables, reports, logs, benchmark timing files, provenance, or QC files. A stale
Snakemake lock must be investigated separately; do not remove it merely to
force cleanup while jobs may still be active. It also prunes empty
`.snakemake/tmp.*` directories, but does not touch Conda or metadata content.

Workflow dry-run:

```bash
snakemake -s workflow/rules/read_graph_benchmark.smk \
  --configfile config/read_graph_benchmark.yaml \
  --profile ~/.config/snakemake/zslurm/ \
  --use-conda \
  -n read_graph_benchmark_all
```

Graph/sample preflight only (writes reports and does not invoke an expensive mapping gate):

```bash
snakemake -s workflow/rules/read_graph_benchmark.smk \
  --configfile config/read_graph_benchmark.yaml \
  --profile ~/.config/snakemake/zslurm/ \
  --use-conda \
  read_graph_preflight
```

Recall and convert the two selected local-V1 member controls before the mapping pilot. This target uses the configured archive resource and can wait for tape recall:

```bash
snakemake -s workflow/rules/read_graph_benchmark.smk \
  --configfile config/read_graph_benchmark.yaml \
  --profile ~/.config/snakemake/zslurm/ \
  --use-conda \
  prepare_long_read_member_controls
```

Small normalized long-read pilot:

```bash
snakemake -s workflow/rules/read_graph_benchmark.smk \
  --configfile config/read_graph_benchmark.yaml \
  --profile ~/.config/snakemake/zslurm/ \
  --use-conda \
  long_read_pilot
```

The long-read target also pulls in the preparation rules automatically if either selected member FASTQ is absent. Running the preparation target separately makes tape-recall failures easier to diagnose.

V2 rescue validation (gated only on V1 and V2):

```bash
snakemake -s workflow/rules/read_graph_benchmark.smk \
  --configfile config/read_graph_benchmark.yaml \
  --profile ~/.config/snakemake/zslurm/ \
  --use-conda \
  v2_rescue_validation
```

Build only the version-matched TMEM Giraffe indexes (no reads are mapped):

```bash
snakemake -s workflow/rules/read_graph_benchmark.smk \
  --configfile config/read_graph_benchmark.yaml \
  --profile ~/.config/snakemake/zslurm/ \
  --use-conda \
  prepare_tmem_giraffe_indexes
```

TMEM pilot (also builds the indexes automatically when absent):

```bash
snakemake -s workflow/rules/read_graph_benchmark.smk \
  --configfile config/read_graph_benchmark.yaml \
  --profile ~/.config/snakemake/zslurm/ \
  --use-conda \
  tmem_pilot
```

Optional catalog-only stage:

```bash
snakemake -s workflow/rules/read_graph_benchmark.smk \
  --configfile config/read_graph_benchmark.yaml \
  --profile ~/.config/snakemake/zslurm/ \
  --use-conda \
  unified_sv_catalog
```

For a full-coverage sensitivity run, set the selected samples' `selected_coverage: full`, choose a separate `results` directory so the normalized pilot is not overwritten, run preflight, then run `long_read_pilot`. To add more samples, add registry entries with explicit roles and inputs; do not relabel controls as held out.

## Required server validation

Before interpreting results:

1. review `graph_preflight.tsv` and `graph_preflight.md`;
2. confirm `leakage_audit.tsv` is PASS for every held-out sample/graph pair;
3. inspect tool versions and the actual Minigraph/vg commands in logs;
4. verify that all per-sample metric rows exist before paired contrasts;
5. manually review rescue-event graph paths where traversal naming is absent or repetitive;
6. inspect TMEM candidate, matched-control and blinded rows separately;
7. review the V2 QC warning in both provenance and the final report;
8. only then expand coverage or cohort size.
