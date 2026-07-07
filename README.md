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
- every BLAST database `prefix` and `marker`;
- all three result directories and per-rule resource requests.

`existing_qc.included_assemblies` points to the final inclusion list from the
completed first QC run. The workflow treats this as an external, fixed
selection, so requesting decontamination does not rebuild the full original QC
checkpoint. Missing CHM13/hg38 alignments needed by decontamination can still
be regenerated under `results.qc`. Remove the `existing_qc` section to run the
first QC pass as part of this workflow.

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

Each configured BLAST database must contain only the named non-human group. Do
not use an unfiltered `nt` database because the database name is treated as the
taxonomic label for every hit.

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
snakemake --cores 32 --use-conda all
```

Each stage file is also a standalone entry point. From the project root
(`/gpfs/work3/0/qtholstg/pangenome`), run:

```bash
# Original-assembly QC
snakemake --profile ~/.config/snakemake/zslurm/ \
  --snakefile whole_pans/workflow/rules/QC.smk \
  --use-conda --rerun-incomplete

# Decontamination, reusing existing_qc.included_assemblies
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

With `existing_qc.included_assemblies` configured, `decontamination` starts
from that completed selection. `post_decontamination_qc` automatically runs
missing decontamination work. If `existing_qc` is removed, QC selection uses a
Snakemake checkpoint; run `qc`, then repeat
`snakemake -n post_decontamination_qc` to inspect the later stages in full.

## Main outputs

QC outputs under the configured `results.qc/summary` directory:

- `assembly_qc.tsv`: metrics and PASS/WARN/FAIL reasons per assembly;
- `sample_qc.tsv`: paired-haplotype decisions, including `PARTIAL` samples;
- `graph_included_assemblies.txt`: assemblies accepted for decontamination;
- `graph_excluded_assemblies.tsv`: excluded assemblies and reasons.

Decontamination outputs under the configured `results.decontamination`
directory:

- `summary/contamination_summary.tsv`: contamination counts and bases;
- `summary/contig_actions.tsv`: REMOVE, SPLIT, and REVIEW decisions;
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

Reference alignment is intended for gross outlier and contamination checks.
Real human structural variation and reference-specific sequence should not be
interpreted as assembly error without supporting read evidence.
