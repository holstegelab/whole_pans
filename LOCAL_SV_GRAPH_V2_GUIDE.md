# Local SV graph v2 from read-rescued alleles

September review update: the original bubble-delta coordinates were affected
by a bytes-to-integer parsing error. The parser is now corrected; use
`review_corrected_v2_qc` as described in [REVIEW_FOLLOWUP_GUIDE.md](REVIEW_FOLLOWUP_GUIDE.md)
to regenerate QC in a separate directory. Rescue PASS and nearby graph bubbles
do not by themselves establish independent allele validation or graph capture.

`workflow/rules/local_sv_graph_v2.smk` is the standalone continuation of
`workflow/rules/novel_sv_discovery.smk`. It creates graph 4 of the planned
comparison:

1. the frozen local SV graph (unchanged);
2. the HPRC SV graph (unchanged);
3. the HPRC graph plus the 196 local assemblies (a later workflow);
4. the frozen local graph plus unique HiFi-validated rescued SV alleles
   (this workflow).

The final graph is:

```text
../whole_pangenome/sv_pangenome/local_sv_graph_v2/graphs/
sv_pangenome.local.v2.rescued.minigraph.gfa
```

The HPRC files are deliberately not inputs to this stage, so this graph remains
a clean local-v1-versus-local-v2 comparison.

## What the workflow does

The workflow:

1. selects candidate loci from the completed novel-SV catalog;
2. chooses one excluded or fragmented carrier assembly as the source allele;
3. maps each selected sample to its PacBio BAM on tape;
4. requests the selected tape BAMs with `daget` in chunks of up to 100 and
   waits for every `dals` status to become `DUL` or `REG`;
5. copies each recalled BAM to temporary compute-visible project storage from
   an archive/staging node;
6. checks that the BAM `@SQ` dictionary matches `../hdr.pacbio.txt` and that a
   PacBio CCS read group is present;
7. converts primary BAM records to one reusable, persistent HiFi FASTQ per
   selected sample;
8. recruits spanning reads, runs local `hifiasm`, and applies the existing
   two-sided rescue QC;
9. combines only unique `PASS` sequences;
10. incrementally adds those sequences to the checksum-verified frozen local
   Minigraph rGFA;
11. checks rGFA integrity and the stable V1-to-V2 graph delta;
12. accounts for every existing high-confidence graph-residual evidence row in
   the cohort catalog; and
13. reports whether any such SV evidence outside the included rescue-source
   assemblies is left without a rescue input in V2.

Samples in `../chc_samples_list.txt` use the `centenarian` tape directory. All
other samples in the assembly discovery manifest use `ad`. For example:

```text
/archive/hulsmanm/pacbio_processed/project/holstegelab/Share/pacbio/
data_processed/merged_data_ad_centenarians/ad/1008096/
1008096.merged.hifi.hg38.bam
```

Snellius mounts `/archive` only on interactive and staging/service nodes. After
`daget`, each BAM is therefore copied to `hifi/archive/compute_bams/` before a
compute job reads it. These BAM copies are `temp()` outputs: Snakemake keeps a
copy after failed conversion for retry, then removes it after successful FASTQ
creation. The converted FASTQ is a `protected()` output and remains under
`hifi/fastq/`, so later rescue plans can reuse it without another tape recall
or BAM-to-FASTQ conversion. Header QC tables remain under `hifi/header_qc/`.

During conversion, `samtools fastq` first writes an uncompressed FASTQ into a
hidden temporary directory beside the final FASTQ on project storage. `bgzip`
then reads that file and writes a validated compressed file in the same
directory; only the completed `.fastq.gz` is moved into its final name. Small
header files, selected-read extraction, local hifiasm assembly, the temporary
graph, and bubble calling use `${TMPDIR:-/tmp}`. Every job removes its own
temporary directory. The small tape-staging marker directory is a Snakemake
`temp()` output and is removed after its consumers finish. Ensure the
filesystem containing `local_sv_graph_v2.results` has enough temporary space
for the BAM copies and uncompressed FASTQs that may be in flight concurrently.

The Snakefile contains only the workflow structure. Archive access lives in
`workflow/scripts/archive_hifi_reads.py`, graph construction and provenance in
`workflow/scripts/build_incremental_sv_graph.py`, and rescue biology in
`workflow/scripts/build_clean_sv_rescues.py`. Each helper is also usable as a
standalone command-line program.

Only reusable FASTQs and review/audit outputs remain on project storage. To
measure the FASTQs after conversion:

```bash
find \
  ../whole_pangenome/sv_pangenome/local_sv_graph_v2/hifi/fastq \
  -maxdepth 1 -type f -name '*.hifi.fastq.gz' -printf '%s\n' |
awk '{bytes += $1; files++}
     END {printf "files=%d bytes=%d GiB=%.2f TiB=%.3f\n",
                 files, bytes, bytes / 2^30, bytes / 2^40}'
```

## Recommended first run

Run this on Snellius from the `whole_pans` directory. First create only the
rescue plan:

```bash
snakemake \
  --snakefile workflow/rules/local_sv_graph_v2.smk \
  --profile ~/.config/snakemake/zslurm/ \
  --use-conda \
  ../whole_pangenome/sv_pangenome/local_sv_graph_v2/rescue/plan.tsv
```

Review its statuses before requesting reads:

```bash
awk -F '\t' '
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      if ($i == "rescue_id") id = i
      if ($i == "plan_status") status = i
      if ($i == "sample_id") sample = i
      if ($i == "svtype") type = i
    }
    print "rescue_id", "status", "sample", "svtype"
    next
  }
  {print $id, $status, $sample, $type}
' OFS='\t' \
  ../whole_pangenome/sv_pangenome/local_sv_graph_v2/rescue/plan.tsv |
  column -t
```

The default configuration can select up to 100 events and therefore may stage
many large tape BAMs. For a pilot, set either:

- `local_sv_graph_v2.plan.event_ids` to a reviewed one-ID-per-line file; or
- `local_sv_graph_v2.plan.max_events` to a small number such as 3–5.

Changing the selection configuration requires rebuilding the plan because
Snakemake does not automatically track arbitrary config values as input
timestamps:

```bash
snakemake \
  --snakefile workflow/rules/local_sv_graph_v2.smk \
  --profile ~/.config/snakemake/zslurm/ \
  --use-conda \
  ../whole_pangenome/sv_pangenome/local_sv_graph_v2/rescue/plan.tsv \
  --forcerun plan_local_v2_rescues
```

Do this before launching rescue jobs. If a plan is changed after rescues have
already run, use a new `local_sv_graph_v2.results` directory or first archive
the old result directories. The rescue helper intentionally refuses to
overwrite them.

Once the plan is accepted, inspect the complete dynamic DAG:

```bash
snakemake \
  --snakefile workflow/rules/local_sv_graph_v2.smk \
  --profile ~/.config/snakemake/zslurm/ \
  --use-conda \
  --dry-run
```

Then build graph v2:

```bash
snakemake \
  --snakefile workflow/rules/local_sv_graph_v2.smk \
  --profile ~/.config/snakemake/zslurm/ \
  --use-conda
```

The default target now includes V2 QC. This QC does not remap assemblies or
sequencing reads and does not add biological validation. It audits the frozen
graph, the V2 graph delta, rescue evidence, graph bubbles, and the existing V1
cohort discovery catalog. To request only the final QC report after the graph
already exists:

```bash
snakemake \
  --snakefile workflow/rules/local_sv_graph_v2.smk \
  --profile ~/.config/snakemake/zslurm/ \
  --use-conda \
  ../whole_pangenome/sv_pangenome/local_sv_graph_v2/qc/report/local_sv_graph_v2_qc_report.md
```

The report gives an accounting `PASS` only when every existing
`HIGH_CONFIDENCE` graph-residual evidence row of at least 50 bp in the rest of
the non-graph cohort maps to a catalog event whose rescue input was included in
V2. Otherwise it records `FAIL` and lists/counts the unaccounted catalog events.
Without remapping, this means "no previously detected SV was left
unaccounted for"; it does not prove that V2 absorbs every allele or that no
additional SV exists.

Because rescue jobs are discovered from a checkpoint, a dry-run before the
plan exists shows only jobs through that checkpoint. A second dry-run after
the plan exists shows the sample staging, conversion, and per-rescue jobs.

## Important outputs

```text
local_sv_graph_v2/
├── hifi/
│   ├── hifi_archive_reads.tsv
│   ├── fastq/
│   │   └── <sample>.hifi.fastq.gz
│   └── header_qc/
├── rescue/
│   ├── plan.tsv
│   ├── results/<rescue_id>/
│   ├── rescue_summary.tsv
│   ├── clean_rescues.fa
│   └── clean_rescue_graph_inputs.tsv
├── graphs/
│   └── sv_pangenome.local.v2.rescued.minigraph.gfa
├── metadata/
│   ├── sv_pangenome.local.v2.graph_summary.tsv
│   ├── sv_pangenome.local.v2.provenance.tsv
│   ├── sv_pangenome.local.v2.bubbles.tsv
│   └── tool_versions.tsv
└── qc/
    ├── tables/
    │   ├── graph_integrity.tsv
    │   ├── delta_segments.tsv
    │   ├── delta_links.tsv
    │   ├── bubble_delta.tsv
    │   ├── rescued_sv_graph_qc.tsv
    │   └── cohort_completeness.tsv
    └── report/
        └── local_sv_graph_v2_qc_report.md
```

The graph build stops rather than silently copying graph v1 when there are no
unique `PASS` rescues, when Minigraph adds no new graph structure, or when the
frozen graph checksum no longer matches its inventory. The provenance table
records checksums for the source graph, plan, rescue inputs, BAM-header
template, conda environment, and final graph.

The cleaned source assemblies are ordinary gzip FASTA files rather than BGZF.
When no `.fai` is available, the rescue helper streams the required contig
interval. This is slower than indexed access but avoids making multi-gigabyte
duplicate assemblies solely for random access.
