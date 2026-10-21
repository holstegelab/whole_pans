# Minigraph-Cactus from existing SV graphs

`workflow/rules/minigraph_cactus_reuse.smk` builds two graphs from the SV GFAs
already listed in `config/config.yaml`:

- **Local V1:** `sv_pangenome.minigraph.gfa`
- **HPRC plus local V1:** `hprc-v2.0-mc-chm13.cohort.sv.gfa`

The workflow reads these files without changing them. It checks the FASTA
inputs, rewrites source tags in 16 independent GFA chunks, compresses the
ordered chunks into one seed GFA, then maps each assembly in its own Snakemake
job before chromosome split, batch alignment, and join. Failed seed chunks or
assembly mappings can be retried without repeating completed jobs. It does
**not** rerun `cactus-minigraph`.

Graph mapping uses the shell tools and filters called by Cactus 3.3.0:
`cactus_sanitizeFastaHeaders`, `minigraph`, `gaf2unstable`, `gaffilter`,
`gaf2paf`, `rgfa2paf`, `vg`, and `filter-paf-deletions`. Snakemake owns
one mapping job for each seqfile row. Shared rules prepare the graph FASTA,
merge the PAFs, add the CHM13 rGFA alignment, and apply the deletion filter.
`cactus-graphmap-split` separates reference contigs in the next stage. The
graph mapping rules call no custom Python script. The Conda environment includes
`vg` and `hal2vg` to provide the standalone deletion-filter commands. These
command options follow the pinned
[Cactus 3.3.0 mapping implementation](https://github.com/ComparativeGenomicsToolkit/cactus/blob/v3.3.0/src/cactus/refmap/cactus_graphmap.py)
and [configuration](https://github.com/ComparativeGenomicsToolkit/cactus/blob/v3.3.0/src/cactus/cactus_progressive_config.xml);
review them if the Cactus version or configuration changes.

The join environment also includes `gfaffix=0.2.2` and `vcfbub=0.1.1`, which
Cactus 3.3.0 invokes during graph clipping and VCF generation. If join failed
with `cannot run gfaffix`, rerun the same Snakemake targets after updating the
environment file. The join rule detects its existing Toil job store and passes
`--restart`; keep that store in place. Completed alignment outputs are reused.

The HPRC build uses different reference FASTA files from the local build, so
the final comparison is descriptive rather than an isolated effect of adding
HPRC samples.

## Run on Snellius

Run from the server checkout. First inspect the local build plan, then run
each build and their comparison:

```bash
cd /gpfs/work3/0/qtholstg/pangenome/whole_pans
mc() {
  snakemake --profile ~/.config/snakemake/zslurm/ \
    --snakefile workflow/rules/minigraph_cactus_reuse.smk \
    --use-conda --rerun-incomplete --retries 0 --rerun-triggers mtime \
    --jobs 8 "$@"
}
mc -n minigraph_cactus_reuse_local_v1
mc minigraph_cactus_reuse_local_v1
mc -n /gpfs/work3/0/qtholstg/pangenome/whole_pangenome/minigraph_cactus_reuse/hprc_plus_local_v1_seeded_mc_smk_v2/stages/minigraph/seed-normalization.json
mc /gpfs/work3/0/qtholstg/pangenome/whole_pangenome/minigraph_cactus_reuse/hprc_plus_local_v1_seeded_mc_smk_v2/stages/minigraph/seed-normalization.json
mc -n /gpfs/work3/0/qtholstg/pangenome/whole_pangenome/minigraph_cactus_reuse/hprc_plus_local_v1_seeded_mc_smk_v2/stages/graphmap/assemblies/NA20809.2.paf
mc /gpfs/work3/0/qtholstg/pangenome/whole_pangenome/minigraph_cactus_reuse/hprc_plus_local_v1_seeded_mc_smk_v2/stages/graphmap/assemblies/NA20809.2.paf
mc minigraph_cactus_reuse_hprc_plus_local_v1
mc minigraph_cactus_reuse_compare
```

The new build IDs are `local_v1_seeded_mc_smk_v2` and
`hprc_plus_local_v1_seeded_mc_smk_v2`. They keep the new batch-alignment
outputs separate from any earlier per-bin run. Alignment and join each request
96 cores and 1,300 GB RAM, so they need a fat node. The split, alignment, and
join stages keep their own Toil job stores for restart after failure. The seed
chunk logs are `logs/seed-chunk-00.log` through `logs/seed-chunk-15.log`. To inspect or run a
single stage, target its `stages/minigraph/seed-normalization.json`,
`metadata/graphmap.done`, `split.done`, `align.done`, or `join.done` file under
that build ID. Individual mapping logs are `logs/graphmap-<assembly>.log`, and
their PAFs are `stages/graphmap/assemblies/<assembly>.paf`. The seed compression
uses eight threads. Each assembly mapping requests six threads and 96 GB; the
global PAF filters remain separate jobs. Later Cactus stages run in their
established order.

The October 5 failure was an rGFA source tag such as
`NA20809.2#2#CM094349.1`. Its matching FASTA manifest name is
`NA20809#2#CM094349.1`; Minigraph added the `NA20809.2` seqfile label. The
current normalizer accepts this form and converts it to the manifest PanSN
name. Run the workflow from the updated server checkout so Snakemake stages
the current script, rather than the older source-cache copy shown in that log.

If an input FASTA, source GFA, or Cactus setting changes after a build starts,
choose a fresh build ID in `config/config.yaml`. Completed stages use `.done`
markers and do not automatically detect every change inside a FASTA manifest.

The comparison writes a metrics TSV, per-chromosome TSV, and Markdown summary
under `minigraph_cactus_reuse/comparison/`.
