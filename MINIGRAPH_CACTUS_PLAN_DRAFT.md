# Draft plan: whole-genome Minigraph-Cactus graph

## Scope and existing inputs

Add an explicit `minigraph_cactus` Snakemake target alongside the existing SV-only `pangenome` target. Do not add it to `rule all`: this is a large, separately launched build. Do not overwrite or pass the existing `sv_pangenome.minigraph.gfa` to `--inGFA`. That GFA was built by the project's own Minigraph wrapper with renamed contigs and has not been shown to meet Cactus's seqfile/source naming requirements. Let `cactus-pangenome` build its own SV graph and base-level graph from the accepted FASTAs.

The authoritative haplotype input is `pangenome.assemblies.included_list`, currently the post-decontamination QC `graph_included_assemblies.txt` (196 cleaned haplotypes in the mounted snapshot). The old `sv_pangenome.ordered_assemblies.tsv` has 196 assemblies plus CHM13 and hg38 and is useful as an audit cross-check, but its Mash order must not be imposed on Cactus: Cactus computes its own Minigraph order. Reference FASTAs already exist in `references.CHM13` and `references.hg38`; use CHM13 as the primary chromosome-scale reference and include GRCh38 as a second promoted reference path if its sequence/contig naming passes preflight. `references.hg38` is already a primary chr1–22/X/Y FASTA, avoiding GRCh38 ALT loci. Decide explicitly whether to add chrM to that reference; do not silently mix different contig sets.

## Seqfile preparation

Create `workflow/scripts/prepare_minigraph_cactus_seqfile.py` and a small Snakemake rule. The script reads the QC inclusion list, not the pre-QC cleaned list; accepts the two reference paths as arguments; and writes a two-column whitespace-delimited seqfile plus `input_manifest.tsv` with sample ID, haplotype, original FASTA path, and QC-list line number. Use the server paths recorded in config and the list; never derive server paths from this local mount or call `realpath`.

Map filenames of the documented form `SAMPLE.hifi.hifiasm.bp.hap1.p_ctg.clean.fa.gz` and `...hap2...` to `SAMPLE.1` and `SAMPLE.2`; use `CHM13` and `GRCh38` for the two reference rows. Preserve a passing haplotype of a `PARTIAL` sample under its actual `.1` or `.2` label. Reject unrecognized names, duplicated sample/haplotype labels, repeated paths, missing or empty inputs, whitespace in paths, and a missing/empty inclusion list. Check that each input FASTA's contig IDs are unique **within that FASTA** and compatible with Cactus's `#` handling on the execution server (a lightweight streaming preflight, not a local run). Record the exact input list and Cactus version for provenance. The cleaned FASTAs already prefix contig IDs with donor/assembly IDs, which may make output path names longer; do not rewrite them without an explicit coordinate and provenance policy.

## Snakemake structure and outputs

Add `workflow/rules/minigraph_cactus.smk`, include it in `Snakefile`, and give it a dedicated `minigraph_cactus` target. Keep a standalone entry point usable with `--snakefile workflow/rules/minigraph_cactus.smk` if it can be done without duplicating rules. Add a distinct `minigraph_cactus:` section to `config/config.yaml` containing input-list/reference keys, server-side shared result directory, `jobStore`, `outName`, Cactus release, output switches, Toil batch system, work directory, and Cactus/Toil resource knobs. Paths must be configurable; the stage should not depend on the older SV graph or change the default target.

Use one `cactus-pangenome` invocation in one Snakemake build rule so Toil manages Cactus's own internal stages and restart state. The command skeleton is:

```bash
cactus-pangenome JOBSTORE SEQFILE \
  --outDir OUTDIR --outName local_mc_chm13 \
  --reference CHM13 GRCh38 --vcfReference CHM13 \
  --gfa --gbz --vcf --giraffe --validate \
  --maxCores MAX_CORES --mgCores MG_CORES \
  --mapCores MAP_CORES --consCores CONS_CORES \
  --indexCores INDEX_CORES
```

Make the Toil batch system and memory/disk settings configurable rather than embedding a scheduler choice in workflow code. For the first build, request the clipped `<outName>.gfa.gz`, `<outName>.gbz`, CHM13 `<outName>.vcf.gz`, the default `.d2.gbz/.dist/.min` Giraffe indexes, `<outName>.sv.gfa.gz`, and `<outName>.stats/` reports. The `.d2` indexes fit the repository's existing read benchmark naming, but they represent a frequency-filtered mapping graph, while the ordinary GFA/GBZ is clipped. Expose `--filter` if the benchmark later needs a different threshold. Do **not** combine `--giraffe` with `--haplo`; a later dedicated mapping profile may use `--haplo` with its extra k-mer preparation instead. HAL, ODGI and a GRCh38 VCF can be added only when required, because they materially increase storage and run time.

Keep Cactus's `--outDir`, Toil job store, work directory, logs and Snakemake benchmark separate. The Snakemake rule should publish a small success marker/manifest only after Cactus exits successfully and required files are present and nonempty; the marker should include seqfile digest, command, version and exact output paths. Preserve the job store and partial `outDir` on failure for restart. Do not use a Snakemake `directory()` output that causes cleanup of the whole Cactus output tree after an interrupted run. A subsequent run with the same inputs can pass Toil's `--restart` when the job store exists; a changed seqfile/config should use a new versioned job store and output directory, rather than restarting incompatible state. Avoid running two Snakemake jobs against the same store.

## Environment and resources

Add a dedicated `workflow/envs/minigraph_cactus.yaml` with pinned `cactus=3.3.0` from Bioconda/conda-forge (Bioconda lists this release) and only helper dependencies actually needed. Confirm on the **server** that the environment supplies `cactus-pangenome` and its required pangenome binaries; if the cluster deployment instead uses the official release binaries or container, expose an executable/runtime setting and pin its version. Do not install or run anything locally. The existing `workflow/envs/pangenome.yaml` is for the SV-only Minigraph stage and should stay separate.

Treat Snakemake's rule resource request as the Toil driver allocation, not the sum of Toil worker requests. Cactus's `--mgCores`, `--mgMemory`, `--mapCores`, `--consCores`, `--consMemory`, `--indexCores`, `--indexMemory`, plus `--maxCores`, `--maxMemory` and `--maxDisk`, need server-calibrated config values. The official human example cites roughly 64 cores, 512 GB RAM and 3 TB working disk for a single-machine-style build; do not assume those capacities exist here. Use shared storage for the job store and output, and a worker-accessible work directory with sufficient scratch. If Toil uses the cluster scheduler internally, ensure the Snakemake profile does not constrain its spawned worker jobs; keep scheduler flags in the server profile/config, not in this rule.

## Server verification and open decisions

On the server, first inspect the generated seqfile and manifest, count rows and reference entries, confirm all selected FASTAs exist, and check representative and duplicate contig names. Run a tiny pilot with two or three accepted haplotypes and both references in a distinct output/job store before committing to all 196. Then use Snakemake's `-n` for the production target and launch with the server profile, `--use-conda` and `--rerun-incomplete`. After completion, inspect Cactus `<outName>.WARNING` if present, `<outName>.stats/clipped-by-genome.tsv`, `graph-stats.tsv`, input-contig accounting, reference walk names in GFA, and whether GBZ/VCF/index outputs match the declared target. Do not treat a missing `.WARNING` alone as a sufficient QC result.

Open before production: (1) confirm Cactus 3.3.0 and its dependency deployment on Snellius; (2) choose Toil batch system and realistic memory/disk/cores for this site; (3) confirm CHM13 and GRCh38 contig sets/builds, especially chrM; (4) decide whether only CHM13 VCF or both coordinate systems are required; (5) decide whether the benchmark should use default `.d2` mapping indexes or a chosen filter/`--haplo` profile; and (6) check whether any cleaned FASTA headers contain `#`, duplicate IDs, or names Cactus will transform.

Sources: [official Cactus pangenome documentation](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md), [official Cactus installation README](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/README.md), [Bioconda Cactus recipe](https://bioconda.github.io/recipes/cactus/README.html), and local `pangenome-papers` construction/tool notes.
