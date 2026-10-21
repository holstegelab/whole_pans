# Staged Minigraph-Cactus on Snellius

`minigraph_cactus` is an explicit target, separate from `rule all` and the SV-only Minigraph graph. It reads the post-decontamination inclusion list, CHM13 and the configured GRCh38 primary FASTA. Preflight keeps surviving `.1` or `.2` haplotypes exactly as selected, hashes decoded FASTA content, records contig lengths, and lists both reference contig sets. GRCh38's configured FASTA lacks chrM.

The build has five [Cactus 3.3.0](https://github.com/ComparativeGenomicsToolkit/cactus/blob/v3.3.0/doc/pangenome.md) stages after preflight: Minigraph, graph mapping, chromosome split, batch alignment and graph join. Each Cactus command runs directly from its Snakemake rule. Stage products are under `<results_root>/<build_id>/stages/`, the final graph under `graph/`, and logs and benchmarks in their own directories. Every stage has a distinct Toil store under `<jobstore_root>/<build_id>/<stage>`. Use a new build ID for this layout; the previous JSON stage markers are incompatible with the current `.done` markers.

## Before a server run

From `/gpfs/work3/0/qtholstg/pangenome/whole_pans` **on Snellius**, inspect the `minigraph_cactus` section of `config/config.yaml`, the inclusion list, CHM13/GRCh38 build and contig names, writable GPFS results/store locations and the available local `TMPDIR` space/inodes. The rule Conda environment pins Cactus 3.3.0; inspect `--help` and the installed package version for `cactus-minigraph`, `cactus-graphmap`, `cactus-graphmap-split`, `cactus-align` and `cactus-graphmap-join`. Confirm the site profile maps `mem_mb`, `disk_mb`, `runtime_min` and small/fat node selection as intended. The baseline uses Toil `single_machine` within each Snakemake allocation and inherits `TMPDIR`; it does not submit separate Toil workers. Do not point worker scratch at a leader-only directory.

Use a **new build ID** and separate inclusion list for each pilot: first 2–3 passing haplotypes for CLI/output compatibility, then 12–20 diverse haplotypes for resource measurement. Those pilots do not establish the full-cohort join peak. The initial reservations are 64 cores/288 GB for Minigraph, 96/288 for mapping, 32/192 for split on small nodes, and 96/1,300 for align and join on fat nodes. Measure real memory, local scratch, GPFS occupancy and walltime before full production. If a small stage exceeds its node's memory, move that stage to a fat-node allocation through the server profile; do not raise its request above 336 GB on the small node.

Server-only commands, with the site's Snakemake profile and partition options added where required:

```bash
snakemake --snakefile Snakefile --use-conda -n minigraph_cactus
snakemake --snakefile Snakefile --use-conda --profile <site-profile> prepare_minigraph_cactus
snakemake --snakefile Snakefile --use-conda --profile <site-profile> minigraph_cactus
```

The preparation result is `metadata/seqfile.txt`, `input_manifest.tsv`, `contig_map.tsv.gz` and `preflight.json`. Cactus graphmap edits its input seqfile, so the rule gives it a private `metadata/seqfile.graphmap.txt` copy. After Cactus succeeds, a small input checker verifies the original genome paths and tree leaves, and the single `_MINIGRAPH_` entry pointing to `graph.fa.gz`. The rule then copies it to `metadata/seqfile.with-minigraph.txt` for split. The original seqfile remains untouched. The split stage's three-column `chromfile.txt` determines the HAL/raw VG pairs passed to join. Before marking split complete, the rule checks `minigraph.split.log`, `contig_sizes.tsv`, and each listed chromosome seqfile and PAF. Review the split log and `_AMBIGUOUS_` directory before accepting chromosome representation; interchromosomal and ambiguous sequence may be lost during assignment.

The whole-genome SV GFA is retained in `stages/minigraph/<out_name>.sv.gfa.gz`. Standalone join does not accept this one GFA as `--sv-gfa` for many chromosome VGs. Final `graph/` contains clipped GFA/GBZ, CHM13 VCF and raw VCF with indexes, full HAL, `.d2` Giraffe indexes, and graph statistics. The `.d2` index is frequency filtered. Review `.WARNING`, graph stats, CHM13/GRCh38 paths, sex chromosomes and representative loci. A five-column Cactus `--inputContigSizes` baseline is deferred until a server pilot confirms sanitized names; no join clipping/exclusion report is claimed for this version.

## Recovery

Inspect `logs/<stage>.log`, the stage products and its Toil store after a failure. Request the same target again: if the store exists, the rule passes `--restart`; if `metadata/<stage>.cactus_finished` exists, it skips Cactus and repeats the output checks before writing `<stage>.done`. A Cactus failure leaves its store and partial products; a later output-check failure may occur after `--clean=onSuccess` removed the store. If products exist but neither store nor `cactus_finished` exists, the rule refuses a fresh run: inspect the log and use a new build ID if the outcome is uncertain. A small checker compares input file metadata, the original seqfile, and locked preflight/settings before each stage. Use a new build ID whenever inputs, settings or workflow commands change. The no-output `minigraph_cactus` target rechecks the retained graph, chromosome files and alignments each time it is requested; `.done` files are completion markers, not content hashes.

Do not run local Snakemake, tests, syntax checks or validation on the SURF mount, as required by `whole_pans/AGENTS.md`.
