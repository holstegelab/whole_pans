# Implementation plan: whole-genome Minigraph-Cactus

This plan supersedes `MINIGRAPH_CACTUS_PLAN_DRAFT.md` and incorporates the review in `MINIGRAPH_CACTUS_PLAN_REVIEW.md`. It describes one separately launched graph build. No workflow or FASTA validation was run on this mount.

## Scope and fixed choices

- Add an explicit `minigraph_cactus` Snakemake target, outside `rule all`. Keep the existing SV-only `pangenome` rule and its GFA independent. Do not pass that GFA to `--inGFA`; let Cactus construct its own Minigraph graph.
- Use the post-decontamination `graph_included_assemblies.txt` as the sole haplotype selection source. It has 196 paths in the mounted snapshot, but the code must accept a changed count. Do not use the older Mash-ordered TSV to order Cactus inputs. Keep a surviving `PARTIAL` sample's actual haplotype number.
- Put `CHM13` first and `GRCh38` second in the seqfile and `--reference CHM13 GRCh38`; use `--vcfReference CHM13`. The current CHM13 reference has 25 contigs, `chr1`–`chr22`, `chrX`, `chrY`, `chrM`; the current GRCh38 primary FASTA has the 24 nuclear contigs and no chrM or ALT loci. Report this intentional asymmetry in preflight; do not manufacture a GRCh38 mitochondrial contig. Verify reference versions/sources on the server before production.
- Build clipped GFA and GBZ, CHM13 VCF, and default `.d2` Giraffe indexes with `--gfa --gbz --vcf --giraffe --validate`. The default `.d2` graph is frequency-filtered and is distinct from the clipped `.gfa.gz`/`.gbz`. Do not combine `--haplo` with `--giraffe`. Leave `--filter`, GRCh38 VCF, ODGI, and extra mapping profiles for a later build ID.

## Paths and config

Use one new `minigraph_cactus:` section in `config/config.yaml`. `included_list` should point to the server-absolute post-QC list. Read reference FASTAs from the **top-level** `references.CHM13` and `references.hg38`, which are already server-absolute; do not use the relative `pangenome.references` values. The Snakemake launch directory is the server `workflow_root` (`/gpfs/work3/0/qtholstg/pangenome/whole_pans` in this config). Resolve only project-relative metadata/log paths against that configured root with lexical `pathlib` joining; never derive a server path from this mount or use `realpath`. Require absolute paths for every FASTA emitted to the seqfile, Cactus result directory, Toil job store, and `--workDir`. Check readability/writability on the server.

Suggested keys and types (values below are shapes, not calibrated site allocations):

```yaml
minigraph_cactus:
  build_id: chm13_196hap_v1       # unique for each changed input/setting set
  included_list: /gpfs/.../graph_included_assemblies.txt
  results_root: /gpfs/.../whole_pangenome/minigraph_cactus
  jobstore_root: /gpfs/.../whole_pangenome/minigraph_cactus_jobstores
  work_dir: /gpfs/.../toil_work   # absolute, visible to each Toil worker
  out_name: local_mc_chm13
  cactus_executable: cactus-pangenome
  cactus_version: '3.3.0'
  batch_system: single_machine    # or a server-supported Toil batch system
  run_mode: fresh                  # fresh | restart; operational, excluded from fingerprint
  driver_threads: 64
  driver_mem_mb: 524288
  driver_disk_mb: 3145728
  runtime_min: 10080
  max_cores: 64
  max_memory: 512GiB
  max_disk: 3TiB
  mg_cores: 64
  mg_memory: 256GiB
  map_cores: 6
  cons_cores: 64
  cons_memory: 256GiB
  index_cores: 63
  index_memory: 256GiB
```

Do not assume those allocations fit Snellius. Set them after a server pilot. For `single_machine`, Snakemake must reserve enough cores, RAM, and disk for **all concurrent** Toil jobs, so `driver_*` must cover `max_*` and worker peaks. With scheduler-backed Toil, `driver_*` reserves only a long-lived leader able to submit jobs, while Toil requests workers separately; keep scheduler flags in the server profile, not the Snakefile. `work_dir` capacity is separate from the shared job store and final results. The chosen execution mode and all resource knobs belong in the build fingerprint; `run_mode` does not, because changing `fresh` to `restart` must not change the build identity. Use `--clean=onSuccess`: a failed store remains restartable, while a successful run may remove the store before the marker is published.

Derive `BUILD = results_root/build_id`, `OUTDIR = BUILD/graph`, `METADATA = BUILD/metadata`, `JOBSTORE = jobstore_root/build_id`; keep logs and benchmark under `BUILD`. Require a new `build_id` (and hence fresh `OUTDIR`/`JOBSTORE`) when inputs, version, or settings change. Disallow path reuse with a different fingerprint.

## Files to implement and rule interfaces

1. Add `workflow/scripts/prepare_minigraph_cactus_seqfile.py` with CLI `--included-list`, `--chm13`, `--grch38`, `--seqfile`, `--input-manifest`, `--contig-map`, `--preflight-json`. It writes all outputs atomically after successful validation. `seqfile.txt` has two whitespace-delimited columns: `CHM13 <absolute path>`, `GRCh38 <absolute path>`, then one row per accepted haplotype. Derive `SAMPLE.1`/`SAMPLE.2` from filenames of form `SAMPLE.hifi.hifiasm.bp.hap[12].p_ctg.clean.fa[.gz]`; reject unrecognized names, duplicate labels/paths, whitespace in paths, missing/empty input, and an empty inclusion list. `input_manifest.tsv` records seqfile name, sample, haplotype, source path, source list line, byte size, nanosecond modification time, and SHA-256 of uncompressed FASTA content. `contig_map.tsv.gz` records original header and Cactus-effective contig ID for each input. `preflight.json` records the inclusion-list hash plus input hashes, seqfile hash, counts, reference contig sets, and preflight version. Stream each FASTA once on the server, hashing its decoded bytes while checking headers. Apply the documented name handling before collision checks: strip through the last `#`, ignore text after the first space, preserve a valid `:start-end` subrange, replace other `:` with `_`; reject empty or duplicate effective IDs within one FASTA. Also reject malformed FASTA/compression. No contig rewrite is in scope. If this build ID already has state or a success marker, compare newly computed preflight identity before replacing any metadata; refuse a changed input set and tell the user to choose a new build ID. The full scan is a one-time cached stage, not part of every restart.
2. Add `workflow/scripts/run_minigraph_cactus.py` with CLI `--preflight-json`, `--seqfile`, `--settings-json` (canonical JSON **text** of current config values), `--outdir`, `--out-name`, `--jobstore`, `--work-dir`, `--marker`, `--log`, `--mode fresh|restart`. It records the exact argument vector (subprocess list, no shell interpolation), Cactus version, SHA-256 of canonical settings and preflight, and resulting build fingerprint in `settings.json`/`state.json` under `METADATA`. Exclude `run_mode` from the canonical settings. Its sole Snakemake output is `METADATA/success.json`; use atomic rename only after Cactus success and output checks. Write state before launch, update it after the process exits, and retain it on failure. Verify that the executable reports the configured version before launching. Use an exclusive build lock under `METADATA` so two leaders cannot use one store.
3. Add `workflow/scripts/verify_minigraph_cactus.py` with CLI `--marker`, `--preflight-json`, `--settings-json` (the **current** canonical JSON text from Snakemake config), `--outdir`, `--out-name`, `--check-input-stat`. It recomputes the settings/build fingerprint, compares it with the marker, checks listed output files are nonempty, and checks path/size/mtime for all inputs against preflight (fast metadata checks, not a full FASTA reread). It reports changed inputs and asks for a new build ID and preflight if anything differs. This script is called by a no-output target on **every** requested `minigraph_cactus` invocation, including when `success.json` is already current.
4. Add `workflow/rules/minigraph_cactus.smk`, include it once from the root `Snakefile`, and keep it usable directly via `--snakefile workflow/rules/minigraph_cactus.smk`. Follow this repo's `configfile: workflow.source_path("../../config/config.yaml")` and server `workdir` pattern, with distinct `MC_*` constants and rule names. No duplicated rule definitions. `prepare_minigraph_cactus` tracks the inclusion list and both reference FASTAs as inputs and declares `seqfile.txt`, `input_manifest.tsv`, `contig_map.tsv.gz`, `preflight.json` as outputs. `build_minigraph_cactus` takes those outputs plus the same list/references as inputs, passes a `json.dumps(..., sort_keys=True, separators=(",", ":"))` settings string through `params`, declares only `success.json` as output, and sets `threads`, `resources` (`mem_mb`, `disk_mb`, `runtime_min`), `log`, and `benchmark`. The build script persists that settings string as `METADATA/settings.json`. `rule minigraph_cactus` takes `success.json` as input, has no output, and passes a newly rendered settings string from the current config to the verifier every time. If the server Snakemake version does not support standalone `workdir` exactly as expected, retain the single source file and document root Snakefile invocation as the supported entry point.
5. Add `workflow/envs/minigraph_cactus.yaml` with a Cactus 3.3.0 pin and minimal Python dependencies. The Cactus launcher may instead be an official release/container command configured as `cactus_executable`; verify binaries and version on the server. Keep this separate from `workflow/envs/pangenome.yaml`.

The build command is one `cactus-pangenome` invocation, with no scheduler-specific shell code:

```text
cactus-pangenome JOBSTORE SEQFILE --outDir OUTDIR --outName local_mc_chm13
  --reference CHM13 GRCh38 --vcfReference CHM13
  --gfa --gbz --vcf --giraffe --validate
  --batchSystem BATCH_SYSTEM --workDir WORK_DIR --clean=onSuccess
  --maxCores N --maxMemory SIZE --maxDisk SIZE
  --mgCores N --mgMemory SIZE --mapCores N
  --consCores N --consMemory SIZE --indexCores N --indexMemory SIZE
```

Use version-checked options from the installed Cactus/Toil before the pilot. Avoid a generic arbitrary flags string in config: it could change outputs without changing the declared checks. For an alternate container launcher, document its exact executable and version in settings.

## Expected artifacts and output checks

For the selected switches, require nonempty `OUTDIR/local_mc_chm13.gfa.gz`, `.gbz`, `.vcf.gz`, `.vcf.gz.tbi`, `.raw.vcf.gz`, `.raw.vcf.gz.tbi`, `.d2.gbz`, `.d2.dist`, `.d2.min`, `.sv.gfa.gz`, `.full.hal`, `.input-contig-sizes.tsv.gz`, `.stats/clipped-by-genome.tsv`, and `.stats/graph-stats.tsv`. The actual `out_name` substitutes for `local_mc_chm13`. Keep the Cactus output tree outside declared Snakemake outputs so an incomplete Snakemake job cannot delete a large partial graph. The marker lists exact paths, sizes, and SHA-256 of the small metadata/settings files; checking huge graph files for nonzero size is sufficient for routine verification. A `.WARNING` file is optional and must be reviewed if present. The server post-build QC also examines input-contig accounting, clipping statistics, reference `W` walks, and VCF/index usability; a successful exit and missing warning alone are not sufficient QC.

## Restart and adoption state machine

- **Fresh:** require no existing marker, unused/empty `OUTDIR` and nonexistent `JOBSTORE`. Create a locked state file with fingerprint before launching. Refuse occupied paths from another fingerprint. `--clean=onSuccess` removes the Toil store after successful Cactus execution; write a `cactus_completed` state immediately on exit zero, verify outputs, then atomically publish `success.json`.
- **Failed/interrupted:** retain partial `OUTDIR`, state and job store. Set `run_mode: restart` only after inspecting the server log. The wrapper requires a matching fingerprint, a prior `running`/`failed` state, an existing store, and no success marker, then invokes the **same** command with `--restart`. Let Toil reject a store that is not actually an unfinished workflow. Never decide to restart solely because a directory exists or because Snakemake used `--rerun-incomplete`.
- **Cactus completed, marker absent:** if `state.json` says `cactus_completed` and all required files pass verification, adopt the result and publish the marker without relaunching Cactus, regardless of whether `--clean=onSuccess` removed the store. Also allow adoption after an interrupted marker write if the log records successful Cactus exit and matching fingerprint; do not infer success solely from files. If any required file is missing, stop. A new build ID is needed if there is no restartable store.
- **Marker present:** the no-output target verifies outputs and fingerprint on every request. A changed config, changed input metadata, or missing graph file is an error; it must not silently reuse the marker or overwrite the old graph. All state and output writes stay under the configured build directory.

## Server-only sequence

1. On the server, confirm Cactus 3.3.0 command availability and `cactus-pangenome --help` option names, reference provenance, FASTA access, shared store/work directory capacity, and execution mode. Set site-calibrated resources in config.
2. Run `prepare_minigraph_cactus` for a **separate pilot build ID and a small explicit pilot inclusion list** containing two or three passing haplotypes plus both references. Review seqfile labels, 25/24 reference contig report, normalized-name collisions, and manifest hashes. Keep pilot `results_root`/`jobstore_root` paths isolated from production.
3. Run the pilot build and verifier; inspect exact output names, GFA `W` lines, `.WARNING` if present, clipped-by-genome, graph-stats, and input-contig-sizes. Then use a production build ID, restore the authoritative inclusion list, run `snakemake -n minigraph_cactus` on the server, and launch with its server profile and `--use-conda`. The production job should be a deliberate explicit target, never `rule all`.
4. If interrupted, inspect state/log/store on the server, set `run_mode: restart`, and rerun the same explicit target. After success, set the mode back to `fresh` only for a *new* build ID. Do not delete an existing job store or partial output automatically.

Sources: [Cactus pangenome interface, naming and output documentation](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md); [Toil job store, work directory, cleanup and restart options](https://toil.readthedocs.io/en/releases-9.0.0/running/cliOptions.html); local `pangenome-papers` construction/tool notes.
