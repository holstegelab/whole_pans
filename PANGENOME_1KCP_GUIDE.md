# Chunked Cactus graphs from the existing Minigraph SV graphs

Updated 19 October 2026 after the user's clarification: **Cactus is the
base-resolution alignment engine**. The earlier seqwish/smoothxg implementation
was the wrong construction route for this goal and is replaced in this entrypoint.
The original saved plan's input inventory remains useful; its engine choice is
superseded. No local tests, syntax checks, Conda solves, Snakemake commands or
real-data feasibility scans have been run, as required by `AGENTS.md`. The
20 October failure investigation inspected existing server logs/reports and the
offending seed record; the code changes below still need server verification.

## Goal and source data

Reuse the completed Minigraph SV graphs and their assembly-to-node mappings for:

| New Cactus build | Frozen source build | Assemblies including references |
|---|---|---:|
| `local_1kcp_cactus_v1` | `local_v1_seeded_mc_smk_v2` | 196 local haplotypes + CHM13 + GRCh38 = 198 |
| `hprc_local_1kcp_cactus_v1` | `hprc_plus_local_v1_seeded_mc_smk_v2` | 660 cohort/HPRC haplotypes + 2 references = 662 |

The frozen manifests determine membership; these counts refer to FASTA/haplotype
assemblies, not diploid individuals. Source builds are beneath
`/gpfs/work3/0/qtholstg/pangenome/whole_pangenome/minigraph_cactus_reuse`.
New IDs keep any earlier seqwish results and old Cactus stores separate.

The 1KCP final graph used seqwish/smoothxg after partitioning. This implementation
borrows the partitioned-construction idea but uses Cactus inside each partition.
It is a hybrid method, not an exact reproduction of the paper's final graph.

## Entrypoint and resource choices

Run `Snakefile.1kcp`, loading `config/pangenome_1kcp.yaml` and
`workflow/rules/pangenome_1kcp.smk`. The existing `minigraph_cactus_reuse.smk`
remains available for its older chromosome-batch route; its original outputs and
job stores are untouched.

| Setting | Local-only | HPRC + local |
|---|---:|---:|
| Target reference span | 20 Mb | 5 Mb |
| Maximum component reference span | 30 Mb | 8 Mb |
| Maximum combined seed + assembly input | 4.5 billion bases | 4.5 billion bases |
| Boundary search slack | 1 Mb | 1 Mb |
| Threads | 16 | 16 |
| Outer Snakemake RAM request | 300,000 MB | 300,000 MB |
| Toil maximum memory | 290 billion bytes | 290 billion bytes |
| Cactus consolidated request | 280 billion bytes | 280 billion bytes |
| Heavy-stage scratch allowance | 500,000 MB | 500,000 MB |

The previous 40 Mb default was chosen for a different engine. Cactus's shipped
memory estimator is already approximately 256 GB at 4 GB combined sequence and
512 GB at 10 GB. A 40 Mb reference span across 662 assemblies alone is about
26.5 billion input bases. Therefore larger chunks cannot be assumed to fit
ordinary Genoa with Cactus. These smaller cohort-specific defaults and input
limits are starting points; representative pilots still have to establish fit.
The actual graph can contain more sequence than reference span times cohort size.

`--consMemory` overrides a resource estimate; it does not make Cactus use that
amount of RAM. `--maxMemory` and `--maxDisk` constrain Toil resource admission.
The outer compute allocation provides real enforcement. There are no automatic
RAM increases: both Snakemake retries and Toil retries are disabled in the launch
examples/rules. `--maxLen 10000` controls Cactus's internal recovery/alignment
intervals, not the size of a reference partition.

## Execution, step by step

1. **Inventory and freeze inputs.** Read the old `seqfile.txt`, input manifest,
   contig crosswalk and preflight. Record exact mapping graph and per-assembly
   PAF identities; check membership, lengths, names and frozen FASTA identities.
   Metadata have content hashes; large graph/PAF files use size and mtime.
   Original FASTA hashes are verified during extraction. Scientific settings,
   including the Cactus engine/version, are locked to the new build ID.

2. **Partition the exact existing mapping graph.** Read all segments and links,
   including non-reference topology. Use CHM13 rank-zero tags for positions.
   Accept an existing reference edge near a target only if it is a bridge of the
   complete undirected graph. Oriented parallel links are retained correctly.
   Node IDs and coordinates never change. Compute components after selected cuts;
   retain off-reference components and write ownership, boundaries and sizes.
   Accept Cactus-supported ambiguity symbols in the frozen graph. Derived seed
   FASTAs are uppercased and `UWSMKRYBDHV` becomes `N`, matching the pinned Cactus
   sanitizer; source GFA bytes and node lengths remain unchanged. The report
   records the base policy and counts affected seed bases/nodes.

3. **Gate infeasible partitions.** `partition/report.json` survives even when
   `feasible` is false. Multiple reference contigs in a component, excessive span
   or excessive seed sequence prevent readiness. Consecutive bubbles may have no
   suitable bridge: a smaller target cannot create a missing boundary. This case
   needs coordinate-aware node splitting or remapping and is explicitly stopped.

4. **Prepare assembly intervals.** Snakemake sorts each existing pre-deletion-
   filter PAF. Check query/target identities, lengths and CIGAR consumption.
   Group consecutive same-chunk visits, retaining gaps up to 1 Mb. Repeated visits
   after another chunk remain separate; conflicting cross-chunk overlaps stop.
   Query coordinates shift onto extracted intervals; CIGAR, target coordinates
   and strand stay intact. Preserve invariant CIGAR/filter annotations.

5. **Extract original sequence.** Stream each original FASTA once, checking
   contig lengths and its frozen uncompressed SHA256. Use original first-token
   identifiers through the crosswalk. Apply the same Cactus base normalization
   as for seed FASTAs; unsupported symbols still stop with an error. Check the
   frozen SHA256 on original bytes, and hash extracted intervals after
   normalization so final path checks use the sequences supplied to Cactus.
   Extraction statistics count masked ambiguity and uppercased bases across
   the entire original FASTA, including excluded regions.
   Originals and indexes are untouched. Excluded ends
   and long gaps are counted. Snakemake compresses interval FASTAs and PAF pieces.

6. **Rank and collect full-cohort chunk inputs.** `input_sizes.tsv` reports actual
   combined input bases, alignments and contributing haplotypes. It flags
   oversized and seed-only chunks. For an admitted chunk, write one FASTA per
   biological seqfile event plus the `_MINIGRAPH_` FASTA and a `seqfile.txt`.
   Translate PAF names into Cactus `id=EVENT|CONTIG` names using explicit
   crosswalks, without changing coordinates/CIGAR. Short internal contig IDs avoid
   sanitizer/subrange ambiguity. Keep a combined expected-sequence FASTA and
   `path_names.tsv` for downstream verification. Only present references enter
   `references.txt`, with CHM13 first. Seed-only chunks stop explicitly rather
   than creating a fictitious Cactus alignment.

7. **Run Cactus once per chunk.** The rule runs the native four-positional-argument
   interface: `cactus-align JOBSTORE SEQFILE PAF OUTPUT.hal --pangenome --outVG
   --outGFA`. Present references are passed dynamically; no reference argument is
   forced onto an off-reference component. Cactus builds its usual star tree;
   `_MINIGRAPH_` is a leaf, not `--root`. Each chunk has its own persistent Toil
   store, log and benchmark. Unique scratch is used via `--workDir`. A failed
   matching store is restarted; successful stores are cleaned by Toil.
   Each `{build}/{chunk}` is a separate Snakemake job and scheduler submission.
   Cactus's `single_machine` setting applies inside that chunk's allocation;
   it does not combine different chunks into one job. Postprocessing rules also
   run separately per chunk. Only shared preparation and final aggregation run
   once per build.

8. **Export and normalize the Cactus graph.** Cactus produces `cactus.hal`,
   `cactus.vg` and `cactus.gfa.gz` together from its alignment. Snakemake decompresses
   the GFA, then runs GFAffix separately with transformation checking. Native
   reference paths are protected from collapse. An explicit Cactus config keeps
   Minigraph paths for completeness checks and excludes ancestral paths.

9. **Verify and restore biological paths.** Translate native Cactus/hal2vg path
   names through the explicit alias table. Check every expected path's oriented
   links, length and sequence checksum. Remove artificial seed path records while
   retaining graph segments and links. Restore biological paths as GFA 1.1 `W`
   records with original sample, haplotype, contig and half-open interval
   coordinates. Densely renumber each chunk.

10. **Assign global IDs and report.** Once all chunks complete, offset their IDs
    cumulatively within each build. Write final per-chunk `graphs/ID.gfa`, graph
    and HAL/VG manifests, coverage accounting, benchmarks and scratch summaries.
    Combine the pieces into `graph/pangenome.gfa.gz` with one header and unique
    node IDs. This file is a collection of disconnected chunks; no boundary
    edges are invented by concatenation.
    A pilot completes through `chunks/ID/normalized.gfa`; global offsets require
    the complete build.

All coordinates in extraction and final GFA walks are **0-based, half-open**.
Per-chunk HALs contain extracted-fragment identifiers/coordinates; interval
sidecars provide the relation to full original contigs. No whole-genome HAL or
boundary connectivity is reconstructed. Graph products are disconnected pieces.
No read-k-mer filtering, whole-genome indexing, VCF or Giraffe products are added.
This is not the full clipping/indexing route of `cactus-graphmap-join`.

Tool execution, sorting, compression, scheduling and summaries stay in Snakemake.
Python helpers handle graph topology, frozen input/coordinate processing and path
formatting. No external orchestration wrapper runs the pipeline.

## Server commands

Use the established server Snakemake environment and your existing `zslurm`
profile. Snakemake submits each chunk rule as an independent job. No separate
Python commands are required; Snakemake runs the necessary helpers itself.

```bash
cd /gpfs/work3/0/qtholstg/pangenome/whole_pans
snakemake -s Snakefile.1kcp --profile ~/.config/snakemake/zslurm \
  --use-conda --jobs 8 --retries 0 --rerun-incomplete \
  pangenome_1kcp_prepare
```

Reuse the same command, changing its final target as follows:

| Target | Purpose |
|---|---|
| `pangenome_1kcp_prepare` | Inventory inputs and write partition feasibility reports. |
| `pangenome_1kcp_inputs` | Check feasibility automatically, prepare sequences/mappings and write `input_sizes.tsv`. |
| `pangenome_1kcp_pilot` | Build the `pilot_chunks` selected in the configuration. |
| `pangenome_1kcp_all` | Complete both graphs after `release_all_chunks: true` is enabled. |

Outputs are beneath the configured results root, separately for
`local_1kcp_cactus_v1` and `hprc_local_1kcp_cactus_v1`. Select pilot IDs from the
partition/input-size reports and measure typical, large and difficult chunks
before releasing all chunks. These are full-cohort pilots. A single chunk can
also be requested by its concrete `chunks/ID/normalized.gfa` path.

Snakemake resolves dependencies automatically: there is no need to run its
internal rules or Python helpers manually. Add `-n` for an optional dry run.

`--jobs 8` allows up to eight submitted workflow jobs at a time; it does not
group eight chunks into one job. Each Cactus chunk separately requests 16 threads
and 300,000 MB. Increase concurrency only to match the available worker pool.
Use the same profile command for the later targets. Adding `--executor local`
would override the profile and execute jobs in the current allocation.

Rules explicitly supply `resources.n`, integer-minute `time` and
`limit_auto_threads`. The profile's `compute` pool does not itself select Genoa
hardware; the manager/worker allocation must use the intended nodes.

Reuse the same targets with `--rerun-incomplete` after failure. A running
Snakemake process does not adopt updated rules. Do not resume old seqwish or
chromosome-batch Cactus stores with this entrypoint. Scratch defaults to shared
storage; free-space checks are not atomic reservations. Usage samples can miss
peaks, SIGKILL can leave temporary files, and retained HAL/VG/GFA inputs consume
substantial storage.

## 20 October partition failure

The local job `197240` stopped on node `s318398`, whose sequence is
`WGTCCATTCG` (`SN:Z:id=GRCh38|chr22`, `SO:i:16320468`, `LN:i:10`).
The previous ACGTN-only check rejected `W`. Seed and assembly extraction now
share the Cactus 3.3.0 ambiguity policy above; the derived sequence is
`NGTCCATTCG`, with the same length and coordinates.

The HPRC+local job `197241` completed its diagnostic report but was **infeasible**:
1 bridge, 0 selected cuts, 25 components, and a largest reference span of
248,387,328 bp against the 8 Mb limit. Its largest component has 321,829,372 seed
bases against the 120 Mb limit. Fixing ambiguity handling does not resolve this
topology limitation. Readiness remains blocked; a revised partition method
(with coordinate-aware mapping handling where needed) is required before
construction. Do not interpret successful `pangenome_1kcp_prepare` completion as
permission to construct chromosome-sized components. New partition logs print
`PASS` or `INFEASIBLE` explicitly.

Run the regression fixtures on the server, then rerun the preparation command
above to obtain the local build's report. No local tests were run.

```bash
cd /gpfs/work3/0/qtholstg/pangenome/whole_pans
python -m unittest discover -s tests -p 'test_one_kcp_sequences.py' -v
```

## Sources and verification boundary

- [Minigraph-Cactus v3.3.0 documentation](https://github.com/ComparativeGenomicsToolkit/cactus/blob/v3.3.0/doc/pangenome.md).
- [Cactus alignment interface and HAL/VG/GFA export](https://github.com/ComparativeGenomicsToolkit/cactus/blob/v3.3.0/src/cactus/setup/cactus_align.py).
- [Cactus resource estimates/configuration](https://github.com/ComparativeGenomicsToolkit/cactus/blob/v3.3.0/src/cactus/cactus_progressive_config.xml).
- [Cactus header sanitation](https://github.com/ComparativeGenomicsToolkit/cactus/blob/v3.3.0/src/cactus/preprocessor/checkUniqueHeaders.py).
- [Cactus 3.3.0 base masking table](https://github.com/ComparativeGenomicsToolkit/cactus/blob/v3.3.0/preprocessor/cactus_sanitizeFastaHeaders.c).
- [1KCP published construction workflow](https://github.com/JianYang-Lab/1KCP-analysis/blob/main/02_Pangenome/rules/graph_construction.smk).
- [GFAffix](https://github.com/codialab/GFAffix).

Pinned environment: Cactus 3.3.0, hal2vg 1.1.8, vg 1.76.1 and GFAffix 0.2.2.
The server Conda solve, input compatibility, real boundary feasibility, path
export behavior and resource peaks still need server verification. Static review
alone establishes neither successful execution nor guaranteed ordinary-node fit.
