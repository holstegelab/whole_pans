# Running the novel-SV discovery workflow

This guide turns the analysis described in
[`NOVEL_SV_DISCOVERY_PLAN.md`](NOVEL_SV_DISCOVERY_PLAN.md) into a practical
server runbook for the implemented Snakemake workflow. The first pass freezes
the existing Minigraph graph, screens all cleaned haplotype assemblies against
it, calls assembly-to-reference SVs, builds a provisional assembly-SV catalog,
and ranks samples for HiFi-read retrieval.

The workflow entry point is
[`workflow/rules/novel_sv_discovery.smk`](workflow/rules/novel_sv_discovery.smk).
Its settings are under `novel_sv_discovery` in
[`config/config.yaml`](config/config.yaml). Graph screening, reference calling,
catalog processing, and optional PAV execution use separate Conda environments
under [`workflow/envs/`](workflow/envs/).

## Important scope

The implemented target is the assembly-only first pass from the plan:

```text
frozen graph inventory
        |
        +--> cleaned-assembly manifest --> graph residual screen --+
        |                                                        |
        +--> CHM13/hg38 assembly callers ------------------------+--> catalog
                                                                  |
                                                                  +--> HiFi ranking
```

It does **not** yet perform read-based discovery, cohort-wide read genotyping,
or graph-v2 construction. Those stages require the HiFi data and a reviewed
set of alleles. The generated catalog therefore remains provisional and uses
`ASSEMBLY_ONLY_REVIEW` or `PENDING_HIFI` validation states.

PAV is implemented as a third assembly caller but is disabled in the default
configuration until its source workflow is checked out on the server. With the
current defaults, the reference route runs diploid SVIM-asm and dipcall for
samples that have both haplotypes. Unpaired samples are still screened against
the graph.

## 1. Work from the server project root

Run the commands below from the directory that contains `whole_pans/`,
`whole_pangenome/`, `refs/`, and `tools/`:

```bash
cd /path/on/the/server/to/pangenome
test -d whole_pans
test -d whole_pangenome
```

Do not replace configuration paths with `realpath` output from the mounted
machine. Snakemake sets its working directory to `whole_pans/`, so paths in
`whole_pans/config/config.yaml` are intentionally relative to
`whole_pans/Snakefile`. For example:

```yaml
novel_sv_discovery:
  results: ../whole_pangenome/sv_pangenome/novel_sv_discovery
  frozen_graph: ../whole_pangenome/sv_pangenome/results/graphs/sv_pangenome.minigraph.gfa
  cleaned_manifest: ../whole_pangenome/assembly_qc_decontaminated/results/resources/all_cleaned_assemblies.tsv
  references:
    CHM13: ../refs/CHM13.fa.gz
    hg38: ../refs/GRCh38_primary_chr1-22_XY.fa
```

These resolve on the execution server even when the same repository has a
different absolute path through this mount.

## 2. Check the required inputs

The first pass depends on the frozen graph, its build provenance, the cleaned
assembly manifest, both QC tables, the contamination summary, the previous
feasibility table, and both references. Check them from `whole_pans/` so the
same relative paths are used as at runtime:

```bash
cd whole_pans

for path in \
  ../whole_pangenome/sv_pangenome/results/graphs/sv_pangenome.minigraph.gfa \
  ../whole_pangenome/sv_pangenome/results/metadata/sv_pangenome.ordered_assemblies.tsv \
  ../whole_pangenome/sv_pangenome/results/metadata/tool_versions.tsv \
  ../whole_pangenome/sv_pangenome/results/metadata/sv_pangenome.graph_summary.tsv \
  ../whole_pangenome/sv_pangenome/results/logs/build_sv_pangenome_graph.log \
  ../whole_pangenome/sv_pangenome/qc_analysis/report/pangenome_qc_report.md \
  ../whole_pangenome/assembly_qc_decontaminated/results/resources/all_cleaned_assemblies.tsv \
  ../whole_pangenome/assembly_qc/results/summary/assembly_qc.tsv \
  ../whole_pangenome/assembly_qc_decontaminated/results/summary/assembly_qc.tsv \
  ../whole_pangenome/assembly_decontamination/results/summary/contamination_summary.tsv \
  ../whole_pangenome/sv_pangenome/qc_analysis/tables/excluded_haplotype_feasibility.csv \
  ../refs/CHM13.fa.gz \
  ../refs/GRCh38_primary_chr1-22_XY.fa
do
  test -s "$path" || printf 'MISSING OR EMPTY\t%s\n' "$path"
done

cd ..
```

No output means every listed file exists and is non-empty. The workflow also
checks that every cleaned assembly name follows the current convention:

```text
SAMPLE.hifi.hifiasm.bp.hap1.p_ctg
SAMPLE.hifi.hifiasm.bp.hap2.p_ctg
```

Diploid SVIM-asm, dipcall, and PAV require both haplotypes. Samples missing one
haplotype are intentionally omitted from those callers but remain eligible for
the frozen-graph screen; they do not make the workflow fail.

## 3. Review the discovery configuration

The current configuration is ready for the known server layout. The main
settings to review before a production run are:

```yaml
novel_sv_discovery:
  sample_metadata: ""
  read_manifest: ""

  graph_screen:
    batch_size: 2
    threads: 16
    mem_mb: 16000
    runtime: 720
    primary_min_chain_score: 5000
    primary_secondary: 5
    rescue_callable_fraction: 0.85
    sensitivity_min_chain_score: 1000
    sensitivity_secondary: 20
    min_sv_size: 50
    raw_min_size: 30
    adjacent_indel_gap: 50
    min_alignment: 5000
    min_anchor: 2000
    min_mapq: 5
    min_identity: 0.90

  reference_calling:
    threads: 16
    alignment_preset: asm5
    sort_mem_mb_per_thread: 2048
    alignment_mem_mb: 64000
    alignment_runtime: 360
    caller_mem_mb: 64000
    caller_runtime: 720
    sex_chromosomes:
      mode: autosomes_only
      exclude_contig_regex: "(?:chr)?[XY]"
      par_bed:
        CHM13: ""
        hg38: ""
      chrY_par_hard_masked:
        CHM13: false
        hg38: false
    callers:
      svim_asm: true
      dipcall: true
      pav: false

  clustering:
    breakpoint_distance: 500
    length_similarity: 0.70

  hifi_recommendations:
    callable_threshold: 0.85
    control_count: 20
```

`sample_metadata` and `read_manifest` are optional. Leave them empty for the
assembly-only pass. If supplied later, their paths must also remain relative to
`whole_pans/`. The workflow rejects absolute paths and paths that escape the
pangenome project tree. Expected optional columns are documented next to those
keys in the config file.

`autosomes_only` is the safe default until sex, paternal haplotype origin, and
reference-specific pseudoautosomal-region BED files are available. Raw X/Y
records remain in caller VCFs, but the merger excludes them from the catalog.
To use `sex_aware`, provide `sample_id`, `sex`, and, for male samples,
`paternal_haplotype` (`hap1` or `hap2`) plus both PAR BED paths. Dipcall then
receives paternal sequence first and the reference-specific PAR intervals.
Dipcall also requires PAR sequence on reference chrY to be hard-masked: point
the configured references at verified masked FASTAs and set the corresponding
`chrY_par_hard_masked` acknowledgements to `true`. In `sex_aware` mode, the
merger retains X/Y records only from PAR-aware dipcall; it conservatively drops
X/Y records from SVIM-asm and PAV.

The graph screen uses two Minigraph passes. Every assembly receives the primary
pass. An assembly with a primary callable fraction below `0.85` receives the
more sensitive pass. The initial event parser retains operations from 30 bp so
that nearby indels can be joined, but only events of at least 50 bp enter the
SV outputs.

## 4. Prepare Conda environments

Use the existing Snakemake driver environment, then let `--use-conda` create
the four purpose-specific rule environments without running jobs:

```bash
snakemake \
  --snakefile whole_pans/workflow/rules/novel_sv_discovery.smk \
  --use-conda \
  --conda-create-envs-only \
  --cores 1 \
  novel_sv_discovery
```

The environments are intentionally separated:

- `sv_graph.yaml`: Python, Minigraph, and gfatools;
- `sv_reference.yaml`: minimap2, samtools, bcftools/htslib, SVIM-asm 1.0.3,
  dipcall 0.3, and make;
- `sv_catalog.yaml`: lightweight Python-only manifest/catalog processing;
- `sv_discovery.yaml`: PAV's broader native dependencies.

The last environment does not install PAV itself because PAV is a source
Snakemake workflow.

For the full three-caller union specified in the plan, check out the pinned PAV
source into the path expected by the config:

```bash
git clone --recursive --branch v2.4.6 \
  https://github.com/EichlerLab/pav.git \
  tools/pav
```

Then enable it in `whole_pans/config/config.yaml`:

```yaml
reference_calling:
  callers:
    svim_asm: true
    dipcall: true
    pav: true
  pav_snakefile: ../tools/pav/Snakefile
```

Keep `pav: false` for a staging run if the checkout is not available. Do not
describe the resulting two-caller catalog as the complete planned multi-caller
union.

## 5. Run lightweight checks before real data

Run the Python tests from the project root:

```bash
python -m unittest discover \
  -s whole_pans/tests \
  -p 'test_*.py'
```

This command uses Python's standard `unittest` discovery mode. `-s` selects the
test directory and `-p` limits discovery to files named `test_*.py`, matching
the repository's current test layout. Run it before any large Snakemake target
because it catches parser, manifest, and catalog logic regressions without
touching the large assembly inputs.

The novel-SV tests cover rGFA indexing, GAF CIGAR parsing, adjacent-indel
handling, anchor checks, split alignments, deterministic task creation, and
catalog clustering.

Next, ask Snakemake to build the production DAG without executing it:

```bash
snakemake \
  --snakefile whole_pans/workflow/rules/novel_sv_discovery.smk \
  --use-conda \
  --cores 1 \
  --dry-run \
  --printshellcmds \
  novel_sv_discovery
```

A dry-run validates workflow structure and config paths visible on the server;
it does not validate the biological content of large inputs.

## 6. Freeze and inventory the current graph

This stage does not modify the GFA. It records the graph size and checksum,
counts segments and links, captures tool versions, extracts existing bubbles,
and copies the ordered input list and build/QC provenance into the discovery
results.

```bash
snakemake \
  --snakefile whole_pans/workflow/rules/novel_sv_discovery.smk \
  --use-conda \
  --cores 1 \
  inventory_frozen_graph_for_novel_sv \
  index_frozen_graph_for_novel_sv
```

Inspect the recorded identity before continuing:

```bash
RESULTS=whole_pangenome/sv_pangenome/novel_sv_discovery

column -t -s $'\t' "$RESULTS/frozen_graph/graph_inventory.tsv"
head "$RESULTS/frozen_graph/tool_versions.tsv"
sha256sum whole_pangenome/sv_pangenome/results/graphs/sv_pangenome.minigraph.gfa
```

The checksum printed by `sha256sum` should match the `sha256` row in
`graph_inventory.tsv`. Keeping this graph frozen is essential: a residual
allele is only meaningful relative to a fixed target.

## 7. Build and inspect the discovery manifest

The manifest builder joins every cleaned assembly to graph membership,
original and post-decontamination QC, contamination results, prior feasibility,
and optional sample/read metadata. It also refreshes rescue tiers and writes a
provisional HiFi request list based on QC and graph membership alone.

```bash
snakemake \
  --snakefile whole_pans/workflow/rules/novel_sv_discovery.smk \
  --use-conda \
  --cores 1 \
  build_novel_sv_discovery_manifest
```

Summarize the manifest without loading assembly FASTAs:

```bash
python - <<'PY'
import csv
from collections import Counter
from pathlib import Path

path = Path("whole_pangenome/sv_pangenome/novel_sv_discovery/manifest/assembly_discovery_manifest.tsv")
with path.open() as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

print("assemblies", len(rows))
print("samples", len({row["sample_id"] for row in rows}))
print("graph_member", Counter(row["graph_member"] for row in rows))
print("rescue_tier", Counter(row["rescue_tier"] for row in rows))
print("priority", Counter(row["discovery_priority"] for row in rows))
PY
```

The `python - <<'PY'` form runs the following inline Python program from
standard input. The quoted `PY` delimiter keeps the shell from expanding values
inside the snippet before Python sees them. This snippet reads the
tab-separated manifest with `csv.DictReader`, treats each row as a dictionary
keyed by the header names, then prints cohort-level counts for assemblies,
unique samples, graph membership, rescue tier, and discovery priority. It is a
metadata check only; it does not open any assembly FASTA files.

For the current cohort, the expected high-level checks are 982 assemblies from
491 samples, including 196 graph members and 786 assemblies outside the graph.
Treat a mismatch as a reason to review the upstream manifests before starting
hundreds of mappings.

The key files from this stage are:

```text
manifest/assembly_discovery_manifest.tsv
manifest/refreshed_excluded_haplotype_feasibility.tsv
hifi/provisional_hifi_samples.tsv
```

`discovery_priority` controls pilot ordering, not eligibility. All assemblies
are eventually screened. Assemblies already used in the graph are retained as
negative controls for mapping ambiguity and representation failures.

## 8. Create deterministic graph-screen tasks

Task creation groups assemblies according to `graph_screen.batch_size`. The
default is two assemblies per task to amortize each graph load while keeping
failures and reruns reasonably isolated.

```bash
snakemake \
  --snakefile whole_pans/workflow/rules/novel_sv_discovery.smk \
  --use-conda \
  --cores 1 \
  create_novel_sv_graph_tasks

RESULTS=whole_pangenome/sv_pangenome/novel_sv_discovery
head "$RESULTS/graph_screen/tasks.tsv"
```

The task identifiers are stable four-digit values such as `0001`. The
checkpoint reads these IDs from `tasks.tsv`; downstream expansion does not
independently calculate a task count. Increasing `batch_size` reduces repeated
graph loads but puts more work and more failure scope into each job. Change it
only after a pilot has measured memory, runtime, and filesystem pressure.

## 9. Calibrate a graph-screen pilot

Before launching the full screen, select a small pilot containing:

- 5--10 graph-member assemblies;
- 10--20 missing-mate or `best_rescue` assemblies;
- several `fragmented_rescue` or `not_recommended` assemblies;
- positive controls or leave-one-out alleles when available.

Snakemake creates the inputs and task table reproducibly. To run one selected
task manually for calibration, change into `whole_pans/` and use the same
script and paths as the rule:

```bash
cd whole_pans
RESULTS=../whole_pangenome/sv_pangenome/novel_sv_discovery
TASK_ID=0001

python workflow/scripts/screen_novel_graph_svs.py run \
  --graph ../whole_pangenome/sv_pangenome/results/graphs/sv_pangenome.minigraph.gfa \
  --segment-index "$RESULTS/frozen_graph/graph_segments.tsv.gz" \
  --manifest "$RESULTS/manifest/assembly_discovery_manifest.tsv" \
  --graph-assemblies ../whole_pangenome/sv_pangenome/results/metadata/sv_pangenome.ordered_assemblies.tsv \
  --tasks "$RESULTS/graph_screen/tasks.tsv" \
  --output-dir "$RESULTS/graph_screen/pilot_tasks/$TASK_ID" \
  --task-id "$TASK_ID" \
  --threads 16

cd ..
```

This Python command invokes the graph-screening helper directly for one task.
The `run` subcommand reads the frozen graph, the precomputed segment index, the
assembly discovery manifest, the original graph-assembly list, and the
deterministic `tasks.tsv` file. `--task-id` selects a single batch from
`tasks.tsv`, `--output-dir` points at a pilot-only destination, and `--threads`
sets the Minigraph thread count. Use this direct command for calibration and
debugging; use the Snakemake target for production so logs, resources,
benchmarks, Conda activation, and completion markers remain managed by the
workflow.

For scheduler execution, prefer the Snakemake task output so resources,
logging, Conda activation, benchmarking, and completion markers are retained:

```bash
snakemake \
  --snakefile whole_pans/workflow/rules/novel_sv_discovery.smk \
  --use-conda \
  --cores 16 \
  ../whole_pangenome/sv_pangenome/novel_sv_discovery/graph_screen/tasks/0001
```

Each task is published atomically as one declared directory. The directory is
created only after all expected non-empty outputs have been checked and it
contains `.complete`, `task_outputs.tsv`, mappings, parsed tables, summaries,
and logs for every assembly in the batch. A failed temporary directory is kept
with a `.failed` suffix for diagnosis. A manual pilot output directory must not
already exist; use Snakemake to manage production reruns.

Review the pilot’s per-assembly summaries, logs, and Snakemake benchmarks:

```bash
RESULTS=whole_pangenome/sv_pangenome/novel_sv_discovery

find "$RESULTS/graph_screen/tasks/0001" -maxdepth 3 -type f | sort | head -40
column -t -s $'\t' "$RESULTS/benchmarks/graph_screen/0001.tsv"
```

The graph is several gigabytes and each Minigraph process loads its own copy.
Use observed peak memory to set concurrency; do not infer safe concurrency from
the number of CPUs alone.

## 10. Run and interpret the frozen-graph screen

Run all graph tasks and combine their results:

```bash
snakemake \
  --snakefile whole_pans/workflow/rules/novel_sv_discovery.smk \
  --use-conda \
  --cores 32 \
  --rerun-incomplete \
  --printshellcmds \
  summarize_novel_sv_graph_screen
```

On a cluster, add the site’s profile rather than embedding scheduler directives
in the workflow:

```bash
snakemake \
  --profile /path/to/the/server/snakemake-profile \
  --snakefile whole_pans/workflow/rules/novel_sv_discovery.smk \
  --use-conda \
  --rerun-incomplete \
  summarize_novel_sv_graph_screen
```

The high-confidence tier requires a primary alignment of at least 5 kb, at
least 90% identity, MAPQ at least 5, and at least 2 kb of sequence anchoring
both sides of the event. Other residuals are retained for review rather than
silently discarded. Split alignments, orientation changes, and terminal events
without two-sided support are written to the complex/review outputs.

The aggregate outputs are:

```text
graph_screen/summary/all_assembly_novel_sv_summary.tsv
graph_screen/summary/all_residual_sv_candidates.tsv.gz
graph_screen/summary/all_complex_alignments.tsv.gz
graph_screen/summary/all_contig_alignments.tsv.gz
```

The aggregator streams these tables instead of loading cohort-wide output into
memory. It verifies that its checkpoint-derived task-directory set exactly
matches `tasks.tsv` and fails on any missing or empty per-assembly file.

Before catalog construction, the workflow reannotates the saved residual table
using the rGFA segment index. This converts Minigraph path labels such as
`chr15:start-end`, as well as bare rank-0 paths such as `chr5`, into stable
source coordinates without rerunning Minigraph:

```bash
snakemake \
  --snakefile whole_pans/workflow/rules/novel_sv_discovery.smk \
  --use-conda --cores 1 \
  reannotate_novel_sv_graph_coordinates
```

The corrected candidates are written to
`graph_screen/coordinate_qc/all_residual_sv_candidates.annotated.tsv.gz`.
`coordinate_qc.tsv` records the resolved fraction, and the rule fails if more
than `graph_screen.max_unresolved_coordinate_fraction` of >=50 bp candidates
lack a graph segment and offset.

Summarize callability and evidence tiers:

```bash
python - <<'PY'
import csv
import gzip
from collections import Counter
from pathlib import Path

root = Path("whole_pangenome/sv_pangenome/novel_sv_discovery/graph_screen/summary")
with (root / "all_assembly_novel_sv_summary.tsv").open() as handle:
    summary = list(csv.DictReader(handle, delimiter="\t"))
with gzip.open(root / "all_residual_sv_candidates.tsv.gz", "rt") as handle:
    candidates = list(csv.DictReader(handle, delimiter="\t"))

callable_fraction = [
    max(float(row.get("primary_callable_fraction") or 0),
        float(row.get("sensitivity_callable_fraction") or 0))
    for row in summary
]
print("assemblies", len(summary))
print("below_0.85_callable", sum(value < 0.85 for value in callable_fraction))
print("candidate_tiers", Counter(row["confidence_tier"] for row in candidates))
print("candidate_types", Counter(row["svtype"] for row in candidates))
PY
```

This inline Python check loads the per-assembly graph-screen summary as TSV and
the residual-candidate table through `gzip.open` because it is compressed. For
each assembly it compares the primary and sensitivity callable fractions and
uses the larger value as the effective callable fraction. The printed counters
show how many assemblies remain below the configured callability threshold and
how graph residual candidates are distributed by confidence tier and SV type.

A clean end-to-end alignment with no residual event supports graph
representability only in callable sequence. A no-call in repetitive,
low-identity, low-MAPQ, fragmented, or otherwise uncallable sequence is not
evidence that the allele is absent.

## 11. Run the assembly-to-reference callers

The independent reference route catches events that may not be expressed as a
single residual indel against the rGFA. It runs against both CHM13 and GRCh38:

- SVIM-asm consumes the two haplotype BAMs in documented diploid mode;
- dipcall compares each phased haplotype pair and emits a VCF plus callable
  BED;
- PAV runs each phased pair through its nested source workflow when enabled.

All three callers skip samples without a complete haplotype pair. Before
calling, each configured reference is copied to a writable, uncompressed
`reference.fa`, then indexed once as `.fai`, `.dict`, and `.mmi`. Minimap2 and
dipcall reuse the same prepared `.mmi`; the original reference is never
modified. `samtools sort` receives the configured per-thread memory explicitly.

Run all enabled callers and write their input manifest:

```bash
snakemake \
  --snakefile whole_pans/workflow/rules/novel_sv_discovery.smk \
  --use-conda \
  --cores 32 \
  --rerun-incomplete \
  write_novel_sv_reference_call_manifest
```

Reference alignments and callers request up to 64 GB by default. A cluster
profile should translate the rule’s `threads`, `mem_mb`, and `runtime`
resources into the local scheduler settings.

Outputs are separated by coordinate system and caller:

```text
reference_calls/CHM13/svim_asm/
reference_calls/CHM13/dipcall/
reference_calls/CHM13/pav/
reference_calls/CHM13/reference.fa{,.fai}
reference_calls/CHM13/reference.dict
reference_calls/CHM13/reference.mmi
reference_calls/hg38/svim_asm/
reference_calls/hg38/dipcall/
reference_calls/hg38/pav/
reference_calls/hg38/reference.fa{,.fai}
reference_calls/hg38/reference.dict
reference_calls/hg38/reference.mmi
reference_calls/reference_call_manifest.tsv
```

CHM13, GRCh38, and graph coordinates remain explicitly separate. Reference
calls are normalized evidence, but they cannot create catalog events by
themselves. CHM13 calls may validate a coordinate-compatible graph residual;
GRCh38 calls remain provenance evidence until an explicit liftover/reconciliation
step is added. This provisional table does not report allele frequencies; those
require a reconciled, genotyped catalog in one coordinate system.

With the default `autosomes_only` mode, the raw caller VCFs still contain their
X/Y records but known X/Y records are filtered from the merged catalog.
Unprojected graph events whose chromosome of origin cannot be established are
retained with `sex_chromosome_status=UNRESOLVED_GRAPH_ORIGIN`; do not assume
those rows are autosomal. Do not switch to `sex_aware` until complete sex,
paternal-haplotype, reference-specific PAR metadata, and verified chrY PAR
hard-masking pass workflow validation.

## 12. Merge the provisional assembly-SV catalog

The merger normalizes graph residuals and enabled reference-caller VCFs into a
common evidence table, sorts them through a temporary SQLite database, and
clusters only compatible records in the same coordinate system and chromosome.
It emits a catalog row only when a nonmember assembly has graph-residual
evidence. Reference-only clusters and graph-member-only self-alignment
residuals are omitted. A graph-member residual overlapping a nonmember cluster
is retained as control evidence and forces review rather than increasing
confidence.
The default cluster criteria are:

- same SV type;
- breakpoints within 500 bp;
- length ratio of at least 0.70;
- exact sequence agreement when both insertion alleles are sequence-resolved.

Point `TMPDIR` to sufficiently large server or node-local scratch before this
stage. The temporary SQLite database is deleted after a successful or failed
run, while the normalized evidence output is retained.

```bash
export TMPDIR=/path/to/large/scratch/for/novel-sv
mkdir -p "$TMPDIR"

snakemake \
  --snakefile whole_pans/workflow/rules/novel_sv_discovery.smk \
  --use-conda \
  --cores 1 \
  --rerun-incomplete \
  merge_assembly_only_novel_sv_catalog
```

The outputs are:

```text
catalog/provisional_graph_residual_sv_catalog.tsv
catalog/all_normalized_assembly_evidence.tsv.gz
```

Inspect confidence, validation state, method support, and coordinate systems:

```bash
python - <<'PY'
import csv
from collections import Counter
from pathlib import Path

path = Path("whole_pangenome/sv_pangenome/novel_sv_discovery/catalog/provisional_graph_residual_sv_catalog.tsv")
with path.open() as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

for field in (
    "coordinate_system",
    "sex_chromosome_status",
    "svtype",
    "confidence",
    "validation_status",
    "graph_representation_status",
):
    print(field, Counter(row[field] for row in rows))
PY
```

This final inline Python summary reads the provisional catalog as a
tab-separated table and prints value counts for the fields that drive review:
coordinate frame, sex-chromosome handling, SV type, confidence, validation
state, and graph-representation status. It is intended to expose obvious
catalog-shape problems, such as unexpected coordinate systems or a large shift
into low-confidence or unresolved graph-origin categories, before the catalog is
used for HiFi prioritization.

`RESIDUAL_TO_FROZEN_GRAPH` is a screening result, not proof that an event is
biologically valid. `HIGH` confidence requires high-quality nonmember graph
evidence plus either a compatible passing linear call or support from at least
two nonmember samples. Candidates supported only by fragmented or
`not_recommended` assemblies remain `PENDING_HIFI`; graph-member control overlap
is labelled `GRAPH_MEMBER_CONTROL_OVERLAP` and remains uncertain.

## 13. Rank samples for HiFi retrieval

The final first-pass step combines assembly QC, graph callability, candidate
types, confidence, validation state, and current read-access metadata. It keeps
all discovery blind spots, selects at most `hifi_recommendations.validation_count`
additional validation samples, and reserves
`hifi_recommendations.control_count` callable graph-member controls.

```bash
snakemake \
  --snakefile whole_pans/workflow/rules/novel_sv_discovery.smk \
  --use-conda \
  --cores 1 \
  rank_recommended_hifi_samples

RESULTS=whole_pangenome/sv_pangenome/novel_sv_discovery
column -t -s $'\t' "$RESULTS/hifi/recommended_hifi_samples.tsv" | head -30
```

The ranking prioritizes samples whose assemblies are below the 0.85 callable
threshold even when no candidate was found. This preserves the distinction
between “no event detected” and “sequence could not be assessed reliably.”

## 14. Run the complete first pass

After the pilot and configuration review, the aggregate target runs every
implemented stage:

```bash
snakemake \
  --snakefile whole_pans/workflow/rules/novel_sv_discovery.smk \
  --use-conda \
  --cores 32 \
  --rerun-incomplete \
  --printshellcmds \
  novel_sv_discovery
```

For cluster execution, add the server profile:

```bash
snakemake \
  --profile /path/to/the/server/snakemake-profile \
  --snakefile whole_pans/workflow/rules/novel_sv_discovery.smk \
  --use-conda \
  --rerun-incomplete \
  --printshellcmds \
  novel_sv_discovery
```

The workflow itself contains no SLURM or other scheduler-specific directives.

## 15. Monitor and resume safely

Snakemake writes separate logs and benchmarks under the result directory:

```bash
RESULTS=whole_pangenome/sv_pangenome/novel_sv_discovery

find "$RESULTS/logs" -type f | sort | tail
find "$RESULTS/benchmarks" -type f | sort | tail
```

Inspect planned and completed files without rerunning work:

```bash
snakemake \
  --snakefile whole_pans/workflow/rules/novel_sv_discovery.smk \
  --summary \
  novel_sv_discovery
```

After a preemption or recoverable job failure, fix the underlying resource or
input problem and rerun the same target with `--rerun-incomplete`. Snakemake
will retain valid completed outputs and schedule missing or incomplete jobs.

Useful failure checks are:

```bash
rg -n -i 'error|failed|killed|out of memory|no space' \
  whole_pangenome/sv_pangenome/novel_sv_discovery/logs \
  whole_pangenome/sv_pangenome/novel_sv_discovery/graph_screen/task_logs
```

If many graph tasks fail with memory pressure, reduce simultaneous jobs through
the profile or increase `graph_screen.mem_mb`. If the SQLite merge fills the
filesystem, move `TMPDIR` to larger scratch and rerun only
`merge_assembly_only_novel_sv_catalog`.

## 16. Handoff to the read-dependent second pass

The assembly-only run is complete when these files exist and have been
reviewed:

```text
frozen_graph/graph_inventory.tsv
manifest/assembly_discovery_manifest.tsv
graph_screen/summary/all_residual_sv_candidates.tsv.gz
graph_screen/coordinate_qc/all_residual_sv_candidates.annotated.tsv.gz
graph_screen/coordinate_qc/coordinate_qc.tsv
reference_calls/reference_call_manifest.tsv
catalog/provisional_graph_residual_sv_catalog.tsv
catalog/all_normalized_assembly_evidence.tsv.gz
hifi/recommended_hifi_samples.tsv
```

Use `recommended_hifi_samples.tsv` to request sample-level reads and populate
the optional read manifest. The later workflow extension should then perform
read-based discovery for poorly callable assemblies, breakpoint validation,
cohort-wide genotyping, and only after review construct graph version 2.

Do not add an entire low-quality assembly to the graph merely because it
supports one candidate. The plan’s central separation remains: an assembly may
provide discovery evidence without being suitable as a graph-construction
input.
