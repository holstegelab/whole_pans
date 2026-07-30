# Plan to Detect Population SVs Missing from the Whole-Genome Graph

## Reviewer revisions for maximum recall (2026-07-17)

This revision keeps the original plan's structure and its evidence-vs-modification
discipline, and adds the changes that most increase how many real SVs the cohort
yields. Summary of deltas, detailed in the phases below:

1. **Assembly-to-reference discovery is promoted to a co-primary route, run as a
   multi-caller union** (PAV + SVIM-asm + dipcall) against *both* CHM13 and
   GRCh38, rather than one caller against one reference. Minigraph rGFA is
   SV-only and reference-collapsed; the residual-graph screen alone
   systematically misses balanced inversions, translocations, and events the
   aligner expresses as split alignments. See Phase 6.
2. **A new read-based discovery route (Phase 6B)** recovers SVs from the 527
   `not_recommended` and 144 `fragmented_rescue` haplotype assemblies whose
   assembly-based calls are least reliable. These 671 haplotypes belong to 341
   unique samples and comprise about 85% of the excluded haplotypes. This is the
   single largest recall gap in the original plan, which used raw reads only for
   validation. Because HiFi reads are not available for the first run-through,
   the assembly-only pass now produces a ranked, sample-level HiFi access list
   before read-based discovery is run.
3. **Cohort-wide genotyping is added as an explicit recall step (Phase 8B):**
   after building the merged candidate catalog, genotype every one of the 982
   haplotypes (graph members included) against the augmented graph to recover
   carriers not seen at discovery and to produce allele frequencies.
4. **A base-level pangenome graph (minigraph–Cactus or PGGB) is flagged as an
   alternative, higher-completeness substrate** (Phase 9 note). If maximum SV
   catalog completeness — not just an updated minigraph rGFA — is the true
   objective, a base-level graph captures nested and smaller variation the
   SV-level rGFA collapses.
5. **Factual corrections:** the screening script
   `workflow/scripts/screen_novel_graph_svs.py` and its tests **do not exist
   yet** — only the README documents the intended interface (there is no
   `tests/test_screen_novel_graph_svs.py`; the only test file present is
   `tests/test_critical_qc_fixes.py`). The frozen graph is **~3.47 GB** on disk
   (3,475,361,198 bytes), not 3.3 GB. `vg` v1.75.1 is available; `minigraph`,
   `gfatools`, and `mash` live in the `whole-pans-pangenome` conda env, but
   `PAV`, `SVIM-asm`, `dipcall`, `truvari`, `jasmine`, and read-based SV callers
   are **not installed** and need a dedicated environment.
6. **Recoverable excluded haplotypes:** of the 786 excluded, `required_action`
   splits into **259 flagged `rerun_decontamination_and_post_qc_then_rebuild`**
   (potentially recoverable by reprocessing, not just by discovery) and **527
   `do_not_add_without_new_assembly_or_manual_review`**. Reprocessing the 259
   before the screen may itself add usable assemblies.
7. **HiFi access is explicitly staged:** the assembly-only portions of Phases
   1--7 can run while Phase 6B is deferred. They must emit
   `recommended_hifi_samples.tsv`, ranked by the expected recall or validation
   gain from retrieving the original sample-level reads. Samples with poorly
   callable assemblies remain priorities even when the first pass finds no
   candidate, because no call in uncallable sequence is not negative evidence.
   Read-dependent validation and genotyping remain pending until the selected
   reads are retrieved.

## Objective

Determine whether haplotype assemblies excluded from the current Minigraph
SV graph contain structural-variant alleles that the graph does not represent,
and construct an updated graph that captures as many high-confidence SVs from
this population as possible.

The current dataset contains:

- 982 cleaned haplotype assemblies from 491 samples in
  `../whole_pangenome/assembly_qc_decontaminated/results/resources/all_cleaned_assemblies.tsv`;
- 196 haplotype assemblies used to construct the current graph;
- 786 cleaned haplotype assemblies outside the graph, belonging to 400 unique
  samples;
- 37 `best_rescue`, 78 `reasonable_rescue`, 144 `fragmented_rescue`, and
  527 `not_recommended` assemblies in the existing feasibility analysis;
- the 671 `fragmented_rescue`/`not_recommended` haplotypes belong to 341 unique
  samples (the natural maximum-recall starting set for sample-level HiFi
  retrieval); 330 of these samples have two low-tier excluded haplotypes and 11
  have one;
- 14 excluded haplotypes that would restore the missing mate of a sample
  already represented in the graph (all 14 carry `mate_status = WARN` and
  `graph_context = restore_missing_mate_for_graph_sample`);
- within the 786, `required_action` splits into 259
  `rerun_decontamination_and_post_qc_then_rebuild_graph_or_incremental_test`
  (may be recoverable by reprocessing) and 527
  `do_not_add_without_new_assembly_or_manual_review`.

The existing graph QC indicates that the graph is probably not saturated. The
last construction decile still contributed approximately 1.05 Mb of
non-reference sequence per added assembly. However, many excluded assemblies
have completeness, reference-coverage, duplication, or fragmentation problems.
Therefore, using an assembly to *discover or support* an SV must be separated
from allowing that assembly to *modify the graph*.

Relevant existing reports are:

- `../whole_pangenome/sv_pangenome/qc_analysis/report/pangenome_qc_report.md`;
- `../whole_pangenome/sv_pangenome/qc_analysis/report/excluded_haplotype_feasibility.md`;
- `../whole_pangenome/sv_pangenome/qc_analysis/tables/excluded_haplotype_feasibility.csv`.

Before running the analysis, refresh the feasibility table from the current QC
summaries because the existing feasibility report predates the latest complete
post-decontamination QC run.

## Overall Strategy

Use three complementary discovery routes so that no single method's blind spots
cap recall:

1. **Frozen-graph residual screen** — map every cleaned assembly to a frozen
   copy of the current graph and detect residual SV-sized differences
   (Phases 3-5);
2. **Assembly-to-reference multi-caller union** — independently call SVs with
   PAV + SVIM-asm + dipcall against both CHM13 and GRCh38, capturing inversions,
   duplications, translocations, and complex events the single-CIGAR graph
   screen misses (Phase 6);
3. **Read-based discovery** — for the 671 fragmented/`not_recommended`
   haplotype assemblies from 341 unique samples, plus other samples selected for
   candidate validation, call SVs directly from raw long reads after the reads
   have been requested (Phase 6B).

Merge and cluster all three routes' candidates, validate, genotype the full
cohort to recover carriers and allele frequencies (Phase 8B), add only
high-confidence alleles to a new graph version, and re-screen to verify capture.

For the first run-through, execute the two assembly-based routes, create an
assembly-only Phase 7 catalog, and generate the ranked HiFi access list. Do not
block the first pass on read availability and do not finalize a read-dependent
candidate as absent or artifactual merely because its reads have not yet been
retrieved. Resume Phase 6B and the read-dependent parts of Phases 8 and 8B in a
second pass.

## Phase 1: Freeze and Inventory the Current Graph

Keep the existing graph unchanged as version 1:

`../whole_pangenome/sv_pangenome/results/graphs/sv_pangenome.minigraph.gfa`

Record:

- GFA checksum and file size;
- Minigraph version and construction parameters;
- ordered reference and assembly inputs;
- graph QC outputs and build logs;
- a catalog of existing graph bubbles, extracted with `gfatools bubble`;
- graph node coordinates and the rGFA `SN`, `SO`, and `SR` tags.

The frozen graph provides a stable target for deciding whether an allele is
already represented. Do not add assemblies during the discovery screen.

## Phase 2: Build the Discovery Manifest

Create one table containing every cleaned assembly and join:

- assembly ID, sample ID, and haplotype;
- cleaned FASTA path;
- current graph membership;
- original and post-decontamination QC metrics;
- contamination results;
- Mash distance to CHM13 and to the graph input cohort;
- rescue tier;
- whether the assembly restores a missing mate;
- raw-read access status (`not_checked`, `recommended`, `requested`,
  `available`, or `unavailable`) and HiFi/ONT paths, platform, and coverage when
  known;
- population, phenotype, or sequencing-batch metadata, if available.

Use the following priority order for pilots and graph inclusion:

1. missing mates of already represented samples;
2. `best_rescue` assemblies;
3. `reasonable_rescue` assemblies;
4. `fragmented_rescue` assemblies;
5. `not_recommended` assemblies.

This is an evidence ranking, not an exclusion rule. Screen all 786 assemblies
eventually because the goal is maximum within-cohort SV recall.

At the end of Phase 2, emit a provisional sample-level HiFi request list from
QC tier and graph membership alone so data-access work can start in parallel.
After the assembly-only discovery pass, add callable fractions and candidate
evidence and rerank it as the final `recommended_hifi_samples.tsv` described in
Phase 6B. Keep the provisional and final versions for provenance.

## Phase 3: Implement and Calibrate Frozen-Graph Screening

The expected interface for
`workflow/scripts/screen_novel_graph_svs.py` is documented in the README
(the `index-graph`, `tasks`, `run`, and `summarize` subcommands), but neither
the script nor a test for it currently exists — the only test file present is
`tests/test_critical_qc_fixes.py`. Implement the script **and** its unit tests
before the production screen, and connect it to the workflow. Add focused tests
using small synthetic rGFA/GAF fixtures so CIGAR parsing, anchor checking,
adjacent-indel handling, and split-alignment logic are verified before the code
is trusted on 786 assemblies; wire the four subcommands into `workflow/rules/`
so the screen is reproducible rather than a set of manual shell calls.

Provision the tooling before running: create a dedicated conda environment for
discovery/merging tools (`PAV`, `SVIM-asm`, `dipcall`, `truvari`, `jasmine` or
`SVanalyzer`, `samtools`, `bcftools`) alongside the existing
`whole-pans-pangenome` env (which already provides `minigraph`, `gfatools`,
`mash`). `vg` v1.75.1 is available at `../vg`. None of the assembly-to-reference
or merge tools are currently installed.

The script should provide separate commands to:

1. index rGFA segments and stable coordinates;
2. create a deterministic task table;
3. map and analyse one assembly or assembly batch;
4. combine all per-assembly results.

Run a calibration pilot containing:

- 5--10 assemblies used to build the graph;
- 10--20 top `best_rescue` or missing-mate assemblies;
- several fragmented or otherwise poor assemblies.

Graph-member assemblies are essential controls. Residual events found in them
measure mapping ambiguity and incomplete graph representation, and help set
filters without relying only on theoretical thresholds.

Also include positive-recall controls. Negative graph-member controls alone can
make stringent filters look good while hiding false negatives. Use a small
truth set spanning insertions, deletions, inversions, duplications, complex
events, size bins, and repeat contexts. Suitable controls include known alleles
temporarily omitted from a small test graph, leave-one-assembly-out graph
builds, and synthetic spike-ins with known breakpoints. Report recall by SV
class/size/context and callable fraction, not only the residual false-positive
rate. Tune the sensitive and high-confidence tiers separately.

Use the pilot to measure runtime, peak memory, temporary storage, and suitable
batch size. The graph is approximately 3.47 GB on disk (3,475,361,198 bytes;
677,995 segments, 983,547 links per the current QC report), so repeatedly
loading it in hundreds of simultaneous jobs could stress the shared filesystem.
When practical, map several assemblies per job while keeping per-assembly
outputs. Because Minigraph loads the whole graph into memory per invocation,
size the job memory request from the pilot's observed peak (expect several GB of
resident memory beyond the on-disk size) and cap concurrency accordingly rather
than launching one job per assembly.

## Phase 4: Map Complete Assemblies to the Frozen Graph

Use Minigraph assembly-to-graph alignment and retain CIGAR operations and
secondary mappings. A starting command is:

```bash
minigraph -cxasm -l5000 -t 16 \
  ../whole_pangenome/sv_pangenome/results/graphs/sv_pangenome.minigraph.gfa \
  assembly.clean.fa.gz \
  > assembly.gaf
```

Confirm the exact options against the installed Minigraph version before the
production run.

Do not let a single `-l5000` mapping pass define whether an assembly is
assessable. For assemblies or contigs with low aligned/callable fraction, run a
second, calibrated sensitivity pass with settings that recover shorter anchors
and additional secondary/split alignments. Keep the primary pass for
comparability and tag which pass produced each event; use the positive controls
from Phase 3 to choose the rescue settings rather than lowering thresholds
blindly.

For every query contig, record:

- query length and aligned fraction;
- primary and secondary alignments;
- alignment identity and MAPQ;
- alignment and query coordinates;
- strand and graph traversal;
- CIGAR operations;
- split or multiple-primary alignments;
- terminal and internal unaligned sequence.

Report both the callable fraction of each assembly and candidate events. A lack
of residual SVs in poorly aligned sequence is not evidence that the graph
represents that sequence.

## Phase 5: Detect Candidate Graph-Missing SVs

Use 50 bp as the minimum SV size. Starting high-confidence criteria are:

- primary alignment;
- alignment length of at least 5 kb;
- identity of at least 90%;
- MAPQ of at least 5;
- at least 2 kb of aligned sequence on both sides of the event;
- an internal insertion or deletion of at least 50 bp.

Retain candidates that fail these criteria in a review tier instead of
discarding them. This preserves recall in repetitive regions.

For the raw sensitivity tier, also retain 30--49 bp indels near the SV boundary
and clusters of nearby CIGAR indels whose combined or locally realigned allele
length is at least 50 bp. Normalize and locally realign these before deciding
whether they meet the final 50 bp definition; alignment fragmentation and
microhomology can otherwise turn one true SV into several sub-threshold edits.

Identify and classify:

- internal insertions and deletions from `cg:Z` CIGAR operations;
- internal unaligned sequence with two-sided graph anchors;
- split alignments and unexpected orientation, which may indicate inversions or
  complex SVs;
- alignments joining distant loci or chromosomes;
- possible duplications and CNVs;
- repeat-associated events with multiple similar graph placements.

Do not call terminal clipping as an SV unless an independent alignment provides
the second breakpoint. Keep ambiguous, low-MAPQ, repeat-rich, and contig-end
events as `uncertain` rather than treating them as absent or false.

Suggested per-assembly outputs are:

- compressed raw GAF;
- one row per contig alignment;
- one row per raw residual event;
- assembly-level callable and candidate summary;
- explicit lists of complex and unresolved alignments.

## Phase 6: Independent Assembly-to-Reference Discovery (co-primary route)

This route is co-primary with the frozen-graph screen, not a supplement. An
SV-level Minigraph rGFA is reference-collapsed and expresses many real variants
(balanced inversions, translocations, dispersed duplications, nested/complex
events) as split or multi-primary alignments rather than as single CIGAR
indels, so the residual-graph screen alone under-recalls exactly the SV classes
that are hardest to represent. Run assembly-to-reference discovery on every
assembly good enough to call from, in parallel with Phase 4-5.

Run a **multi-caller union**, not a single caller, because callers disagree
substantially on inversions, duplications, and breakpoint placement:

- **PAV** (alignment-based, strong on inversions and complex loci);
- **SVIM-asm** (fast, good indel/duplication sensitivity);
- **dipcall** (produces a phased, well-normalized VCF and a callable-region BED;
  its BED defines where absence of a call means "callable and reference-like"
  versus "not assessed").

Run each caller against **both CHM13 and GRCh38**. Use CHM13 as the primary
coordinate system because it is the graph backbone, but GRCh38 recovers alleles
in regions where the two references disagree and makes results comparable with
external SV catalogs (HPRC, gnomAD-SV, 1000G). Project the union to a single
coordinate system in Phase 7; retain both liftovers.

Pair the two haplotypes of each sample where both are available (dipcall and PAV
both consume a diploid pair) so heterozygous and hap-specific SVs are called
correctly rather than as two independent haploid calls.

This route should detect or clarify:

- inversions;
- tandem and interspersed duplications;
- CNVs;
- translocations;
- nested and complex SVs;
- events that graph mapping expresses as multiple alignments rather than one
  insertion or deletion.

Take the union of frozen-graph residual candidates and assembly-to-reference
calls. For each reference-based SV, construct its alternate allele with
sufficient flanking sequence and test whether that allele is representable in
the frozen graph.

`minigraph --call` may be used to genotype bubbles already present in the graph,
but it must not be used as the only missing-allele discovery method.

## Phase 6B: Read-Based Discovery for Low-Quality Assemblies

The largest recall gap in the original plan is that it uses raw reads only for
*validation*. But 144 `fragmented_rescue` and 527 `not_recommended` haplotype
assemblies — 671 of 786 excluded haplotypes (about 85%), from 341 unique samples
— are precisely the ones for which assembly-based calling (Phases 4--6) is
least reliable. Relegating their reads to validation-only would discard much of
their discoverable SV signal. Where raw reads are retrieved, call SVs directly
from reads and treat those calls as first-class discovery evidence.

**First-pass constraint — reads are not initially accessible.** Do not block
the assembly-only screen while locating reads. Complete Phases 4, 5, and 6,
defer Phase 6B, and run an assembly-only Phase 7 merge, carrying `PENDING_HIFI`
as an explicit validation state where appropriate. Absence of a call from a
low-callability assembly must not be interpreted as reference genotype or as
evidence that its sample does not warrant read retrieval.

**Required first-pass output — ranked HiFi access list.** Generate one row per
sample (not per haplotype) in `recommended_hifi_samples.tsv`. Original HiFi
reads are sample-level and usually contain both haplotypes; do not request the
same read set twice because both haplotype assemblies are listed. Include every
recommended sample and a transparent priority tier rather than imposing an
arbitrary top-N cutoff.

For maximum recall, the default discovery request set is all 341 unique samples
with at least one `fragmented_rescue` or `not_recommended` excluded haplotype.
Use the ranking to retrieve them in batches if access or compute is limited, not
to remove lower-ranked samples from the final sensitivity analysis. Add samples
outside this set when their first-pass calls need validation or they provide
controls. In the current table, 330 samples have two such low-tier haplotypes
and form the first QC-based batch; 11 have one. Refresh these counts after
regenerating the feasibility table.

Assign each sample all applicable reason codes and use the highest applicable
priority:

1. `P1_DISCOVERY_BLIND_SPOT`: one or both haplotypes have very low callable
   fraction, failed mapping, or a `fragmented_rescue`/`not_recommended`
   assembly. Put samples with both haplotypes poorly callable first. A sample
   with zero first-pass candidates remains in this tier because the assembly
   screen had little opportunity to find an SV.
2. `P1_VALIDATE_CANDIDATE`: the sample carries a singleton, a candidate found
   only in a poor assembly, a proposed graph-addition allele, a complex or
   balanced event, a method-discordant event, or an event in a repeat,
   segmental duplication, low-MAPQ region, or contig end.
3. `P2_GENOTYPE_OR_PHASE`: read evidence would resolve an uncertain carrier,
   genotype, breakpoint, allele sequence, or haplotype assignment for an
   otherwise credible event.
4. `P3_CONTROL`: a small population- and batch-aware set of graph-member and
   high-quality assembly samples, including both positive-candidate and
   no-candidate controls, for estimating read/assembly concordance and the
   read-calling false-positive rate.

Within a tier, rank first by number of uncallable haplotypes and increasing
callable fraction, then by number and importance of unresolved candidates. Use
underrepresented population/sequencing batches and missing-mate status as
tie-breakers. Do not rank only by raw candidate count: misassembled samples can
produce many artifacts, while the most incomplete samples can produce no calls.

The list should contain at least:

- `sample_id`, both assembly IDs, graph membership, missing-mate status, rescue
  tiers, and the callable fraction of each haplotype;
- `priority_tier`, semicolon-separated `reason_codes`, and a short
  human-readable recommendation;
- candidate counts by discovery method, confidence, and SV type, plus the most
  important candidate IDs;
- population and sequencing batch where available;
- intended use (`discovery`, `validation`, `genotyping`, `phasing`, or
  `control`) and requested read type;
- read-access status, path/accession when later known, platform, yield or
  estimated coverage, and checksums.

After access is granted, build the read manifest from this list and join it to
the discovery manifest. If reads cannot be obtained for a recommended sample,
mark it `unavailable` and keep its affected regions/candidates explicitly
read-unassessable; do not silently drop the sample.

**Method.** For each sample with reads:

1. Verify sample identity, read yield, length distribution, and estimated
   genomic coverage. Retain the full read set for the first sensitive call; if
   coverage is low or uneven, propagate read-callable regions rather than
   interpreting a no-call as a reference genotype.
2. Map long reads to CHM13 (and GRCh38) with `minimap2 -ax map-hifi` /
   `map-ont`, and to the frozen graph with `vg giraffe --preset hifi` (when the
   required indexes have been built) or `GraphAligner` where a graph-relative
   callset is wanted;
3. Call SVs with both `sniffles2` and `cuteSV` (and optionally `pbsv` for HiFi)
   against the linear references. Retain the multi-caller **union** as the
   sensitivity set and annotate caller support as a confidence feature; do not
   discard a single-caller event solely because it is outside the intersection;
4. Feed these calls into the same Phase 7 clustering, tagged with discovery
   method `reads_<caller>` so their independence from assembly calls is tracked.

Read-based calls that recur across independent samples, or that corroborate an
assembly-based candidate, are strong evidence even when the sample's assembly is
`not_recommended`. This is the mechanism by which poor-assembly samples still
contribute confirmed SVs to the catalog without their assemblies ever entering
the graph.

## Phase 7: Normalize and Cluster Population SV Alleles

Project candidates to CHM13 through reference-ranked rGFA nodes using the `SN`,
`SO`, and `SR` tags. Retain graph node/offset coordinates and flanking sequences
for events that cannot be projected unambiguously.

Cluster candidates across assemblies using:

- compatible SV type;
- breakpoint proximity;
- length similarity;
- reciprocal overlap for reference-spanning events;
- inserted-sequence identity for insertions;
- orientation and both breakpoint pairs for inversions and complex SVs.

Truvari or Jasmine are suitable starting tools, but clustering parameters
should be calibrated separately by SV type, size, and repeat context. Avoid
merging nearby but biologically distinct alleles solely by breakpoint distance.

The clustered table should contain one row per population SV allele with:

- stable event ID;
- CHM13 coordinates and optional GRCh38 coordinates;
- graph coordinates;
- SV type and size;
- reference and alternate allele sequence where resolved;
- carrier assemblies, haplotypes, and independent samples;
- discovery methods;
- assembly QC tiers;
- graph-representation status;
- validation and confidence status.

## Phase 8: Validate Missing-Allele Candidates

Evaluate candidates using independent evidence. A strong candidate should have:

- clean two-sided sequence anchors;
- consistent breakpoints and allele sequence;
- recurrence in independent samples, or agreement between graph-based and
  reference-based discovery;
- no contamination evidence;
- plausible repeat and sequence-complexity context;
- raw HiFi read support where data are available.

For singleton SVs, require breakpoint-spanning long reads, allele-specific
k-mer support, local reassembly, or another strong orthogonal signal. Two
haplotypes from the same sample should not count as two independent carriers.

Candidates found only in fragmented or `not_recommended` assemblies should not
enter the graph without raw-read validation or recurrence in a better assembly.
During the read-free first pass, label such candidates `PENDING_HIFI` rather
than downgrading them solely because the requested reads are not yet accessible.
Give additional scrutiny to:

- centromeres and telomeres;
- segmental duplications;
- long tandem repeats;
- contig ends;
- regions with multiple high-scoring alignments;
- inserted sequence with non-human contamination evidence.

Use confidence levels such as:

- `HIGH`: strong anchors plus independent assembly/method/read support;
- `MEDIUM`: strong sequence representation but only one evidence source;
- `UNCERTAIN`: repeat-associated, weakly anchored, fragmented, or unsupported;
- `ARTIFACT`: evidence supports contamination, misassembly, or mapping error.

## Phase 8B: Cohort-Wide Genotyping and Allele Frequencies

Discovery finds an allele in whichever assemblies happened to carry it well
enough to call; it does not tell you how many of the 982 haplotypes actually
carry it. Recover the remaining carriers — and produce the allele-frequency
table the cohort ultimately needs — by genotyping every haplotype against the
candidate set, including the 196 graph members.

Two complementary genotypers:

- **Assembly-path genotyping:** once the augmented graph exists (or an
  incremental test graph carrying the candidate bubbles), re-map all 982
  assemblies and record which traverse each candidate bubble. This reuses the
  Phase 4 machinery and needs no reads.
- **Read-based genotyping:** for samples with reads, genotype the candidate SVs
  with `vg call`/`vg giraffe` against the augmented graph, or with a long-read
  genotyper against the linear references. This confirms carriers whose
  assemblies are too poor to traverse the bubble cleanly. Read genotypes are
  diploid sample-level results; assign them to a specific haplotype only when
  phasing evidence supports that assignment.

Report per-allele carrier counts, independent-sample counts, and allele
frequency, separated by QC tier and by discovery-vs-genotyping evidence. This
step both raises recall (carriers missed at discovery) and converts the catalog
from a presence list into a frequency-aware population resource.

## Phase 9: Construct Graph Version 2

Select representative high-confidence source assemblies rather than adding
every carrier of every allele. Suggested inclusion order is:

1. validated missing mates;
2. validated `best_rescue` assemblies;
3. validated `reasonable_rescue` assemblies;
4. selected fragmented assemblies only when their alleles have strong
   independent support.

First build an incremental test graph and inspect how much sequence and
topological complexity each source adds. For the final graph, perform a
reproducible rebuild from:

1. CHM13;
2. hg38 compatibility sequence;
3. the original 196 graph assemblies;
4. selected new assemblies in a recorded deterministic order.

For a validated allele found only in an unusable whole-genome assembly, either
keep it in a separate population SV catalog or experimentally add a curated
allele-plus-flanks sequence. Do not add the complete poor assembly merely to
capture one allele.

**Substrate choice — consider a base-level graph.** The current graph is an
SV-level Minigraph rGFA (`-l5000` scale), which by design collapses reference
sequence and omits sub-SV and nested variation. If the true objective is a
maximally complete SV catalog rather than only an updated rGFA, evaluate
building a **base-level graph** with minigraph–Cactus (the standard "add
base-level detail to the minigraph backbone" pipeline that yields a
`vg`-genotypable graph plus a decomposed VCF) or PGGB. A base-level graph:

- represents nested and smaller variants the SV rGFA cannot hold;
- is directly `vg`-genotypable (supports Phase 8B without a separate step);
- emits a normalized, decomposed VCF suitable for merge with the reference-based
  callsets.

Weigh this against cost: minigraph–Cactus over ~200 high-quality haplotypes is
substantially heavier than an rGFA rebuild. A reasonable path is to keep the
rGFA as the discovery/QC backbone and build a base-level graph over the
*validated high-confidence source assemblies only*, as the genotyping and
final-catalog substrate. Record which substrate each deliverable is built on.

Repeat the existing graph QC and reject a proposed build if it introduces:

- missing or inconsistent rGFA tags;
- disconnected reference chromosomes;
- unexpected rank or source-name problems;
- extreme graph hubs or implausible growth;
- large contributions dominated by contamination or unanchored sequence.

## Phase 10: Verify Capture and Measure Saturation

Repeat the complete assembly-to-graph screen against graph version 2. For every
validated missing allele, verify that:

- it now maps without the original residual SV-sized difference;
- both flanks remain correctly anchored;
- the locus did not become excessively ambiguous;
- graph integrity remains acceptable.

Produce saturation summaries showing:

- cumulative validated population SV alleles versus assemblies screened;
- new SV alleles and new SV base pairs per additional haplotype;
- results by QC tier and population or phenotype subgroup;
- callable assembly fraction;
- SV type and size distributions;
- singleton and recurrent allele counts;
- unresolved candidate count;
- fraction of validated candidates captured by graph version 2.

Use several randomized assembly orders, plus population-stratified orders, to
place uncertainty around the discovery curve. Because the objective is maximum
recall within this cohort, screening should finish only after all 786 excluded
assemblies have been analysed. The saturation curve determines which validated
alleles and source assemblies must enter the graph, not which assemblies are
allowed to contribute discovery evidence.

## Main Deliverables

1. Complete assembly discovery manifest with raw-read access status; after
   retrieval, join it to a raw-read manifest (path/platform/coverage/checksum per
   sample; unavailable and read-unassessable samples flagged).
2. Per-assembly graph-mapping and callable-region summaries.
3. Raw and filtered graph-residual SV candidates.
4. Independent assembly-to-reference SV callset (PAV + SVIM-asm + dipcall,
   CHM13 and GRCh38), with per-caller and consensus tiers.
5. Ranked, sample-level `recommended_hifi_samples.tsv` with priority tiers,
   reason codes, callable fractions, linked candidate IDs, and intended use.
6. Read-based SV callset for retrieved low-quality-assembly and validation
   samples (Phase 6B).
7. Deduplicated population SV allele catalog with carriers, discovery methods,
   and confidence.
8. Cohort-wide genotype matrix and allele-frequency table over all 982
   haplotypes (Phase 8B).
9. Validation evidence for singleton and difficult candidates.
10. Reproducible graph version 2 input list and build provenance (and, if built,
   the base-level graph and its decomposed VCF).
11. Graph version 2 QC report.
12. Before/after capture and population saturation report.
13. `screen_novel_graph_svs.py` implementation plus its unit tests and workflow
    rules (currently absent).

## Success Criteria

The first assembly-only run-through is complete when:

- all 982 cleaned assemblies have a graph-screen result or a documented reason
  why they are unassessable;
- the assembly-to-reference and graph-residual calls have been merged into a
  provisional catalog with read-dependent evidence marked `PENDING_HIFI`;
- `recommended_hifi_samples.tsv` has been generated at sample level, includes
  poorly callable samples even when they have no candidate calls, and links
  validation recommendations to specific candidate IDs.

The full analysis is complete when:

- every sample whose reads were retrieved has a read-based
  discovery/genotyping result, and every recommended sample whose reads could
  not be retrieved remains explicitly marked unavailable/read-unassessable;
- the discovery routes are reconciled: candidates unique to each of the three
  routes are quantified, so the marginal recall contribution of the
  reference-based and read-based routes over the graph screen is reported rather
  than assumed;
- all high-confidence graph-missing candidates are either represented in graph
  version 2 or explicitly retained in the external SV catalog;
- graph-member controls show an acceptably low residual false-positive rate;
- every graph-added allele has traceable assembly and validation evidence;
- the updated graph passes the existing integrity checks;
- re-screening confirms capture without introducing major mapping ambiguity;
- the remaining new-allele discovery rate and unresolved burden are reported,
  rather than assuming that the graph is closed.
