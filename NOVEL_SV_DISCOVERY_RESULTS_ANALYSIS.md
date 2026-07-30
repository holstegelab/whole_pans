# Analysis of the novel-SV discovery results

**Analysis date:** 2026-07-29  
**Workflow:** `workflow/rules/novel_sv_discovery.smk`  
**Results:** `../whole_pangenome/sv_pangenome/novel_sv_discovery/`

## Executive summary

The assembly-only workflow completed its intended first pass for all 982
haplotype assemblies from 491 samples. The frozen Minigraph rGFA was
successfully inventoried, all assemblies were screened, both enabled linear
callers ran against CHM13 and GRCh38, graph coordinates passed QC, a provisional
catalog was built, and 420 samples were ranked for HiFi retrieval.

The results support the conclusion that the frozen graph is not saturated, but
they do **not** establish 895,941 confirmed novel SVs:

- 3,070,340 graph-residual observations of at least 50 bp were retained before
  clustering.
- These observations formed 895,941 provisional catalog loci.
- Only 39,976 loci (4.46%) are labelled `HIGH`; 812,134 (90.65%) remain
  `UNCERTAIN`.
- 676,596 loci (75.52%) are supported only by assemblies in the two lowest QC
  tiers and are therefore `PENDING_HIFI`.
- 110,451 loci (12.33%) overlap residual evidence from graph-member controls,
  showing that alignment/representation background is substantial.
- 700,153 loci (78.15%) occur in only one excluded sample, and 745,877 (83.25%)
  are shorter than 500 bp. Both features increase sensitivity to assembly and
  alignment artefacts.

The most defensible immediate review set is the 39,976 `HIGH` loci, especially
the 17,622 that also have a `LINEAR_CALLER_SUPPORTED` validation label. These
are still candidates rather than truth-set-confirmed alleles.

A major workflow-scope issue must be kept in mind: the final catalog only emits
clusters that contain a Minigraph residual from an assembly outside the frozen
graph. SVIM-asm/dipcall-only calls are written to the normalized evidence file
but do not become catalog rows. Consequently, the catalog is an annotated
**graph-residual catalog**, not a full union of the assembly-to-reference
callers. PAV and read-based discovery were not run.

## 1. Run identity and cohort coverage

### Frozen graph

| Metric | Result |
|---|---:|
| Graph type | Minigraph SV-level rGFA |
| Graph size | 3,475,361,198 bytes |
| SHA-256 | `9a6a8a0044bb1e32b53ff5b6e4aa45ab28ca8bfa8d887b4f6a48f10c142ead4f` |
| Segments | 677,995 |
| Links | 983,547 |
| Graph input haplotypes | 196 |
| Minigraph | 0.21-r606 |
| gfatools | 0.5-r235-dirty |

The copied graph QC report records 271.9 Mb of non-reference sequence and a
mean contribution of approximately 1.05 Mb per added assembly in the final
construction decile. That continued late growth is independent evidence that
the graph remained open when construction stopped.

### Discovery cohort

| Category | Haplotypes |
|---|---:|
| Total cleaned assemblies | 982 |
| Frozen-graph members | 196 |
| Assemblies outside the graph | 786 |
| `best_rescue` | 40 |
| `reasonable_rescue` | 75 |
| `fragmented_rescue` | 144 |
| `not_recommended` | 527 |

All 491 samples have both haplotypes in the cleaned manifest. The current
refreshed tier counts differ slightly from older planning numbers: there are 40
`best_rescue` and 75 `reasonable_rescue` haplotypes in this run.

Among the excluded haplotypes:

- 14 restore a missing mate of a graph-represented sample;
- 245 are marked for decontamination/post-QC rerun before possible graph use;
- 527 require a new assembly or manual review before graph inclusion.

No read paths, population labels, or sequencing-batch labels were supplied.
All 982 assemblies have `read_access_status=not_checked`.

## 2. Frozen-graph screen

### Callability

| Metric | Result |
|---|---:|
| Mean primary callable fraction | 0.953679 |
| Mean, graph members | 0.954633 |
| Mean, excluded assemblies | 0.953442 |
| Minimum callable fraction | 0.936802 |
| Assemblies below 0.85 | 0 |
| Assemblies below 0.90 | 0 |
| Assemblies below 0.95 | 131 |
| Sensitivity-pass runs | 0 |

The primary Minigraph pass was broadly callable for every assembly, so the
configured rescue pass was never triggered. The similarity of graph-member and
excluded-assembly callability suggests that gross mappability is not the main
source of the difference in residual counts.

Callability here measures how much assembly sequence met the screen's alignment
criteria. It is not an SV-recall estimate, particularly in repeats, duplicated
sequence, and complex graph regions.

### Residual observations

| Screen category | All assemblies | Excluded assemblies | Graph-member controls |
|---|---:|---:|---:|
| Raw events from 30 bp | 10,559,290 | 8,597,768 | 1,961,522 |
| High-confidence, at least 50 bp | 1,986,693 | 1,673,214 | 313,479 |
| Review, at least 50 bp | 1,083,647 | 906,388 | 177,259 |
| Below final 50 bp threshold | 7,488,950 | 6,018,166 | 1,470,784 |
| High-confidence + review | 3,070,340 | 2,579,602 | 490,738 |
| Split/complex contigs | 241,143 | 216,055 | 25,088 |

Only 29.08% of raw events passed the final size/confidence route into the
coordinate-QC input; 70.92% remained subthreshold.

Excluded assemblies produced approximately 3,282 retained observations per
haplotype, compared with 2,504 per graph-member control. The excess is
consistent with missing graph sequence, but the large control baseline means
the difference cannot be interpreted as an equivalent number of novel alleles.

Every assembly received `REVIEW_SPLIT_OR_COMPLEX_ALIGNMENT`. This happens
because every assembly has at least one contig classified as split/complex.
The status is therefore not discriminative at the assembly level in this run;
the underlying counts and loci are more informative.

### Coordinate QC

All 3,070,340 retained residual observations received a stable graph coordinate:

| Coordinate-QC metric | Result |
|---|---:|
| Resolved | 3,070,340 |
| Unresolved | 0 |
| Resolved fraction | 1.000000 |
| Allowed unresolved fraction | 0.01 |
| Status | `PASS` |

This means coordinate recovery worked technically. It does not mean that every
event could be projected to the CHM13 backbone. In the final catalog, 179,017
loci remain in an off-reference `GRAPH` coordinate frame.

## 3. Provisional catalog

The clustering settings were a maximum breakpoint distance of 500 bp and a
minimum length similarity of 0.70. Insertions with two available sequences
were only merged when their sequences agreed exactly.

The catalog contains **895,941 loci**.

### Coordinate frames

| Coordinate system | Loci | Percent |
|---|---:|---:|
| CHM13 | 716,924 | 80.02% |
| Off-reference graph (`GRAPH`) | 179,017 | 19.98% |

`GRAPH` loci have stable graph coordinates, but they are not ordinary CHM13
positions. They require graph-aware review or a separate off-reference
coordinate/liftover strategy.

All `position_0` and `end_0` fields are **0-based**. They should not be copied
directly into a 1-based VCF position without conversion.

### SV classes

| SV type | Loci | Percent |
|---|---:|---:|
| DEL | 350,556 | 39.13% |
| INS | 331,400 | 36.99% |
| COMPLEX_INDEL | 213,985 | 23.88% |

No inversion, duplication, translocation, BND, or CNV rows occur in the final
catalog. This reflects the graph-residual gating described below, not evidence
that those classes are absent from the cohort.

### Size distribution

| Absolute SV length | Loci | Percent |
|---|---:|---:|
| 50–99 bp | 431,151 | 48.12% |
| 100–499 bp | 314,726 | 35.13% |
| 500–999 bp | 49,017 | 5.47% |
| 1–9.9 kb | 91,920 | 10.26% |
| 10–49.9 kb | 8,399 | 0.94% |
| 50–999.9 kb | 728 | 0.08% |

There are no catalog representatives of at least 1 Mb. The strong enrichment
near the 50 bp boundary makes repeat/context annotation and breakpoint-level
validation especially important.

### Independent sample support

| Excluded samples supporting a locus | Loci | Percent |
|---|---:|---:|
| 1 | 700,153 | 78.15% |
| 2–4 | 130,854 | 14.60% |
| 5–9 | 37,676 | 4.21% |
| 10–49 | 24,091 | 2.69% |
| At least 50 | 3,167 | 0.35% |

Mean support is 2.27 excluded samples per locus and the maximum is 397. These
are discovery-support counts, **not allele frequencies**: the workflow has not
yet genotyped every catalog locus in every sample.

### Evidence concordance

| Methods represented in a catalog cluster | Loci | Percent |
|---|---:|---:|
| Minigraph residual only | 702,461 | 78.40% |
| Minigraph + SVIM-asm only | 55,644 | 6.21% |
| Minigraph + dipcall only | 29,352 | 3.28% |
| Minigraph + both linear callers | 108,484 | 12.11% |
| **Any enabled linear-caller support** | **193,480** | **21.60%** |

Linear-caller support is present for about one fifth of the graph-residual
catalog. Concordance is useful evidence, but it is not completely independent:
all routes use the same assemblies, and both reference callers depend on
assembly-to-reference alignment.

### Confidence and validation labels

| Confidence | Loci | Percent |
|---|---:|---:|
| `HIGH` | 39,976 | 4.46% |
| `MEDIUM` | 43,831 | 4.89% |
| `UNCERTAIN` | 812,134 | 90.65% |

| Validation status | Loci | Percent |
|---|---:|---:|
| `PENDING_HIFI` | 676,596 | 75.52% |
| `GRAPH_MEMBER_CONTROL_OVERLAP` | 110,451 | 12.33% |
| `ASSEMBLY_ONLY_REVIEW` | 86,796 | 9.69% |
| `LINEAR_CALLER_SUPPORTED` | 22,098 | 2.47% |

The validation label is hierarchical. A locus supported by a linear caller can
still be labelled `PENDING_HIFI` if all discovery assemblies are low-QC, or
`GRAPH_MEMBER_CONTROL_OVERLAP` if a graph-member residual is in the same
cluster. Therefore, `discovery_methods` identifies all 193,480 loci with linear
support, while only 22,098 receive the explicit
`LINEAR_CALLER_SUPPORTED` validation status.

The `HIGH` set consists of:

- 17,622 `LINEAR_CALLER_SUPPORTED` loci;
- 22,354 `ASSEMBLY_ONLY_REVIEW` loci with high graph evidence and support from
  at least two excluded samples;
- 27,452 CHM13-coordinate and 12,524 off-reference graph-coordinate loci.

The confidence labels are deterministic workflow heuristics. They have not
been calibrated against a positive truth set and should not be renamed
“validated.”

## 4. Examples from the `HIGH` candidate set

### Highest independent-sample support

| Event | Coordinate (`position_0`) | Type | Length | Samples | Validation |
|---|---|---:|---:|---:|---|
| `PSV_f0c2bea4f61e3f4a` | `GRAPH:h1tg000831l\|SR110:342044` | INS | 52 | 82 | Assembly-only review |
| `PSV_b85a8d90ca0f5f7b` | `CHM13:chr17:83843738` | INS | 56 | 61 | Linear-caller supported |
| `PSV_42cc459eeed19347` | `GRAPH:h1tg000162l\|SR17:4665828` | INS | 235 | 61 | Assembly-only review |
| `PSV_1fd4c09605f9edeb` | `CHM13:chr17:83852838` | DEL | 308 | 59 | Linear-caller supported |
| `PSV_4de77b4530a1cfe9` | `CHM13:chr10:42202336` | INS | 71 | 57 | Assembly-only review |

Several of the highest-support events are short and occur close to other
events, so repeat-driven or multi-representation loci should be checked before
interpreting the support as population prevalence.

### Largest `HIGH` events

| Event | Coordinate (`position_0`) | Type | Length | Samples | Validation |
|---|---|---:|---:|---:|---|
| `PSV_e80198f305d55977` | `CHM13:chr7:112730158` | DEL | 140,020 | 2 | Assembly-only review |
| `PSV_3497f8fd547eafc0` | `CHM13:chr4:48379366` | INS | 130,411 | 2 | Assembly-only review |
| `PSV_19d2d45f1e54cfe2` | `CHM13:chr4:48379366` | COMPLEX_INDEL | 130,403 | 2 | Assembly-only review |
| `PSV_10382796676f487b` | `CHM13:chr16:34873911` | COMPLEX_INDEL | 126,768 | 2 | Assembly-only review |
| `PSV_b3357e6eaa55bf8f` | `CHM13:chr15:26240728` | DEL | 122,040 | 5 | Assembly-only review |

The two different chr4 representations at the same start position illustrate
why breakpoint inspection and representation harmonization are needed even in
the `HIGH` set.

## 5. HiFi retrieval recommendations

The final ranking contains 420 samples:

| Priority | Samples | Meaning |
|---|---:|---|
| `P1_DISCOVERY_BLIND_SPOT` | 355 | 341 samples with at least one low-tier excluded haplotype, plus 14 missing-mate samples |
| `P1_VALIDATE_CANDIDATE` | 45 | Non-blind-spot samples with prioritized catalog evidence |
| `P3_CONTROL` | 20 | Fully represented, callable graph-member controls with zero catalog carrier loci |

All 420 have `read_access_status=not_checked`, and no HiFi path or accession is
present. All also have `uncallable_haplotype_count=0`; their blind-spot
classification is driven by assembly QC tier or missing-mate status rather than
failure of the graph-callability threshold.

The first six ranked low-QC discovery/validation samples are:

| Rank | Sample | QC tiers | Catalog loci containing sample | `HIGH` loci |
|---:|---|---|---:|---:|
| 1 | `19R0915` | fragmented / fragmented | 7,059 | 246 |
| 2 | `6922020` | not recommended / not recommended | 6,147 | 395 |
| 3 | `9634972` | not recommended / not recommended | 5,162 | 276 |
| 4 | `19R0755` | fragmented / fragmented | 5,715 | 310 |
| 5 | `6053710` | fragmented / fragmented | 5,703 | 303 |
| 6 | `5648514` | not recommended / not recommended | 4,954 | 494 |

Rank 7, sample `2648705`, is the first missing-graph-mate case in the list and
has 2,343 catalog loci, including 336 `HIGH` loci.

Candidate counts are catalog memberships, not confirmed per-sample genotypes.
The ordering also considers QC/blind-spot status and callability, so it is not
strictly descending by candidate count.

## 6. Interpretation

### What the results support

1. **The frozen graph is probably incomplete for this cohort.** Excluded
   assemblies produce more residual observations than graph-member controls,
   many residuals recur across samples, 193,480 catalog loci have support from
   at least one linear caller, and the original graph growth curve was still
   open.

2. **A tractable high-priority subset exists.** The 17,622 `HIGH` loci with an
   explicit linear-caller-supported validation label are the strongest first
   set for breakpoint review and HiFi validation. The remaining 22,354 `HIGH`
   loci provide a second, repeated graph-evidence tier.

3. **Off-reference graph sequence matters.** Nearly 20% of catalog loci use
   graph rather than CHM13 coordinates. Restricting follow-up to linear
   coordinates would discard a substantial part of the discovery signal.

4. **Assembly quality is the largest unresolved evidence problem.** Three
   quarters of loci are supported only by `fragmented_rescue` and/or
   `not_recommended` assemblies and are correctly held as `PENDING_HIFI`.

### What the results do not support

1. **895,941 is not a confirmed novel-SV count.** The catalog includes
   singletons, mapping background, alternate representations of complex loci,
   and evidence from poor assemblies.

2. **Sample support is not allele frequency.** Cohort-wide genotyping has not
   been performed, and absence from an assembly's discovery calls is not a
   reference genotype.

3. **The catalog is not a complete multi-caller SV union.** The merger requires
   at least one non-graph-member Minigraph residual in every emitted cluster.
   Calls found only by SVIM-asm or dipcall are omitted from the final catalog.

4. **The GRCh38 calls do not contribute catalog loci directly.** The reference
   manifest contains all 491 samples for both callers in both reference frames
   (1,964 rows), but final catalog coordinates are only `CHM13` or `GRAPH`.
   There is no cross-reference liftover/harmonization step that can merge an
   hg38-only cluster with a CHM13/graph residual.

5. **Complex SV discovery is incomplete.** Split/complex alignments are
   retained in a separate table, but the catalog contains only DEL, INS, and
   COMPLEX_INDEL. Caller-only INV, DUP, TRA, BND, and CNV evidence is excluded
   by graph-residual gating.

6. **Sex-chromosome SVs are out of scope.** The run used `autosomes_only` and
   excluded X/Y records from the merged catalog.

## 7. Recommended next actions

### Immediate review set

1. Start with the 17,622 `HIGH` + `LINEAR_CALLER_SUPPORTED` loci.
2. Annotate repeat class, segmental duplications, centromeric/telomeric context,
   genes, and known SV catalogs.
3. Inspect the largest events and recurrent short events in the source
   alignments/VCFs; explicitly collapse alternative representations at the same
   complex locus.
4. Keep CHM13 and off-reference graph-coordinate review tracks separate.

### Calibrate the graph-residual background

1. Use the 196 graph-member assemblies as a formal negative/background set.
2. Build per-region and per-event-type background rates rather than relying
   only on the current `GRAPH_MEMBER_CONTROL_OVERLAP` flag.
3. Review the 110,451 overlapping clusters and recurrent control hotspots
   before adding any allele to a new graph.
4. Add positive controls, such as leave-one-assembly-out alleles or synthetic
   spike-ins, so stricter background filtering does not silently remove true
   SVs.

### Correct or clarify catalog scope

If the intended product is a graph-residual catalog, rename/document it
accordingly and keep the current conservative gate.

If the intended product is the planned assembly-only multi-caller union:

1. retain caller-only SVIM-asm/dipcall clusters;
2. harmonize CHM13 and GRCh38 coordinates before cross-reference clustering;
3. add the disabled PAV route after server validation;
4. normalize and reconcile complex SV representations;
5. expose separate evidence tiers for graph-only, linear-only, and concordant
   events.

### Read-based validation and genotyping

1. Resolve data access for the 420 ranked samples.
2. Begin with the 14 missing-mate samples, the 45 dedicated validation samples,
   a QC-diverse subset of the 341 low-tier samples, and all 20 controls.
3. Perform breakpoint-spanning HiFi validation before graph modification.
4. Genotype the reviewed catalog across all 491 samples before estimating
   allele frequencies.
5. Re-screen against a graph-v2 candidate and require the selected residual
   signal to decrease without reducing positive-control recall.

### Reporting improvements

The per-assembly `screen_status` should be refined because it is identical for
all 982 assemblies. A more useful summary would report rates per assembled Gb,
numbers of high/review complex contigs, graph-member-background percentiles,
and explicit threshold-based severity levels.

## 8. Source files used

- `../whole_pangenome/sv_pangenome/novel_sv_discovery/frozen_graph/graph_inventory.tsv`
- `../whole_pangenome/sv_pangenome/novel_sv_discovery/frozen_graph/pangenome_qc_report.md`
- `../whole_pangenome/sv_pangenome/novel_sv_discovery/provenance/discovery_tool_versions.tsv`
- `../whole_pangenome/sv_pangenome/novel_sv_discovery/manifest/assembly_discovery_manifest.tsv`
- `../whole_pangenome/sv_pangenome/novel_sv_discovery/graph_screen/summary/all_assembly_novel_sv_summary.tsv`
- `../whole_pangenome/sv_pangenome/novel_sv_discovery/graph_screen/coordinate_qc/coordinate_qc.tsv`
- `../whole_pangenome/sv_pangenome/novel_sv_discovery/catalog/provisional_graph_residual_sv_catalog.tsv`
- `../whole_pangenome/sv_pangenome/novel_sv_discovery/catalog/all_normalized_assembly_evidence.tsv.gz`
- `../whole_pangenome/sv_pangenome/novel_sv_discovery/reference_calls/reference_call_manifest.tsv`
- `../whole_pangenome/sv_pangenome/novel_sv_discovery/hifi/recommended_hifi_samples.tsv`
- `workflow/rules/novel_sv_discovery.smk`
- `workflow/scripts/screen_novel_graph_svs.py`
- `workflow/scripts/merge_novel_sv_catalog.py`
- `workflow/scripts/recommend_hifi_samples.py`

The exact catalog statistics above were calculated by streaming the complete
895,941-row catalog. The graph summary, discovery manifest, call manifest,
coordinate-QC table, provenance tables, and HiFi ranking were also scanned in
full. The 7.6 GB compressed normalized-evidence archive was inspected for its
schema and example records but was not fully decompressed for an independent
row count; method concordance was calculated from the complete final catalog.

