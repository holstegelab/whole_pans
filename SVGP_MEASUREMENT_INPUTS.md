# Measurement inputs and limits

The native mapping target is `svgp_pilot_mapping`. It maps the configured pilot
people to exactly two graphs. `svgp_pilot` also extracts junction evidence and
applies any calibrated carrier locks. Header-only lock/molecule templates are
valid preparation inputs; they yield uncalibrated evidence or no-calls. They do
not produce reference genotypes. Set the two development extraction settings
explicitly before extracting support.
Full pilot aggregates require 12–20 distinct people and exactly the two frozen
graphs. For the initial single-person check, request that person's explicit
mapping/support output after listing them in `pilot_people`. The inventory
checkpoint resolves aliases first, then makes the actual paired FASTQs mapping
inputs so Snakemake tracks missing or changed files.

`svgp_pilot_supporting_calls` optionally runs `vg pack`/`vg call`. This requires
an approved direct-genotyping review, a stable allele/call crosswalk and verified
reference coordinate paths. The resulting VCF is a supporting callset. Pack depth
counts reads, including duplicates, and is not the molecule evidence used by the
carrier assay. No `vg call -v` or cover-restricted calling is assumed.

Each graph configuration requires `gbz`, `dist`, `minimizer`, `zipcodes`,
`manifest`, `targets`, `graph_sha256`, `node_space_id`, `reference_build`, and
`reference_paths`. Use an empty reference-path list for a carrier-only assay on
a graph with artificial covers. Covers never supply biological phase. The
manifest has one `VERIFIED` row per index component and all four actual
checksums are checked once. `graph_sha256` identifies the exact GBZ used for
mapping, not its ancestral GFA.
The shared distance index must already have no write permission bits: Giraffe
can update writable distance-index mappings. The rule checks this condition and
does not change the source index or its permissions.

For coordinate-based genotyping, also supply `reference_review` and
`genotyping_review`. Every configured reference path must have a sequence
checksum equal to its frozen FASTA chromosome, with a review evidence file.
This workflow consumes an externally prepared verified bundle. It does not
rebuild missing biological paths or invent a reference from a greedy cover.

Use `svgp_assay_targets.example.tsv` for the post-conversion node-space
crosswalk. It has one row per feature and allele. Features `ref`, `left` and
`right` specify an oriented edge with `edge_from` and `edge_to`, for example
`>123` and `<456`; reverse-complement traversals are also counted. These edges
must discriminate the exact complete allele/locus configuration, with genomic
competitors considered during review. Features `flank` and `internal` specify
node intervals (`node_id`, `node_length`, `start0`, `end0`), in forward-node
0-based half-open coordinates. Overlapping intervals within a feature are
merged for both the target-length denominator and covered-base counts.
Unknown features fail. Unverified/missing junctions remain no-calls.

`catalogue_graph_sha256` in both assay targets and measurement locks must match
the panel allele's raw discovery-graph checksum. It links that allele to the
exact GBZ checksum and node space; node IDs are never equated between graphs.
The target `evidence` file must document translation/chopping, orientation,
edge specificity and interval lengths for that frozen bundle.

The extractor streams `vg view -aj`, keeping only panel-overlapping
observations. Consecutive exact alignment matches must extend across each
junction by `anchor_matches` bases on both sides. Low MAPQ and secondary
observations contribute to the ambiguity summary, not allele support. Depth
uses matched graph bases, with overlapping mates/duplicate observations from
the same assigned molecule counted once at each node. Internal depth cannot
replace junction support.

Use `svgp_molecules.example.tsv` to supply reviewed physical-molecule assignments
from a duplicate-handling method suitable for the library/alignment. `read_name`
is the FASTQ fragment name, with only a terminal `/1` or `/2` removed.
`molecule_id` must group both mates and PCR/optical duplicates within the same
`person_id,library_id`. Independent libraries are kept separate. A raw read
name is not proof of a unique molecule. A `VERIFIED` assignment needs `method`
and an existing `evidence` file. Unassigned panel observations make that
sample–allele a no-call. This explicit input is necessary because the old TMEM
helper does not implement general duplicate handling. The workflow does not
pretend to have implemented a general graph-aware PCR duplicate marker.

Carrier thresholds come only from `svgp_measurement_lock.example.tsv`. A locus
advances with `status=LOCKED`, `technical_status=PASS`, an identical graph/node
space/extraction definition and existing calibration, threshold sensitivity,
depth/batch and technical-review files. Set `locked_before_phenotypes=yes`
only after a phenotype-blind freeze. Positive thresholds are mandatory. The
current assay returns `ALT_SUPPORTED`, `REF_SUPPORTED`, or `NO_CALL`. Even
two-sided ALT support is not a calibrated heterozygous genotype. Dosage remains
missing and allele frequency is not calculated. The carrier matrix uses 1/0
only for resolved carrier states and `.` otherwise.

For ONT, list all available people in `measurement.ont.people` and provide a
sample entry with a whole-genome `bam`+`bai`, or `fastqs` and a chemistry-matched
`preset` (`map-ont` or `lr:hq`). Every entry also needs `reference_build`,
`read_scope=whole_genome` and an `alignment_review` evidence file. The known-site
VCF needs the matching reference/build review template. Force calls are generated
by Sniffles, but read inspection is supplied separately in the adjudication TSV;
force calling alone cannot establish truth. Include source/discovery overlap,
development/held-out role, REF and each junction, full-allele spanning status,
depth, sequence agreement and alternatives. SNP-fingerprint confirmation is
required to score ONT/Illumina agreement. An unprovided cell is
`NOT_EVALUABLE`, including a missing Sniffles record.
The fingerprint evidence file must exist and be nonempty; a confirmation label
with only an unavailable evidence path cannot pass identity confirmation.

`svgp_cohort` requires a reviewed cohort release plus locked, callable
main-graph pilot alleles with resolved ONT cells at those alleles. The release
cites actual pilot, ONT, calibration and compute-budget files, a reviewer and
a retention policy.
The pilot and ONT evidence paths must identify the exact files consumed by the
release rule; normalizing relative path spelling does not resolve symlinks or
translate between local and server mount paths.
Configure only technically ready people in `cohort_people`; complete manifests and all paired
lanes must pass before mapping. Discovery donors and relatives can be included
as labelled development/sensitivity measurements. Primary association
eligibility is applied separately. Long evidence and frequency summaries retain
`analysis_group`, keeping independent-primary and development denominators
separate. Main-graph pilot mappings are reused. New cohort
people are mapped to one main graph only. Samples/graphs can run in parallel
within the existing profile's CPU/memory budget; libraries within one person
run sequentially to limit resources and simplify read accounting.

Current releases checked from official release pages on 20 October 2026:
[vg 1.77.0](https://github.com/vgteam/vg/releases/tag/v1.77.0),
[minimap2 2.31](https://github.com/lh3/minimap2/releases/tag/v2.31),
[Sniffles 2.8.1](https://github.com/fritzsedlazeck/Sniffles/releases/tag/v2.8.1),
[SAMtools 1.24](https://github.com/samtools/samtools/releases/tag/1.24),
[BCFtools 1.24](https://github.com/samtools/bcftools/releases/tag/1.24).
The dedicated environment does not upgrade the existing environments.
Consult the [Giraffe guide](https://github.com/vgteam/vg/wiki/Mapping-short-reads-with-Giraffe),
[graph genotyping guide](https://github.com/vgteam/vg/wiki/SV-Genotyping-and-variant-calling),
[minimap2 presets](https://github.com/lh3/minimap2#map-long-genomic) and
[Sniffles force-genotyping instructions](https://github.com/fritzsedlazeck/Sniffles#d-genotyping-a-known-set-of-svs-force-calling).

No tools, tests, syntax checks or dry runs were executed locally. All execution
and environment solving belong on the server under the repository instructions.
