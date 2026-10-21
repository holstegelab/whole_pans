# September review: immediate follow-up

This implements the immediate next phase of
`../whole_pangenome_review_2026-09-15.txt`. It prepares corrected QC, donor
provenance, a frozen TMEM comparison and a small sequence-resolved validation
panel. It does not run cohort genotyping or association analysis.

## Run on the server

Start in `/gpfs/work3/0/qtholstg/pangenome/whole_pans`. The standalone entry
point avoids loading the assembly/decontamination workflow:

```bash
snakemake -s workflow/rules/review_followup.smk -n -p review_followup
snakemake -s workflow/rules/review_followup.smk --use-conda --cores 8 review_followup
```

Use your existing execution profile for cluster submission. No new scheduler
settings are embedded in these rules. If the dry run schedules unfinished
whole-genome TMEM mapping/indexing, those jobs retain their existing resource
requirements; use the appropriate cluster profile rather than a login node.
The same `review_followup` target is available from the main Snakefile; its
default assembly-QC target remains unchanged.

Individual stages can be requested as follows:

```bash
snakemake -s workflow/rules/review_followup.smk -n -p review_corrected_v2_qc
snakemake -s workflow/rules/review_followup.smk -n -p review_hprc_donor_provenance
snakemake -s workflow/rules/review_followup.smk -n -p review_freeze_tmem
snakemake -s workflow/rules/review_followup.smk -n -p review_sequence_resolved_panel
```

Paths and panel settings are in `config/review_followup.yaml`. Existing graph,
reference, read and TMEM paths are inherited from `config/config.yaml`.
No large datasets are downloaded, moved or removed by the new rules. Graph
construction and HiFi rescue are treated as existing inputs. Completing the
existing TMEM analysis can run its missing upstream jobs.

## Outputs and interpretation

New results go to `../whole_pangenome/review_followup_2026-09-18/`:

| Output | Purpose |
| --- | --- |
| `v2_qc/` | Regenerated bubble coordinates and rescue/cohort QC with the bytes parser fixed. Existing graph-integrity results are reused. |
| `hprc_provenance/` | Corrected copies of all four example tables, actual sample-grouped input order, and insertion sequences. Old donor labels remain in `previous_*` columns. |
| `tmem_frozen/` | Completed configured sample/graph evidence, report, membership, mapping summaries, targets, scripts, configuration and checksums. |
| `panel/panel.tsv` | Up to 30 loci by default, including TMEM as one multiallelic locus. |
| `panel/*.queries.fa` | REF and ALT queries with flanks, separated by reference build. |
| `panel/placement/` | Whole-genome minimap2 PAFs with secondary hits, CIGARs, cs tags and tool versions. |
| `panel/validation_worksheet.tsv` | Blank read/source-assembly validation worksheet; copy it before editing. |
| `report.md`, `checksums.sha256` | Scope and interpretation notes, plus checksums of the new preparation outputs. |

V2 correction reads the existing discovery catalog and normalized evidence;
this can take time. HPRC provenance scans the combined graph and bubble catalog.
Panel preparation streams CHM13 once to retrieve small reference windows.
Placement then indexes/maps against CHM13 and GRCh38 separately. These are
server analyses, not lightweight local checks.

The first local source rank defaults to **464**, as recorded for the graph
reviewed on September 15. Donors are assigned in first-seen sample order, with
that sample's assemblies grouped exactly as in `hprc_plus_cohort.smk`. Each
assignment must also match the graph's source contig, rank and offset. A mismatch
stops the rule instead of guessing a donor. This is introducing-donor
provenance, not a carrier list or allele frequency.

Simple insertion bubbles are checked against their GFA nodes. The FASTA uses
the bubble traversal orientation, which can be the reverse complement of the
node sequence. The exact insertion and zero-length reference interval are
preserved; distinct alleles are not merged by breakpoint proximity.

Automatic selection reserves four difficult/repeat-review candidates and up
to 15 additional genic candidates, then fills the remaining slots from Tier A/B
candidates. With no BEDs, low 5-mer diversity is only a repeat-review proxy.
Optional CHM13 `repeat_beds` and `segdup_beds` annotate reference flanks;
unprovided annotations are `NOT_ASSESSED`. They do not annotate novel inserted
sequence. TMEM retains GRCh38 coordinates. No liftover is performed.

To choose particular candidates, set `panel.selection` to a TSV with columns
`candidate_id` and `panel_role`; use `candidate`, `genic_candidate` or
`repeat_review_control`. Select at most `max_loci - 1` HPRC candidate IDs; TMEM
is added separately. The complete corrected candidate table is kept regardless
of panel selection. The default panel is an insertion-focused starting set,
not a representative sample of every residual SV class.

Freezing requires exactly one evidence row for every configured sample/graph
pair (currently 20 × 4), graph-checksum agreement and callable corrected targets.
Uninformative sample-level evidence remains an evidence category, never a
reference genotype. Prior high-absence samples are listed separately as
negative assay controls with an explicit non-truth label. The snapshot is a
protected Snakemake output. Use a new `review_followup.results` directory when
creating a revised snapshot. Graph checksums come from the existing inventory;
the freeze stage hashes copied artifacts, not the large graphs again.

Placement status is a review aid. A query-end-spanning alignment can still
contain a large insertion/deletion; inspect the CIGAR, alternative hits,
junctions and repeats. An unmapped ALT query is unresolved, not proof of
absence. None of these stages declares independent HiFi validation, V2 allele
capture, calibrated genotypes or external-resource novelty. Those require
the source-assembly/read review recorded in the worksheet.

The original `sv_examples/` donor labels and V2 `qc/tables/bubble_delta.tsv`
are superseded for interpretation by the new corrected copies. Old direct TMEM
outputs outside `analysis_v2/` and the July novel-SV report remain historical.
The old files are retained for provenance.

## Server verification

Local tests/dry runs were not performed, in accordance with `AGENTS.md`.
Small regression fixtures cover byte coordinates, donor grouping, allele
orientation, streamed reference windows and incomplete snapshot rejection:

```bash
python -m unittest discover -s tests -p 'test_review_followup.py'
snakemake -s workflow/rules/review_followup.smk -n -p review_followup
```

The implementation follows the
[Minigraph bubble format](https://github.com/lh3/minigraph#calling-structural-variations)
and [minimap2 assembly-alignment options](https://github.com/lh3/minimap2#full-genomeassembly-alignment).
