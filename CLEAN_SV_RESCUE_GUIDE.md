# Clean HiFi rescue sequences for novel SVs

Use `workflow/scripts/build_clean_sv_rescues.py` to replace a questionable
assembly interval with a local sequence reconstructed from the original HiFi
reads. The script does not modify the frozen graph. It creates auditable
candidate sequences that can be screened before building an incremental test
graph.

The intended use is narrow:

- rescue an SV allele carried by a fragmented or otherwise excluded assembly;
- use one representative local sequence per validated allele;
- require the read-derived contig to span both sides of the SV;
- never add `candidate_rescue.fa` to a graph when its `qc.tsv` status is
  `FAIL`.

Do not replace good whole assemblies with local patches. Keep complete,
well-QC'd assemblies as whole-assembly graph inputs. Use these rescue
sequences only when the assembly containing the allele is not safe to insert
as a whole.

## 1. Environment

From the `whole_pans` directory:

```bash
conda env create -f workflow/envs/sv_rescue.yaml
conda activate whole-pans-sv-rescue
```

The workflow uses `minimap2`, `hifiasm`, and `samtools`. It expects local HiFi
reads in ordinary FASTA or four-line FASTQ, optionally gzip-compressed. If the
download is a PacBio BAM, convert it once:

```bash
samtools fastq -@ 16 sample.hifi.bam |
  bgzip -@ 16 > sample.hifi.fastq.gz
```

The assembly FASTA files must already have `samtools faidx` indexes:

```bash
samtools faidx /path/to/assembly.clean.fa.gz
```

## 2. Read manifest

Create a sample-level TSV, for example `config/hifi_reads.tsv`:

```text
sample_id	hifi_path
17R1158	/data/hifi/17R1158.hifi.fastq.gz
9058208	/data/hifi/9058208.hifi.fastq.gz
```

Use exactly one row per sample. Multiple FASTA/FASTQ files for one sample can
be placed in one `hifi_path` cell separated by semicolons. Paths must refer to
the machine where the rescue command will run.

## 3. Build a conservative plan

The default planner selects at most 100 recurrent `UNCERTAIN`,
`PENDING_HIFI` catalog events supported by graph residuals and both linear
assembly callers. It seeks source alleles specifically in
`fragmented_rescue` and `not_recommended` assemblies.

```bash
RESCUE_ROOT=../whole_pangenome/sv_pangenome/novel_sv_discovery/rescue

python workflow/scripts/build_clean_sv_rescues.py plan \
  --catalog ../whole_pangenome/sv_pangenome/novel_sv_discovery/catalog/provisional_graph_residual_sv_catalog.tsv \
  --graph-candidates ../whole_pangenome/sv_pangenome/novel_sv_discovery/graph_screen/coordinate_qc/all_residual_sv_candidates.annotated.tsv.gz \
  --assembly-manifest ../whole_pangenome/sv_pangenome/novel_sv_discovery/manifest/assembly_discovery_manifest.tsv \
  --read-manifest config/hifi_reads.tsv \
  --output "$RESCUE_ROOT/plan.tsv"
```

Planning streams the large catalog and graph-candidate files, so run it on the
server where those files are available. It does not load them completely into
memory.

To rescue reviewed catalog IDs rather than use the default filters:

```bash
python workflow/scripts/build_clean_sv_rescues.py plan \
  --catalog ../whole_pangenome/sv_pangenome/novel_sv_discovery/catalog/provisional_graph_residual_sv_catalog.tsv \
  --graph-candidates ../whole_pangenome/sv_pangenome/novel_sv_discovery/graph_screen/coordinate_qc/all_residual_sv_candidates.annotated.tsv.gz \
  --assembly-manifest ../whole_pangenome/sv_pangenome/novel_sv_discovery/manifest/assembly_discovery_manifest.tsv \
  --read-manifest config/hifi_reads.tsv \
  --event-ids reviewed_psv_ids.txt \
  --source-qc-tiers all \
  --output "$RESCUE_ROOT/plan.reviewed.tsv"
```

Review `plan_status` before running. Only `READY` rows can run. Common other
states are `WAITING_FOR_HIFI`, `SOURCE_NEAR_CONTIG_END`, and
`NO_MATCHING_SOURCE`. Coordinates in the plan are 0-based, half-open.

## 4. Run local assemblies

Run one rescue:

```bash
python workflow/scripts/build_clean_sv_rescues.py run \
  --plan "$RESCUE_ROOT/plan.tsv" \
  --rescue-id PSV_example \
  --output-dir "$RESCUE_ROOT/results" \
  --threads 16
```

Or run every `READY` row sequentially:

```bash
awk -F '\t' '
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      if ($i == "rescue_id") id = i
      if ($i == "plan_status") status = i
    }
    next
  }
  $status == "READY" {print $id}
' "$RESCUE_ROOT/plan.tsv" |
while read -r rescue_id; do
  python workflow/scripts/build_clean_sv_rescues.py run \
    --plan "$RESCUE_ROOT/plan.tsv" \
    --rescue-id "$rescue_id" \
    --output-dir "$RESCUE_ROOT/results" \
    --threads 16
done
```

Submit the same one-rescue command as an array through the cluster profile if
parallel execution is preferable. The script itself contains no
scheduler-specific settings.

Each result directory contains:

- `seed.fa`: the source assembly interval, normally with 50 kb flanks;
- `reads_to_seed.paf` and `selected_reads.fa.gz`: recruitment evidence;
- `local_contigs.fa` and `local_contigs_to_seed.paf`: local assembly evidence;
- `candidate_rescue.fa`: best two-sided local contig, even if final QC fails;
- `clean_rescue.fa`: written only for a `PASS`;
- `qc.tsv`, `run_status.json`, tool versions, and logs.

If read support is insufficient or no two-sided local contig is recovered,
`qc.tsv` records `FAIL` and no candidate FASTA may be present. Tool or malformed
input failures are instead recorded as `ERROR` in `run_status.json`.

A default run requires at least three reads spanning 1 kb on both sides and
having no more than 30 bp of local indel disagreement with the source allele.
A final `PASS` also requires an assembled contig spanning 10 kb on both sides,
MAPQ at least 20, alignment identity at least 0.99, no `N` bases, at least
20 kb of rescued sequence, and no more than 20 bp of local indel disagreement
with the source allele. These are technical filters, not proof that the allele
is biologically correct.

Result directories are never overwritten. If an interrupted or erroneous run
must be repeated, use a new output directory or archive the earlier result
first.

## 5. Combine passing rescues

```bash
python workflow/scripts/build_clean_sv_rescues.py summarize \
  --plan "$RESCUE_ROOT/plan.tsv" \
  --results-dir "$RESCUE_ROOT/results" \
  --summary "$RESCUE_ROOT/rescue_summary.tsv" \
  --fasta "$RESCUE_ROOT/clean_rescues.fa" \
  --graph-inputs "$RESCUE_ROOT/clean_rescue_graph_inputs.tsv"
```

The combined FASTA contains only unique sequences whose `qc.tsv` status is
`PASS`. Exact sequence duplicates are recorded but included only once.

Before adding `clean_rescues.fa` to a graph:

1. confirm the event with read alignments and, where possible, a second
   carrier;
2. remap the rescue to CHM13 and the frozen graph and inspect its placement;
3. build an incremental test graph and check that the expected residual SV is
   absorbed without unrelated graph growth;
4. retain `clean_rescue_graph_inputs.tsv` as the provenance record.

These are curated allele-plus-flank inputs, not substitutes for complete
assemblies.
