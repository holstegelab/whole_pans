# Flat reference mapping comparison

Run from `/gpfs/work3/0/qtholstg/pangenome/whole_pans` on Snellius. This is a
separate target in `workflow/rules/personalized_sv_comparison.smk`:

```bash
snakemake -s workflow/rules/personalized_sv_comparison.smk \
  flat_reference_comparison --use-conda --cores 48 \
  --resources mem_mb=288000 disk_mb=500000 \
  --rerun-triggers mtime --dry-run
```

After checking the planned jobs, run the same command without `--dry-run`.
The target uses the graph workflow's prepared DemenTree sample manifest and raw
paired FASTQs. It builds DRAGMAP hashes for the exact CHM13 and
`GRCh38_primary_chr1-22_XY.fa` inputs listed under `pangenome.references` in
`config/config.yaml`. DRAGMAP maps each sample to both references. The hg38
FASTA intentionally has chromosomes 1–22, X and Y and no mitochondrial contig.

Results are under `personalized_sv_comparison.results/flat_reference/`:

- `mapping/{chm13,hg38}/{sample}.bam` and `.bam.bai`: coordinate-sorted flat
  alignments, kept for inspection or downstream analysis.
- `statistics/flat_mapping_metrics.tsv`: one row per sample and flat reference.
- `statistics/graph_minus_flat_per_sample.tsv`: paired graph-minus-flat values
  for each sample, graph, flat reference and metric.
- `statistics/graph_minus_flat_summary.tsv`: mean and median paired differences
  across samples.
- `report.md`: compact summary and interpretation notes.

The comparison reads the existing per-sample graph mapping statistics. If those
are absent, Snakemake schedules graph indexing and mapping as dependencies; it
does not require the separate TMEM insertion-evidence analysis. All mapping
fractions use the same input FASTQ read or base counts. BAMs are retained, so
check storage before running all samples.
The `--rerun-triggers mtime` option avoids rebuilding existing graph products
solely because the rule file changed.

The reference checksum records the configured source file bytes (compressed
for CHM13); hash directory size includes its staged uncompressed FASTA.

MAPQ scores from Giraffe and DRAGMAP have different calibration. Identity uses
GAF matches/block length for graphs and BAM `NM`/CIGAR for flat alignments.
Reference coverage is shown separately because graph BAMs are surjected to the
CHM13 backbone and the hg38 reference has a different sequence and length.
