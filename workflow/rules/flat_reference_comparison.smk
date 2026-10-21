"""Map the personalized-SV pilot reads to the two graph input FASTAs."""

# Included at the end of personalized_sv_comparison.smk so the manifest,
# samples, FASTQ counts, and per-graph mapping statistics are shared.
FRC = PSC["flat_reference"]
FRC_OUT = f"{PSC_OUT}/flat_reference"
FRC_ENV = str(PSC_ROOT / "workflow/envs/flat_reference_comparison.yaml")
FRC_SCRIPT = str(PSC_ROOT / "workflow/scripts/flat_reference_comparison.py")
FRC_REFERENCES = {
    "chm13": config["pangenome"]["references"]["CHM13"],
    "hg38": config["pangenome"]["references"]["hg38"],
}
FRC_REFERENCE_RE = "(?:chm13|hg38)"
FRC_FLAT_METRICS = f"{FRC_OUT}/statistics/flat_mapping_metrics.tsv"
FRC_PAIRS = f"{FRC_OUT}/statistics/graph_minus_flat_per_sample.tsv"
FRC_SUMMARY = f"{FRC_OUT}/statistics/graph_minus_flat_summary.tsv"
FRC_REPORT = f"{FRC_OUT}/report.md"
if int(FRC.get("mapping_threads", 24)) < 5:
    raise WorkflowError("personalized_sv_comparison.flat_reference.mapping_threads must be at least 5")


def frc_reference(wildcards):
    return FRC_REFERENCES[wildcards.reference]


def frc_mapping_jobs():
    return expand(
        f"{FRC_OUT}/statistics/mapping/{{reference}}/{{sample}}.tsv",
        reference=tuple(FRC_REFERENCES),
        sample=PSC_SAMPLES,
    )


rule flat_reference_comparison:
    input:
        FRC_REPORT


rule build_flat_dragmap_hash:
    input:
        reference=frc_reference
    output:
        hash=directory(f"{FRC_OUT}/indexes/{{reference}}"),
        info=f"{FRC_OUT}/indexes/{{reference}}.reference.tsv"
    log:
        f"{FRC_OUT}/logs/indexes/{{reference}}.log"
    benchmark:
        f"{FRC_OUT}/benchmarks/indexes/{{reference}}.tsv"
    conda:
        FRC_ENV
    threads:
        int(FRC.get("index_threads", 16))
    resources:
        mem_mb=int(FRC.get("index_mem_mb", 96000)),
        runtime_min=int(FRC.get("index_runtime_min", 1440)),
        disk_mb=int(FRC.get("index_disk_mb", 100000))
    wildcard_constraints:
        reference=FRC_REFERENCE_RE
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output.hash:q} "$(dirname {log:q})"
        # DRAGMAP hashes an ordinary FASTA. Keep the staged FASTA with its hash
        # so the hash directory remains self-contained after this job exits.
        case {input.reference:q} in
          *.gz) gzip -cd -- {input.reference:q} > {output.hash:q}/reference.fa ;;
          *) cp -- {input.reference:q} {output.hash:q}/reference.fa ;;
        esac
        test -s {output.hash:q}/reference.fa
        dragen-os --build-hash-table true --ht-num-threads {threads} \
          --ht-reference {output.hash:q}/reference.fa \
          --output-directory {output.hash:q} > {log:q} 2>&1
        test -n "$(find {output.hash:q} -type f ! -name reference.fa -print -quit)"
        source_sha=$(sha256sum {input.reference:q} | cut -d' ' -f1)
        source_size=$(stat -c '%s' {input.reference:q})
        hash_size=$(du -sb {output.hash:q} | cut -f1)
        printf 'target\treference_fasta\treference_sha256\treference_size_bytes\thash_size_bytes\n' > {output.info:q}.tmp
        printf '%s\t%s\t%s\t%s\t%s\n' \
          {wildcards.reference:q} {input.reference:q} "$source_sha" "$source_size" "$hash_size" \
          >> {output.info:q}.tmp
        mv {output.info:q}.tmp {output.info:q}
        """


rule map_sample_to_flat_reference:
    input:
        manifest=psc_read_manifest,
        hash=f"{FRC_OUT}/indexes/{{reference}}",
        info=f"{FRC_OUT}/indexes/{{reference}}.reference.tsv",
        r1=lambda wc: psc_reads(wc, "r1"),
        r2=lambda wc: psc_reads(wc, "r2")
    output:
        bam=f"{FRC_OUT}/mapping/{{reference}}/{{sample}}.bam",
        bai=f"{FRC_OUT}/mapping/{{reference}}/{{sample}}.bam.bai"
    log:
        f"{FRC_OUT}/logs/mapping/{{reference}}/{{sample}}.log"
    benchmark:
        f"{FRC_OUT}/benchmarks/mapping/{{reference}}/{{sample}}.tsv"
    conda:
        FRC_ENV
    threads:
        int(FRC.get("mapping_threads", 24))
    resources:
        mem_mb=int(FRC.get("mapping_mem_mb", 96000)),
        runtime_min=int(FRC.get("mapping_runtime_min", 1440)),
        disk_mb=int(FRC.get("mapping_disk_mb", 150000))
    wildcard_constraints:
        reference=FRC_REFERENCE_RE,
        sample=PSC_SAMPLE_RE
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.bam:q})" "$(dirname {log:q})"
        work=$(mktemp -d -p "${{TMPDIR:-/tmp}}" "flat-{wildcards.reference}-{wildcards.sample}.XXXXXX")
        bam_tmp={output.bam:q}.tmp
        bai_tmp={output.bai:q}.tmp
        : > {log:q}
        trap 'rm -rf -- "$work"; rm -f -- "$bam_tmp" "$bai_tmp"' EXIT
        r1_files=({input.r1:q})
        r2_files=({input.r2:q})
        if [ "${{#r1_files[@]}}" -eq 0 ] || [ "${{#r1_files[@]}}" -ne "${{#r2_files[@]}}" ]; then
          echo "Expected matched FASTQ pairs for {wildcards.sample}" >> {log:q}
          exit 1
        fi
        mapper_threads=$(({threads} - 4))
        lane_bams=()
        for ((i=0; i<${{#r1_files[@]}}; i++)); do
          lane_bam="$work/lane_$i.bam"
          if [ "${{#r1_files[@]}}" -eq 1 ]; then lane_bam="$bam_tmp"; fi
          dragen-os -r {input.hash:q} \
            -1 "${{r1_files[i]}}" -2 "${{r2_files[i]}}" \
            --RGID "{wildcards.sample}.$((i + 1))" --RGSM {wildcards.sample:q} \
            --num-threads "$mapper_threads" 2>> {log:q} \
            | samtools sort -@ 3 -m 1G -T "$work/sort_$i" -o "$lane_bam" - 2>> {log:q}
          lane_bams+=("$lane_bam")
        done
        if [ "${{#lane_bams[@]}}" -gt 1 ]; then
          samtools merge -@ 4 -o "$bam_tmp" "${{lane_bams[@]}}" 2>> {log:q}
        fi
        samtools index -@ 4 "$bam_tmp" "$bai_tmp" 2>> {log:q}
        mv "$bam_tmp" {output.bam:q}
        mv "$bai_tmp" {output.bai:q}
        trap - EXIT
        rm -rf -- "$work"
        """


rule calculate_flat_mapping_metrics:
    input:
        script=FRC_SCRIPT,
        graph_stats_script=PSC_STATS,
        manifest=psc_read_manifest,
        bam=f"{FRC_OUT}/mapping/{{reference}}/{{sample}}.bam",
        bai=f"{FRC_OUT}/mapping/{{reference}}/{{sample}}.bam.bai",
        fastq=f"{PSC_OUT}/fastq_stats/{{sample}}.tsv",
        info=f"{FRC_OUT}/indexes/{{reference}}.reference.tsv",
        mapping_benchmark=f"{FRC_OUT}/benchmarks/mapping/{{reference}}/{{sample}}.tsv",
        index_benchmark=f"{FRC_OUT}/benchmarks/indexes/{{reference}}.tsv"
    output:
        f"{FRC_OUT}/statistics/mapping/{{reference}}/{{sample}}.tsv"
    params:
        group=psc_group,
        batch=psc_batch,
        dragmap_version=str(FRC.get("dragmap_version", "1.3.0"))
    log:
        f"{FRC_OUT}/logs/statistics/mapping/{{reference}}/{{sample}}.log"
    conda:
        FRC_ENV
    threads: 1
    resources:
        mem_mb=8000,
        runtime_min=360
    wildcard_constraints:
        reference=FRC_REFERENCE_RE,
        sample=PSC_SAMPLE_RE
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        work=$(mktemp -d -p "${{TMPDIR:-/tmp}}" "flat-stats.XXXXXX")
        trap 'rm -rf -- "$work"' EXIT
        : > {log:q}
        samtools coverage {input.bam:q} > "$work/coverage.tsv" 2>> {log:q}
        samtools flagstat {input.bam:q} > "$work/flagstat.txt" 2>> {log:q}
        python {input.script:q} metrics --bam {input.bam:q} \
          --fastq-stats {input.fastq:q} --coverage "$work/coverage.tsv" \
          --flagstat "$work/flagstat.txt" --reference-info {input.info:q} \
          --mapping-benchmark {input.mapping_benchmark:q} \
          --index-benchmark {input.index_benchmark:q} \
          --sample {wildcards.sample:q} --sample-group {params.group:q} \
          --batch {params.batch:q} --target {wildcards.reference:q} \
          --mapper-version {params.dragmap_version:q} --output {output:q} \
          >> {log:q} 2>&1
        """


rule compare_flat_and_graph_mapping:
    input:
        script=FRC_SCRIPT,
        graph_stats_script=PSC_STATS,
        samples=PSC_SAMPLES_TSV,
        flat=frc_mapping_jobs(),
        graph=psc_mapping_jobs()
    output:
        flat=FRC_FLAT_METRICS,
        pairs=FRC_PAIRS,
        summary=FRC_SUMMARY,
        report=FRC_REPORT
    params:
        flat=lambda wc, input: " ".join("--flat " + shlex.quote(str(path)) for path in input.flat),
        graph=lambda wc, input: " ".join("--graph " + shlex.quote(str(path)) for path in input.graph)
    log:
        f"{FRC_OUT}/logs/compare.log"
    conda:
        FRC_ENV
    threads: 1
    resources:
        mem_mb=4000,
        runtime_min=60
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.report:q})" "$(dirname {log:q})"
        python {input.script:q} compare --samples {input.samples:q} \
          {params.flat} {params.graph} \
          --flat-output {output.flat:q} --pairs-output {output.pairs:q} \
          --summary-output {output.summary:q} --report-output {output.report:q} \
          > {log:q} 2>&1
        """
