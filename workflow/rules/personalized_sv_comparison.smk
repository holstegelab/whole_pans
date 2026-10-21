"""Direct short-read comparison of four whole-genome SV graphs."""

configfile: workflow.source_path("../../config/config.yaml")


import csv
import re
import shlex
from pathlib import Path


PSC_ROOT = Path(workflow.current_basedir).parents[1]
workdir: str(PSC_ROOT)
shell.executable("/bin/bash")

PSC = config["personalized_sv_comparison"]
PSC_OUT = PSC["results"]
PSC_ANALYSIS = PSC.get("analysis_results", f"{PSC_OUT}/analysis_v2")
PSC_ENV = str(PSC_ROOT / "workflow/envs/personalized_sv_comparison.yaml")
PSC_STATS = str(PSC_ROOT / "workflow/scripts/personalized_sv_statistics.py")
PSC_TMEM = str(PSC_ROOT / "workflow/scripts/tmem_whole_graph_target.py")

PSC_GRAPHS = tuple(PSC["source_graphs"])
PSC_GROUPS = PSC["sample_groups"]
PSC_SAMPLES = tuple(sample for group in PSC_GROUPS.values() for sample in group)
if len(PSC_GRAPHS) != 4:
    raise WorkflowError("personalized_sv_comparison.source_graphs must contain four graphs")
if len(PSC_SAMPLES) != len(set(PSC_SAMPLES)):
    raise WorkflowError("a pilot sample occurs in more than one sample group")

PSC_SAMPLE_GROUP = {
    sample: group for group, samples in PSC_GROUPS.items() for sample in samples
}
PSC_GRAPH_RE = "(?:" + "|".join(map(re.escape, PSC_GRAPHS)) + ")"
PSC_SAMPLE_RE = "(?:" + "|".join(map(re.escape, PSC_SAMPLES)) + ")"

PSC_SAMPLES_TSV = f"{PSC_OUT}/samples.tsv"
PSC_INVENTORY = f"{PSC_OUT}/graph_inventory.tsv"
PSC_MEMBERSHIP = f"{PSC_OUT}/sample_graph_membership.tsv"
PSC_TMEM_TARGETS = f"{PSC_ANALYSIS}/tmem_targets.tsv"
PSC_TMEM_VALIDATION = f"{PSC_ANALYSIS}/tmem_pilot_validation.tsv"
PSC_TMEM_FASTA = f"{PSC_ANALYSIS}/tmem/allele_queries.fa"
PSC_RAW_CATALOG = f"{PSC_OUT}/sv_catalogs.tsv"
PSC_SV_CATALOG = f"{PSC_ANALYSIS}/sv_catalogs.tsv"
PSC_GRAPH_TARGETS = f"{PSC_ANALYSIS}/tmem_graph_targets.tsv"
PSC_MAPPING_METRICS = f"{PSC_ANALYSIS}/mapping_metrics.tsv"
PSC_SV_EVIDENCE = f"{PSC_ANALYSIS}/sv_evidence.tsv"
PSC_TMEM_EVIDENCE = f"{PSC_ANALYSIS}/tmem_19kb_evidence.tsv"
PSC_COMPARISON = f"{PSC_ANALYSIS}/paired_comparison.tsv"
PSC_REPORT = f"{PSC_ANALYSIS}/report.md"


def psc_source_graph(wildcards):
    return PSC["source_graphs"][wildcards.target]


def psc_read_manifest(wildcards=None):
    return str(checkpoints.prepare_dementree_manifest.get().output.manifest)


def psc_sample_row(wildcards):
    with open(psc_read_manifest(), newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["sample_id"] == wildcards.sample:
                return row
    raise WorkflowError(f"sample is absent from prepared manifest: {wildcards.sample}")


def psc_reads(wildcards, mate):
    field = "r1_paths" if mate == "r1" else "r2_paths"
    return [path for path in psc_sample_row(wildcards)[field].split(",") if path]


def psc_group(wildcards):
    return psc_sample_row(wildcards)["sample_group"]


def psc_batch(wildcards):
    return psc_sample_row(wildcards)["batch"]


def psc_inventory_jobs():
    return expand(f"{PSC_OUT}/inventory/{{target}}.tsv", target=PSC_GRAPHS)


def psc_mapping_jobs():
    return expand(
        f"{PSC_OUT}/statistics/mapping/{{target}}/{{sample}}.tsv",
        target=PSC_GRAPHS,
        sample=PSC_SAMPLES,
    )


def psc_evidence_jobs():
    return expand(
        f"{PSC_ANALYSIS}/statistics/sv_evidence/{{target}}/{{sample}}.tsv",
        target=PSC_GRAPHS,
        sample=PSC_SAMPLES,
    )


rule personalized_sv_comparison:
    input:
        PSC_REPORT


checkpoint prepare_dementree_manifest:
    input:
        manifests=PSC["fastq_manifests"]
    output:
        manifest=PSC_SAMPLES_TSV
    params:
        columns=PSC["manifest_columns"],
        groups=";".join(
            f"{sample}={group}" for group, samples in PSC_GROUPS.items() for sample in samples
        ),
        expected=";".join(PSC_SAMPLES),
        translations=";".join(
            f"{item['from']}={item['to']}" for item in PSC.get("path_translations", [])
        )
    log:
        f"{PSC_OUT}/logs/prepare_dementree_manifest.log"
    conda:
        PSC_ENV
    threads: 1
    resources:
        mem_mb=1000,
        runtime_min=30
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.manifest:q})" "$(dirname {log:q})"
        body={output.manifest:q}.body.tmp
        final={output.manifest:q}.tmp
        trap 'rm -f -- "$body" "$final"' EXIT

        awk -F '\t' -v OFS='\t' \
          -v batch_col={params.columns[batch]} \
          -v sample_col={params.columns[sample]} \
          -v type_col={params.columns[read_type]} \
          -v sex_col={params.columns[sex]} \
          -v r1_col={params.columns[r1]} \
          -v r2_col={params.columns[r2]} \
          -v groups={params.groups:q} -v expected={params.expected:q} \
          -v translations={params.translations:q} '
          function translate_paths(value,    n, path, i, j, result) {{
            n = split(value, path, ",")
            for (i = 1; i <= n; i++) {{
              for (j = 1; j <= n_translations; j++) {{
                if (old_prefix[j] != "" && index(path[i], old_prefix[j]) == 1) {{
                  sub("^" old_prefix[j], new_prefix[j], path[i])
                  break
                }}
              }}
              result = result (i == 1 ? "" : ",") path[i]
            }}
            return result
          }}
          BEGIN {{
            n = split(groups, entries, ";")
            for (i = 1; i <= n; i++) {{
              split(entries[i], pair, "=")
              group[pair[1]] = pair[2]
            }}
            n = split(expected, entries, ";")
            for (i = 1; i <= n; i++) wanted[entries[i]] = 1
            n_translations = split(translations, entries, ";")
            for (i = 1; i <= n_translations; i++) {{
              split(entries[i], pair, "=")
              old_prefix[i] = pair[1]
              new_prefix[i] = pair[2]
            }}
          }}
          {{
            sample = $sample_col
            if (!(sample in wanted)) next
            batch = FILENAME
            sub(/^.*NL_DemenTree_/, "", batch)
            sub(/\.tsv$/, "", batch)
            read_type = $type_col
            sex = sex_col > 0 && $sex_col != "" ? $sex_col : "unknown"
            r1 = translate_paths($r1_col)
            r2 = translate_paths($r2_col)
            fingerprint = batch OFS read_type OFS sex OFS r1 OFS r2
            if (sample in seen) {{
              if (seen[sample] != fingerprint) {{
                print "Conflicting duplicate sample: " sample > "/dev/stderr"
                failed = 1
              }}
              next
            }}
            if (read_type != "fastq_paired" || r1 == "" || r2 == "") {{
              print "Invalid paired-read row for " sample > "/dev/stderr"
              failed = 1
              next
            }}
            seen[sample] = fingerprint
            print sample, batch, read_type, sex, r1, r2, group[sample]
          }}
          END {{
            for (sample in wanted) if (!(sample in seen)) {{
              print "Pilot sample is absent from manifests: " sample > "/dev/stderr"
              failed = 1
            }}
            if (failed) exit 1
          }}
        ' {input.manifests:q} | sort -t $'\t' -k7,7 -k1,1 > "$body"

        while IFS=$'\t' read -r sample batch read_type sex r1 r2 group; do
          IFS=',' read -ra r1_files <<< "$r1"
          IFS=',' read -ra r2_files <<< "$r2"
          if [ "${{#r1_files[@]}}" -ne "${{#r2_files[@]}}" ]; then
            echo "$sample has unequal numbers of R1 and R2 files" >&2
            exit 1
          fi
          for reads in "${{r1_files[@]}}" "${{r2_files[@]}}"; do
            if [ ! -r "$reads" ]; then
              echo "FASTQ is not readable for $sample: $reads" >&2
              exit 1
            fi
          done
        done < "$body"

        printf 'sample_id\tbatch\tread_type\tsex\tr1_paths\tr2_paths\tsample_group\n' > "$final"
        cat "$body" >> "$final"
        mv "$final" {output.manifest:q}
        rm -f -- "$body"
        trap - EXIT
        """


rule inventory_source_graph:
    input:
        graph=psc_source_graph
    output:
        inventory=f"{PSC_OUT}/inventory/{{target}}.tsv",
        bubbles=f"{PSC_OUT}/sv_catalogs/{{target}}.bubbles.tsv"
    params:
        graph_type=lambda wc: PSC["graph_metadata"][wc.target]["graph_type"],
        reference=lambda wc: PSC["graph_metadata"][wc.target]["reference_build"]
    log:
        f"{PSC_OUT}/logs/inventory/{{target}}.log"
    benchmark:
        f"{PSC_OUT}/benchmarks/inventory/{{target}}.tsv"
    conda:
        PSC_ENV
    threads: 1
    resources:
        mem_mb=int(PSC["resources"].get("inventory_mem_mb", 32000)),
        runtime_min=int(PSC["resources"].get("inventory_runtime_min", 240)),
        disk_mb=int(PSC["resources"].get("inventory_disk_mb", 10000))
    wildcard_constraints:
        target=PSC_GRAPH_RE
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.inventory:q})" "$(dirname {output.bubbles:q})" "$(dirname {log:q})"
        test -s {input.graph:q} || {{ echo "Missing source graph: {input.graph}" >&2; exit 1; }}
        work=$(mktemp -d -p "${{TMPDIR:-/tmp}}" "psc-inventory-{wildcards.target}.XXXXXX")
        trap 'rm -rf -- "$work"' EXIT
        case {input.graph:q} in
          *.gz) gzip -cd -- {input.graph:q} > "$work/graph.gfa" ;;
          *) cp -- {input.graph:q} "$work/graph.gfa" ;;
        esac
        test -s "$work/graph.gfa"
        read -r segments links paths segment_bp < <(
          awk -F '\t' '
            $1 == "S" {{segments++; bp += length($3)}}
            $1 == "L" {{links++}}
            $1 == "P" || $1 == "W" {{paths++}}
            END {{print segments + 0, links + 0, paths + 0, bp + 0}}
          ' "$work/graph.gfa"
        )
        gfatools bubble "$work/graph.gfa" > {output.bubbles:q} 2> {log:q}
        checksum=$(sha256sum {input.graph:q} | cut -d' ' -f1)
        printf 'target\tsource_graph\tgraph_checksum\treference_build\tgraph_type\tsize_bytes\tsegments\tlinks\tembedded_paths\tsegment_bp\tsv_bubbles\n' > {output.inventory:q}
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
          {wildcards.target:q} {input.graph:q} "$checksum" {params.reference:q} {params.graph_type:q} \
          "$(stat -c '%s' {input.graph:q})" "$segments" "$links" "$paths" "$segment_bp" \
          "$(wc -l < {output.bubbles:q})" >> {output.inventory:q}
        """


rule inventory_source_graphs:
    input:
        inventories=psc_inventory_jobs(),
        bubbles=expand(f"{PSC_OUT}/sv_catalogs/{{target}}.bubbles.tsv", target=PSC_GRAPHS)
    output:
        inventory=PSC_INVENTORY,
        catalog=PSC_RAW_CATALOG
    params:
        bubbles=" ".join(
            "--bubble " + shlex.quote(f"{target}={PSC_OUT}/sv_catalogs/{target}.bubbles.tsv")
            for target in PSC_GRAPHS
        )
    log:
        f"{PSC_OUT}/logs/inventory/merge.log"
    conda:
        PSC_ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.inventory:q})" "$(dirname {log:q})"
        head -n 1 {input.inventories[0]:q} > {output.inventory:q}.tmp
        for table in {input.inventories:q}; do tail -n +2 "$table" >> {output.inventory:q}.tmp; done
        mv {output.inventory:q}.tmp {output.inventory:q}
        python {PSC_STATS:q} catalog {params.bubbles} --output {output.catalog:q} > {log:q} 2>&1
        """


rule audit_source_graph_membership:
    input:
        samples=PSC_SAMPLES_TSV,
        hprc=PSC["membership_manifests"]["hprc"],
        local=PSC["membership_manifests"]["local_v1"],
        rescue=PSC["membership_manifests"]["local_v2_rescues"]
    output:
        PSC_MEMBERSHIP
    log:
        f"{PSC_OUT}/logs/sample_graph_membership.log"
    conda:
        PSC_ENV
    threads: 1
    resources:
        mem_mb=1000,
        runtime_min=30
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        work=$(mktemp -d -p "${{TMPDIR:-/tmp}}" "psc-membership.XXXXXX")
        trap 'rm -rf -- "$work"' EXIT

        awk -F '\t' '
          NR == 1 {{
            for (i = 1; i <= NF; i++) if ($i == "Sample ID" || $i == "sample_id") sample_col = i
            if (!sample_col) {{print "No Sample ID column in HPRC membership manifest" > "/dev/stderr"; exit 1}}
            next
          }}
          $sample_col != "" && toupper($sample_col) !~ /^(CHM13|GRCH38)$/ {{print $sample_col}}
        ' {input.hprc:q} | sort -u > "$work/hprc"

        awk -F '\t' '
          NR == 1 {{for (i = 1; i <= NF; i++) header[$i] = i; next}}
          $(header["role"]) == "assembly" {{
            sample = $(header["assembly_id"])
            sub(/\.hifi\.hifiasm\.bp\.hap[12]\.p_ctg(\.clean)?$/, "", sample)
            print sample
          }}
        ' {input.local:q} | sort -u > "$work/local_v1"

        awk -F '\t' '
          NR == 1 {{for (i = 1; i <= NF; i++) header[$i] = i; next}}
          $(header["status"]) == "PASS" && $(header["sample_id"]) != "" {{print $(header["sample_id"])}}
        ' {input.rescue:q} | sort -u > "$work/rescues"
        sort -u "$work/local_v1" "$work/rescues" > "$work/local_v2"
        sort -u "$work/hprc" "$work/local_v1" > "$work/hprc_plus_local"

        printf 'sample_id\tsample_group\ttarget\tis_graph_member\tmembership_status\n' > {output:q}.tmp
        tail -n +2 {input.samples:q} | while IFS=$'\t' read -r sample batch type sex r1 r2 group; do
          for target in {PSC_GRAPHS}; do
            if grep -Fxq -- "$sample" "$work/$target"; then member=true; else member=false; fi
            printf '%s\t%s\t%s\t%s\tknown\n' "$sample" "$group" "$target" "$member" >> {output:q}.tmp
          done
        done
        mv {output:q}.tmp {output:q}
        if awk -F '\t' 'NR > 1 && $4 == "true" {{print; found=1}} END {{exit !found}}' {output:q} > {log:q}; then
          echo "Pilot contains a source-graph member; see {output}" >&2
          exit 1
        fi
        """


rule prepare_tmem_targets:
    input:
        script=PSC_TMEM,
        samples=PSC_SAMPLES_TSV,
        vcf=PSC["tmem"]["vcf"],
        targets=PSC["tmem"]["target_insertions"],
        reference=PSC["tmem"]["regional_reference"],
        candidates=PSC["tmem"]["candidate_table"],
        kmers=PSC["tmem"]["kmer_validation"]
    output:
        targets=PSC_TMEM_TARGETS,
        fasta=PSC_TMEM_FASTA,
        validation=PSC_TMEM_VALIDATION
    params:
        offset=int(PSC["tmem"]["vcf_position_offset"]),
        genome_chrom=PSC["tmem"]["grch38_chrom"],
        genome_start=int(PSC["tmem"]["grch38_region_start"]),
        flank=int(PSC["tmem"].get("flank_bp", 1000)),
        high=",".join(PSC_GROUPS["high_insertion_probability"]),
        low=",".join(PSC_GROUPS["lower_insertion_probability"]),
        absent=",".join(PSC_GROUPS["high_insertion_absence_probability"]),
        untested=",".join(PSC_GROUPS["not_tested_on_tmem_graph"])
    log:
        f"{PSC_ANALYSIS}/logs/prepare_tmem_targets.log"
    conda:
        PSC_ENV
    threads: 1
    resources:
        mem_mb=4000,
        runtime_min=60
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.fasta:q})" "$(dirname {log:q})"
        raw_targets={output.targets:q}.raw.tmp
        trap 'rm -f -- "$raw_targets" {output.validation:q}.tmp' EXIT
        python {PSC_TMEM:q} queries \
          --vcf {input.vcf:q} --targets {input.targets:q} --reference {input.reference:q} \
          --vcf-position-offset {params.offset} --flank-bp {params.flank} \
          --output-fasta {output.fasta:q} --output-manifest "$raw_targets" \
          > {log:q} 2>&1
        awk -F '\t' -v OFS='\t' -v chrom={params.genome_chrom:q} -v start={params.genome_start} '
          NR == 1 {{
            for (i=1; i<=NF; i++) header[$i] = i
            print $0, "whole_genome_chrom", "whole_genome_pos_1based", "source_regional_coordinate"
            next
          }}
          {{
            pos = $(header["vcf_pos"])
            print $0, chrom, start + pos - 1, $(header["vcf_chrom"]) ":" pos
          }}
        ' "$raw_targets" > {output.targets:q}.tmp
        mv {output.targets:q}.tmp {output.targets:q}

        awk -F '\t' -v OFS='\t' -v high={params.high:q} -v low={params.low:q} \
          -v absent={params.absent:q} -v untested={params.untested:q} '
          BEGIN {{
            print "sample_id", "sample_group", "runtime_validation"
            n = split(high, a, ","); for (i=1; i<=n; i++) wanted[a[i]] = "high_insertion_probability"
            n = split(low, a, ","); for (i=1; i<=n; i++) wanted[a[i]] = "lower_insertion_probability"
            n = split(absent, a, ","); for (i=1; i<=n; i++) wanted[a[i]] = "high_insertion_absence_probability"
            n = split(untested, a, ","); for (i=1; i<=n; i++) wanted[a[i]] = "not_tested_on_tmem_graph"
          }}
          FNR == 1 {{
            delete header
            for (i=1; i<=NF; i++) header[$i] = i
            next
          }}
          FILENAME == ARGV[1] {{
            sample = $(header["sample"])
            if (!(sample in wanted)) next
            if (wanted[sample] == "high_insertion_probability" && $(header["candidate_tier"]) == "very_high") candidate_ok[sample] = 1
            if (wanted[sample] == "lower_insertion_probability" && $(header["candidate_tier"]) == "moderate_candidate") candidate_ok[sample] = 1
            next
          }}
          FILENAME == ARGV[2] {{
            sample = $(header["sample"])
            if (!(sample in wanted)) next
            if (wanted[sample] == "high_insertion_probability" && $(header["kmer_call"]) == "kmer_supported_alt_edges") kmer_ok[sample] = 1
            if (wanted[sample] == "high_insertion_absence_probability" && $(header["kmer_call"]) == "no_alt_edge_kmer_signal") kmer_ok[sample] = 1
            if (wanted[sample] == "not_tested_on_tmem_graph") unexpectedly_tested[sample] = 1
          }}
          END {{
            for (sample in wanted) {{
              group = wanted[sample]
              ok = (group == "high_insertion_probability" ? candidate_ok[sample] && kmer_ok[sample] :
                    group == "lower_insertion_probability" ? candidate_ok[sample] :
                    group == "high_insertion_absence_probability" ? kmer_ok[sample] :
                    !unexpectedly_tested[sample])
              print sample, group, (ok ? "PASS" : "FAIL")
              if (!ok) failed = 1
            }}
            if (failed) exit 1
          }}
        ' {input.candidates:q} {input.kmers:q} > {output.validation:q}.tmp
        mv {output.validation:q}.tmp {output.validation:q}
        rm -f -- "$raw_targets"
        trap - EXIT
        """


rule build_source_giraffe_index:
    input:
        graph=psc_source_graph,
        inventory=PSC_INVENTORY,
        membership=PSC_MEMBERSHIP
    output:
        gbz=f"{PSC_OUT}/indexes/{{target}}/source.giraffe.gbz",
        dist=f"{PSC_OUT}/indexes/{{target}}/source.dist",
        minimizer=f"{PSC_OUT}/indexes/{{target}}/source.shortread.withzip.min",
        zipcodes=f"{PSC_OUT}/indexes/{{target}}/source.shortread.zipcodes",
        backbone=f"{PSC_OUT}/indexes/{{target}}/backbone.paths"
    params:
        reference=PSC.get("common_backbone_sample", "CHM13")
    log:
        f"{PSC_OUT}/logs/indexes/{{target}}.log"
    benchmark:
        f"{PSC_OUT}/benchmarks/indexes/{{target}}.tsv"
    conda:
        PSC_ENV
    threads:
        int(PSC["resources"].get("index_threads", 32))
    resources:
        mem_mb=int(PSC["resources"].get("index_mem_mb", 256000)),
        runtime_min=int(PSC["resources"].get("index_runtime_min", 2880)),
        disk_mb=int(PSC["resources"].get("index_disk_mb", 500000))
    wildcard_constraints:
        target=PSC_GRAPH_RE
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.gbz:q})" "$(dirname {log:q})"
        work=$(mktemp -d -p "${{TMPDIR:-/tmp}}" "psc-index-{wildcards.target}.XXXXXX")
        trap 'rm -rf -- "$work"' EXIT
        mkdir -p "$work/tmp"
        case {input.graph:q} in
          *.gz) gzip -cd -- {input.graph:q} > "$work/source.gfa" ;;
          *) cp -- {input.graph:q} "$work/source.gfa" ;;
        esac
        vg autoindex --workflow sr-giraffe --prefix "$work/source" \
          --gfa "$work/source.gfa" --threads {threads} \
          --target-mem {resources.mem_mb}M --tmp-dir "$work/tmp" > {log:q} 2>&1
        for file in "$work/source.giraffe.gbz" "$work/source.dist" \
                    "$work/source.shortread.withzip.min" "$work/source.shortread.zipcodes"; do
          test -s "$file"
        done
        vg validate "$work/source.giraffe.gbz" >> {log:q} 2>&1
        # vg memory-maps the distance index. Keep the shared copy read-only so
        # concurrent mapping jobs cannot open it for writing.
        chmod a-w -- "$work/source.dist"
        vg paths --list -x "$work/source.giraffe.gbz" \
          | awk -v sample={params.reference:q} '$0 == sample || index($0, sample "#") == 1 || index($0, sample ".") == 1' \
          > "$work/backbone.paths"
        test -s "$work/backbone.paths" || {{ echo "No {params.reference} backbone paths in {input.graph}" >&2; exit 1; }}
        mv "$work/source.giraffe.gbz" {output.gbz:q}
        mv "$work/source.dist" {output.dist:q}
        mv "$work/source.shortread.withzip.min" {output.minimizer:q}
        mv "$work/source.shortread.zipcodes" {output.zipcodes:q}
        mv "$work/backbone.paths" {output.backbone:q}
        """


rule count_sample_fastq:
    input:
        manifest=psc_read_manifest,
        r1=lambda wc: psc_reads(wc, "r1"),
        r2=lambda wc: psc_reads(wc, "r2")
    output:
        f"{PSC_OUT}/fastq_stats/{{sample}}.tsv"
    log:
        f"{PSC_OUT}/logs/fastq_stats/{{sample}}.log"
    conda:
        PSC_ENV
    threads: 1
    resources:
        mem_mb=2000,
        runtime_min=120
    wildcard_constraints:
        sample=PSC_SAMPLE_RE
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        seqkit stats -T {input.r1:q} {input.r2:q} > {output:q} 2> {log:q}
        """


rule export_mapping_node_translation:
    input:
        gbz=f"{PSC_OUT}/indexes/{{target}}/source.giraffe.gbz"
    output:
        f"{PSC_ANALYSIS}/tmem/translations/{{target}}.tsv"
    log:
        f"{PSC_ANALYSIS}/logs/tmem/translations/{{target}}.log"
    conda:
        PSC_ENV
    threads: 1
    resources:
        mem_mb=16000,
        runtime_min=120,
        disk_mb=10000
    wildcard_constraints:
        target=PSC_GRAPH_RE
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        trap 'rm -f -- {output:q}.tmp' EXIT
        vg gbwt -Z --translation {output:q}.tmp {input.gbz:q} > {log:q} 2>&1
        test -s {output:q}.tmp
        mv {output:q}.tmp {output:q}
        trap - EXIT
        """


rule locate_tmem_in_source_graph:
    input:
        script=PSC_TMEM,
        graph=psc_source_graph,
        queries=PSC_TMEM_FASTA,
        manifest=PSC_TMEM_TARGETS,
        translation=f"{PSC_ANALYSIS}/tmem/translations/{{target}}.tsv"
    output:
        gaf=f"{PSC_ANALYSIS}/tmem/targets/{{target}}.gaf",
        table=f"{PSC_ANALYSIS}/tmem/targets/{{target}}.tsv"
    params:
        min_coverage=float(PSC["tmem"].get("target_min_query_coverage", 0.95)),
        min_identity=float(PSC["tmem"].get("target_min_identity", 0.90)),
        min_mapq=int(PSC["tmem"].get("target_min_mapq", 20)),
        min_flank=int(PSC["tmem"].get("target_min_flank_aligned_bp", 900))
    log:
        f"{PSC_ANALYSIS}/logs/tmem/targets/{{target}}.log"
    benchmark:
        f"{PSC_ANALYSIS}/benchmarks/tmem/targets/{{target}}.tsv"
    conda:
        PSC_ENV
    threads:
        int(PSC["resources"].get("target_threads", 8))
    resources:
        mem_mb=int(PSC["resources"].get("target_mem_mb", 32000)),
        runtime_min=int(PSC["resources"].get("target_runtime_min", 120)),
        disk_mb=10000
    wildcard_constraints:
        target=PSC_GRAPH_RE
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.gaf:q})" "$(dirname {log:q})"
        minigraph -c -xasm --vc --secondary=yes -t {threads} {input.graph:q} {input.queries:q} \
          > {output.gaf:q}.tmp 2> {log:q}
        mv {output.gaf:q}.tmp {output.gaf:q}
        python {PSC_TMEM:q} audit --graph {wildcards.target:q} --gaf {output.gaf:q} \
          --translation {input.translation:q} \
          --query-manifest {input.manifest:q} --min-coverage {params.min_coverage} \
          --min-identity {params.min_identity} --min-mapq {params.min_mapq} \
          --min-flank-aligned {params.min_flank} --output {output.table:q} >> {log:q} 2>&1
        """


rule map_sample_to_source_graph:
    input:
        manifest=psc_read_manifest,
        gbz=f"{PSC_OUT}/indexes/{{target}}/source.giraffe.gbz",
        dist=f"{PSC_OUT}/indexes/{{target}}/source.dist",
        minimizer=f"{PSC_OUT}/indexes/{{target}}/source.shortread.withzip.min",
        zipcodes=f"{PSC_OUT}/indexes/{{target}}/source.shortread.zipcodes",
        backbone=f"{PSC_OUT}/indexes/{{target}}/backbone.paths",
        r1=lambda wc: psc_reads(wc, "r1"),
        r2=lambda wc: psc_reads(wc, "r2")
    output:
        gaf=f"{PSC_OUT}/mapping/{{target}}/{{sample}}.gaf.gz",
        bam=temp(f"{PSC_OUT}/mapping/{{target}}/{{sample}}.backbone.bam"),
        bai=temp(f"{PSC_OUT}/mapping/{{target}}/{{sample}}.backbone.bam.bai")
    log:
        f"{PSC_OUT}/logs/mapping/{{target}}/{{sample}}.log"
    benchmark:
        f"{PSC_OUT}/benchmarks/mapping/{{target}}/{{sample}}.tsv"
    conda:
        PSC_ENV
    threads:
        int(PSC["resources"].get("mapping_threads", 24))
    resources:
        mem_mb=int(PSC["resources"].get("mapping_mem_mb", 96000)),
        runtime_min=int(PSC["resources"].get("mapping_runtime_min", 1440)),
        disk_mb=int(PSC["resources"].get("mapping_disk_mb", 150000)),
        gaf_writers=1
    wildcard_constraints:
        target=PSC_GRAPH_RE,
        sample=PSC_SAMPLE_RE
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.gaf:q})" "$(dirname {log:q})"
        gaf_tmp={output.gaf:q}.tmp
        bam_tmp={output.bam:q}.tmp
        bai_tmp={output.bai:q}.tmp
        trap 'rm -f -- "$gaf_tmp" "$bam_tmp" "$bai_tmp"' EXIT
        rm -f -- "$gaf_tmp" "$bam_tmp" "$bai_tmp"
        stream_fastq() {{
          for reads in "$@"; do
            case "$reads" in *.gz) gzip -cd -- "$reads" ;; *) cat -- "$reads" ;; esac
          done
        }}

        # Shared distance indexes must be read-only when several Giraffe jobs
        # use the same index on a cluster node.
        chmod a-w -- {input.dist:q}

        # Write GAF directly instead of first storing a whole-genome GAM in
        # node-local scratch. The same GAF is then surjected to the backbone.
        vg giraffe -p -t {threads} -b default -o GAF \
          -Z {input.gbz:q} -d {input.dist:q} \
          -m {input.minimizer:q} -z {input.zipcodes:q} \
          -f <(stream_fastq {input.r1:q}) -f <(stream_fastq {input.r2:q}) \
          2> {log:q} | bgzip -c > "$gaf_tmp"
        test -s "$gaf_tmp"
        bgzip -t "$gaf_tmp"

        bgzip -cd "$gaf_tmp" \
          | vg surject -x {input.gbz:q} -G -b -i -t {threads} \
              -F {input.backbone:q} - 2>> {log:q} \
          | samtools sort -@ 4 -o "$bam_tmp" -
        samtools index -@ 4 "$bam_tmp" "$bai_tmp"
        mv "$gaf_tmp" {output.gaf:q}
        mv "$bam_tmp" {output.bam:q}
        mv "$bai_tmp" {output.bai:q}
        trap - EXIT
        """


rule calculate_mapping_metrics:
    input:
        manifest=psc_read_manifest,
        gaf=f"{PSC_OUT}/mapping/{{target}}/{{sample}}.gaf.gz",
        bam=f"{PSC_OUT}/mapping/{{target}}/{{sample}}.backbone.bam",
        bai=f"{PSC_OUT}/mapping/{{target}}/{{sample}}.backbone.bam.bai",
        fastq=f"{PSC_OUT}/fastq_stats/{{sample}}.tsv",
        graph=f"{PSC_OUT}/indexes/{{target}}/source.giraffe.gbz",
        inventory=PSC_INVENTORY,
        benchmark=f"{PSC_OUT}/benchmarks/mapping/{{target}}/{{sample}}.tsv",
        index_benchmark=f"{PSC_OUT}/benchmarks/indexes/{{target}}.tsv",
        dist=f"{PSC_OUT}/indexes/{{target}}/source.dist",
        minimizer=f"{PSC_OUT}/indexes/{{target}}/source.shortread.withzip.min",
        zipcodes=f"{PSC_OUT}/indexes/{{target}}/source.shortread.zipcodes"
    output:
        f"{PSC_OUT}/statistics/mapping/{{target}}/{{sample}}.tsv"
    params:
        group=psc_group,
        batch=psc_batch,
        vg_version=str(PSC["vg_version"])
    log:
        f"{PSC_OUT}/logs/statistics/mapping/{{target}}/{{sample}}.log"
    conda:
        PSC_ENV
    threads: 1
    resources:
        mem_mb=4000,
        runtime_min=120
    wildcard_constraints:
        target=PSC_GRAPH_RE,
        sample=PSC_SAMPLE_RE
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        work=$(mktemp -d -p "${{TMPDIR:-/tmp}}" "psc-stats.XXXXXX")
        trap 'rm -rf -- "$work"' EXIT
        samtools coverage {input.bam:q} > "$work/coverage.tsv"
        samtools flagstat {input.bam:q} > "$work/flagstat.txt"
        python {PSC_STATS:q} mapping --gaf {input.gaf:q} --fastq-stats {input.fastq:q} \
          --coverage "$work/coverage.tsv" --flagstat "$work/flagstat.txt" \
          --mapping-graph {input.graph:q} \
          --index {input.dist:q} --index {input.minimizer:q} --index {input.zipcodes:q} \
          --inventory {input.inventory:q} --benchmark {input.benchmark:q} \
          --index-benchmark {input.index_benchmark:q} \
          --sample {wildcards.sample:q} --sample-group {params.group:q} --batch {params.batch:q} \
          --target {wildcards.target:q} --mapper-version {params.vg_version:q} \
          --output {output:q} > {log:q} 2>&1
        """


rule calculate_sv_evidence:
    input:
        script=PSC_STATS,
        manifest=psc_read_manifest,
        gaf=f"{PSC_OUT}/mapping/{{target}}/{{sample}}.gaf.gz",
        target=f"{PSC_ANALYSIS}/tmem/targets/{{target}}.tsv",
        inventory=PSC_INVENTORY
    output:
        f"{PSC_ANALYSIS}/statistics/sv_evidence/{{target}}/{{sample}}.tsv"
    params:
        group=psc_group,
        batch=psc_batch,
        vg_version=str(PSC["vg_version"]),
        min_mapq=int(PSC["evidence_thresholds"].get("min_mapq", 5)),
        min_junction=int(PSC["evidence_thresholds"].get("min_junction_reads_per_side", 3)),
        min_reference=int(PSC["evidence_thresholds"].get("min_reference_reads", 3))
    log:
        f"{PSC_ANALYSIS}/logs/statistics/sv_evidence/{{target}}/{{sample}}.log"
    conda:
        PSC_ENV
    threads: 1
    resources:
        mem_mb=4000,
        runtime_min=360
    wildcard_constraints:
        target=PSC_GRAPH_RE,
        sample=PSC_SAMPLE_RE
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {PSC_STATS:q} evidence --gaf {input.gaf:q} --target-table {input.target:q} \
          --inventory {input.inventory:q} --sample {wildcards.sample:q} \
          --sample-group {params.group:q} --batch {params.batch:q} --target {wildcards.target:q} \
          --mapper-version {params.vg_version:q} --min-mapq {params.min_mapq} \
          --min-junction-reads {params.min_junction} --min-reference-reads {params.min_reference} \
          --output {output:q} > {log:q} 2>&1
        """


rule summarize_19kb_insertion:
    input:
        psc_evidence_jobs()
    output:
        PSC_TMEM_EVIDENCE
    params:
        inputs=lambda wc, input: " ".join("--input " + shlex.quote(str(path)) for path in input)
    log:
        f"{PSC_ANALYSIS}/logs/statistics/merge_tmem_evidence.log"
    conda:
        PSC_ENV
    shell:
        r"""
        set -euo pipefail
        python {PSC_STATS:q} merge {params.inputs} --output {output:q} > {log:q} 2>&1
        """


rule compare_source_graphs:
    input:
        script=PSC_STATS,
        metrics=psc_mapping_jobs(),
        tmem=PSC_TMEM_EVIDENCE
    output:
        metrics=PSC_MAPPING_METRICS,
        evidence=PSC_SV_EVIDENCE,
        comparison=PSC_COMPARISON
    params:
        metrics=lambda wc, input: " ".join("--input " + shlex.quote(str(path)) for path in input.metrics),
        permutations=int(PSC.get("paired_permutations", 10000)),
        seed=int(PSC.get("statistics_seed", 20260902))
    log:
        f"{PSC_ANALYSIS}/logs/statistics/compare.log"
    conda:
        PSC_ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.metrics:q})" "$(dirname {log:q})"
        python {PSC_STATS:q} merge {params.metrics} --output {output.metrics:q} > {log:q} 2>&1
        cp -- {input.tmem:q} {output.evidence:q}
        python {PSC_STATS:q} compare --metrics {output.metrics:q} \
          --permutations {params.permutations} --seed {params.seed} \
          --output {output.comparison:q} >> {log:q} 2>&1
        """


rule corrected_sv_catalog:
    input:
        script=PSC_STATS,
        bubbles=expand(f"{PSC_OUT}/sv_catalogs/{{target}}.bubbles.tsv", target=PSC_GRAPHS)
    output:
        PSC_SV_CATALOG
    params:
        reference=PSC.get("common_backbone_sample", "CHM13"),
        bubbles=" ".join(
            "--bubble " + shlex.quote(f"{target}={PSC_OUT}/sv_catalogs/{target}.bubbles.tsv")
            for target in PSC_GRAPHS
        )
    log:
        f"{PSC_ANALYSIS}/logs/sv_catalog.log"
    conda:
        PSC_ENV
    threads: 1
    resources:
        mem_mb=4000,
        runtime_min=60
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {PSC_STATS:q} catalog {params.bubbles} --reference-sample {params.reference:q} \
          --output {output:q} > {log:q} 2>&1
        """


rule summarize_tmem_graph_targets:
    input:
        tables=expand(f"{PSC_ANALYSIS}/tmem/targets/{{target}}.tsv", target=PSC_GRAPHS)
    output:
        PSC_GRAPH_TARGETS
    params:
        tables=lambda wc, input: " ".join("--input " + shlex.quote(str(path)) for path in input.tables)
    log:
        f"{PSC_ANALYSIS}/logs/tmem/merge_targets.log"
    conda:
        PSC_ENV
    threads: 1
    resources:
        mem_mb=1000,
        runtime_min=30
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log:q})"
        python {PSC_TMEM:q} merge {params.tables} --output {output:q} > {log:q} 2>&1
        """


rule direct_sv_graph_comparison_report:
    input:
        script=PSC_STATS,
        samples=PSC_SAMPLES_TSV,
        inventory=PSC_INVENTORY,
        membership=PSC_MEMBERSHIP,
        catalog=PSC_SV_CATALOG,
        metrics=PSC_MAPPING_METRICS,
        evidence=PSC_SV_EVIDENCE,
        comparison=PSC_COMPARISON,
        targets=PSC_GRAPH_TARGETS,
        tmem_validation=PSC_TMEM_VALIDATION
    output:
        PSC_REPORT
    log:
        f"{PSC_ANALYSIS}/logs/report.log"
    conda:
        PSC_ENV
    threads: 1
    resources:
        mem_mb=4000,
        runtime_min=60
    shell:
        r"""
        set -euo pipefail
        python {PSC_STATS:q} report --samples {input.samples:q} --inventory {input.inventory:q} \
          --membership {input.membership:q} --catalog {input.catalog:q} \
          --metrics {input.metrics:q} --evidence {input.evidence:q} \
          --comparison {input.comparison:q} --targets {input.targets:q} \
          --validation {input.tmem_validation:q} --output {output:q} > {log:q} 2>&1
        """


include: "flat_reference_comparison.smk"
