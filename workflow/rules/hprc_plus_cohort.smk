"""Add the selected local assemblies to the HPRC Minigraph SV graphs."""

configfile: workflow.source_path("../../config/config.yaml")


import csv
from pathlib import Path


HPRC_COHORT_PROJECT_ROOT = Path(workflow.current_basedir).parents[1]
workdir: str(HPRC_COHORT_PROJECT_ROOT)

HPRC_COHORT = config["hprc_plus_cohort"]
HPRC_COHORT_OUTDIR = HPRC_COHORT["results"]
HPRC_COHORT_GRAPHS = HPRC_COHORT["graphs"]
HPRC_COHORT_ASSEMBLIES = HPRC_COHORT["assemblies"]
HPRC_COHORT_PANGENOME = config.get("pangenome", {})
HPRC_COHORT_ENV = str(workflow.source_path("../envs/pangenome.yaml"))
HPRC_COHORT_SAMPLE_MANIFEST = (
    f"{HPRC_COHORT_OUTDIR}/metadata/assembly_samples.tsv"
)


def hprc_cohort_sample_graph(reference, sample):
    return (
        f"{HPRC_COHORT_OUTDIR}/graphs/checkpoints/"
        f"hprc-v2.0-mc-{reference}.sample-{sample}.cohort.sv.gfa"
    )


def hprc_cohort_sample_manifest(wildcards):
    return checkpoints.prepare_hprc_cohort_samples.get().output.manifest


def hprc_cohort_samples():
    manifest = checkpoints.prepare_hprc_cohort_samples.get().output.manifest
    manifest_path = Path(str(manifest))
    if not manifest_path.is_absolute():
        manifest_path = HPRC_COHORT_PROJECT_ROOT / manifest_path

    with manifest_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    samples = list(dict.fromkeys(row["sample"] for row in rows))
    if not samples:
        raise WorkflowError(f"No cohort samples found in {manifest}")
    return samples


def hprc_cohort_sample_backbone(wildcards):
    samples = hprc_cohort_samples()
    try:
        index = samples.index(wildcards.sample)
    except ValueError:
        raise WorkflowError(
            f"Sample {wildcards.sample} is not present in {HPRC_COHORT_SAMPLE_MANIFEST}"
        )
    if index == 0:
        return HPRC_COHORT_GRAPHS[wildcards.reference]
    return hprc_cohort_sample_graph(wildcards.reference, samples[index - 1])


def hprc_cohort_last_sample_graph(wildcards):
    return hprc_cohort_sample_graph(wildcards.reference, hprc_cohort_samples()[-1])


rule hprc_plus_cohort:
    input:
        expand(
            f"{HPRC_COHORT_OUTDIR}/graphs/hprc-v2.0-mc-{{reference}}.cohort.sv.gfa",
            reference=HPRC_COHORT_GRAPHS,
        )


checkpoint prepare_hprc_cohort_samples:
    input:
        assemblies=HPRC_COHORT_ASSEMBLIES
    output:
        manifest=HPRC_COHORT_SAMPLE_MANIFEST
    log:
        f"{HPRC_COHORT_OUTDIR}/logs/prepare_assembly_samples.log"
    threads: 1
    resources:
        mem_mb=1000,
        runtime_min=10
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {output.manifest:q})" "$(dirname {log:q})"
        temporary={output.manifest:q}.tmp
        trap 'rm -f -- "$temporary"' EXIT

        awk -F '\t' '
            BEGIN {{
                OFS = "\t"
                print "sample", "assembly_id", "path"
            }}
            NR > 1 && $2 == "assembly" {{
                sample = $3
                if (!sub(/\.hifi\.hifiasm\.bp\.hap[12]\.p_ctg(\.clean)?$/, "", sample)) {{
                    print "Cannot derive sample ID from assembly: " $3 > "/dev/stderr"
                    invalid = 1
                    next
                }}
                if (!(sample in seen)) {{
                    count++
                    seen[sample] = 1
                }}
                print sample, $3, $4
            }}
            END {{
                if (invalid || count == 0) exit 1
            }}
        ' {input.assemblies:q} > "$temporary" 2> {log:q}

        test "$(wc -l < "$temporary")" -gt 1
        mv "$temporary" {output.manifest:q}
        trap - EXIT
        """


rule build_hprc_plus_cohort_sv_graph_sample:
    input:
        backbone=hprc_cohort_sample_backbone,
        manifest=hprc_cohort_sample_manifest
    output:
        gfa=(
            f"{HPRC_COHORT_OUTDIR}/graphs/checkpoints/"
            "hprc-v2.0-mc-{reference}.sample-{sample}.cohort.sv.gfa"
        )
    params:
        preset=HPRC_COHORT.get(
            "minigraph_preset",
            HPRC_COHORT_PANGENOME.get("minigraph_preset", "ggs"),
        ),
        benchmark_dir=f"{HPRC_COHORT_OUTDIR}/benchmarks/samples"
    log:
        f"{HPRC_COHORT_OUTDIR}/logs/{{reference}}.sample-{{sample}}.minigraph.log"
    benchmark:
        f"{HPRC_COHORT_OUTDIR}/benchmarks/samples/{{reference}}.sample-{{sample}}.minigraph.tsv"
    conda:
        HPRC_COHORT_ENV
    threads:
        int(HPRC_COHORT.get("threads", HPRC_COHORT_PANGENOME.get("threads", 48)))
    resources:
        mem_mb=int(
            HPRC_COHORT.get("mem_mb", HPRC_COHORT_PANGENOME.get("mem_mb", 120000))
        ),
        disk_mb=int(
            HPRC_COHORT.get("disk_mb", HPRC_COHORT_PANGENOME.get("disk_mb", 20000))
        ),
        runtime_min=int(
            HPRC_COHORT.get(
                "runtime_min", HPRC_COHORT_PANGENOME.get("runtime_min", 2880)
            )
        )
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {output.gfa:q})" \
                 "$(dirname {log:q})" \
                 {params.benchmark_dir:q}

        mapfile -t assemblies < <(
            awk -F '\t' -v sample={wildcards.sample:q} \
                'NR > 1 && $1 == sample {{print $3}}' {input.manifest:q}
        )
        if [ "${{#assemblies[@]}}" -eq 0 ]; then
            echo "No cohort assemblies found for sample {wildcards.sample}" >&2
            exit 1
        fi
        for assembly in "${{assemblies[@]}}"; do
            if [ ! -s "$assembly" ]; then
                echo "Assembly is missing or empty: $assembly" >&2
                exit 1
            fi
        done

        temporary={output.gfa:q}.tmp
        trap 'rm -f -- "$temporary"' EXIT
        minigraph -cx{params.preset} -t {threads} \
          {input.backbone:q} "${{assemblies[@]}}" \
          > "$temporary" 2> {log:q}
        test -s "$temporary"
        mv "$temporary" {output.gfa:q}
        trap - EXIT
        """


rule finalize_hprc_plus_cohort_sv_graph:
    input:
        gfa=hprc_cohort_last_sample_graph
    output:
        gfa=f"{HPRC_COHORT_OUTDIR}/graphs/hprc-v2.0-mc-{{reference}}.cohort.sv.gfa"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {output.gfa:q})"
        temporary={output.gfa:q}.tmp
        trap 'rm -f -- "$temporary"' EXIT
        cp --reflink=auto {input.gfa:q} "$temporary"
        test -s "$temporary"
        mv "$temporary" {output.gfa:q}
        trap - EXIT
        """


wildcard_constraints:
    reference="grch38|chm13",
    sample="[A-Za-z0-9][A-Za-z0-9_.-]*"
