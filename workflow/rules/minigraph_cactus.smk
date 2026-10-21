"""Explicit five-stage Minigraph-Cactus build; excluded from rule all."""

configfile: workflow.source_path("../../config/config.yaml")

import json
from pathlib import Path


MC = config["minigraph_cactus"]
MC_ROOT = Path(config["workflow_root"])
if not MC_ROOT.is_absolute():
    raise WorkflowError("workflow_root must be an absolute server path")
workdir: str(MC_ROOT)

MC_BUILD_ID = str(MC["build_id"])
if not MC_BUILD_ID or MC_BUILD_ID in {".", ".."} or "/" in MC_BUILD_ID or "\\" in MC_BUILD_ID:
    raise WorkflowError("minigraph_cactus.build_id must be one directory name")
if not MC["out_name"] or "/" in MC["out_name"] or "\\" in MC["out_name"]:
    raise WorkflowError("minigraph_cactus.out_name must be one filename prefix")
if MC["batch_system"] != "single_machine":
    raise WorkflowError("minigraph_cactus requires single_machine Toil")
for key in ("results_root", "jobstore_root", "included_list"):
    if not Path(MC[key]).is_absolute():
        raise WorkflowError(f"minigraph_cactus.{key} must be an absolute server path")

MC_OUT = str(Path(MC["results_root"]) / MC_BUILD_ID)
MC_METADATA = f"{MC_OUT}/metadata"
MC_STAGES = f"{MC_OUT}/stages"
MC_GRAPH = f"{MC_OUT}/graph"
MC_JOBSTORES = str(Path(MC["jobstore_root"]) / MC_BUILD_ID)
MC_CHM13 = str(Path(config["references"]["CHM13"]))
MC_GRCH38 = str(Path(config["references"]["hg38"]))
for path in (MC_CHM13, MC_GRCH38):
    if not Path(path).is_absolute():
        raise WorkflowError("Minigraph-Cactus references must be absolute server paths")

MC_SEQFILE = f"{MC_METADATA}/seqfile.txt"
MC_MANIFEST = f"{MC_METADATA}/input_manifest.tsv"
MC_CONTIG_MAP = f"{MC_METADATA}/contig_map.tsv.gz"
MC_PREFLIGHT = f"{MC_METADATA}/preflight.json"
MC_GFA = f"{MC_STAGES}/minigraph/{MC['out_name']}.sv.gfa.gz"
MC_PAF = f"{MC_STAGES}/graphmap/mappings.paf"
MC_GRAPH_FASTA = f"{MC_STAGES}/graphmap/graph.fa.gz"
MC_GRAPHMAP_SEQFILE = f"{MC_METADATA}/seqfile.graphmap.txt"
MC_AUGMENTED_SEQFILE = f"{MC_METADATA}/seqfile.with-minigraph.txt"
MC_SPLIT_DIR = f"{MC_STAGES}/split"
MC_CHROMFILE = f"{MC_SPLIT_DIR}/chromfile.txt"
MC_ALIGN_DIR = f"{MC_STAGES}/align"
MC_ENV = str(workflow.source_path("../envs/minigraph_cactus.yaml"))
MC_PREPARE = str(workflow.source_path("../scripts/prepare_minigraph_cactus_seqfile.py"))
MC_GUARD = str(workflow.source_path("../scripts/check_minigraph_cactus_inputs.py"))
MC_SETTINGS_JSON = json.dumps({"minigraph_cactus": MC,
                               "references": {"CHM13": MC_CHM13, "GRCh38": MC_GRCH38}},
                              sort_keys=True, separators=(",", ":"))

shell.executable("/bin/bash")


def mc_done(stage):
    return f"{MC_METADATA}/{stage}.done"


def mc_finished(stage):
    return f"{MC_METADATA}/{stage}.cactus_finished"


def mc_store(stage):
    return f"{MC_JOBSTORES}/{stage}"


rule prepare_minigraph_cactus:
    input:
        included=MC["included_list"], chm13=MC_CHM13, grch38=MC_GRCH38
    output:
        seqfile=MC_SEQFILE, manifest=MC_MANIFEST,
        contig_map=MC_CONTIG_MAP, preflight=MC_PREFLIGHT
    params:
        script=MC_PREPARE, log_dir=f"{MC_OUT}/logs"
    log:
        f"{MC_OUT}/logs/prepare_minigraph_cactus.log"
    conda: MC_ENV
    threads: 1
    resources: mem_mb=8192, disk_mb=10000, runtime_min=1440
    shell:
        r"""
        mkdir -p {params.log_dir:q}
        python {params.script:q} --included-list {input.included:q} \
          --chm13 {input.chm13:q} --grch38 {input.grch38:q} \
          --seqfile {output.seqfile:q} --input-manifest {output.manifest:q} \
          --contig-map {output.contig_map:q} --preflight-json {output.preflight:q} \
          > {log:q} 2>&1
        """


rule mc_minigraph:
    input:
        seqfile=MC_SEQFILE, manifest=MC_MANIFEST,
        contig_map=MC_CONTIG_MAP, preflight=MC_PREFLIGHT
    output: done=mc_done("minigraph")
    params:
        store=mc_store("minigraph"), finished=mc_finished("minigraph"),
        guard=MC_GUARD, settings=MC_SETTINGS_JSON, metadata=MC_METADATA, gfa=MC_GFA,
        mg_cores=MC["stages"]["minigraph"]["mg_cores"],
        mg_memory=MC["stages"]["minigraph"]["mg_memory"],
        max_cores=MC["stages"]["minigraph"]["max_cores"],
        max_memory=MC["stages"]["minigraph"]["max_memory"],
        benchmark_dir=f"{MC_OUT}/benchmarks"
    log: f"{MC_OUT}/logs/minigraph.log"
    benchmark: f"{MC_OUT}/benchmarks/minigraph.tsv"
    conda: MC_ENV
    threads: int(MC["stages"]["minigraph"]["threads"])
    resources:
        mem_mb=int(MC["stages"]["minigraph"]["mem_mb"]),
        disk_mb=int(MC["stages"]["minigraph"]["disk_mb"]),
        runtime_min=int(MC["stages"]["minigraph"]["runtime_min"])
    shell:
        r"""
        mkdir -p "$(dirname {params.store:q})" "$(dirname {params.gfa:q})" \
          "$(dirname {log:q})" {params.benchmark_dir:q}
        python {params.guard:q} --preflight-json {input.preflight:q} \
          --seqfile {input.seqfile:q} --settings-json {params.settings:q} \
          --metadata-dir {params.metadata:q} --initialize
        if [[ ! -e {params.finished:q} ]]; then
          restart=()
          if [[ -e {params.store:q} ]]; then
            restart=(--restart)
          elif [[ -e {params.gfa:q} ]]; then
            echo 'Minigraph output exists without its Toil store; inspect it or use a new build ID' >&2
            exit 1
          fi
          cactus-minigraph {params.store:q} {input.seqfile:q} {params.gfa:q} \
            --reference CHM13 GRCh38 --mgCores {params.mg_cores} \
            --mgMemory {params.mg_memory:q} --batchSystem single_machine \
            --clean=onSuccess --maxCores {params.max_cores} \
            --maxMemory {params.max_memory:q} "${{restart[@]}}" > {log:q} 2>&1
          touch {params.finished:q}
        fi
        test -s {params.gfa:q}
        touch {output.done:q}
        """


rule mc_graphmap:
    input:
        seqfile=MC_SEQFILE, preflight=MC_PREFLIGHT,
        previous=mc_done("minigraph")
    output: done=mc_done("graphmap")
    params:
        store=mc_store("graphmap"), finished=mc_finished("graphmap"),
        guard=MC_GUARD, settings=MC_SETTINGS_JSON, metadata=MC_METADATA,
        gfa=MC_GFA, seq_copy=MC_GRAPHMAP_SEQFILE,
        augmented=MC_AUGMENTED_SEQFILE, paf=MC_PAF, fasta=MC_GRAPH_FASTA,
        map_cores=MC["stages"]["graphmap"]["map_cores"],
        max_cores=MC["stages"]["graphmap"]["max_cores"],
        max_memory=MC["stages"]["graphmap"]["max_memory"],
        benchmark_dir=f"{MC_OUT}/benchmarks"
    log: f"{MC_OUT}/logs/graphmap.log"
    benchmark: f"{MC_OUT}/benchmarks/graphmap.tsv"
    conda: MC_ENV
    threads: int(MC["stages"]["graphmap"]["threads"])
    resources:
        mem_mb=int(MC["stages"]["graphmap"]["mem_mb"]),
        disk_mb=int(MC["stages"]["graphmap"]["disk_mb"]),
        runtime_min=int(MC["stages"]["graphmap"]["runtime_min"])
    shell:
        r"""
        python {params.guard:q} --preflight-json {input.preflight:q} \
          --seqfile {input.seqfile:q} --settings-json {params.settings:q} \
          --metadata-dir {params.metadata:q}
        test -s {params.gfa:q}
        mkdir -p "$(dirname {params.store:q})" "$(dirname {params.paf:q})" \
          "$(dirname {log:q})" {params.benchmark_dir:q}
        if [[ ! -e {params.finished:q} ]]; then
          restart=()
          if [[ -e {params.store:q} ]]; then
            test -s {params.seq_copy:q}
            restart=(--restart)
          else
            if [[ -e {params.seq_copy:q} || -e {params.paf:q} || -e {params.fasta:q} ]]; then
              echo 'Graphmap output exists without its Toil store; inspect it or use a new build ID' >&2
              exit 1
            fi
            cp {input.seqfile:q} {params.seq_copy:q}
          fi
          cactus-graphmap {params.store:q} {params.seq_copy:q} {params.gfa:q} {params.paf:q} \
            --outputFasta {params.fasta:q} --reference CHM13 GRCh38 \
            --mapCores {params.map_cores} --batchSystem single_machine \
            --clean=onSuccess --maxCores {params.max_cores} \
            --maxMemory {params.max_memory:q} "${{restart[@]}}" > {log:q} 2>&1
          touch {params.finished:q}
        fi
        test -s {params.paf:q}
        test -s {params.fasta:q}
        python {params.guard:q} --preflight-json {input.preflight:q} \
          --seqfile {input.seqfile:q} --settings-json {params.settings:q} \
          --metadata-dir {params.metadata:q} --graphmap-copy {params.seq_copy:q} \
          --graph-fasta {params.fasta:q}
        cp {params.seq_copy:q} {params.augmented:q}
        touch {output.done:q}
        """


rule mc_split:
    input:
        seqfile=MC_SEQFILE, preflight=MC_PREFLIGHT, previous=mc_done("graphmap")
    output: done=mc_done("split")
    params:
        store=mc_store("split"), finished=mc_finished("split"),
        guard=MC_GUARD, settings=MC_SETTINGS_JSON, metadata=MC_METADATA,
        seqfile=MC_AUGMENTED_SEQFILE,
        gfa=MC_GFA, paf=MC_PAF, outdir=MC_SPLIT_DIR,
        max_cores=MC["stages"]["split"]["max_cores"],
        max_memory=MC["stages"]["split"]["max_memory"],
        benchmark_dir=f"{MC_OUT}/benchmarks"
    log: f"{MC_OUT}/logs/split.log"
    benchmark: f"{MC_OUT}/benchmarks/split.tsv"
    conda: MC_ENV
    threads: int(MC["stages"]["split"]["threads"])
    resources:
        mem_mb=int(MC["stages"]["split"]["mem_mb"]),
        disk_mb=int(MC["stages"]["split"]["disk_mb"]),
        runtime_min=int(MC["stages"]["split"]["runtime_min"])
    shell:
        r"""
        python {params.guard:q} --preflight-json {input.preflight:q} \
          --seqfile {input.seqfile:q} --settings-json {params.settings:q} \
          --metadata-dir {params.metadata:q}
        test -s {params.seqfile:q} && test -s {params.gfa:q} && test -s {params.paf:q}
        mkdir -p "$(dirname {params.store:q})" {params.outdir:q} \
          "$(dirname {log:q})" {params.benchmark_dir:q}
        if [[ ! -e {params.finished:q} ]]; then
          restart=()
          if [[ -e {params.store:q} ]]; then
            restart=(--restart)
          elif [[ -n "$(find {params.outdir:q} -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
            echo 'Split output exists without its Toil store; inspect it or use a new build ID' >&2
            exit 1
          fi
          cactus-graphmap-split {params.store:q} {params.seqfile:q} {params.gfa:q} {params.paf:q} \
            --outDir {params.outdir:q} --reference CHM13 GRCh38 \
            --batchSystem single_machine --clean=onSuccess \
            --maxCores {params.max_cores} --maxMemory {params.max_memory:q} \
            "${{restart[@]}}" > {log:q} 2>&1
          touch {params.finished:q}
        fi
        test -s {params.outdir:q}/chromfile.txt
        test -s {params.outdir:q}/minigraph.split.log
        test -s {params.outdir:q}/contig_sizes.tsv
        awk 'NF != 3 || $1 == "." || $1 == ".." || $1 ~ /\// || seen[$1]++ {{bad=1}} END {{exit NR == 0 || bad}}' {params.outdir:q}/chromfile.txt
        while read -r name seqfile paf; do
          test -s "$seqfile" && test -s "$paf"
        done < {params.outdir:q}/chromfile.txt
        touch {output.done:q}
        """


rule mc_align:
    input:
        seqfile=MC_SEQFILE, preflight=MC_PREFLIGHT, previous=mc_done("split")
    output: done=mc_done("align")
    params:
        store=mc_store("align"), finished=mc_finished("align"),
        guard=MC_GUARD, settings=MC_SETTINGS_JSON, metadata=MC_METADATA,
        chromfile=MC_CHROMFILE, outdir=MC_ALIGN_DIR,
        cons_cores=MC["stages"]["align"]["cons_cores"],
        cons_memory=MC["stages"]["align"]["cons_memory"],
        max_cores=MC["stages"]["align"]["max_cores"],
        max_memory=MC["stages"]["align"]["max_memory"],
        benchmark_dir=f"{MC_OUT}/benchmarks"
    log: f"{MC_OUT}/logs/align.log"
    benchmark: f"{MC_OUT}/benchmarks/align.tsv"
    conda: MC_ENV
    threads: int(MC["stages"]["align"]["threads"])
    resources:
        mem_mb=int(MC["stages"]["align"]["mem_mb"]),
        disk_mb=int(MC["stages"]["align"]["disk_mb"]),
        runtime_min=int(MC["stages"]["align"]["runtime_min"])
    shell:
        r"""
        python {params.guard:q} --preflight-json {input.preflight:q} \
          --seqfile {input.seqfile:q} --settings-json {params.settings:q} \
          --metadata-dir {params.metadata:q}
        test -s {params.chromfile:q}
        mkdir -p "$(dirname {params.store:q})" {params.outdir:q} \
          "$(dirname {log:q})" {params.benchmark_dir:q}
        if [[ ! -e {params.finished:q} ]]; then
          restart=()
          if [[ -e {params.store:q} ]]; then
            restart=(--restart)
          elif [[ -n "$(find {params.outdir:q} -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
            echo 'Alignment output exists without its Toil store; inspect it or use a new build ID' >&2
            exit 1
          fi
          cactus-align {params.store:q} {params.chromfile:q} {params.outdir:q} \
            --batch --pangenome --outVG --reference CHM13 GRCh38 --maxLen 10000 \
            --consCores {params.cons_cores} --consMemory {params.cons_memory:q} \
            --batchSystem single_machine --clean=onSuccess \
            --maxCores {params.max_cores} --maxMemory {params.max_memory:q} \
            "${{restart[@]}}" > {log:q} 2>&1
          touch {params.finished:q}
        fi
        while read -r name seqfile paf; do
          test -s {params.outdir:q}/"$name".raw.vg
          test -s {params.outdir:q}/"$name".hal
        done < {params.chromfile:q}
        touch {output.done:q}
        """


rule mc_join:
    input:
        seqfile=MC_SEQFILE, preflight=MC_PREFLIGHT, previous=mc_done("align")
    output: done=mc_done("join")
    params:
        store=mc_store("join"), finished=mc_finished("join"),
        guard=MC_GUARD, settings=MC_SETTINGS_JSON, metadata=MC_METADATA,
        chromfile=MC_CHROMFILE, align_dir=MC_ALIGN_DIR,
        outdir=MC_GRAPH, out_name=MC["out_name"],
        index_cores=MC["stages"]["join"]["index_cores"],
        index_memory=MC["stages"]["join"]["index_memory"],
        max_cores=MC["stages"]["join"]["max_cores"],
        max_memory=MC["stages"]["join"]["max_memory"],
        benchmark_dir=f"{MC_OUT}/benchmarks"
    log: f"{MC_OUT}/logs/join.log"
    benchmark: f"{MC_OUT}/benchmarks/join.tsv"
    conda: MC_ENV
    threads: int(MC["stages"]["join"]["threads"])
    resources:
        mem_mb=int(MC["stages"]["join"]["mem_mb"]),
        disk_mb=int(MC["stages"]["join"]["disk_mb"]),
        runtime_min=int(MC["stages"]["join"]["runtime_min"])
    shell:
        r"""
        python {params.guard:q} --preflight-json {input.preflight:q} \
          --seqfile {input.seqfile:q} --settings-json {params.settings:q} \
          --metadata-dir {params.metadata:q}
        mkdir -p "$(dirname {params.store:q})" {params.outdir:q} \
          "$(dirname {log:q})" {params.benchmark_dir:q}
        vgs=()
        hals=()
        while read -r name seqfile paf; do
          test -s {params.align_dir:q}/"$name".raw.vg
          test -s {params.align_dir:q}/"$name".hal
          vgs+=({params.align_dir:q}/"$name".raw.vg)
          hals+=({params.align_dir:q}/"$name".hal)
        done < {params.chromfile:q}
        test "${{#vgs[@]}}" -gt 0
        if [[ ! -e {params.finished:q} ]]; then
          restart=()
          if [[ -e {params.store:q} ]]; then
            restart=(--restart)
          elif [[ -n "$(find {params.outdir:q} -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
            echo 'Join output exists without its Toil store; inspect it or use a new build ID' >&2
            exit 1
          fi
          cactus-graphmap-join {params.store:q} --vg "${{vgs[@]}}" --hal "${{hals[@]}}" \
            --outDir {params.outdir:q} --outName {params.out_name:q} \
            --reference CHM13 GRCh38 --vcfReference CHM13 \
            --gfa --gbz --vcf --giraffe --indexCores {params.index_cores} \
            --indexMemory {params.index_memory:q} --batchSystem single_machine \
            --clean=onSuccess --maxCores {params.max_cores} \
            --maxMemory {params.max_memory:q} "${{restart[@]}}" > {log:q} 2>&1
          touch {params.finished:q}
        fi
        for suffix in .gfa.gz .gbz .vcf.gz .raw.vcf.gz .d2.gbz .d2.dist \
                      .d2.shortread.withzip.min .d2.shortread.zipcodes .full.hal \
                      .stats/graph-stats.tsv; do
          test -s {params.outdir:q}/{params.out_name}"$suffix"
        done
        for vcf in .vcf.gz .raw.vcf.gz; do
          if [[ ! -s {params.outdir:q}/{params.out_name}"$vcf".tbi ]]; then
            test -s {params.outdir:q}/{params.out_name}"$vcf".csi
          fi
        done
        touch {output.done:q}
        """


rule minigraph_cactus:
    input:
        seqfile=MC_SEQFILE, preflight=MC_PREFLIGHT, done=mc_done("join")
    params:
        guard=MC_GUARD, settings=MC_SETTINGS_JSON, metadata=MC_METADATA,
        gfa=MC_GFA, paf=MC_PAF, fasta=MC_GRAPH_FASTA,
        augmented=MC_AUGMENTED_SEQFILE, chromfile=MC_CHROMFILE,
        align_dir=MC_ALIGN_DIR, outdir=MC_GRAPH, out_name=MC["out_name"]
    conda: MC_ENV
    threads: 1
    resources: mem_mb=1024, runtime_min=30
    shell:
        r"""
        python {params.guard:q} --preflight-json {input.preflight:q} \
          --seqfile {input.seqfile:q} --settings-json {params.settings:q} \
          --metadata-dir {params.metadata:q}
        test -s {params.gfa:q} && test -s {params.paf:q} && test -s {params.fasta:q}
        test -s {params.augmented:q} && test -s {params.chromfile:q}
        while read -r name seqfile paf; do
          test -s "$seqfile" && test -s "$paf"
          test -s {params.align_dir:q}/"$name".raw.vg
          test -s {params.align_dir:q}/"$name".hal
        done < {params.chromfile:q}
        for suffix in .gfa.gz .gbz .vcf.gz .raw.vcf.gz .d2.gbz .d2.dist \
                      .d2.shortread.withzip.min .d2.shortread.zipcodes .full.hal \
                      .stats/graph-stats.tsv; do
          test -s {params.outdir:q}/{params.out_name}"$suffix"
        done
        for vcf in .vcf.gz .raw.vcf.gz; do
          if [[ ! -s {params.outdir:q}/{params.out_name}"$vcf".tbi ]]; then
            test -s {params.outdir:q}/{params.out_name}"$vcf".csi
          fi
        done
        """
