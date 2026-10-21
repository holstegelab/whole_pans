"""Reuse frozen SV GFAs for two staged Minigraph-Cactus builds."""

configfile: workflow.source_path("../../config/config.yaml")

import json
from pathlib import Path
from shlex import quote


MR = config["minigraph_cactus_reuse"]
MR_BASE = config["minigraph_cactus"]
workdir: str(config["workflow_root"])

MR_RESULTS = str(Path(MR["results_root"]))
MR_STORES = str(Path(MR["jobstore_root"]))
MR_IDS = {entry["build_id"]: entry for entry in MR["builds"].values()}
MR_LOCAL = MR["builds"]["local_v1"]["build_id"]
MR_HPRC = MR["builds"]["hprc_plus_local_v1"]["build_id"]
MR_ENV = str(workflow.source_path("../envs/minigraph_cactus.yaml"))
MR_PREPARE = str(workflow.source_path("../scripts/prepare_minigraph_cactus_reuse.py"))
MR_NORMALIZE = str(workflow.source_path("../scripts/normalize_minigraph_cactus_seed.py"))
MR_COMPARE = str(workflow.source_path("../scripts/compare_minigraph_cactus_reuse.py"))
MR_COMMON = str(workflow.source_path("../scripts/minigraph_cactus_common.py"))
MR_SEQFILE = str(workflow.source_path("../scripts/prepare_minigraph_cactus_seqfile.py"))
MR_ROOT = f"{MR_RESULTS}/{{build_id}}"
MR_META = f"{MR_ROOT}/metadata"
MR_STAGE = f"{MR_ROOT}/stages"
MR_SEED_CHUNKS = int(MR["seed"]["chunks"])
if not 2 <= MR_SEED_CHUNKS <= 99:
    raise WorkflowError("reuse seed chunks must be between 2 and 99")
MR_SEED_IDS = [f"{i:02d}" for i in range(MR_SEED_CHUNKS)]

shell.executable("/bin/bash")

if len(MR_IDS) != 2 or MR["batch_system"] != "single_machine":
    raise WorkflowError("reuse needs two different build IDs and single_machine Toil")


def mr_build(wildcards):
    return MR_IDS[wildcards.build_id]


def mr_profile(wildcards, stage):
    base = MR["seed"] if stage == "seed" else MR_BASE["stages"][stage]
    return {**base, **mr_build(wildcards).get("stage_overrides", {}).get(stage, {})}


def mr_store(wildcards, stage):
    return f"{MR_STORES}/{wildcards.build_id}/{stage}"


def mr_settings(build_id):
    return json.dumps({"results_root": MR_RESULTS, "build_id": build_id,
                       "out_name": MR_IDS[build_id]["out_name"]})


def mr_graphmap_seqfile(wildcards):
    return str(checkpoints.prepare_minigraph_cactus_reuse.get(
        build_id=wildcards.build_id).output.seqfile)


def mr_graphmap_rows(wildcards):
    with open(mr_graphmap_seqfile(wildcards), encoding="utf-8") as handle:
        rows = [line.split() for line in handle
                if line.strip() and not line.startswith("#")]
    if not rows or any(len(row) != 2 for row in rows):
        raise WorkflowError("graphmap seqfile must have one sample and FASTA path per line")
    return rows


def mr_graphmap_fasta(wildcards):
    paths = dict(mr_graphmap_rows(wildcards))
    if wildcards.sample not in paths:
        raise WorkflowError(f"sample absent from graphmap seqfile: {wildcards.sample}")
    return paths[wildcards.sample]


def mr_graphmap_pafs(wildcards):
    return expand(f"{MR_STAGE}/graphmap/assemblies/{{sample}}.paf",
                  build_id=wildcards.build_id,
                  sample=[row[0] for row in mr_graphmap_rows(wildcards)])


wildcard_constraints:
    build_id="[^/]+",
    chunk="[0-9][0-9]"


rule minigraph_cactus_reuse_local_v1:
    input: f"{MR_RESULTS}/{MR_LOCAL}/metadata/join.done"


checkpoint prepare_minigraph_cactus_reuse:
    input:
        source=lambda wc: mr_build(wc)["source_gfa"],
        local=MR["local_manifest"],
        script=MR_PREPARE,
        common=MR_COMMON,
        seqfile_script=MR_SEQFILE,
        extra=lambda wc: ([mr_build(wc)["hprc_manifest"]] if "hprc_manifest" in mr_build(wc)
                          else [mr_build(wc)["chm13"], mr_build(wc)["grch38"]])
    output:
        seqfile=f"{MR_META}/seqfile.txt",
        manifest=f"{MR_META}/input_manifest.tsv",
        contig_map=f"{MR_META}/contig_map.tsv.gz",
        preflight=f"{MR_META}/preflight.json"
    params:
        extra=lambda wc: ("--hprc-manifest " + quote(mr_build(wc)["hprc_manifest"]) +
                          " --hprc-downloads " + quote(mr_build(wc)["hprc_downloads"])
                          if "hprc_manifest" in mr_build(wc) else
                          "--chm13 " + quote(mr_build(wc)["chm13"]) +
                          " --grch38 " + quote(mr_build(wc)["grch38"]))
    log: f"{MR_ROOT}/logs/preflight.log"
    conda: MR_ENV
    threads: 1
    resources: mem_mb=16000, disk_mb=20000, runtime_min=2880
    shell:
        "mkdir -p $(dirname {log:q}) && "
        "python {input.script:q} --seed-gfa {input.source:q} "
        "--local-manifest {input.local:q} {params.extra} "
        "--seqfile {output.seqfile:q} --input-manifest {output.manifest:q} "
        "--contig-map {output.contig_map:q} --preflight-json {output.preflight:q} "
        "> {log:q} 2>&1"


rule mc_reuse_lock_preflight:
    input: f"{MR_META}/preflight.json"
    output: f"{MR_META}/preflight.locked.json"
    threads: 1
    resources: mem_mb=1000, runtime_min=10
    shell:
        r"""
        if [[ -e {output:q} ]]; then
          cmp -s {input:q} {output:q} || \
            (echo 'Preflight changed; choose a new build ID' >&2; exit 1)
        else
          cp {input:q} {output:q}
        fi
        """


rule mc_reuse_seed_chunk:
    input:
        source=lambda wc: mr_build(wc)["source_gfa"],
        preflight=f"{MR_META}/preflight.json",
        locked=f"{MR_META}/preflight.locked.json",
        contig_map=f"{MR_META}/contig_map.tsv.gz",
        script=MR_NORMALIZE,
        common=MR_COMMON
    output:
        gfa=f"{MR_STAGE}/minigraph/chunks/{{chunk}}.gfa",
        report=f"{MR_STAGE}/minigraph/chunks/{{chunk}}.json"
    log: f"{MR_ROOT}/logs/seed-chunk-{{chunk}}.log"
    conda: MR_ENV
    threads: 1
    resources: mem_mb=4000, disk_mb=5000, runtime_min=720
    shell:
        r"""
        mkdir -p "$(dirname {log:q})"
        cmp -s {input.preflight:q} {input.locked:q} || \
          (echo 'Preflight changed; choose a new build ID' >&2; exit 1)
        python {input.script:q} chunk --source-gfa {input.source:q} \
          --output-gfa {output.gfa:q} --contig-map {input.contig_map:q} \
          --preflight-json {input.preflight:q} --report {output.report:q} \
          --chunk-index {wildcards.chunk} --chunks {MR_SEED_CHUNKS} \
          > {log:q} 2>&1
        """


rule mc_reuse_seed:
    input:
        chunks=lambda wc: expand(f"{MR_STAGE}/minigraph/chunks/{{chunk}}.gfa",
                                 build_id=wc.build_id, chunk=MR_SEED_IDS),
        reports=lambda wc: expand(f"{MR_STAGE}/minigraph/chunks/{{chunk}}.json",
                                  build_id=wc.build_id, chunk=MR_SEED_IDS),
        source=lambda wc: mr_build(wc)["source_gfa"],
        preflight=f"{MR_META}/preflight.json",
        locked=f"{MR_META}/preflight.locked.json",
        script=MR_NORMALIZE,
        common=MR_COMMON
    output:
        gfa=f"{MR_STAGE}/minigraph/seed.sv.gfa.gz",
        report=f"{MR_STAGE}/minigraph/seed-normalization.json"
    log: f"{MR_ROOT}/logs/seed.log"
    conda: MR_ENV
    threads: lambda wc: mr_profile(wc, "seed")["threads"]
    resources:
        mem_mb=lambda wc: mr_profile(wc, "seed")["mem_mb"],
        disk_mb=lambda wc: mr_profile(wc, "seed")["disk_mb"],
        runtime_min=lambda wc: mr_profile(wc, "seed")["runtime_min"]
    shell:
        r"""
        mkdir -p "$(dirname {log:q})"
        cmp -s {input.preflight:q} {input.locked:q} || \
          (echo 'Preflight changed; choose a new build ID' >&2; exit 1)
        cat {input.chunks:q} | pigz -1 -p {threads} -c > {output.gfa:q}.partial 2> {log:q}
        mv {output.gfa:q}.partial {output.gfa:q}
        python {input.script:q} finalize --source-gfa {input.source:q} \
          --output-gfa {output.gfa:q} --preflight-json {input.preflight:q} \
          --chunk-reports {input.reports:q} --report {output.report:q} \
          >> {log:q} 2>&1
        """


rule mc_reuse_graphmap_prepare:
    input:
        gfa=f"{MR_STAGE}/minigraph/seed.sv.gfa.gz",
        seqfile=f"{MR_META}/seqfile.txt",
        preflight=f"{MR_META}/preflight.json",
        locked=f"{MR_META}/preflight.locked.json"
    output:
        gfa=f"{MR_STAGE}/graphmap/cactus.gfa",
        fasta=f"{MR_STAGE}/graphmap/graph.fa.gz",
        augmented=f"{MR_META}/seqfile.with-minigraph.txt"
    log: f"{MR_ROOT}/logs/graphmap-prepare.log"
    conda: MR_ENV
    threads: 8
    resources: mem_mb=32000, disk_mb=30000, runtime_min=720
    shell:
        r"""
        mkdir -p "$(dirname {output.gfa:q})" "$(dirname {log:q})"
        exec > {log:q} 2>&1
        cmp -s {input.preflight:q} {input.locked:q} || \
          (echo 'Preflight changed; choose a new build ID' >&2; exit 1)
        gzip -dc {input.gfa:q} | awk -F '\t' '
          BEGIN {{OFS=FS}}
          FNR==NR {{
            if (split($0, row, /[[:space:]]+/) != 2) {{
              print "Malformed graphmap seqfile row " FNR > "/dev/stderr"; exit 1
            }}
            names[row[1]]=1; next
          }}
          /^S\t/ {{
            for (i=4; i<=NF; i++) if (index($i,"SN:Z:")==1) {{
              source=substr($i,6); first=index(source,"#")
              if (!first) break
              rest=substr(source,first+1); second=index(rest,"#")
              if (!second) exit 1
              sample=substr(source,1,first-1); hap=substr(rest,1,second-1)
              event=(sample in names ? sample : sample "." hap)
              if (!(event in names)) {{print "Unknown GFA source: " source > "/dev/stderr"; exit 1}}
              $i="SN:Z:id=" event "|" substr(rest,second+1)
              break
            }}
          }}
          {{print}}
        ' {input.seqfile:q} - > {output.gfa:q}.partial
        mv {output.gfa:q}.partial {output.gfa:q}
        gfatools gfa2fa {output.gfa:q} | \
          sed -e 's/^>\(.\)/>id=_MINIGRAPH_|\1/g' | \
          bgzip --threads {threads} -c > {output.fasta:q}.partial
        mv {output.fasta:q}.partial {output.fasta:q}
        cp {input.seqfile:q} {output.augmented:q}
        printf '_MINIGRAPH_ %s\n' {output.fasta:q} >> {output.augmented:q}
        test -s {output.gfa:q} && test -s {output.fasta:q}
        """


rule mc_reuse_graphmap_one:
    input:
        gfa=f"{MR_STAGE}/graphmap/cactus.gfa",
        fasta=mr_graphmap_fasta,
        seqfile=f"{MR_META}/seqfile.txt"
    output:
        paf=f"{MR_STAGE}/graphmap/assemblies/{{sample}}.paf",
        gaf=f"{MR_STAGE}/graphmap/assemblies/{{sample}}.gaf.gz"
    params: benchmark_dir=f"{MR_ROOT}/benchmarks"
    log: f"{MR_ROOT}/logs/graphmap-{{sample}}.log"
    benchmark: f"{MR_ROOT}/benchmarks/graphmap-{{sample}}.tsv"
    conda: MR_ENV
    threads: lambda wc: mr_profile(wc, "graphmap")["map_cores"]
    resources: mem_mb=96000, disk_mb=100000, runtime_min=1440
    shell:
        r"""
        mkdir -p "$(dirname {output.paf:q})" "$(dirname {log:q})" \
          {params.benchmark_dir:q}
        exec > {log:q} 2>&1
        temp_dir=$(mktemp -d "${{TMPDIR:-/tmp}}/mc-map.XXXXXX")
        trap 'rm -rf "$temp_dir"' EXIT
        if [[ {input.fasta:q} == *.gz ]]; then
          gzip -dc {input.fasta:q}
        else
          cat {input.fasta:q}
        fi | cactus_sanitizeFastaHeaders - {wildcards.sample:q} -p > "$temp_dir/query.fa"
        minigraph {input.gfa:q} "$temp_dir/query.fa" -c -xasm -t {threads} \
          -o "$temp_dir/mapping.gaf"
        gaf2unstable "$temp_dir/mapping.gaf" -g {input.gfa:q} \
          -o "$temp_dir/node_lengths.tsv" | \
          gaffilter - -r 5 -m 0.25 -q 5 -b 250000 -o 0 -i 0.5 \
          > "$temp_dir/unstable.gaf"
        gaf2paf "$temp_dir/unstable.gaf" -l "$temp_dir/node_lengths.tsv" | \
          awk 'BEGIN {{OFS="\t"}} {{$6="id=_MINIGRAPH_|"$6; print}}' > {output.paf:q}
        bgzip --threads {threads} -c "$temp_dir/mapping.gaf" > {output.gaf:q}
        test -s {output.gaf:q}
        """


rule mc_reuse_graphmap_merge:
    input:
        pafs=mr_graphmap_pafs,
        gfa=f"{MR_STAGE}/graphmap/cactus.gfa"
    output:
        mapped=f"{MR_STAGE}/graphmap/mapped.paf",
        reference=f"{MR_STAGE}/graphmap/reference.paf",
        merged=f"{MR_STAGE}/graphmap/merged.paf"
    log: f"{MR_ROOT}/logs/graphmap-merge.log"
    conda: MR_ENV
    threads: 1
    resources: mem_mb=32000, disk_mb=100000, runtime_min=720
    shell:
        r"""
        mkdir -p "$(dirname {output.merged:q})" "$(dirname {log:q})"
        exec > {log:q} 2>&1
        cat {input.pafs:q} > {output.mapped:q}
        rgfa2paf {input.gfa:q} -T 'id=_MINIGRAPH_|' -P 'id=CHM13|' \
          -i {output.mapped:q} > {output.reference:q}
        cat {output.mapped:q} {output.reference:q} > {output.merged:q}
        test -s {output.merged:q}
        """


rule mc_reuse_graphmap_vg:
    input: f"{MR_STAGE}/graphmap/cactus.gfa"
    output:
        vg=f"{MR_STAGE}/graphmap/cactus.vg",
        translation=f"{MR_STAGE}/graphmap/cactus.trans"
    log: f"{MR_ROOT}/logs/graphmap-vg.log"
    conda: MR_ENV
    threads: 1
    resources: mem_mb=288000, disk_mb=1000000, runtime_min=2880
    shell:
        r"""
        mkdir -p "$(dirname {output.vg:q})" "$(dirname {log:q})"
        vg convert -r 0 -g {input:q} -p -T {output.translation:q} \
          > {output.vg:q} 2> {log:q}
        test -s {output.vg:q} && test -s {output.translation:q}
        """


rule mc_reuse_graphmap:
    input:
        paf=f"{MR_STAGE}/graphmap/merged.paf",
        vg=f"{MR_STAGE}/graphmap/cactus.vg",
        translation=f"{MR_STAGE}/graphmap/cactus.trans"
    output:
        paf=f"{MR_STAGE}/graphmap/mappings.paf",
        done=f"{MR_META}/graphmap.done"
    log: f"{MR_ROOT}/logs/graphmap-deletion-filter.log"
    conda: MR_ENV
    threads: lambda wc: mr_profile(wc, "graphmap")["map_cores"]
    resources: mem_mb=288000, disk_mb=1000000, runtime_min=2880
    shell:
        r"""
        mkdir -p "$(dirname {output.done:q})" "$(dirname {log:q})"
        filter-paf-deletions {input.vg:q} {input.translation:q} {input.paf:q} \
          -d 10000000 -v -p -t {threads} -m 0.01 -s 2 \
          > {output.paf:q} 2> {log:q}
        test -s {output.paf:q}
        touch {output.done:q}
        """


rule mc_reuse_split:
    input:
        done=f"{MR_META}/graphmap.done",
        seqfile=f"{MR_META}/seqfile.with-minigraph.txt",
        gfa=f"{MR_STAGE}/minigraph/seed.sv.gfa.gz",
        paf=f"{MR_STAGE}/graphmap/mappings.paf",
        fasta=f"{MR_STAGE}/graphmap/graph.fa.gz"
    output: f"{MR_META}/split.done"
    params:
        store=lambda wc: mr_store(wc, "split"),
        seqfile=f"{MR_META}/seqfile.with-minigraph.txt",
        gfa=f"{MR_STAGE}/minigraph/seed.sv.gfa.gz",
        paf=f"{MR_STAGE}/graphmap/mappings.paf",
        outdir=f"{MR_STAGE}/split",
        max_cores=lambda wc: mr_profile(wc, "split")["max_cores"],
        max_memory=lambda wc: mr_profile(wc, "split")["max_memory"]
    log: f"{MR_ROOT}/logs/split.log"
    benchmark: f"{MR_ROOT}/benchmarks/split.tsv"
    conda: MR_ENV
    threads: lambda wc: mr_profile(wc, "split")["threads"]
    resources:
        mem_mb=lambda wc: mr_profile(wc, "split")["mem_mb"],
        disk_mb=lambda wc: mr_profile(wc, "split")["disk_mb"],
        runtime_min=lambda wc: mr_profile(wc, "split")["runtime_min"]
    shell:
        r"""
        mkdir -p "$(dirname {params.store:q})" {params.outdir:q} "$(dirname {log:q})"
        restart=()
        if [[ -e {params.store:q} ]]; then restart=(--restart); fi
        cactus-graphmap-split {params.store:q} {params.seqfile:q} {params.gfa:q} {params.paf:q} \
          --outDir {params.outdir:q} --reference CHM13 GRCh38 \
          --batchSystem single_machine --clean=onSuccess \
          --maxCores {params.max_cores} --maxMemory {params.max_memory:q} \
          "${{restart[@]}}" > {log:q} 2>&1
        test -s {params.outdir:q}/chromfile.txt
        test -s {params.outdir:q}/contig_sizes.tsv
        awk 'NF != 3 || $1 == "." || $1 == ".." || $1 ~ /\// || seen[$1]++ {{bad=1}} END {{exit NR == 0 || bad}}' \
          {params.outdir:q}/chromfile.txt
        while read -r name seqfile paf; do
          test -s "$seqfile" && test -s "$paf"
        done < {params.outdir:q}/chromfile.txt
        touch {output:q}
        """


rule mc_reuse_align:
    input: f"{MR_META}/split.done"
    output: f"{MR_META}/align.done"
    params:
        chromfile=f"{MR_STAGE}/split/chromfile.txt",
        store=lambda wc: mr_store(wc, "align"),
        outdir=f"{MR_STAGE}/align",
        max_cores=lambda wc: mr_profile(wc, "align")["max_cores"],
        max_memory=lambda wc: mr_profile(wc, "align")["max_memory"],
        cons_cores=lambda wc: mr_profile(wc, "align")["cons_cores"],
        cons_memory=lambda wc: mr_profile(wc, "align")["cons_memory"]
    log: f"{MR_ROOT}/logs/align.log"
    benchmark: f"{MR_ROOT}/benchmarks/align.tsv"
    conda: MR_ENV
    threads: lambda wc: mr_profile(wc, "align")["threads"]
    resources:
        mem_mb=lambda wc: mr_profile(wc, "align")["mem_mb"],
        disk_mb=lambda wc: mr_profile(wc, "align")["disk_mb"],
        runtime_min=lambda wc: mr_profile(wc, "align")["runtime_min"]
    shell:
        r"""
        mkdir -p "$(dirname {params.store:q})" {params.outdir:q} "$(dirname {log:q})"
        restart=()
        if [[ -e {params.store:q} ]]; then restart=(--restart); fi
        cactus-align {params.store:q} {params.chromfile:q} {params.outdir:q} \
          --batch --pangenome --outVG --reference CHM13 GRCh38 --maxLen 10000 \
          --consCores {params.cons_cores} --consMemory {params.cons_memory:q} \
          --batchSystem single_machine --clean=onSuccess \
          --maxCores {params.max_cores} --maxMemory {params.max_memory:q} \
          "${{restart[@]}}" > {log:q} 2>&1
        while read -r name seqfile paf; do
          test -s {params.outdir:q}/"$name".raw.vg
          test -s {params.outdir:q}/"$name".hal
        done < {params.chromfile:q}
        touch {output:q}
        """


rule mc_reuse_join:
    input: f"{MR_META}/align.done"
    output: f"{MR_META}/join.done"
    params:
        chromfile=f"{MR_STAGE}/split/chromfile.txt",
        align_root=f"{MR_STAGE}/align",
        store=lambda wc: mr_store(wc, "join"),
        outdir=f"{MR_ROOT}/graph",
        out_name=lambda wc: mr_build(wc)["out_name"],
        max_cores=lambda wc: mr_profile(wc, "join")["max_cores"],
        max_memory=lambda wc: mr_profile(wc, "join")["max_memory"],
        index_cores=lambda wc: mr_profile(wc, "join")["index_cores"],
        index_memory=lambda wc: mr_profile(wc, "join")["index_memory"]
    log: f"{MR_ROOT}/logs/join.log"
    benchmark: f"{MR_ROOT}/benchmarks/join.tsv"
    conda: MR_ENV
    threads: lambda wc: mr_profile(wc, "join")["threads"]
    resources:
        mem_mb=lambda wc: mr_profile(wc, "join")["mem_mb"],
        disk_mb=lambda wc: mr_profile(wc, "join")["disk_mb"],
        runtime_min=lambda wc: mr_profile(wc, "join")["runtime_min"]
    shell:
        r"""
        mkdir -p "$(dirname {params.store:q})" {params.outdir:q} "$(dirname {log:q})"
        vgs=()
        hals=()
        while read -r name seqfile paf; do
          test -s {params.align_root:q}/"$name".raw.vg
          test -s {params.align_root:q}/"$name".hal
          vgs+=({params.align_root:q}/"$name".raw.vg)
          hals+=({params.align_root:q}/"$name".hal)
        done < {params.chromfile:q}
        test "${{#vgs[@]}}" -gt 0
        for tool in gfaffix vcfbub; do
          command -v "$tool" >/dev/null || {{ echo "Missing Cactus join tool: $tool" >&2; exit 1; }}
        done
        restart=()
        if [[ -e {params.store:q} ]]; then restart=(--restart); fi
        cactus-graphmap-join {params.store:q} --vg "${{vgs[@]}}" --hal "${{hals[@]}}" \
          --outDir {params.outdir:q} --outName {params.out_name:q} \
          --reference CHM13 GRCh38 --vcfReference CHM13 \
          --gfa --gbz --vcf --giraffe --indexCores {params.index_cores} \
          --indexMemory {params.index_memory:q} --batchSystem single_machine \
          --clean=onSuccess --maxCores {params.max_cores} --maxMemory {params.max_memory:q} \
          "${{restart[@]}}" > {log:q} 2>&1
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
        touch {output:q}
        """


rule minigraph_cactus_reuse_hprc_plus_local_v1:
    input: f"{MR_RESULTS}/{MR_HPRC}/metadata/join.done"


rule minigraph_cactus_reuse_compare:
    input:
        local=f"{MR_RESULTS}/{MR_LOCAL}/metadata/join.done",
        hprc=f"{MR_RESULTS}/{MR_HPRC}/metadata/join.done",
        script=MR_COMPARE,
        common=MR_COMMON
    output:
        tsv=f"{MR_RESULTS}/comparison/local_v1_vs_hprc_plus_local_v1.tsv",
        chrom=f"{MR_RESULTS}/comparison/local_v1_vs_hprc_plus_local_v1.chromosomes.tsv",
        report=f"{MR_RESULTS}/comparison/local_v1_vs_hprc_plus_local_v1.md"
    params:
        local=mr_settings(MR_LOCAL),
        hprc=mr_settings(MR_HPRC),
        prefix=f"{MR_RESULTS}/comparison/local_v1_vs_hprc_plus_local_v1"
    log: f"{MR_RESULTS}/comparison/compare.log"
    conda: MR_ENV
    threads: 1
    resources: mem_mb=4096, runtime_min=120
    shell:
        "mkdir -p $(dirname {log:q}) && "
        "python {input.script:q} --local-settings {params.local:q} "
        "--hprc-settings {params.hprc:q} --out-prefix {params.prefix:q} "
        "> {log:q} 2>&1"
