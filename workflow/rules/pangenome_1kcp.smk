"""Chunked Cactus alignment from frozen Minigraph mappings. Use Snakefile.1kcp."""
import csv
import json
import re
from itertools import chain
from pathlib import Path
from textwrap import dedent

K = config['pangenome_1kcp']
KROOT = str(Path(K['results_root']))
KB = K['builds']
KR = K['resources']
KE = str(workflow.source_path('../envs/pangenome_1kcp.yaml'))
KC = str(workflow.source_path('../scripts/one_kcp_common.py'))
KP = str(workflow.source_path('../scripts/prepare_1kcp_inputs.py'))
KPART = str(workflow.source_path('../scripts/partition_1kcp_graph.py'))
KF = str(workflow.source_path('../scripts/format_1kcp_graph.py'))
KTMP = str(workflow.source_path('../scripts/one_kcp_scratch.sh'))
B = KROOT + '/{build}'
C = B + '/chunks/{chunk}'
A = B + '/assemblies/{sample}'
shell.executable('/bin/bash')

if not KB or any(not re.fullmatch(r'[A-Za-z0-9_.-]+', name) for name in KB):
    raise WorkflowError('Provide distinct filesystem-safe new build IDs')
if not 1 <= int(KR['threads']) <= 192 or not 0 < int(KR['mem_mb']) <= 300000:
    raise WorkflowError('Construction must remain within 192 CPUs / 300000 MB')
if not 0 < int(KR['prepare_mem_mb']) <= 300000:
    raise WorkflowError('Preparation memory exceeds ordinary-node policy')
if K['engine'] != 'cactus' or str(K['cactus_version']) != '3.3.0':
    raise WorkflowError('This entrypoint requires the pinned Cactus 3.3.0 environment')
if not 0 < int(KR['cons_memory']) <= int(KR['toil_max_memory']) < int(KR['mem_mb'])*1000000:
    raise WorkflowError('Cactus/Toil memory must leave headroom within the Snakemake reservation')
if any(int(KR[key]) <= 0 for key in ('time','disk_mb','prepare_time','prepare_disk_mb')):
    raise WorkflowError('Resource requests must be positive')
for build, data in KB.items():
    if int(data.get('max_gap_bp', K['max_gap_bp'])) < 0:
        raise WorkflowError('max_gap_bp must be nonnegative')
    if any(int(data.get(key, K[key])) <= 0 for key in (
            'max_interval_bp', 'max_chunk_input_bp', 'max_chunk_alignments')):
        raise WorkflowError('Interval and chunk input limits must be positive')
    source, dest = Path(data['source']), Path(KROOT)/build
    if not source.is_absolute() or not dest.is_absolute():
        raise WorkflowError('Use absolute SERVER paths for source and result roots')
    if source == dest or source in dest.parents or dest in source.parents:
        raise WorkflowError('New result directories must be separate from legacy inputs')


def ksetting(wc, key):
    return KB[wc.build].get(key, K[key])


def ksource(wc, suffix):
    return str(Path(KB[wc.build]['source']) / suffix)


def kinventory(wc):
    return checkpoints.one_kcp_inventory.get(build=wc.build).output.inventory


def ksamples(wc):
    return json.loads(Path(kinventory(wc)).read_text())['samples']


def ksample(wc):
    return next(row for row in ksamples(wc) if row['seqfile_name'] == wc.sample)


def kpartition(wc):
    return checkpoints.one_kcp_partition.get(build=wc.build).output.partition


def kchunks(wc):
    folder = Path(kpartition(wc))
    partition_info = json.loads((folder/'report.json').read_text())
    if not partition_info['feasible']:
        raise WorkflowError(f'Partition infeasible; inspect {folder}/report.json and chunks.tsv')
    with open(folder/'chunks.tsv') as handle:
        return [row['chunk'] for row in csv.DictReader(handle, delimiter='\t')]


def kadmit(wc):
    if wc.chunk not in kchunks(wc):
        raise WorkflowError('Chunk absent from approved partition')
    if not ksetting(wc, 'release_all_chunks') and wc.chunk not in ksetting(wc, 'pilot_chunks'):
        raise WorkflowError('Choose reviewed pilot_chunks in config, or release_all_chunks after measured pilots')
    return str(Path(kpartition(wc)))


def kassembly_products(wc):
    return [f'{KROOT}/{wc.build}/assemblies/{row["seqfile_name"]}/{kind}'
            for row in ksamples(wc) for kind in ('plan','sequences')]


def kstats(wc):
    return [f'{KROOT}/{wc.build}/chunks/{cid}/normalized.stats.json' for cid in kchunks(wc)]


def kgraphs(wc):
    return [f'{KROOT}/{wc.build}/graphs/{cid}.gfa' for cid in kchunks(wc)]


def kpilots(wc):
    paths = []
    for build in KB:
        pilots = KB[build].get('pilot_chunks', K['pilot_chunks'])
        if not pilots:
            raise WorkflowError(f'Select pilot_chunks for {build} from its preparation report first')
        paths.extend(f'{KROOT}/{build}/chunks/{cid}/normalized.gfa' for cid in pilots)
    return paths


def ksettings_json(wc):
    # These change the scientific result and are frozen. Resource/pilot/release
    # settings may change without repurposing the scientific build identity.
    keys = ('engine','cactus_version','cactus_max_len',
            'reference_prefix','target_prefix','chunk_bp','boundary_slack_bp',
            'max_reference_span_bp','max_seed_bp','max_interval_bp','max_gap_bp',
            'max_chunk_input_bp','max_chunk_alignments')
    return json.dumps({key: ksetting(wc, key) for key in keys}, sort_keys=True)


def ktable(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def kwrite_table(path, rows, fields=None):
    fields = fields or list(dict.fromkeys(chain.from_iterable(row for row in rows)))
    with open(path, 'w') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


wildcard_constraints:
    build='|'.join(re.escape(x) for x in KB),
    chunk='[0-9]{6}',
    sample='[A-Za-z0-9_.-]+'


# The default target stops at diagnostic reports, including infeasible results.
rule pangenome_1kcp_prepare:
    input: expand(KROOT + '/{build}/partition', build=list(KB))

rule pangenome_1kcp_ready:
    input: expand(KROOT + '/{build}/metadata/partition.ready', build=list(KB))

rule pangenome_1kcp_pilot:
    input: kpilots

rule pangenome_1kcp_inputs:
    input: expand(KROOT + '/{build}/input_sizes.tsv', build=list(KB))

rule pangenome_1kcp_all:
    input: expand(KROOT + '/{build}/report.json', build=list(KB))


checkpoint one_kcp_inventory:
    input:
        manifest=lambda wc: ksource(wc, 'metadata/input_manifest.tsv'),
        seqfile=lambda wc: ksource(wc, 'metadata/seqfile.txt'),
        crosswalk=lambda wc: ksource(wc, 'metadata/contig_map.tsv.gz'),
        frozen=lambda wc: ksource(wc, 'metadata/preflight.json'),
        graph=lambda wc: ksource(wc, 'stages/graphmap/cactus.gfa'),
        script=KP, common=KC
    output: inventory=B + '/metadata/inventory.json'
    params: source=lambda wc: KB[wc.build]['source'], settings=ksettings_json
    threads: 1
    resources:
        n=1, mem_mb=KR['prepare_mem_mb'], time=KR['prepare_time'],
        disk_mb=1000, limit_auto_threads=1
    log: B + '/logs/inventory.log'
    benchmark: B + '/benchmarks/inventory.tsv'
    conda: KE
    shell:
        'python {input.script:q} inventory --source {params.source:q} '
        '--settings {params.settings:q} --output {output.inventory:q} > {log:q} 2>&1'


checkpoint one_kcp_partition:
    input:
        inventory=kinventory,
        graph=lambda wc: ksource(wc, 'stages/graphmap/cactus.gfa'),
        script=KPART, common=KC
    output: partition=directory(B + '/partition')
    params:
        reference=lambda wc: ksetting(wc, 'reference_prefix'),
        target=lambda wc: ksetting(wc, 'chunk_bp'),
        slack=lambda wc: ksetting(wc, 'boundary_slack_bp'),
        span=lambda wc: ksetting(wc, 'max_reference_span_bp'),
        seed=lambda wc: ksetting(wc, 'max_seed_bp')
    threads: 1
    resources:
        n=1, mem_mb=KR['prepare_mem_mb'], time=KR['prepare_time'],
        disk_mb=KR['prepare_disk_mb'], limit_auto_threads=1
    log: B + '/logs/partition.log'
    benchmark: B + '/benchmarks/partition.tsv'
    conda: KE
    shell:
        'python {input.script:q} --gfa {input.graph:q} --inventory {input.inventory:q} '
        '--outdir {output.partition:q} --reference-prefix {params.reference:q} '
        '--chunk-bp {params.target} --slack-bp {params.slack} '
        '--max-span-bp {params.span} --max-seed-bp {params.seed} > {log:q} 2>&1'


rule one_kcp_ready:
    input: partition=kpartition
    output: B + '/metadata/partition.ready'
    threads: 1
    resources: n=1, mem_mb=1000, time=10, disk_mb=100, limit_auto_threads=1
    log: B + '/logs/ready.log'
    run:
        partition_info = json.loads((Path(input.partition)/'report.json').read_text())
        if not partition_info['feasible']:
            raise WorkflowError(f'Partition infeasible: {input.partition}/report.json. ' + partition_info['remedy'])
        Path(output[0]).write_text('Existing-node partition passed configured component limits.\n')
        Path(log[0]).write_text(f"Approved {partition_info['chunks']} components.\n")


rule one_kcp_sort_paf:
    input:
        ready=B + '/metadata/partition.ready', inventory=kinventory,
        paf=lambda wc: ksample(wc)['paf_identity']['path'], scratch=KTMP
    output: paf=temp(A + '/sorted.paf')
    params:
        scratch=K['scratch_root'], usage=A + '/sort.scratch.tsv',
        outdir=A
    threads: 1
    resources:
        n=1, mem_mb=KR['prepare_mem_mb'], time=KR['prepare_time'],
        disk_mb=KR['prepare_disk_mb'], limit_auto_threads=1
    log: A + '/sort.log'
    benchmark: B + '/benchmarks/assemblies/{sample}.sort.tsv'
    conda: KE
    shell:
        r'''
        exec > {log:q} 2>&1
        mkdir -p {params.outdir:q}
        source {input.scratch:q}
        one_kcp_scratch {params.scratch:q} {resources.disk_mb} {params.usage:q}
        LC_ALL=C sort --parallel={threads} -s -S 4G -T "$ONE_KCP_TMP" -k1,1 -k3,3n -k4,4n {input.paf:q} > {output.paf:q}
        '''


rule one_kcp_plan_intervals:
    input:
        paf=A + '/sorted.paf', inventory=kinventory, partition=kpartition,
        ready=B + '/metadata/partition.ready', script=KP, common=KC
    output: plan=directory(A + '/plan')
    params:
        prefix=lambda wc: ksetting(wc, 'target_prefix'),
        max_interval=lambda wc: ksetting(wc, 'max_interval_bp'),
        max_gap=lambda wc: ksetting(wc, 'max_gap_bp')
    threads: 1
    resources:
        n=1, mem_mb=KR['prepare_mem_mb'], time=KR['prepare_time'],
        disk_mb=KR['prepare_disk_mb'], limit_auto_threads=1
    log: A + '/plan.log'
    benchmark: B + '/benchmarks/assemblies/{sample}.plan.tsv'
    conda: KE
    shell:
        'python {input.script:q} plan --inventory {input.inventory:q} --sample {wildcards.sample:q} '
        '--nodes {input.partition:q}/nodes.tsv --paf {input.paf:q} --target-prefix {params.prefix:q} '
        '--max-interval-bp {params.max_interval} --max-gap-bp {params.max_gap} '
        '--outdir {output.plan:q} > {log:q} 2>&1; '
        'find {output.plan:q} -maxdepth 1 -name "*.paf" -exec gzip -n {{}} +'


rule one_kcp_extract_sequences:
    input:
        plan=A + '/plan', inventory=kinventory,
        fasta=lambda wc: ksample(wc)['source_path'], script=KP, common=KC
    output: sequences=directory(A + '/sequences')
    threads: 1
    resources:
        n=1, mem_mb=KR['prepare_mem_mb'], time=KR['prepare_time'],
        disk_mb=KR['prepare_disk_mb'], limit_auto_threads=1
    log: A + '/extract.log'
    benchmark: B + '/benchmarks/assemblies/{sample}.extract.tsv'
    conda: KE
    shell:
        'python {input.script:q} extract --inventory {input.inventory:q} --sample {wildcards.sample:q} '
        '--plan {input.plan:q} --outdir {output.sequences:q} > {log:q} 2>&1; '
        'find {output.sequences:q} -maxdepth 1 -name "*.fa" -exec gzip -n {{}} +'


rule one_kcp_collect:
    input:
        inventory=kinventory, partition=kadmit, assemblies=kassembly_products,
        sizes=B + '/input_sizes.tsv',
        ready=B + '/metadata/partition.ready', script=KP, common=KC
    output: inputs=directory(C + '/inputs')
    params:
        root=B,
        max_bases=lambda wc: ksetting(wc, 'max_chunk_input_bp'),
        max_alignments=lambda wc: ksetting(wc, 'max_chunk_alignments')
    threads: 1
    resources:
        n=1, mem_mb=KR['prepare_mem_mb'], time=KR['prepare_time'],
        disk_mb=KR['prepare_disk_mb'], limit_auto_threads=1
    log: C + '/collect.log'
    benchmark: B + '/benchmarks/chunks/{chunk}.collect.tsv'
    conda: KE
    shell:
        'python {input.script:q} collect --inventory {input.inventory:q} --build-root {params.root:q} '
        '--partition {input.partition:q} --chunk {wildcards.chunk:q} --outdir {output.inputs:q} '
        '--max-input-bp {params.max_bases} --max-alignments {params.max_alignments} > {log:q} 2>&1'


rule one_kcp_input_sizes:
    input:
        inventory=kinventory, partition=kpartition, assemblies=kassembly_products
    output: B + '/input_sizes.tsv'
    params:
        root=B,
        max_bases=lambda wc: ksetting(wc, 'max_chunk_input_bp'),
        max_alignments=lambda wc: ksetting(wc, 'max_chunk_alignments')
    threads: 1
    resources: n=1, mem_mb=4000, time=60, disk_mb=1000, limit_auto_threads=1
    log: B + '/logs/input_sizes.log'
    run:
        root = Path(params.root)
        inventory = json.loads(Path(input.inventory).read_text())
        chunks = {r['chunk']: dict(chunk=r['chunk'], reference=r['reference'],
                  reference_span=int(r['max_reference_span']), seed_bp=int(r['seed_bp']),
                  input_bases=int(r['seed_bp']), alignments=0, intervals=0, assembly_haplotypes=0)
                  for r in ktable(Path(input.partition)/'chunks.tsv')}
        for sample in inventory['samples']:
            stats = json.loads((root/'assemblies'/sample['seqfile_name']/'sequences/stats.json').read_text())
            for cid, values in stats['chunks'].items():
                row = chunks[cid]
                row['input_bases'] += values['bases']
                row['alignments'] += values['alignments']
                row['intervals'] += values['intervals']
                row['assembly_haplotypes'] += 1
        for row in chunks.values():
            row['status'] = ('SEED_ONLY' if row['assembly_haplotypes'] == 0 else
                             'PASS' if row['input_bases'] <= params.max_bases and
                             row['alignments'] <= params.max_alignments else 'OVERSIZED')
        ordered = sorted(chunks.values(), key=lambda r: (-r['input_bases'], -r['alignments'], r['chunk']))
        kwrite_table(output[0], ordered)
        Path(log[0]).write_text(f'Ranked {len(ordered)} full-cohort chunks.\n')


rule one_kcp_cactus_config:
    output: B + '/metadata/cactus.xml'
    threads: 1
    resources: n=1, mem_mb=1000, time=10, disk_mb=100, limit_auto_threads=1
    log: B + '/logs/cactus_config.log'
    conda: KE
    shell:
        dedent(r'''
        python - {output:q} > {log:q} 2>&1 <<'PY'
        import sys
        import xml.etree.ElementTree as ET
        from pathlib import Path
        from cactus.shared.common import cactusRootPath
        tree = ET.parse(Path(cactusRootPath()) / 'cactus_progressive_config.xml')
        node = tree.getroot().find('hal2vg')
        node.set('includeMinigraph', '1')
        node.set('includeAncestor', '0')
        tree.write(sys.argv[1], encoding='utf-8', xml_declaration=True)
        PY
        ''')


rule one_kcp_cactus_align:
    input:
        inputs=C + '/inputs', config=B + '/metadata/cactus.xml', scratch=KTMP
    output:
        hal=C + '/cactus.hal', vg=C + '/cactus.vg', gfa=C + '/cactus.gfa.gz'
    params:
        scratch=K['scratch_root'], usage=C + '/cactus.scratch.tsv',
        store=B + '/jobstores/{chunk}',
        max_memory=KR['toil_max_memory'], cons_memory=KR['cons_memory'],
        max_len=lambda wc: ksetting(wc, 'cactus_max_len')
    threads: KR['threads']
    resources:
        n=lambda wc, threads: threads, mem_mb=KR['mem_mb'], time=KR['time'],
        disk_mb=KR['disk_mb'], limit_auto_threads=lambda wc, threads: threads
    log: C + '/cactus.log'
    benchmark: B + '/benchmarks/chunks/{chunk}.cactus.tsv'
    conda: KE
    shell:
        r'''
        exec > {log:q} 2>&1
        source {input.scratch:q}
        one_kcp_scratch {params.scratch:q} {resources.disk_mb} {params.usage:q}
        export OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
        mkdir -p "$(dirname {params.store:q})"
        refs=()
        while IFS= read -r ref; do refs+=("$ref"); done < {input.inputs:q}/references.txt
        ref_args=()
        if (( ${{#refs[@]}} )); then ref_args=(--reference "${{refs[@]}}"); fi
        restart=()
        if [[ -d {params.store:q} ]]; then restart=(--restart); fi
        cactus-align {params.store:q} {input.inputs:q}/seqfile.txt \
            {input.inputs:q}/input.paf {output.hal:q} \
            --pangenome --outVG --outGFA "${{ref_args[@]}}" \
            --configFile {input.config:q} --maxLen {params.max_len} \
            --consCores {threads} --maxCores {threads} \
            --consMemory {params.cons_memory} --maxMemory {params.max_memory} \
            --maxDisk "$(({resources.disk_mb} * 1000000))" --workDir "$ONE_KCP_TMP" \
            --batchSystem single_machine --binariesMode local \
            --retryCount 0 --clean onSuccess "${{restart[@]}}"
        test -s {output.hal:q} && test -s {output.vg:q} && test -s {output.gfa:q}
        '''


rule one_kcp_cactus_gfa:
    input: C + '/cactus.gfa.gz'
    output: C + '/cactus.gfa'
    threads: 1
    resources: n=1, mem_mb=1000, time=120, disk_mb=KR['prepare_disk_mb'], limit_auto_threads=1
    log: C + '/cactus_gfa.log'
    conda: KE
    shell: 'gzip -dc {input:q} > {output:q} 2> {log:q}'


rule one_kcp_gfaffix:
    input: gfa=C + '/cactus.gfa', inputs=C + '/inputs', scratch=KTMP
    output: gfa=C + '/gfaffix.gfa'
    params: scratch=K['scratch_root'], usage=C + '/gfaffix.scratch.tsv'
    threads: KR['threads']
    resources:
        n=lambda wc, threads: threads, mem_mb=KR['mem_mb'], time=KR['time'],
        disk_mb=KR['disk_mb'], limit_auto_threads=lambda wc, threads: threads
    log: C + '/gfaffix.log'
    benchmark: B + '/benchmarks/chunks/{chunk}.gfaffix.tsv'
    conda: KE
    shell:
        r'''
        exec > {log:q} 2>&1
        source {input.scratch:q}
        one_kcp_scratch {params.scratch:q} {resources.disk_mb} {params.usage:q}
        export OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
        gfaffix --version
        protect=$(cat {input.inputs:q}/reference_regex.txt)
        gfaffix -c -p {threads} -x "$protect" {input.gfa:q} -o {output.gfa:q}
        '''


rule one_kcp_normalize:
    input: gfa=C + '/gfaffix.gfa', inputs=C + '/inputs', script=KF, common=KC
    output: gfa=C + '/normalized.gfa', stats=C + '/normalized.stats.json'
    threads: 1
    resources:
        n=1, mem_mb=KR['mem_mb'], time=KR['prepare_time'],
        disk_mb=KR['prepare_disk_mb'], limit_auto_threads=1
    log: C + '/normalize.log'
    benchmark: B + '/benchmarks/chunks/{chunk}.normalize.tsv'
    conda: KE
    shell:
        'python {input.script:q} normalize --gfa {input.gfa:q} --fasta {input.inputs:q}/input.fa '
        '--intervals {input.inputs:q}/intervals.tsv --path-names {input.inputs:q}/path_names.tsv '
        '--output {output.gfa:q} --stats {output.stats:q} > {log:q} 2>&1'


rule one_kcp_offsets:
    input: stats=kstats, partition=kpartition, script=KF, common=KC
    output: B + '/offsets.tsv'
    params: root=B
    threads: 1
    resources: n=1, mem_mb=1000, time=10, disk_mb=100, limit_auto_threads=1
    log: B + '/logs/offsets.log'
    conda: KE
    shell:
        'python {input.script:q} offsets --partition {input.partition:q} --root {params.root:q} '
        '--output {output:q} > {log:q} 2>&1'


rule one_kcp_global_ids:
    input: gfa=C + '/normalized.gfa', offsets=B + '/offsets.tsv', script=KF, common=KC
    output: B + '/graphs/{chunk}.gfa'
    threads: 1
    resources:
        n=1, mem_mb=KR['prepare_mem_mb'], time=KR['prepare_time'],
        disk_mb=KR['prepare_disk_mb'], limit_auto_threads=1
    log: C + '/global_ids.log'
    benchmark: B + '/benchmarks/chunks/{chunk}.global_ids.tsv'
    conda: KE
    shell:
        'python {input.script:q} offset --gfa {input.gfa:q} --offsets {input.offsets:q} '
        '--chunk {wildcards.chunk:q} --output {output:q} > {log:q} 2>&1'


rule one_kcp_combined_gfa:
    input: graphs=kgraphs
    output: B + '/graph/pangenome.gfa.gz'
    threads: 1
    resources:
        n=1, mem_mb=1000, time=KR['time'], disk_mb=KR['disk_mb'], limit_auto_threads=1
    log: B + '/logs/combined_gfa.log'
    benchmark: B + '/benchmarks/combined_gfa.tsv'
    conda: KE
    shell:
        r'''
        exec 2> {log:q}
        awk -F '\t' 'BEGIN {{print "H\tVN:Z:1.1\tRS:Z:CHM13 GRCh38"}}
          $1 != "H" {{print}}' {input.graphs:q} | gzip -n -1 -c > {output:q}
        '''


rule one_kcp_report:
    input:
        graphs=kgraphs, inventory=kinventory, assemblies=kassembly_products,
        offsets=B + '/offsets.tsv', combined=B + '/graph/pangenome.gfa.gz'
    output:
        summary=B + '/report.json', coverage=B + '/coverage.tsv',
        graphs=B + '/graphs.tsv', benchmarks=B + '/benchmarks.tsv',
        scratch=B + '/scratch.tsv'
    params: root=B
    threads: 1
    resources: n=1, mem_mb=4000, time=60, disk_mb=1000, limit_auto_threads=1
    log: B + '/logs/report.log'
    run:
        root = Path(params.root)
        inventory = json.loads(Path(input.inventory).read_text())
        coverage = []
        for sample in inventory['samples']:
            stats = json.loads((root/'assemblies'/sample['seqfile_name']/'sequences/stats.json').read_text())
            coverage.append({key: stats[key] for key in ('sample','input_bases','represented_bases',
                             'excluded_unaligned_bases','alignments','intervals','ambiguous_assignments')})
        kwrite_table(output.coverage, coverage)
        offsets = {r['chunk']: r for r in ktable(input.offsets)}
        chunks = []
        for row in ktable(root/'partition/chunks.tsv'):
            cid = row['chunk']
            stats = json.loads((root/'chunks'/cid/'normalized.stats.json').read_text())
            inputs = json.loads((root/'chunks'/cid/'inputs/stats.json').read_text())
            chunks.append({**row, **stats, 'input_bases': inputs['input_bases'],
                           'assembly_haplotypes': inputs['assembly_haplotypes'],
                           'id_offset': offsets[cid]['offset'], 'gfa': str(root/'graphs'/f'{cid}.gfa'),
                           'hal': str(root/'chunks'/cid/'cactus.hal'),
                           'vg': str(root/'chunks'/cid/'cactus.vg')})
        kwrite_table(output.graphs, chunks)
        benchmarks = []
        for path in sorted((root/'benchmarks').rglob('*.tsv')):
            benchmarks.extend(dict(job=str(path.relative_to(root/'benchmarks')), **row) for row in ktable(path))
        kwrite_table(output.benchmarks, benchmarks, fields=None if benchmarks else ['job'])
        scratch = []
        for path in sorted(root.glob('chunks/*/*.scratch.tsv')):
            scratch.extend(dict(job=str(path.relative_to(root)), **row) for row in ktable(path))
        kwrite_table(output.scratch, scratch, ['job','sampled_peak_kib'])
        summary_data = dict(engine='cactus', cactus_version=K['cactus_version'],
                      chunks=len(chunks), graph_nodes=sum(r['nodes'] for r in chunks),
                      graph_bases=sum(r['graph_bases'] for r in chunks),
                      input_bases=sum(r['input_bases'] for r in coverage),
                      represented_bases=sum(r['represented_bases'] for r in coverage),
                      excluded_unaligned_bases=sum(r['excluded_unaligned_bases'] for r in coverage),
                      graph_manifest=str(output.graphs), coverage=str(output.coverage),
                      combined_gfa=str(input.combined),
                      benchmarks=str(output.benchmarks), scratch=str(output.scratch),
                      interpretation='Cactus bp-resolution chunk graphs and per-chunk HALs; no inter-chunk edges, whole-genome HAL, VCF, GBZ or Giraffe indexes. No read-k-mer filtering.',
                      limitations='Existing-node boundary adaptation; inherited mapping provenance. Resource ceilings require measured pilots; scratch peaks are sampled.')
        Path(output.summary).write_text(json.dumps(summary_data, indent=2, sort_keys=True)+'\n')
        Path(log[0]).write_text(f'Reported {len(chunks)} completed graph chunks.\n')
