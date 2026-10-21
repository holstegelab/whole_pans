#!/usr/bin/env python3
"""Freeze reuse inventory, check feasibility, plan disjoint intervals, stream FASTAs."""
import argparse
import csv
import gzip
import hashlib
import itertools
import json
import re
import shutil
from collections import defaultdict
from pathlib import Path
from one_kcp_common import (table, rows, open_text, write_table, write_json, identity,
                            unchanged, check_paf, Writers, new_output_dir,
                            cactus_bases, CACTUS_BASE_POLICY)


def inventory(a):
    root, out = Path(a.source), Path(a.output)
    manifest = table(root/'metadata/input_manifest.tsv')
    seqfile = dict(line.split() for line in (root/'metadata/seqfile.txt').read_text().splitlines() if line.strip())
    if not manifest or len(seqfile) != len(manifest):
        raise ValueError('Empty or inconsistent frozen manifest/seqfile')
    contigs = table(root/'metadata/contig_map.tsv.gz')
    names, paths, haplotypes = set(), set(), set()
    samples = []
    for row in manifest:
        label, path = row['seqfile_name'], row['source_path']
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', label) or label in names or path in paths:
            raise ValueError('Duplicate or unsafe assembly label/path')
        hap = (row['sample'], row['haplotype'])
        if hap in haplotypes or not hap[1].isdigit() or re.search(r'[\s#]', hap[0]) or hap[0] == '_MINIGRAPH_':
            raise ValueError('Duplicate or invalid sample/haplotype identity')
        names.add(label); paths.add(path); haplotypes.add(hap)
        if seqfile.get(label) != path:
            raise ValueError(f'Seqfile differs from frozen input manifest: {label}')
        source = identity(path)
        if source['size'] != int(row['size']) or source['mtime_ns'] != int(row['mtime_ns']):
            raise ValueError(f'Original FASTA identity changed: {path}')
        paf = identity(root/'stages/graphmap/assemblies'/f'{label}.paf')
        # Empty assembly PAFs are allowed and explicitly accounted as unmapped.
        samples.append(dict(**row, fasta_identity=source, paf_identity=paf))
    by_sample = defaultdict(list)
    seen_queries, seen_originals, seen_pansn = set(), set(), set()
    for row in contigs:
        label = row['seqfile_name']
        if label not in names or int(row['length_bp']) <= 0:
            raise ValueError('Invalid crosswalk row')
        query = f"id={label}|{row['cactus_contig_id']}"
        original = (label, row['original_header'].split()[0])
        if query in seen_queries or original in seen_originals or row['pansn_source_name'] in seen_pansn:
            raise ValueError('Duplicate contig crosswalk identity')
        seen_queries.add(query); seen_originals.add(original); seen_pansn.add(row['pansn_source_name'])
        by_sample[label].append(row)
    for sample in samples:
        cs = by_sample[sample['seqfile_name']]
        if len(cs) != int(sample['contig_count']) or sum(int(c['length_bp']) for c in cs) != int(sample['total_bases']):
            raise ValueError('Crosswalk totals disagree with frozen manifest')
        prefix = f"{sample['sample']}#{sample['haplotype']}#"
        if any(not c['pansn_source_name'].startswith(prefix) for c in cs):
            raise ValueError('Crosswalk sample/haplotype differs from manifest')
    preflight = json.loads((root/'metadata/preflight.json').read_text())
    if preflight['seqfile_sha256'] != hashlib.sha256((root/'metadata/seqfile.txt').read_bytes()).hexdigest():
        raise ValueError('Frozen seqfile checksum mismatch')
    frozen = {r['name']: r for r in preflight['inputs']}
    if set(frozen) != names or len(frozen) != len(preflight['inputs']):
        raise ValueError('Preflight and manifest assembly membership differ')
    for sample in samples:
        prior = frozen[sample['seqfile_name']]
        for key, old_key in (('source_path','path'), ('sample','sample'), ('haplotype','haplotype'),
                             ('size','size'), ('mtime_ns','mtime_ns'),
                             ('sha256_uncompressed','sha256'), ('contig_count','contig_count'),
                             ('total_bases','total_bases')):
            if str(sample[key]) != str(prior[old_key]):
                raise ValueError(f"Frozen preflight disagrees with manifest: {sample['seqfile_name']} {key}")
    metadata = [identity(root/'metadata'/name, True) for name in
                ('seqfile.txt','input_manifest.tsv','contig_map.tsv.gz','preflight.json')]
    graph = identity(root/'stages/graphmap/cactus.gfa')
    if not graph['size']:
        raise ValueError('Empty mapping graph')
    settings = json.loads(a.settings)
    reference_lengths = {f"id={r['seqfile_name']}|{r['cactus_contig_id']}": int(r['length_bp'])
                         for r in contigs if f"id={r['seqfile_name']}|{r['cactus_contig_id']}".startswith(settings['reference_prefix'])}
    if not reference_lengths:
        raise ValueError('Reference prefix does not match the frozen contig crosswalk')
    result = dict(source=str(root), samples=samples, graph=graph, metadata=metadata,
                  settings=settings, reference_lengths=reference_lengths,
                  identity_policy='metadata SHA256; large-file size+mtime; FASTA uncompressed SHA256 checked during extraction')
    # A persistent lock is outside Snakemake-managed outputs. Never repurpose a build ID.
    lock = out.with_name('inventory.lock.json')
    out.parent.mkdir(parents=True, exist_ok=True)
    if lock.exists():
        if json.loads(lock.read_text()) != result:
            raise ValueError('Build ID is locked to different inputs/settings; use a new build ID')
    else:
        with lock.open('x') as handle:
            json.dump(result, handle, indent=2, sort_keys=True)
    write_json(out, result)


def sample_data(a):
    inv = json.loads(Path(a.inventory).read_text())
    for record in inv['metadata']:
        unchanged(record)
    sample = next(s for s in inv['samples'] if s['seqfile_name'] == a.sample)
    contigs = [r for r in table(Path(inv['source'])/'metadata/contig_map.tsv.gz') if r['seqfile_name'] == a.sample]
    return sample, contigs


def plan(a):
    if a.max_gap_bp < 0 or a.max_interval_bp <= 0:
        raise ValueError('max-gap-bp must be nonnegative and max-interval-bp positive')
    sample, contigs = sample_data(a)
    unchanged(sample['paf_identity'])
    crosswalk = {f"id={a.sample}|{c['cactus_contig_id']}": c for c in contigs}
    nodes = {r['node']: (r['chunk'], int(r['length'])) for r in table(a.nodes)}
    out = new_output_dir(a.outdir)
    writers, intervals = Writers(), []
    chunk_stats = defaultdict(lambda: dict(bases=0, alignments=0, intervals=0))
    seen, count, last_query = set(), 0, None
    try:
        for query, records in itertools.groupby(rows(a.paf), key=lambda r: r[0]):
            if query not in crosswalk or query in seen or (last_query is not None and query < last_query):
                raise ValueError(f'Unknown/unsorted PAF query: {query}')
            seen.add(query); last_query = query
            meta = crosswalk[query]
            group, previous_start = None, -1

            def emit(group):
                if group is None:
                    return
                chunk, begin, end, alignments = group
                if end-begin > a.max_interval_bp:
                    raise ValueError(f'Interval exceeds limit: {query}:{begin}-{end}')
                name = f"{meta['pansn_source_name']}[{begin}-{end}]"
                if re.search(r'\s', name) or meta['pansn_source_name'].count('#') != 2:
                    raise ValueError(f'Invalid PanSN name: {name}')
                intervals.append(dict(query=query, original=meta['original_header'].split()[0],
                                      name=name, chunk=chunk, start=begin, end=end,
                                      sample=sample['sample'], haplotype=sample['haplotype']))
                chunk_stats[chunk]['bases'] += end-begin
                chunk_stats[chunk]['alignments'] += len(alignments)
                chunk_stats[chunk]['intervals'] += 1
                for row in alignments:
                    row = row.copy()
                    row[0], row[1] = name, str(end-begin)
                    row[2], row[3] = str(int(row[2])-begin), str(int(row[3])-begin)
                    row[5] = '_MINIGRAPH_#0#' + row[5][len(a.target_prefix):]
                    # Keep CIGAR and parent-alignment evidence used by Cactus's PAF
                    # filter. Coordinate-bearing annotations would be stale after slicing.
                    row = row[:12] + [tag for tag in row[12:]
                                     if tag.startswith(('cg:Z:', 'gl:i:', 'gi:f:', 'AS:i:'))]
                    writers.get(out/f'{chunk}.paf').write('\t'.join(row)+'\n')

            for row in records:
                start, end = check_paf(row)
                if start < previous_start or int(row[1]) != int(meta['length_bp']):
                    raise ValueError(f'Unsorted coordinates or query-length mismatch: {query}')
                previous_start = start
                if not row[5].startswith(a.target_prefix):
                    raise ValueError(f'Unexpected PAF target namespace: {row[5]}')
                node = row[5][len(a.target_prefix):]
                if node not in nodes or int(row[6]) != nodes[node][1]:
                    raise ValueError(f'PAF target absent or length mismatch: {node}')
                chunk = nodes[node][0]
                count += 1
                if group and start < group[2] and chunk != group[0]:
                    conflict = dict(query=query, previous_chunk=group[0], chunk=chunk,
                                    overlap_start=start, overlap_end=min(end, group[2]))
                    write_json(out/'ambiguity.json', conflict)
                    raise ValueError(f'Ambiguous overlapping chunk assignments: {conflict}')
                if group and chunk == group[0] and start <= group[2] + a.max_gap_bp:
                    group[2] = max(group[2], end)
                    group[3].append(row)
                else:
                    emit(group)
                    group = [chunk, start, end, [row]]
            emit(group)
    finally:
        writers.close()
    represented = sum(r['end']-r['start'] for r in intervals)
    write_table(out/'intervals.tsv', ['query','original','name','chunk','start','end','sample','haplotype'], intervals)
    write_json(out/'stats.json', dict(sample=a.sample, input_bases=int(sample['total_bases']),
               represented_bases=represented, excluded_unaligned_bases=int(sample['total_bases'])-represented,
               alignments=count, intervals=len(intervals), ambiguous_assignments=0,
               chunks=dict(chunk_stats),
               gap_policy=f'Keep consecutive same-chunk visits across gaps <= {a.max_gap_bp} bp; '
                          'split at other chunks and larger gaps; no end padding.'))
    unchanged(sample['paf_identity'])


def extract(a):
    sample, contigs = sample_data(a)
    unchanged(sample['fasta_identity'])
    expected = {c['original_header'].split()[0]: c for c in contigs}
    requests = defaultdict(list)
    for row in table(Path(a.plan)/'intervals.tsv'):
        row['start'], row['end'] = int(row['start']), int(row['end'])
        requests[row['original']].append(row)
    for iv in requests.values():
        iv.sort(key=lambda r: r['start'])
    out = new_output_dir(a.outdir)
    writers, seen, sequence_hashes = Writers(), set(), {}
    digest, current, pos, iv_index = hashlib.sha256(), None, 0, 0
    lower_bases, ambiguity_bases = 0, 0
    opener = gzip.open if sample['source_path'].endswith('.gz') else open

    def finish_contig():
        if current is not None and pos != int(expected[current]['length_bp']):
            raise ValueError(f'FASTA/crosswalk length mismatch: {current}')

    try:
        with opener(sample['source_path'], 'rb') as handle:
            for raw in handle:
                digest.update(raw)
                if raw.startswith(b'>'):
                    finish_contig()
                    current = raw[1:].split()[0].decode()
                    if current not in expected or current in seen:
                        raise ValueError(f'Unknown/duplicate FASTA identifier: {current}')
                    seen.add(current); pos = 0; iv_index = 0
                    continue
                bases = raw.rstrip(b'\r\n')
                if current is None or not bases:
                    raise ValueError('Unsupported FASTA layout: sequence before header or empty sequence line')
                lower_bases += len(bases) - len(bases.translate(None, b'acgtnuwsmkrybdhv'))
                bases, masked = cactus_bases(bases, f'{a.sample}:{current} at offset {pos}')
                ambiguity_bases += masked
                iv = requests[current]
                line_end = pos + len(bases)
                while iv_index < len(iv) and iv[iv_index]['start'] < line_end:
                    r = iv[iv_index]
                    left, right = max(pos, r['start']), min(line_end, r['end'])
                    if left < right:
                        writer = writers.get(out/f"{r['chunk']}.fa")
                        if left == r['start']:
                            writer.write('>'+r['name']+'\n')
                            sequence_hashes[r['name']] = hashlib.sha256()
                        piece = bases[left-pos:right-pos]
                        sequence_hashes[r['name']].update(piece)
                        writer.write(piece.decode()+'\n')
                    if r['end'] <= line_end:
                        iv_index += 1
                    else:
                        break
                pos = line_end
            finish_contig()
    finally:
        writers.close()
    if seen != set(expected) or digest.hexdigest() != sample['sha256_uncompressed']:
        raise ValueError('Original FASTA differs from frozen contig map or uncompressed SHA256')
    unchanged(sample['fasta_identity'])
    intervals = table(Path(a.plan)/'intervals.tsv')
    for r in intervals:
        r['sha256'] = sequence_hashes[r['name']].hexdigest()
    write_table(out/'intervals.tsv', ['query','original','name','chunk','start','end','sample','haplotype','sha256'], intervals)
    stats = json.loads((Path(a.plan)/'stats.json').read_text())
    stats.update(fasta_sha256=digest.hexdigest(), lowercase_bases_uppercased=lower_bases,
                 ambiguity_bases_masked=ambiguity_bases, base_normalization=CACTUS_BASE_POLICY,
                 sequence_access='one sequential integrity-checked pass; originals never indexed or rewritten')
    write_json(out/'stats.json', stats)


def fasta_lengths(path):
    name, length = None, 0
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    yield name, length
                name, length = line[1:].strip(), 0
            else:
                if name is None:
                    raise ValueError('Sequence before FASTA header')
                length += len(line.strip())
        if name is not None:
            yield name, length


def collect(a):
    """Write a Cactus seqfile, per-event FASTAs and renamed, sliced mappings."""
    inv = json.loads(Path(a.inventory).read_text())
    root, out = Path(a.build_root), new_output_dir(a.outdir)
    fasta_files = [('_MINIGRAPH_', '_MINIGRAPH_', '0',
                    Path(a.partition)/'seed'/f'{a.chunk}.fa')]
    paf_files, interval_rows, haps = [], [], set()
    for sample in inv['samples']:
        label = sample['seqfile_name']
        extracted = root/'assemblies'/label/'sequences'
        selected = [r for r in table(extracted/'intervals.tsv') if r['chunk'] == a.chunk]
        if selected:
            interval_rows.extend(selected)
            haps.add((sample['sample'], sample['haplotype']))
            fasta_files.append((label, sample['sample'], sample['haplotype'],
                                extracted/f'{a.chunk}.fa.gz'))
            paf_files.append(root/'assemblies'/label/'plan'/f'{a.chunk}.paf.gz')
    if not haps:
        seed_bp = sum(length for _, length in fasta_lengths(fasta_files[0][3]))
        raise ValueError(f'Chunk {a.chunk} has {seed_bp} seed bases but no mapped assembly intervals; '
                         'Cactus requires at least one biological event alongside _MINIGRAPH_. '
                         'Review this component before releasing the build.')
    lengths, cactus_names, path_rows = {}, {}, []
    for event, sample, haplotype, path in fasta_files:
        parts = event.rsplit('.', 1)
        event_identity = (parts[0], parts[1]) if len(parts) == 2 and parts[1].isdigit() else (event, '0')
        if event_identity != (sample, haplotype):
            raise ValueError(f'Cactus event name does not encode manifest sample/haplotype: {event}')
        for number, (name, length) in enumerate(fasta_lengths(path), 1):
            if name in lengths or length <= 0:
                raise ValueError(f'Duplicate/empty sequence: {name}')
            lengths[name] = length
            if event == '_MINIGRAPH_':
                if not name.startswith('_MINIGRAPH_#0#'):
                    raise ValueError(f'Unexpected seed FASTA name: {name}')
                contig = name[len('_MINIGRAPH_#0#'):]
                if not re.fullmatch(r'[A-Za-z0-9_.-]+', contig) or '_sub_' in contig:
                    raise ValueError(f'Unsafe seed node name for Cactus: {contig}')
            else:
                # Short local IDs avoid Cactus/hal2vg interpreting original brackets,
                # colons or _sub_ suffixes as native path subranges.
                contig = f'kcp{number:09d}'
            cactus_name = f'id={event}|{contig}'
            cactus_names[name] = cactus_name
            path_rows.append(dict(name=name, event=event, contig=contig,
                                  cactus_name=cactus_name,
                                  pansn_name=f'{sample}#{haplotype}#{contig}',
                                  length=length,
                                  kind='seed' if event == '_MINIGRAPH_' else 'assembly'))
    bases = sum(lengths.values())
    if bases > a.max_input_bp:
        raise ValueError(f'Chunk {a.chunk}: {bases} input bases exceed limit {a.max_input_bp}; review chunking')
    for row in interval_rows:
        if lengths.get(row['name']) != int(row['end'])-int(row['start']):
            raise ValueError(f"Extracted sequence length mismatch: {row['name']}")
    required_bytes = 3*bases + 20*sum(path.stat().st_size for path in paf_files)
    if shutil.disk_usage(out).free < required_bytes:
        raise ValueError(f'Insufficient output-filesystem capacity: need {required_bytes} free bytes')
    fasta_dir = out/'fastas'
    fasta_dir.mkdir()
    seqfile_rows = []
    with open(out/'input.fa', 'w') as combined:
        for event, sample, haplotype, path in fasta_files:
            fasta_path = fasta_dir/f'{event}.fa'
            seqfile_rows.append((event, str(fasta_path.absolute())))
            with open_text(path) as source, fasta_path.open('w') as dest:
                for line in source:
                    combined.write(line)
                    if line.startswith('>'):
                        name = line[1:].strip()
                        dest.write('>' + cactus_names[name].split('|', 1)[1] + '\n')
                    else:
                        dest.write(line)
    with (out/'seqfile.txt').open('w') as handle:
        for event, path in seqfile_rows:
            handle.write(f'{event}\t{path}\n')
    present_events = {event for event, _ in seqfile_rows}
    references = [event for event in ('CHM13', 'GRCh38') if event in present_events]
    (out/'references.txt').write_text(''.join(event+'\n' for event in references))
    protected = []
    if references:
        for row in path_rows:
            if row['event'] == references[0]:
                protected.extend(re.escape(alias) for alias in
                                 (row['pansn_name'], row['pansn_name']+'#0',
                                  row['cactus_name'], row['name']))
    (out/'reference_regex.txt').write_text('^(?:'+'|'.join(protected)+')$\n' if protected else '^$\n')
    alignments = 0
    with (out/'input.paf').open('w') as dest:
        for path in paf_files:
            for row in rows(path):
                check_paf(row)
                if lengths.get(row[0]) != int(row[1]) or lengths.get(row[5]) != int(row[6]):
                    raise ValueError(f'Collected FASTA/PAF namespace or length mismatch: {row[0]}')
                alignments += 1
                if alignments > a.max_alignments:
                    raise ValueError('Chunk exceeds alignment-count limit')
                row[0], row[5] = cactus_names[row[0]], cactus_names[row[5]]
                dest.write('\t'.join(row)+'\n')
    write_table(out/'intervals.tsv', ['query','original','name','chunk','start','end','sample','haplotype','sha256'], interval_rows)
    write_table(out/'path_names.tsv',
                ['name','event','contig','cactus_name','pansn_name','length','kind'], path_rows)
    write_json(out/'stats.json', dict(chunk=a.chunk, input_bases=bases,
               sequences=len(lengths), seed_sequences=sum(n.startswith('_MINIGRAPH_#') for n in lengths),
               alignments=alignments, assembly_haplotypes=len(haps),
               haplotypes=sorted(haps), seed_only=not haps,
               cactus_events=[event for event, _ in seqfile_rows],
               cactus_references=references,
               cactus_seqfile='seqfile.txt', cactus_path_names='path_names.tsv'))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    for command, func, keys in [
        ('inventory', inventory, ('source','settings','output')),
        ('plan', plan, ('inventory','sample','nodes','paf','target-prefix','outdir')),
        ('extract', extract, ('inventory','sample','plan','outdir')),
        ('collect', collect, ('inventory','build-root','partition','chunk','outdir'))]:
        s = sub.add_parser(command)
        for key in keys:
            s.add_argument('--'+key, required=True)
        if command == 'plan':
            s.add_argument('--max-interval-bp', type=int, required=True)
            s.add_argument('--max-gap-bp', type=int, default=1000000)
        if command == 'collect':
            s.add_argument('--max-input-bp', type=int, required=True)
            s.add_argument('--max-alignments', type=int, required=True)
        s.set_defaults(func=func)
    args = p.parse_args()
    args.func(args)
