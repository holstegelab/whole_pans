#!/usr/bin/env python3
"""Check exact path spelling, remove only artificial paths, and assign dense IDs."""
import argparse
import hashlib
import json
import re
from pathlib import Path
from one_kcp_common import rows, table, write_json, write_table


COMPLEMENT = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def expected_sequences(fasta):
    result, name, digest, length = {}, None, None, 0
    with open(fasta) as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    result[name] = (length, digest.hexdigest())
                name = line[1:].strip()
                if name in result:
                    raise ValueError(f'Duplicate input sequence: {name}')
                digest, length = hashlib.sha256(), 0
            else:
                seq = line.strip().upper().encode()
                digest.update(seq); length += len(seq)
        if name is not None:
            result[name] = (length, digest.hexdigest())
    return result


def cactus_path_aliases(path_names, expected):
    """Map Cactus's explicit native names back to the extracted interval names."""
    aliases, canonical_names = {}, set()
    for row in table(path_names):
        name = row['name']
        if name in canonical_names or name not in expected or int(row['length']) != expected[name][0]:
            raise ValueError(f'Cactus path crosswalk differs from collected FASTA: {name}')
        canonical_names.add(name)
        # hal2vg emits reference-sense paths without a phase block, and
        # haplotype-sense paths with phase block 0. Both spellings are explicit.
        for alias in (name, row['pansn_name'], row['pansn_name']+'#0', row['cactus_name']):
            if alias in aliases and aliases[alias] != name:
                raise ValueError(f'Ambiguous Cactus path alias: {alias}')
            aliases[alias] = name
    if canonical_names != set(expected):
        raise ValueError('Cactus path crosswalk does not cover every collected input sequence')
    return aliases


def normalize(a):
    expected = expected_sequences(a.fasta)
    aliases = cactus_path_aliases(a.path_names, expected)
    intervals = {r['name']: r for r in table(a.intervals)}
    seqs, links = {}, set()
    for row in rows(a.gfa):
        if row[0] == 'S':
            if row[1] in seqs or row[2] == '*' or any(c not in 'ACGTNacgtn' for c in row[2]):
                raise ValueError('Invalid or duplicate GFA segment')
            seqs[row[1]] = row[2].upper()
        elif row[0] == 'L':
            if row[2] not in ('+','-') or row[4] not in ('+','-') or row[5] != '0M':
                raise ValueError('Expected oriented blunt-ended GFA links')
            links.add((row[1]+row[2], row[3]+row[4]))
            flip = {'+':'-', '-':'+'}
            links.add((row[3]+flip[row[4]], row[1]+flip[row[2]]))
        elif row[0] not in ('H','P'):
            raise ValueError(f'Expected S/L/P graph, found {row[0]}')
    # Keep ALL segments and edges, including seed-only content. Only seed P lines
    # are removed, so off-reference graph sequence cannot silently disappear.
    mapping = {node: str(i+1) for i, node in enumerate(sorted(seqs, key=lambda n: (len(n), n)))}
    seen, assembly_nodes, link_count, path_count = set(), set(), 0, 0
    with open(a.output, 'w') as out:
        references = [sample for sample in ('CHM13','GRCh38')
                      if any(meta['sample'] == sample for meta in intervals.values())]
        out.write('H\tVN:Z:1.1' + ('\tRS:Z:'+' '.join(references) if references else '') + '\n')
        for row in rows(a.gfa):
            if row[0] == 'S':
                out.write(f'S\t{mapping[row[1]]}\t{seqs[row[1]]}\n')
            elif row[0] == 'L':
                if row[5] != '0M':
                    raise ValueError('Cannot use W walks with overlapping GFA links')
                out.write(f'L\t{mapping[row[1]]}\t{row[2]}\t{mapping[row[3]]}\t{row[4]}\t0M\n')
                link_count += 1
            elif row[0] == 'P':
                raw_name = row[1]
                name = aliases.get(raw_name)
                if name is None or name in seen or (row[3] != '*' and any(x != '0M' for x in row[3].split(','))):
                    raise ValueError(f'Unexpected/duplicate path or path overlaps: {raw_name}')
                seen.add(name)
                digest, length, walk = hashlib.sha256(), 0, []
                steps = row[2].split(',')
                if any((left, right) not in links for left, right in zip(steps, steps[1:])):
                    raise ValueError(f'Path traverses a missing oriented link: {name}')
                for step in steps:
                    if step[-1] not in '+-':
                        raise ValueError('Invalid path orientation')
                    node = step[:-1]
                    seq = seqs[node]
                    if step[-1] == '-':
                        seq = seq.translate(COMPLEMENT)[::-1]
                    digest.update(seq.encode()); length += len(seq)
                    walk.append(('>' if step[-1]=='+' else '<') + mapping[node])
                if (length, digest.hexdigest()) != expected[name]:
                    raise ValueError(f'Path sequence changed during construction: {name}')
                if name.startswith('_MINIGRAPH_#0#'):
                    continue
                meta = intervals[name]
                if length != int(meta['end'])-int(meta['start']) or digest.hexdigest() != meta['sha256']:
                    raise ValueError(f'Path differs from frozen extracted interval: {name}')
                assembly_nodes.update(s[:-1] for s in steps)
                # Use original FASTA first-token contig IDs; preserve the full
                # sanitized/PanSN crosswalk in intervals.tsv for downstream use.
                out.write('\t'.join(['W',meta['sample'],meta['haplotype'],meta['original'],
                                     meta['start'],meta['end'],''.join(walk)])+'\n')
                path_count += 1
    if seen != set(expected):
        raise ValueError(f'Construction lost {len(set(expected)-seen)} input paths')
    write_json(a.stats, dict(nodes=len(mapping), edges=link_count, assembly_walks=path_count,
               graph_bases=sum(map(len, seqs.values())),
               nodes_without_assembly_walk=len(set(seqs)-assembly_nodes),
               cactus_path_names_verified=True,
               all_input_path_sequences_verified=True))


def offsets(a):
    chunks = table(Path(a.partition)/'chunks.tsv')
    result, total = [], 0
    for chunk in chunks:
        stats = json.loads((Path(a.root)/'chunks'/chunk['chunk']/'normalized.stats.json').read_text())
        result.append(dict(chunk=chunk['chunk'], offset=total, nodes=stats['nodes']))
        total += stats['nodes']
    write_table(a.output, ['chunk','offset','nodes'], result)


def offset(a):
    entry = next(r for r in table(a.offsets) if r['chunk'] == a.chunk)
    shift = int(entry['offset'])
    count = 0
    with open(a.output, 'w') as out:
        for row in rows(a.gfa):
            if row[0] == 'S':
                row[1] = str(int(row[1])+shift); count += 1
            elif row[0] == 'L':
                row[1], row[3] = str(int(row[1])+shift), str(int(row[3])+shift)
            elif row[0] == 'W':
                row[6] = re.sub(r'[0-9]+', lambda m: str(int(m[0])+shift), row[6])
            out.write('\t'.join(row)+'\n')
    if count != int(entry['nodes']):
        raise ValueError('Dense node count differs from offset manifest')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(required=True)
    for name, func, keys in [
        ('normalize', normalize, ('gfa','fasta','intervals','output','stats')),
        ('offsets', offsets, ('partition','root','output')),
        ('offset', offset, ('gfa','offsets','chunk','output'))]:
        s = sub.add_parser(name)
        for key in keys:
            s.add_argument('--'+key, required=True)
        if name == 'normalize':
            s.add_argument('--path-names', required=True)
        s.set_defaults(func=func)
    args = p.parse_args(); args.func(args)
