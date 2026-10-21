#!/usr/bin/env python3
"""Conservative, existing-node bridge partitioning of the ORIGINAL cactus.gfa.

Links are canonicalized by oriented node ends. Reciprocal duplicate records
represent one physical link; different orientations remain parallel edges.
An undirected bridge is a sufficient (not necessary) bubble-safe boundary.
Source nodes, sequences, IDs and mapping coordinates are never edited.
Derived seed FASTAs use Cactus-compatible base normalization.
"""
import argparse
import bisect
from collections import defaultdict
from pathlib import Path
from one_kcp_common import (rows, write_table, write_json, Writers, unchanged,
                            new_output_dir, cactus_bases, CACTUS_BASE_POLICY)
import json


def link_key(a, ao, b, bo):
    # L a + b + joins a's right end to b's left end.
    ends = ((a, ao == '+'), (b, bo == '-'))
    return tuple(sorted(ends))


def bridges(adjacency):
    """Iterative Tarjan DFS with edge IDs (handles parallel edges and self loops)."""
    discovery, low, result = {}, {}, set()
    tick = 0
    for start in adjacency:
        if start in discovery:
            continue
        discovery[start] = low[start] = tick
        tick += 1
        stack = [(start, None, None, iter(adjacency[start]))]
        while stack:
            node, parent, incoming, edges = stack[-1]
            try:
                other, eid = next(edges)
            except StopIteration:
                stack.pop()
                if parent is not None:
                    low[parent] = min(low[parent], low[node])
                    if low[node] > discovery[parent]:
                        result.add(incoming)
                continue
            if eid == incoming:
                continue
            if other in discovery:
                low[node] = min(low[node], discovery[other])
            else:
                discovery[other] = low[other] = tick
                tick += 1
                stack.append((other, node, eid, iter(adjacency[other])))
    return result


def main(args):
    inventory = json.loads(Path(args.inventory).read_text())
    unchanged(inventory['graph'])
    out = new_output_dir(args.outdir)
    lengths, reference, edges = {}, defaultdict(list), set()
    ambiguity_bases, ambiguity_nodes = 0, 0
    for row in rows(args.gfa):
        if row[0] == 'S':
            node, seq = row[1:3]
            if node in lengths:
                raise ValueError(f'Duplicate seed node: {node}')
            if seq == '*' or not seq:
                raise ValueError(f'Missing seed sequence: {node}')
            _, masked = cactus_bases(seq.encode(), f'seed node {node}')
            ambiguity_bases += masked
            ambiguity_nodes += int(masked > 0)
            lengths[node] = len(seq)
            tags = {x[:2]: x[5:] for x in row[3:]}
            if 'LN' in tags and int(tags['LN']) != len(seq):
                raise ValueError(f'Seed LN mismatch: {node}')
            if tags.get('SR') == '0':
                if not tags.get('SN', '').startswith(args.reference_prefix) or 'SO' not in tags:
                    raise ValueError(f'Unrecognized rank-0 reference tags: {node}')
                start = int(tags['SO'])
                reference[tags['SN']].append((start, start + len(seq), node))
        elif row[0] == 'L':
            if row[5] != '0M' or row[2] not in '+-' or row[4] not in '+-':
                raise ValueError('Only blunt-ended oriented GFA1 links are supported')
            edges.add(link_key(row[1], row[2], row[3], row[4]))
        elif row[0] not in ('H', 'P', 'W'):
            raise ValueError(f'Unsupported GFA record: {row[0]}')
    if not lengths or not reference:
        raise ValueError('Missing seed nodes or rank-0 reference tags')
    if set(reference) != set(inventory['reference_lengths']):
        raise ValueError('Rank-0 reference contigs differ from the frozen crosswalk')
    adjacency = {node: [] for node in lengths}
    edges = sorted(edges)
    edge_ids = {}
    for eid, key in enumerate(edges):
        (a, _), (b, _) = key
        if a not in lengths or b not in lengths:
            raise ValueError(f'Link refers to an undefined segment: {a}, {b}')
        adjacency[a].append((b, eid))
        adjacency[b].append((a, eid))
        edge_ids[key] = eid
    safe = bridges(adjacency)
    cuts, boundaries, ref_nodes, reference_summary = set(), [], {}, []
    for name, intervals in sorted(reference.items()):
        intervals.sort()
        if intervals[-1][1] != inventory['reference_lengths'][name]:
            raise ValueError(f'Reference graph span differs from original FASTA length: {name}')
        if intervals[0][0] != 0 or any(a[1] != b[0] for a, b in zip(intervals, intervals[1:])):
            raise ValueError(f'Reference tags have gaps/overlaps: {name}')
        candidates = []
        for start, end, node in intervals:
            ref_nodes[node] = (name, start, end)
        for left, right in zip(intervals, intervals[1:]):
            eid = edge_ids.get(link_key(left[2], '+', right[2], '+'))
            if eid is None:
                raise ValueError(f'Missing forward reference adjacency: {name}:{left[1]}')
            if eid in safe:
                candidates.append((left[1], eid))
        positions = [x[0] for x in candidates]
        reference_summary.append(dict(reference=name, length=intervals[-1][1],
                                      safe_candidate_count=len(candidates)))
        for target in range(args.chunk_bp, intervals[-1][1], args.chunk_bp):
            index = bisect.bisect_left(positions, target)
            nearby = candidates[max(0, index-1):index+1]
            best = min(nearby, key=lambda x: (abs(x[0]-target), x[0])) if nearby else None
            accepted = best is not None and abs(best[0]-target) <= args.slack_bp
            if accepted:
                cuts.add(best[1])
            boundaries.append(dict(reference=name, target=target,
                                   nearest_safe=best[0] if best else '.',
                                   distance=abs(best[0]-target) if best else '.',
                                   selected=int(accepted),
                                   reason='bridge' if accepted else 'no_bridge_within_slack'))
    owner, components = {}, []
    for start in sorted(lengths):
        if start in owner:
            continue
        cid = f'{len(components):06d}'
        stack, nodes = [start], []
        owner[start] = cid
        while stack:
            node = stack.pop()
            nodes.append(node)
            for other, eid in adjacency[node]:
                if eid not in cuts and other not in owner:
                    owner[other] = cid
                    stack.append(other)
        components.append((cid, nodes))
    summaries, problems = [], []
    for cid, nodes in components:
        refs = defaultdict(list)
        for node in nodes:
            if node in ref_nodes:
                name, start, end = ref_nodes[node]
                refs[name].append((start, end))
        spans = {name: max(x[1] for x in iv)-min(x[0] for x in iv) for name, iv in refs.items()}
        bp = sum(lengths[n] for n in nodes)
        issues = []
        if len(refs) > 1:
            issues.append('multiple_references')
        if any(span > args.max_span_bp for span in spans.values()):
            issues.append('reference_span_limit')
        if bp > args.max_seed_bp:
            issues.append('seed_bases_limit')
        flat = [x for iv in refs.values() for x in iv]
        summaries.append(dict(chunk=cid, reference=';'.join(sorted(refs)) or '.',
                              start=min((x[0] for x in flat), default=0),
                              end=max((x[1] for x in flat), default=0),
                              max_reference_span=max(spans.values(), default=0),
                              seed_bp=bp, nodes=len(nodes), status=','.join(issues) or 'PASS'))
        if issues:
            problems.append(dict(chunk=cid, issues=issues, spans=spans, seed_bp=bp))
    write_table(out/'chunks.tsv', ['chunk','reference','start','end','max_reference_span','seed_bp','nodes','status'], summaries)
    write_table(out/'boundaries.tsv', ['reference','target','nearest_safe','distance','selected','reason'], boundaries)
    write_table(out/'references.tsv', ['reference','length','safe_candidate_count'], reference_summary)
    write_table(out/'nodes.tsv', ['node','chunk','length'],
                (dict(node=n, chunk=owner[n], length=lengths[n]) for n in sorted(lengths)))
    write_table(out/'cut_links.tsv', ['left','left_end','right','right_end'],
                (dict(left=edges[i][0][0], left_end='right' if edges[i][0][1] else 'left',
                      right=edges[i][1][0], right_end='right' if edges[i][1][1] else 'left') for i in sorted(cuts)))
    # A diagnostic report is a successful preparation product even when infeasible.
    # The separate readiness gate refuses downstream sequence extraction/construction.
    write_json(out/'report.json', dict(feasible=not problems, problems=problems,
               nodes=len(lengths), physical_links=len(edges), bridges=len(safe), cuts=len(cuts),
               chunks=len(components), off_reference_chunks=sum(r['reference']=='.' for r in summaries),
               largest_reference_span=max(r['max_reference_span'] for r in summaries),
               largest_seed_bp=max(r['seed_bp'] for r in summaries),
               seed_base_normalization=CACTUS_BASE_POLICY,
               seed_ambiguity_bases=ambiguity_bases, seed_ambiguity_nodes=ambiguity_nodes,
               method='conservative_existing_node_bridges', settings=vars(args),
               remedy='If infeasible: review localized CIGAR lift-over or affected-input remapping; no automatic node splitting.'))
    if not problems:
        writers = Writers()
        try:
            for row in rows(args.gfa):
                if row[0] == 'S':
                    seq, _ = cactus_bases(row[2].encode(), f'seed node {row[1]}')
                    writers.get(out/'seed'/f'{owner[row[1]]}.fa').write(f'>_MINIGRAPH_#0#{row[1]}\n{seq.decode()}\n')
        finally:
            writers.close()
    unchanged(inventory['graph'])
    print(f'Partition {"PASS" if not problems else "INFEASIBLE"}: '
          f'{len(components)} chunks, {len(cuts)} cuts, {len(problems)} oversized/mixed components; '
          f'{ambiguity_bases} ambiguity bases in {ambiguity_nodes} seed nodes. '
          f'See {out / "report.json"}', flush=True)
    if problems:
        print('Diagnostic report written; readiness remains blocked and no seed FASTAs were emitted.',
              flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('gfa','inventory','outdir','reference-prefix'):
        p.add_argument('--'+key, required=True)
    for key in ('chunk-bp','slack-bp','max-span-bp','max-seed-bp'):
        p.add_argument('--'+key, type=int, required=True)
    a = p.parse_args()
    if not 0 < a.chunk_bp <= a.max_span_bp or a.slack_bp < 0 or a.max_seed_bp <= 0:
        p.error('Invalid partition limits')
    main(a)
