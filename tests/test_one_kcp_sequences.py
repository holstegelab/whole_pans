"""Sequence regressions: run on the server, not the mounted checkout."""
import hashlib
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

SCRIPTS = Path(__file__).absolute().parents[1] / 'workflow' / 'scripts'
sys.path.insert(0, str(SCRIPTS))
import one_kcp_common as common
import partition_1kcp_graph as partition
import prepare_1kcp_inputs as prepare


class OneKcpSequenceTests(unittest.TestCase):
    def test_cactus_masking_preserves_length(self):
        for sequence, expected, count in (
            (b'WGTCCATTCG', b'NGTCCATTCG', 1),
            (b'acgtnACGTN', b'ACGTNACGTN', 0),
            (b'UWSMKRYBDHVuwsmkrybdhv', b'N' * 22, 22),
        ):
            with self.subTest(sequence=sequence):
                self.assertEqual(common.cactus_bases(sequence, 'fixture'), (expected, count))
                self.assertEqual(len(sequence), len(expected))
        for sequence in (b'', b'*', b'AC-T', b'AC.T', b'AX', b'AC GT', b'\xff'):
            with self.subTest(sequence=sequence), self.assertRaisesRegex(ValueError, 'fixture'):
                common.cactus_bases(sequence, 'fixture')

    def partition_args(self, root, extra=''):
        graph = root / 'source.gfa'
        graph.write_text(
            'S\ts1\tACGT\tLN:i:4\tSN:Z:id=CHM13|chr1\tSO:i:0\tSR:i:0\n'
            'S\ts318398\tWGTCCATTCG\tLN:i:10\tSN:Z:id=GRCh38|chr22\tSO:i:16320468\tSR:i:1\n'
            'L\ts1\t+\ts318398\t+\t0M\n' + extra)
        inventory = root / 'inventory.json'
        common.write_json(inventory, dict(graph=common.identity(graph),
                                         reference_lengths={'id=CHM13|chr1': 4}))
        return SimpleNamespace(gfa=str(graph), inventory=str(inventory),
                               outdir=str(root / 'partition'), reference_prefix='id=CHM13|',
                               chunk_bp=4, slack_bp=1, max_span_bp=8, max_seed_bp=100)

    def test_ambiguous_seed_keeps_source_and_node_lengths(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            args = self.partition_args(root)
            source = Path(args.gfa).read_bytes()
            with redirect_stdout(io.StringIO()):
                partition.main(args)
            self.assertEqual(Path(args.gfa).read_bytes(), source)
            nodes = {r['node']: r['length'] for r in common.table(root / 'partition/nodes.tsv')}
            self.assertEqual(nodes, {'s1': '4', 's318398': '10'})
            self.assertEqual((root / 'partition/seed/000000.fa').read_text(),
                             '>_MINIGRAPH_#0#s1\nACGT\n>_MINIGRAPH_#0#s318398\nNGTCCATTCG\n')
            report = json.loads((root / 'partition/report.json').read_text())
            self.assertTrue(report['feasible'])
            self.assertEqual((report['seed_ambiguity_bases'], report['seed_ambiguity_nodes']), (1, 1))

    def test_infeasible_partition_still_blocks_seed_emission(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            args = self.partition_args(root)
            args.max_seed_bp = 5
            log = io.StringIO()
            with redirect_stdout(log):
                partition.main(args)
            self.assertIn('INFEASIBLE', log.getvalue())
            self.assertFalse(json.loads((root / 'partition/report.json').read_text())['feasible'])
            self.assertFalse((root / 'partition/seed').exists())

    def test_bad_seed_records_still_fail(self):
        for extra, error in (
            ('S\ts318398\tACGT\n', 'Duplicate seed node: s318398'),
            ('S\tempty\t*\n', 'Missing seed sequence: empty'),
            ('S\tbad\tAC-T\n', 'Unsupported sequence symbols'),
        ):
            with self.subTest(extra=extra), tempfile.TemporaryDirectory() as directory:
                args = self.partition_args(Path(directory), extra)
                with self.assertRaisesRegex(ValueError, error):
                    partition.main(args)

    def test_extraction_masks_before_interval_hash_but_checks_original_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / 'original.fa'
            raw = b'>chr22\nwgtc\nCATTCG\n'
            fasta.write_bytes(raw)
            sample = dict(source_path=str(fasta), fasta_identity=common.identity(fasta),
                          sha256_uncompressed=hashlib.sha256(raw).hexdigest())
            contigs = [dict(original_header='chr22', length_bp='10')]
            plan = root / 'plan'
            plan.mkdir()
            interval = dict(query='id=GRCh38|chr22', original='chr22',
                            name='GRCh38#0#chr22[0-10]', chunk='000000', start=0, end=10,
                            sample='GRCh38', haplotype=0)
            common.write_table(plan / 'intervals.tsv', list(interval), [interval])
            common.write_json(plan / 'stats.json', {})
            args = SimpleNamespace(sample='GRCh38', plan=str(plan), outdir=str(root / 'sequences'))
            with patch.object(prepare, 'sample_data', return_value=(sample, contigs)):
                prepare.extract(args)
                sample['sha256_uncompressed'] = '0' * 64
                args.outdir = str(root / 'bad-checksum')
                with self.assertRaisesRegex(ValueError, 'uncompressed SHA256'):
                    prepare.extract(args)
            output = root / 'sequences'
            self.assertEqual((output / '000000.fa').read_text(),
                             '>GRCh38#0#chr22[0-10]\nNGTC\nCATTCG\n')
            intervals = common.table(output / 'intervals.tsv')
            self.assertEqual(intervals[0]['sha256'], hashlib.sha256(b'NGTCCATTCG').hexdigest())
            stats = json.loads((output / 'stats.json').read_text())
            self.assertEqual(stats['ambiguity_bases_masked'], 1)
            self.assertEqual(stats['lowercase_bases_uppercased'], 4)
            self.assertEqual(stats['fasta_sha256'], hashlib.sha256(raw).hexdigest())
            self.assertEqual(fasta.read_bytes(), raw)


if __name__ == '__main__':
    unittest.main()
