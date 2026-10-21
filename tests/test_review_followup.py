"""Small regression fixtures; run on the server, not the mounted checkout."""

import importlib.util
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


SCRIPTS = Path(__file__).absolute().parents[1] / "workflow" / "scripts"


def load(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


followup = load("review_followup")
qc = load("qc_local_sv_graph_v2")


class ReviewFollowupTests(unittest.TestCase):
    def test_binary_coordinates_and_large_integers(self):
        self.assertEqual(qc.int_value(b"123"), 123)
        self.assertEqual(qc.int_value("9007199254740993"), 9007199254740993)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            frozen, rescued = root / "v1.tsv", root / "v2.tsv"
            frozen.write_bytes(b"")
            rescued.write_bytes(b"chr1\t123\t124\t3\t2\t0\t1\t51\t0\t0\t0\ta,b,c\tA\tACGT\n")
            delta, removed = qc.bubble_delta(frozen, rescued)
            self.assertEqual((delta[0]["start_0"], delta[0]["end_0"], delta[0]["bubble_length_bp"]),
                             (123, 124, 51))
            self.assertEqual(removed, 0)
            rescued.write_bytes(b"chr1\tbroken\t124\t3\t2\t0\t1\t51\t0\t0\t0\ta,b,c\tA\tACGT\n")
            with self.assertRaises(ValueError):
                qc.bubble_delta(frozen, rescued)

    def test_grouped_donor_order_and_reference_allele_orientation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest, graph, bubbles, table = [root / name for name in
                                              ("manifest.tsv", "graph.gfa", "bubbles.tsv", "ranked.tsv")]
            followup.write_table(manifest, [
                dict(sample="A", assembly_id="A.hap1.clean", path="A1.fa"),
                dict(sample="B", assembly_id="B.hap1.clean", path="B1.fa"),
                dict(sample="A", assembly_id="A.hap2.clean", path="A2.fa"),
            ])
            graph.write_text("S\tn2\tCTTA\tSN:Z:A.hap2.contig1\tSO:i:10\tSR:i:11\n")
            bubbles.write_text("CHM13#0#chr1\t1000\t1000\t3\t2\t0\t0\t4\t0\t0\t0\tn1,n2,n3\t*\tTAAG\n")
            followup.write_table(table, [dict(candidate_id="INS1", inserted_node_id="n2",
                source_rank="11", source_name="A.hap2.contig1", source_offset_0based="10",
                introducing_sample="B", introducing_assembly="B.hap1.clean",
                insertion_size_bp="4", reference_path="CHM13#0#chr1", position_0based="1000")])
            args = SimpleNamespace(manifest=manifest, graph=graph, bubbles=bubbles,
                                   tables=[table], first_rank=10, output=root / "out")
            followup.correct_provenance(args)
            corrected = list(followup.rows(root / "out" / "all_candidates.tsv"))[0]
            self.assertEqual(corrected["introducing_sample"], "A")
            self.assertEqual(corrected["introducing_assembly"], "A.hap2.clean")
            self.assertEqual(corrected["previous_introducing_sample"], "B")
            self.assertEqual(followup.read_fasta(root / "out" / "insertion_alleles.fa"), {"INS1": "TAAG"})
            args.first_rank = 11
            with self.assertRaisesRegex(ValueError, "disagrees with source contig"):
                followup.correct_provenance(args)

    def test_reference_windows_cross_lines_and_reject_missing_contigs(self):
        with tempfile.TemporaryDirectory() as directory:
            fasta = Path(directory) / "reference.fa"
            fasta.write_text(">chr1\nAACCGG\nTTAACC\n>chr2\nAAAA\n")
            self.assertEqual(followup.fetch_windows(fasta, {"x": ("chr1", 4, 9)}), {"x": "GGTTA"})
            with self.assertRaisesRegex(ValueError, "missing/truncated"):
                followup.fetch_windows(fasta, {"x": ("1", 4, 9)})

    def test_empty_context_is_not_repeat_free(self):
        candidates = [dict(candidate_id="x", chromosome="chr1", position_0based="100")]
        self.assertEqual(followup.bed_context([], candidates, 10), {"x": "NOT_ASSESSED"})

    def test_freeze_refuses_missing_sample_graph_pairs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            evidence, inventory = root / "evidence.tsv", root / "inventory.tsv"
            followup.write_table(evidence, [dict(sample_id="A", target="g1")])
            followup.write_table(inventory, [dict(target="g1", graph_checksum="abc")])
            args = SimpleNamespace(evidence=evidence, inventory=inventory,
                                   samples=["A", "B"], graphs=["g1"])
            with self.assertRaisesRegex(ValueError, "incomplete/duplicated"):
                followup.freeze_tmem(args)


if __name__ == "__main__":
    unittest.main()
