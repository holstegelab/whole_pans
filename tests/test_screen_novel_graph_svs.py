#!/usr/bin/env python3

import csv
import gzip
import importlib.util
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "workflow" / "scripts" / "screen_novel_graph_svs.py"
SPEC = importlib.util.spec_from_file_location("screen_novel_graph_svs", SCRIPT)
SCREEN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SCREEN)


def task(**updates):
    row = {
        "assembly_id": "sample.hifi.hifiasm.bp.hap1.p_ctg",
        "sample_id": "sample",
        "haplotype": "hap1",
        "graph_member": "false",
        "rescue_tier": "best_rescue",
        "path": "assembly.fa",
    }
    row.update(updates)
    return row


def alignment(cigar, **updates):
    row = {
        "mapping_pass": "primary",
        "query_name": "contig1",
        "query_length": 10000,
        "query_start": 0,
        "query_end": 10000,
        "strand": "+",
        "graph_path": ">s1",
        "path_length": 10000,
        "path_start": 0,
        "path_end": 10000,
        "matches": 9940,
        "block_length": 10000,
        "identity": 0.994,
        "mapq": 60,
        "primary": True,
        "cigar": cigar,
    }
    row.update(updates)
    return row


class GraphIndexTests(unittest.TestCase):
    def test_index_rgfa_reads_sequence_and_ln_lengths(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            gfa = tmp / "graph.gfa"
            output = tmp / "segments.tsv.gz"
            gfa.write_text(
                "H\tVN:Z:1.0\n"
                "S\ts1\tACGT\tSN:Z:CHM13.chr1\tSO:i:100\tSR:i:0\n"
                "S\ts2\t*\tLN:i:7\tSN:Z:sample.ctg\tSO:i:20\tSR:i:2\n"
            )

            SCREEN.index_graph(gfa, output)

            with gzip.open(output, "rt") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0], {"segment_id": "s1", "length": "4", "sn": "CHM13.chr1", "so": "100", "sr": "0"})
            self.assertEqual(rows[1]["length"], "7")

    def test_path_projection_honors_reverse_segment_orientation(self):
        segments = {
            "s1": {"length": 100, "sn": "CHM13.chr1", "so": 1000, "sr": 0},
            "s2": {"length": 50, "sn": "CHM13.chr1", "so": 2000, "sr": 0},
        }
        forward = SCREEN.locate_graph_position(">s1<s2", 25, segments)
        reverse = SCREEN.locate_graph_position(">s1<s2", 110, segments)
        self.assertEqual(forward["stable_position"], 1025)
        self.assertEqual(reverse["segment"], "s2")
        self.assertEqual(reverse["offset"], 39)
        self.assertEqual(reverse["stable_position"], 2039)


class CigarTests(unittest.TestCase):
    def test_internal_insertion_has_two_sided_anchors(self):
        events = SCREEN.extract_cigar_events("2500M60I3000M")
        self.assertEqual(len(events), 1)
        self.assertEqual(events[0]["operation"], "I")
        self.assertEqual(events[0]["query_position"], 2500)
        self.assertEqual(events[0]["target_position"], 2500)
        self.assertEqual(events[0]["left_anchor"], 2500)
        self.assertEqual(events[0]["right_anchor"], 3000)

    def test_adjacent_subthreshold_indels_emit_combined_cluster(self):
        events = SCREEN.extract_cigar_events("2100M30I10M25D2100M", cluster_gap=20)
        atoms = [event for event in events if event["kind"] == "atomic"]
        clusters = [event for event in events if event["kind"] == "cluster"]
        self.assertEqual([event["length"] for event in atoms], [30])
        # Sub-30 bp operations are not retained alone, but do contribute to a
        # nearby cluster that crosses the 50 bp SV definition.
        self.assertEqual(clusters[0]["length"], 55)
        self.assertEqual(clusters[0]["operation"], "ID")

        events = SCREEN.extract_cigar_events("2100M30I10M30D2100M", cluster_gap=20)
        cluster = [event for event in events if event["kind"] == "cluster"][0]
        self.assertEqual(cluster["operation"], "ID")
        self.assertEqual(cluster["length"], 60)
        self.assertEqual(cluster["members"], 2)
        self.assertGreaterEqual(cluster["left_anchor"], 2000)
        self.assertGreaterEqual(cluster["right_anchor"], 2000)

    def test_terminal_event_stays_in_review(self):
        records = [alignment("60I6000M", block_length=6060, matches=6000)]
        segments = {"s1": {"length": 10000, "sn": "CHM13.chr1", "so": 0, "sr": 0}}
        thresholds = {
            "raw_min_size": 30,
            "min_sv_size": 50,
            "cluster_gap": 50,
            "min_alignment": 5000,
            "min_identity": 0.90,
            "min_mapq": 5,
            "min_anchor": 2000,
        }
        event = SCREEN.events_from_alignments(records, task(), segments, thresholds)[0]
        self.assertEqual(event["confidence_tier"], "REVIEW")
        self.assertIn("TERMINAL_EVENT_WITHOUT_SECOND_BREAKPOINT", event["filter_reasons"])

    def test_secondary_event_is_not_high_confidence(self):
        records = [alignment("2500M60D3000M", primary=False)]
        segments = {"s1": {"length": 10000, "sn": "CHM13.chr1", "so": 0, "sr": 0}}
        thresholds = {
            "raw_min_size": 30,
            "min_sv_size": 50,
            "cluster_gap": 50,
            "min_alignment": 5000,
            "min_identity": 0.90,
            "min_mapq": 5,
            "min_anchor": 2000,
        }
        event = SCREEN.events_from_alignments(records, task(), segments, thresholds)[0]
        self.assertEqual(event["confidence_tier"], "REVIEW")
        self.assertIn("SECONDARY_ALIGNMENT", event["filter_reasons"])


class SplitAlignmentTests(unittest.TestCase):
    def test_opposite_orientation_primary_alignments_are_reviewed(self):
        records = [
            alignment("3000M", query_start=0, query_end=3000, graph_path=">s1"),
            alignment(
                "3000M",
                query_start=4000,
                query_end=7000,
                graph_path=">s2",
                strand="-",
            ),
        ]
        rows = SCREEN.classify_split_alignments(records, task())
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["classification"], "POTENTIAL_INVERSION_OR_TRANSLOCATION")
        self.assertIn("inconsistent orientation", rows[0]["notes"])


class TaskTests(unittest.TestCase):
    def test_task_table_is_prioritized_and_batched_deterministically(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            manifest = tmp / "manifest.tsv"
            graph = tmp / "ordered.tsv"
            output = tmp / "tasks.tsv"
            manifest.write_text(
                "assembly_id\tcleaned_fasta_path\tsample_id\thaplotype\trescue_tier\tgraph_context\n"
                "b.hifi.hifiasm.bp.hap1.p_ctg\t/b.fa\tb\thap1\tnot_recommended\toutside_graph\n"
                "a.hifi.hifiasm.bp.hap1.p_ctg\t/a.fa\ta\thap1\tbest_rescue\trestore_missing_mate_for_graph_sample\n"
                "c.hifi.hifiasm.bp.hap1.p_ctg\t/c.fa\tc\thap1\tbest_rescue\toutside_graph\n"
            )
            graph.write_text("order\trole\tassembly_id\tpath\n")

            SCREEN.build_tasks(manifest, graph, output, batch_size=2)

            with output.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["assembly_id"].split(".")[0] for row in rows], ["a", "c", "b"])
            self.assertEqual([row["task_index"] for row in rows], ["1", "1", "2"])
            self.assertEqual([row["task_id"] for row in rows], ["0001", "0001", "0002"])

    def test_run_task_writes_one_atomic_declared_directory(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            assembly = tmp / "assembly.fa"
            assembly.write_text(">contig1\n" + "A" * 6000 + "\n")
            tasks = tmp / "tasks.tsv"
            tasks.write_text(
                "task_id\ttask_index\tassembly_index\tassembly_id\tsample_id\t"
                "haplotype\tpath\tgraph_member\trescue_tier\tgraph_context\n"
                f"0001\t1\t1\tsample.hifi.hifiasm.bp.hap1.p_ctg\tsample\t"
                f"hap1\t{assembly}\tfalse\tbest_rescue\toutside_graph\n"
            )
            segments = tmp / "segments.tsv"
            segments.write_text(
                "segment_id\tlength\tsn\tso\tsr\n"
                "s1\t10000\tCHM13.chr1\t0\t0\n"
            )
            output = tmp / "task-output" / "0001"
            args = SimpleNamespace(
                tasks=str(tasks),
                output_dir=str(output),
                task_id="0001",
                task_index=None,
                segment_index=str(segments),
                graph="graph.gfa",
                graph_assemblies="ordered.tsv",
                manifest="manifest.tsv",
                batch_size=1,
                threads=1,
                primary_min_chain_score=5000,
                primary_secondary=5,
                minigraph_extra="",
                rescue_callable_fraction=0.85,
                sensitivity_min_chain_score=1000,
                sensitivity_secondary=20,
                sensitivity_minigraph_extra="",
                min_sv_size=50,
                raw_min_size=30,
                adjacent_indel_gap=50,
                min_alignment=5000,
                min_anchor=2000,
                min_mapq=5,
                min_identity=0.90,
                completion_marker=None,
            )

            def fake_minigraph(graph, assembly_path, gaf, log, *unused):
                Path(log).parent.mkdir(parents=True, exist_ok=True)
                Path(log).write_text("mock minigraph\n")
                Path(gaf).parent.mkdir(parents=True, exist_ok=True)
                with gzip.open(gaf, "wt") as handle:
                    handle.write(
                        "contig1\t6000\t0\t6000\t+\t>s1\t10000\t0\t5940\t"
                        "5940\t6000\t60\ttp:A:P\tcg:Z:2500M60I3440M\n"
                    )

            with mock.patch.object(SCREEN, "run_minigraph", side_effect=fake_minigraph):
                SCREEN.run_task(args)

            self.assertTrue((output / ".complete").is_file())
            self.assertTrue((output / "task_outputs.tsv").is_file())
            self.assertTrue(
                (output / "candidates" / "sample.hifi.hifiasm.bp.hap1.p_ctg.residual_svs.tsv.gz").is_file()
            )


if __name__ == "__main__":
    unittest.main()
