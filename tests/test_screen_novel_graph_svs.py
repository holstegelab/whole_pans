#!/usr/bin/env python3

import csv
import gzip
import importlib.util
import io
import subprocess
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

    def test_path_projection_accepts_minigraph_stable_interval_labels(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            index = Path(tmpdir) / "segments.tsv"
            index.write_text(
                "segment_id\tlength\tsn\tso\tsr\n"
                "s1\t100\tchr1\t1000\t0\n"
                "s2\t50\tdonor.ctg\t2000\t2\n"
            )
            segments = SCREEN.load_segment_index(index)

            forward = SCREEN.locate_graph_position(
                ">chr1:1000-1100>donor.ctg:2000-2050", 25, segments
            )
            reverse = SCREEN.locate_graph_position(
                ">chr1:1000-1100<donor.ctg:2000-2050", 110, segments
            )

            self.assertEqual(forward["segment"], "chr1:1000-1100")
            self.assertEqual(forward["stable_position"], 1025)
            self.assertEqual(forward["sr"], 0)
            self.assertEqual(reverse["stable_position"], 2039)
            self.assertEqual(reverse["sr"], 2)

    def test_path_projection_accepts_bare_rank_zero_stable_paths(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            index = Path(tmpdir) / "segments.tsv"
            index.write_text(
                "segment_id\tlength\tsn\tso\tsr\n"
                "s1\t50\tchr5\t0\t0\n"
                "s2\t50\tchr5\t50\t0\n"
                "s3\t100\tchr5\t200\t1\n"
            )
            segments = SCREEN.load_segment_index(index)

            location = SCREEN.locate_graph_position("chr5", 75, segments)

            self.assertEqual(location["segment"], "chr5")
            self.assertEqual(location["offset"], 75)
            self.assertEqual(location["stable_position"], 75)
            self.assertEqual(location["sr"], 0)
            self.assertEqual(
                SCREEN.locate_graph_position("chr5", 125, segments),
                SCREEN.empty_graph_location(),
            )


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


class MinigraphTests(unittest.TestCase):
    def test_run_minigraph_streams_stdout_through_gzip_writer(self):
        gaf = (
            b"contig1\t6000\t0\t6000\t+\t>s1\t10000\t0\t6000\t"
            b"6000\t6000\t60\ttp:A:P\tcg:Z:6000M\n"
        )
        process = mock.MagicMock()
        process.stdout = io.BytesIO(gaf)
        process.returncode = 0
        process.__enter__.return_value = process
        process.__exit__.return_value = False

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            output = tmp / "alignment.gaf.gz"
            log = tmp / "minigraph.log"

            with mock.patch.object(SCREEN.subprocess, "Popen", return_value=process) as popen:
                SCREEN.run_minigraph(
                    "graph.gfa",
                    "assembly.fa",
                    output,
                    log,
                    threads=4,
                    min_chain_score=5000,
                    secondary=5,
                    extra="",
                )

            self.assertEqual(output.read_bytes()[:2], b"\x1f\x8b")
            with gzip.open(output, "rb") as handle:
                self.assertEqual(handle.read(), gaf)
            self.assertIs(popen.call_args.kwargs["stdout"], subprocess.PIPE)


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


class ConcatenateTablesTests(unittest.TestCase):
    def test_streams_fields_larger_than_csv_limit(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            first = tmp / "first.tsv.gz"
            second = tmp / "second.tsv.gz"
            output = tmp / "combined.tsv.gz"
            large_cigar = "1M" * 200_000

            with gzip.open(first, "wt") as handle:
                handle.write("query_name\tcigar\n")
                handle.write(f"contig1\t{large_cigar}\n")
            with gzip.open(second, "wt") as handle:
                handle.write("query_name\tcigar\n")
                handle.write("contig2\t100M\n")

            SCREEN.concatenate_tables(
                [first, second], output, ["query_name", "cigar"]
            )

            with gzip.open(output, "rt") as handle:
                self.assertEqual(handle.readline(), "query_name\tcigar\n")
                self.assertEqual(handle.readline(), f"contig1\t{large_cigar}\n")
                self.assertEqual(handle.readline(), "contig2\t100M\n")
                self.assertEqual(handle.read(), "")

    def test_rejects_an_unexpected_header(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            source = tmp / "source.tsv"
            output = tmp / "combined.tsv"
            source.write_text("query_name\twrong_column\ncontig1\t100M\n")

            with self.assertRaisesRegex(ValueError, "Unexpected columns"):
                SCREEN.concatenate_tables(
                    [source], output, ["query_name", "cigar"]
                )


class ReannotateCandidateTests(unittest.TestCase):
    def candidate_table(self, path, graph_path):
        row = {field: "" for field in SCREEN.EVENT_FIELDS}
        row.update(
            {
                "event_id": "GSV_1",
                "assembly_id": "sample.hifi.hifiasm.bp.hap1.p_ctg",
                "sample_id": "sample",
                "haplotype": "hap1",
                "graph_member": "false",
                "graph_path": graph_path,
                "path_position_0": "25",
                "svtype": "INS",
                "svlen": "60",
                "event_size_bp": "60",
                "confidence_tier": "HIGH_CONFIDENCE",
            }
        )
        SCREEN.write_tsv(path, [row], SCREEN.EVENT_FIELDS)

    def test_reannotates_saved_candidates_and_writes_qc(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            index = tmp / "segments.tsv"
            index.write_text(
                "segment_id\tlength\tsn\tso\tsr\n"
                "s1\t100\tchr1\t1000\t0\n"
            )
            source = tmp / "source.tsv.gz"
            output = tmp / "annotated.tsv.gz"
            qc = tmp / "qc.tsv"
            self.candidate_table(source, ">chr1:1000-1100")

            SCREEN.reannotate_candidates(source, index, output, qc)

            with gzip.open(output, "rt") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["graph_segment"], "chr1:1000-1100")
            self.assertEqual(row["stable_source"], "chr1")
            self.assertEqual(row["stable_position_0"], "1025")
            self.assertEqual(row["source_rank"], "0")
            self.assertIn("status\tPASS", qc.read_text())

    def test_rejects_excess_unresolved_coordinates(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            index = tmp / "segments.tsv"
            index.write_text(
                "segment_id\tlength\tsn\tso\tsr\n"
                "s1\t100\tchr1\t1000\t0\n"
            )
            source = tmp / "source.tsv.gz"
            output = tmp / "annotated.tsv.gz"
            qc = tmp / "qc.tsv"
            self.candidate_table(source, ">unknown_path")

            with self.assertRaisesRegex(ValueError, "Unresolved graph coordinates"):
                SCREEN.reannotate_candidates(
                    source,
                    index,
                    output,
                    qc,
                    max_unresolved_fraction=0.0,
                )
            self.assertFalse(output.exists())
            self.assertIn("status\tFAIL", qc.read_text())


if __name__ == "__main__":
    unittest.main()
