#!/usr/bin/env python3
"""Small synthetic tests for the SV-rGFA short-read pilot helpers."""

from __future__ import annotations

import argparse
import csv
import gzip
import subprocess
import sys
import tempfile
import tracemalloc
import unittest
from collections import Counter
from pathlib import Path

import yaml


SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import sv_gaf_metrics as gaf_metrics  # noqa: E402
import sv_short_read_pilot as pilot  # noqa: E402
import sv_short_read_report as report  # noqa: E402
import tmem_whole_graph_target as target  # noqa: E402


def write_fastq(path: Path, names: list[str], mate: int) -> None:
    handle = gzip.open(path, "wt") if path.suffix == ".gz" else path.open("w")
    with handle:
        for name in names:
            handle.write(f"@{name}/{mate}\nACGT\n+\nIIII\n")


def read_one_tsv(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise AssertionError(f"expected one row in {path}, observed {len(rows)}")
    return rows[0]


def gaf_line(
    name: str, path: str, mapq: int, *, secondary: bool = False,
    query_length: int = 100, query_start: int = 0, query_end: int = 100,
) -> str:
    tags = "\ttp:A:S" if secondary else "\ttp:A:P"
    block = query_end - query_start
    return (
        f"{name}\t{query_length}\t{query_start}\t{query_end}\t+\t{path}\t"
        f"1000\t0\t{block}\t{block}\t{block}\t{mapq}{tags}\n"
    )


class FastqTests(unittest.TestCase):
    def test_plain_and_gzip_lanes_and_pair_names(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            r1a, r2a = root / "L1_R1.fastq", root / "L1_R2.fastq"
            r1b, r2b = root / "L2_R1.fastq.gz", root / "L2_R2.fastq.gz"
            write_fastq(r1a, ["a", "b"], 1)
            write_fastq(r2a, ["a", "b"], 2)
            write_fastq(r1b, ["c"], 1)
            write_fastq(r2b, ["c"], 2)
            output = root / "qc.tsv"
            pilot.fastq_qc(
                "sample", [str(r1a), str(r1b)], [str(r2a), str(r2b)], output
            )
            row = read_one_tsv(output)
            self.assertEqual(row["lane_count"], "2")
            self.assertEqual(row["total_input_pairs"], "3")
            self.assertEqual(row["total_input_reads"], "6")

    def test_pair_name_mismatch_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            r1, r2 = root / "R1.fastq", root / "R2.fastq"
            write_fastq(r1, ["one"], 1)
            write_fastq(r2, ["different"], 2)
            with self.assertRaisesRegex(ValueError, "normalized names differ"):
                pilot.fastq_qc("sample", [str(r1)], [str(r2)], root / "qc.tsv")

    def test_deterministic_exact_subsample(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            r1, r2 = root / "R1.fastq", root / "R2.fastq.gz"
            names = [f"read{i}" for i in range(10)]
            write_fastq(r1, names, 1)
            write_fastq(r2, names, 2)
            qc = root / "qc.tsv"
            pilot.fastq_qc("sample", [str(r1)], [str(r2)], qc)

            observed = []
            for suffix in ("a", "b"):
                args = argparse.Namespace(
                    sample="sample", r1=[str(r1)], r2=[str(r2)], qc=qc,
                    max_pairs=4, seed=20260901,
                    output_r1=root / f"{suffix}.R1.fastq.gz",
                    output_r2=root / f"{suffix}.R2.fastq.gz",
                    stats=root / f"{suffix}.tsv",
                )
                pilot.subsample(args)
                with gzip.open(args.output_r1, "rt") as handle:
                    observed.append(handle.read())
                self.assertEqual(read_one_tsv(args.stats)["actual_pairs"], "4")
            self.assertEqual(observed[0], observed[1])


class SelectionTests(unittest.TestCase):
    def test_frozen_default_selection_is_six_plus_four(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = root / "config.yaml"
            manifest = root / "frozen.tsv"
            output = root / "selected.tsv"
            config.write_text(yaml.safe_dump({"sv_short_read_pilot": {
                "enabled": True,
                "sample_roles": ["tmem_candidate", "tmem_random_blinded"],
                "include_matched_controls": False,
                "expected_role_counts": {
                    "tmem_candidate": 6, "tmem_random_blinded": 4,
                    "tmem_matched_control": 6,
                },
            }}))
            source_rows = []
            for role, count in (
                ("tmem_candidate", 6), ("tmem_random_blinded", 4),
                ("tmem_matched_control", 6),
            ):
                source_rows.extend({
                    "sample_id": f"{role}_{index}", "role": role,
                    "r1": f"{role}_{index}.R1.fastq.gz",
                    "r2": f"{role}_{index}.R2.fastq.gz",
                } for index in range(count))
            pilot.write_tsv(manifest, source_rows, ["sample_id", "role", "r1", "r2"])
            pilot.select_samples(config, manifest, output)
            rows = pilot.read_tsv(output)
        self.assertEqual(len(rows), 10)
        self.assertEqual(
            Counter(row["role"] for row in rows),
            Counter({"tmem_candidate": 6, "tmem_random_blinded": 4}),
        )


class IndexReuseTests(unittest.TestCase):
    def test_changed_source_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "graph.gfa"
            source.write_text("S\t1\tACGT\n")
            products = {}
            for label, suffix in (
                ("gbz", ".giraffe.gbz"), ("distance", ".dist"),
                ("minimizer", ".shortread.withzip.min"),
                ("zipcodes", ".shortread.zipcodes"),
            ):
                path = root / f"graph{suffix}"
                path.write_bytes(label.encode())
                products[label] = path
            manifest = root / "manifest.tsv"
            row = {
                "graph_id": "graph", "source_sha256": pilot.sha256_file(source),
                "indexing_mode": "rgfa_artificial_path_cover",
            }
            for label, path in products.items():
                row[f"{label}_path"] = str(path)
                row[f"{label}_size_bytes"] = path.stat().st_size
            pilot.write_tsv(manifest, [row], list(row))
            args = argparse.Namespace(
                graph="graph", source=source, manifest=manifest,
                output=root / "validated",
            )
            pilot.validate_index(args)
            self.assertTrue(args.output.is_file())
            source.write_text("S\t1\tTGCA\n")
            with self.assertRaisesRegex(ValueError, "source checksum differs"):
                pilot.validate_index(args)


class GafTests(unittest.TestCase):
    def metadata(self) -> dict[str, object]:
        return {
            "sample_id": "S", "sample_role": "tmem_candidate", "batch": "B",
            "graph_id": "local_v1", "baseline_graph": "hprc_chm13_sv",
            "mapper": "vg_giraffe", "indexing_mode": "rgfa_artificial_path_cover",
            "run_mode": "full",
        }

    def callable_target(self) -> dict[str, object]:
        return {
            "present": "true", "callable": "callable_complete_graph_specific_target",
            "left": {(">1", ">2")}, "right": {(">2", ">3")},
            "reference": {(">1", ">4")}, "insertion_edges": {(">2", ">3")},
            "insertion_nodes": {"2"}, "failure": "", "ambiguity": "none",
        }

    def test_unmapped_mapq0_secondary_and_malformed_rows(self) -> None:
        lines = [
            gaf_line("zero/1", ">1>2>3", 0),
            gaf_line("hit/1", ">1>2>3", 60),
            gaf_line("secondary/1", ">1>2>3", 50, secondary=True),
            gaf_line("unmapped/1", "*", 0, query_end=0),
            "malformed\trow\n",
        ]
        metrics, evidence = gaf_metrics.summarize_stream(
            lines, 8, self.metadata(), self.callable_target(), 5, 1
        )
        self.assertEqual(metrics["gaf_records_observed"], 5)
        self.assertEqual(metrics["primary_mapped_reads"], 2)
        self.assertEqual(metrics["unmapped_reads"], 6)
        self.assertEqual(metrics["gaf_unmapped_records"], 1)
        self.assertEqual(metrics["mapq_0_fraction"], 0.5)
        self.assertEqual(metrics["secondary_multimapping_count"], 1)
        self.assertEqual(metrics["malformed_gaf_rows"], 1)
        self.assertEqual(evidence["left_alt_junction_distinct_primary_reads"], 1)
        self.assertEqual(evidence["right_alt_junction_distinct_primary_reads"], 1)
        self.assertEqual(evidence["both_junctions_meet_threshold"], "true")

    def test_absent_graph_target_is_not_reported_as_zero_support(self) -> None:
        absent = self.callable_target()
        absent.update({
            "present": "false", "callable": "not_callable_graph_target_absent",
            "failure": "no_alternate_query_complete_high_identity_both-flank_mapping",
        })
        _, evidence = gaf_metrics.summarize_stream(
            [gaf_line("read/1", ">1>2>3", 60)], 2,
            self.metadata(), absent, 5, 1,
        )
        self.assertEqual(evidence["evidence_status"], "not_callable_graph_target_absent")
        self.assertEqual(evidence["left_alt_junction_distinct_primary_reads"], "")
        self.assertEqual(evidence["both_junctions_meet_threshold"], "")

    def test_streaming_accumulators_have_bounded_memory(self) -> None:
        absent = self.callable_target()
        absent["callable"] = "not_callable_graph_target_absent"
        line = gaf_line("read", ">1>2>3", 60)
        tracemalloc.start()
        record_count = 25_000
        metrics, _ = gaf_metrics.summarize_stream(
            (line for _ in range(record_count)), record_count,
            self.metadata(), absent, 5, 1,
        )
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        self.assertEqual(metrics["primary_mapped_reads"], record_count)
        self.assertLess(peak, 10 * 1024 * 1024)


class GiraffeLogTests(unittest.TestCase):
    def test_fragment_learning_and_single_end_fallback_are_visible(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            log = Path(directory) / "giraffe.log"
            log.write_text("Using fragment length estimate: 351.5 +/- 42.25\n")
            learned = pilot.parse_giraffe_log(log)
            self.assertEqual(learned["paired_fragment_learning_succeeded"], "true")
            self.assertEqual(learned["single_end_fallback"], "false")
            self.assertEqual(learned["learned_fragment_mean"], "351.5")

            log.write_text(
                "warning: Finalizing fragment length distribution before reaching maximum sample size\n"
                "mean: 340, stdev: 55; mapped123 reads single ended\n"
            )
            fallback = pilot.parse_giraffe_log(log)
            self.assertEqual(fallback["paired_fragment_learning_succeeded"], "false")
            self.assertEqual(fallback["single_end_fallback"], "true")
            self.assertEqual(fallback["single_end_fallback_read_count"], 123)


class TargetAuditTests(unittest.TestCase):
    def query_rows(self) -> list[dict[str, str]]:
        common = {
            "query_length": "2100", "left_flank_bp": "1000",
            "right_flank_bp": "1000", "allele_start": "1000",
            "allele_end": "1100", "derivation_source": "synthetic",
        }
        return [
            {**common, "query_name": "REF", "allele_type": "reference"},
            {**common, "query_name": "ALT", "allele_type": "alternate"},
        ]

    def alignment(self, traversal: list[str]) -> dict[str, object]:
        return {
            "query_length": 2100, "query_start": 0, "query_end": 2100,
            "coverage": 1.0, "identity": 0.99, "mapq": 60,
            "traversal": traversal, "cigar": "2100=", "primary": True,
        }

    def test_unqualified_alternate_does_not_prove_graph_absence(self) -> None:
        result = target.derive_target(
            "hprc_chm13_sv", self.query_rows(),
            {"REF": [self.alignment([">1", ">2", ">5"])]},
            0.95, 0.90, 20, 900,
        )
        self.assertEqual(result["target_present"], "unknown")
        self.assertEqual(result["callable_status"], "not_callable_no_qualified_alternate")

    def test_orientation_aware_junction_derivation(self) -> None:
        result = target.derive_target(
            "local_v1", self.query_rows(), {
                "REF": [self.alignment([">1", ">2", ">5"])],
                "ALT": [self.alignment([">1", ">3", ">4", ">5"])],
            }, 0.95, 0.90, 20, 900,
        )
        self.assertEqual(result["callable_status"], "callable_complete_graph_specific_target")
        self.assertEqual(result["left_junction"], ">1,>3")
        self.assertEqual(result["right_junction"], ">4,>5")
        self.assertEqual(result["insertion_nodes"], "3,4")

    def test_nonunique_reference_mapping_is_ambiguous(self) -> None:
        result = target.derive_target(
            "local_v1", self.query_rows(), {
                "REF": [
                    self.alignment([">1", ">2", ">5"]),
                    self.alignment([">6", ">7", ">8"]),
                ],
                "ALT": [self.alignment([">1", ">3", ">4", ">5"])],
            }, 0.95, 0.90, 20, 900,
        )
        self.assertEqual(result["target_present"], "unknown")
        self.assertEqual(result["callable_status"], "not_callable_target_ambiguous")
        self.assertIn("multiple_complete", result["failure_reason"])


class ReportTests(unittest.TestCase):
    def test_report_has_required_lead_and_limitations(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            def table(name: str, row: dict[str, object]) -> Path:
                path = root / name
                pilot.write_tsv(path, [row], list(row))
                return path

            preflight = table("preflight.tsv", {
                "graph_id": "hprc_chm13_sv", "source_size_bytes": 1024,
                "p_path_lines": 0, "w_walk_lines": 0,
                "indexing_mode": "rgfa_artificial_path_cover",
                "graph_qc_status": "released_input",
            })
            selected = table("selected.tsv", {
                "sample_id": "S", "role": "tmem_candidate", "batch": "B",
                "r1": "R1.fastq.gz", "r2": "R2.fastq.gz",
            })
            membership = table("membership.tsv", {
                "sample_id": "S", "graph_id": "hprc_chm13_sv",
                "audit_status": "PASS", "reason": "absent",
            })
            manifest = table("index.tsv", {
                "graph_id": "hprc_chm13_sv", "vg_version": "1.76.1",
                "gbz_size_bytes": 100, "total_index_size_bytes": 400,
                "artificial_path_cover_created": "true",
                "artificial_path_cover_log_seen": "true", "validation": "PASS",
            })
            benchmark = table("index.benchmark.tsv", {"s": 3600, "max_rss": 1000})
            metrics = table("metrics.tsv", {
                "sample_id": "S", "sample_role": "tmem_candidate",
                "graph_id": "hprc_chm13_sv", "mapped_read_fraction": 0.99,
                "mapq_ge_20_fraction": 0.90, "median_mapq": 60,
                "mean_alignment_identity": 0.99, "mean_query_aligned_fraction": 0.98,
                "mean_clipped_bases_per_read": 2, "single_end_fallback": "false",
                "paired_fragment_learning_succeeded": "true",
                "mapping_wall_time_seconds": 60, "peak_rss_mb": 100,
                "reads_per_second": 1000, "learned_fragment_mean": 350,
                "learned_fragment_sd": 40,
            })
            contrasts = table("contrasts.tsv", {
                "row_type": "summary_overall", "contrast_id": "local_V1_minus_HPRC",
                "summary_group": "overall", "metric": "mapped_read_fraction",
                "n_samples": 1, "median_delta": 0.01, "min_delta": 0.01,
                "max_delta": 0.01,
            })
            targets = table("targets.tsv", {
                "graph_id": "hprc_chm13_sv", "target_present": "false",
                "callable_status": "not_callable_graph_target_absent",
                "failure_reason": "no_alternate_query_complete_high_identity_both-flank_mapping",
            })
            evidence = table("evidence.tsv", {
                "sample_id": "S", "sample_role": "tmem_candidate",
                "graph_id": "hprc_chm13_sv",
                "evidence_status": "not_callable_graph_target_absent",
            })
            tools = table("tools.tsv", {
                "tool": "vg", "version": "1.76.1", "status": "available",
                "required_check": "PASS",
            })
            provenance = table("provenance.tsv", {
                "key": "analysis_label", "value": "SV-graph short-read mapping QC",
            })
            output = root / "report.md"
            result = subprocess.run([
                sys.executable, str(SCRIPTS / "sv_short_read_report.py"),
                "--mode", "full", "--index-preflight", str(preflight),
                "--selected", str(selected), "--membership", str(membership),
                "--index-manifest", str(manifest), "--index-benchmark",
                f"hprc_chm13_sv={benchmark}", "--metrics", str(metrics),
                "--contrasts", str(contrasts), "--target-presence", str(targets),
                "--evidence", str(evidence), "--tools", str(tools),
                "--provenance", str(provenance), "--output", str(output),
            ], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            rendered = output.read_text()
            self.assertEqual(rendered.splitlines()[0], report.LEAD)
            self.assertIn("not mapping accuracy", rendered)
            self.assertIn("not mapper failure", rendered)
            self.assertIn("Local V2 remains experimental", rendered)


if __name__ == "__main__":
    unittest.main()
