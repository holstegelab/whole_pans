#!/usr/bin/env python3

import csv
import gzip
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = PROJECT_ROOT / "workflow" / "scripts"


def load_script(name):
    path = SCRIPTS / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class KrakenReportTests(unittest.TestCase):
    def parse_report(self, contents):
        module = load_script("kraken_report")
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "report.tsv"
            path.write_text(contents)
            return module.parse_kraken_report(path, {"2": "bacteria", "10239": "viruses"})

    def test_six_column_report_inherits_target_group(self):
        groups, names = self.parse_report(
            "100.00\t2\t0\tR\t1\troot\n"
            "50.00\t1\t0\tD\t2\t  Bacteria\n"
            "50.00\t1\t1\tS\t1353891\t    Achromobacter deleyi\n"
        )
        self.assertEqual(groups["1353891"], "bacteria")
        self.assertEqual(names["1353891"], "Achromobacter deleyi")

    def test_eight_column_minimizer_report_inherits_target_group(self):
        groups, names = self.parse_report(
            "100.00\t2\t0\t100\t80\tR\t1\troot\n"
            "50.00\t1\t0\t50\t40\tD\t2\t  Bacteria\n"
            "50.00\t1\t1\t50\t40\tS\t1353891\t    Achromobacter deleyi\n"
        )
        self.assertEqual(groups["1353891"], "bacteria")
        self.assertEqual(names["2"], "Bacteria")

    def test_ambiguous_report_width_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "6 or 8 columns"):
            self.parse_report("100.00\t1\t0\t10\tR\t1\troot\n")


class FilterFastaTests(unittest.TestCase):
    def test_cleaned_headers_are_prefixed_and_map_records_new_name(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            assembly = tmp / "input.fa"
            decisions = tmp / "decisions.tsv"
            split_bed = tmp / "split.bed"
            cleaned = tmp / "clean.fa.gz"
            removed = tmp / "removed.fa.gz"
            review = tmp / "review.fa.gz"
            split_map = tmp / "split_map.tsv"
            assembly.write_text(">h1tg000001l description\nACGT\n>contaminant\nAAAA\n")
            decisions.write_text(
                "contig\tdecision\n"
                "h1tg000001l\tKEEP\n"
                "contaminant\tREMOVE\n"
            )
            split_bed.write_text("")

            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPTS / "filter_fasta.py"),
                    "--assembly",
                    str(assembly),
                    "--decisions",
                    str(decisions),
                    "--split-bed",
                    str(split_bed),
                    "--cleaned",
                    str(cleaned),
                    "--removed",
                    str(removed),
                    "--review",
                    str(review),
                    "--split-map",
                    str(split_map),
                    "--header-prefix",
                    "sample.hap1",
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            with gzip.open(cleaned, "rt") as handle:
                self.assertEqual(handle.readline().strip(), ">sample.hap1.h1tg000001l description")
            with gzip.open(removed, "rt") as handle:
                self.assertEqual(handle.readline().strip(), ">contaminant decision=REMOVE")
            with split_map.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["original_contig"], "h1tg000001l")
            self.assertEqual(rows[0]["output_contig"], "sample.hap1.h1tg000001l")


class PostDecontaminationQcTests(unittest.TestCase):
    def test_post_qc_reuses_cleaned_stats_without_expensive_qc_jobs(self):
        qc_rules = (PROJECT_ROOT / "workflow" / "rules" / "QC.smk").read_text()
        post_rules = (
            PROJECT_ROOT / "workflow" / "rules" / "post_decontamination_QC.smk"
        ).read_text()

        self.assertIn("stats=all_cleaned_stats", post_rules)
        self.assertIn("--sequence-only", post_rules)
        self.assertNotIn("rule post_qc_fasta_stats:", post_rules)
        self.assertNotIn("rule post_qc_align_to_reference:", post_rules)
        self.assertNotIn("rule post_qc_paf_metrics:", post_rules)
        self.assertNotIn("rule post_qc_compleasm:", post_rules)
        self.assertNotIn("all_original_compleasm", post_rules)
        self.assertNotIn("--compleasm-results-dir", post_rules)
        self.assertIn(
            'summary=temp(f"{QC_OUTDIR}/compleasm/{{assembly}}/summary.txt")',
            qc_rules,
        )

    def test_sequence_only_classification_does_not_require_reused_metrics(self):
        module = load_script("summarize_qc")
        thresholds = {
            "fail": {
                "min_total_length_bp": 90,
                "max_total_length_bp": 110,
                "min_contig_n50_bp": 40,
                "max_n_percent": 1.0,
            },
            "warn": {
                "min_total_length_bp": 95,
                "max_total_length_bp": 105,
                "min_contig_n50_bp": 45,
                "max_n_percent": 0.5,
            },
        }
        row = {
            "total_length_bp": 100,
            "contig_n50_bp": 50,
            "n_percent": 0.0,
        }

        status, failures, warnings = module.classify(
            row, thresholds, checks=module.SEQUENCE_CHECKS
        )

        self.assertEqual(status, "PASS")
        self.assertEqual(failures, "")
        self.assertEqual(warnings, "")

    def test_sequence_only_cli_reads_decontamination_stats(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            results = tmp / "decontamination"
            stats = results / "stats"
            stats.mkdir(parents=True)
            assembly_ids = [
                "sample.hifi.hifiasm.bp.hap1.p_ctg",
                "sample.hifi.hifiasm.bp.hap2.p_ctg",
            ]
            manifest = tmp / "manifest.tsv"
            manifest.write_text(
                "assembly_id\tpath\n"
                + "".join(
                    f"{assembly_id}\t/cleaned/{assembly_id}.clean.fa.gz\n"
                    for assembly_id in assembly_ids
                )
            )
            for assembly_id in assembly_ids:
                (stats / f"{assembly_id}.clean.seqkit.tsv").write_text(
                    "num_seqs\tsum_len\tN50\tmax_len\tGC(%)\tsum_n\n"
                    "100\t100\t50\t60\t41.0\t0\n"
                )
            config = tmp / "config.yaml"
            config.write_text(
                "thresholds:\n"
                "  fail:\n"
                "    min_total_length_bp: 90\n"
                "    max_total_length_bp: 110\n"
                "    min_contig_n50_bp: 40\n"
                "    max_n_percent: 1.0\n"
                "  warn:\n"
                "    min_total_length_bp: 95\n"
                "    max_total_length_bp: 105\n"
                "    min_contig_n50_bp: 45\n"
                "    max_n_percent: 0.5\n"
            )
            assembly_output = tmp / "assembly.tsv"
            sample_output = tmp / "sample.tsv"
            included_output = tmp / "included.txt"
            excluded_output = tmp / "excluded.tsv"

            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPTS / "summarize_qc.py"),
                    "--manifest",
                    str(manifest),
                    "--config",
                    str(config),
                    "--results-dir",
                    str(results),
                    "--seqkit-suffix",
                    ".clean.seqkit.tsv",
                    "--sequence-only",
                    "--assembly-output",
                    str(assembly_output),
                    "--sample-output",
                    str(sample_output),
                    "--included-output",
                    str(included_output),
                    "--excluded-output",
                    str(excluded_output),
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            with assembly_output.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            included = included_output.read_text().splitlines()

        self.assertEqual([row["assembly_status"] for row in rows], ["PASS", "PASS"])
        self.assertNotIn("compleasm_complete_percent", rows[0])
        self.assertEqual(len(included), 2)


class PangenomeQcTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = load_script("pangenome_qc_analysis")

    def test_prefixed_reference_names_normalize_to_chromosome(self):
        self.assertEqual(self.module.normalize_chrom("CHM13.chr1"), "1")
        self.assertEqual(self.module.normalize_chrom("GRCh38.chrX"), "X")

    def test_build_log_counts_every_warning_form(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            log = Path(tmpdir) / "build.log"
            log.write_text(
                "Assemblies to add after CHM13 and hg38: 3\n"
                "[W] stable sequence 'x' associated with different ranks on segment 's1': 1 != 2\n"
                '[W::worker_for] stable sequence "x" already present in the graph. This will lead to inconsistent rGFA.\n'
                "[W::gfa_ins_filter] multi-link between >s1 and >s2\n"
            )
            parsed = self.module.parse_build_log(log)
        self.assertEqual(parsed["warning_lines"], 3)
        self.assertEqual(parsed["rgfa_consistency_warnings"], 2)
        self.assertEqual(parsed["warning_categories"]["worker_for"], 1)

    def test_stable_name_conflicts_are_counted(self):
        frame = pd.DataFrame(
            [
                {"sn": "sample.hap1.ctg1", "sr": 1, "ln": 100},
                {"sn": "sample.hap1.ctg1", "sr": 2, "ln": 50},
                {"sn": "sample.hap2.ctg1", "sr": 3, "ln": 75},
            ]
        )
        stats = self.module.stable_name_conflict_stats(frame)
        self.assertEqual(stats["stable_names_multiple_ranks"], 1)
        self.assertEqual(stats["stable_name_conflict_segments"], 2)
        self.assertEqual(stats["stable_name_conflict_bp"], 150)

    def test_mash_outliers_accepts_sequence_only_post_qc_columns(self):
        mash = pd.DataFrame(
            [
                {
                    "assembly_id": "sample.hifi.hifiasm.bp.hap1.p_ctg.clean",
                    "mash_distance": "0.01",
                    "matching_hashes": "50/100",
                }
            ]
        )
        assembly_qc = pd.DataFrame(
            [
                {
                    "assembly_id": "sample.hifi.hifiasm.bp.hap1.p_ctg",
                    "assembly_status": "PASS",
                    "warning_reasons": "",
                    "fail_reasons": "",
                    "contig_n50_bp": "50000000",
                }
            ]
        )
        args = SimpleNamespace(mad_k=3.5, low_matching_fraction=0.9)

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            (output_dir / "tables").mkdir()
            flagged, _paired = self.module.mash_outliers(
                args, mash, assembly_qc, pd.DataFrame(), output_dir
            )

        self.assertIn("contig_n50_bp", flagged.columns)
        self.assertNotIn("best_reference_covered_percent", flagged.columns)


if __name__ == "__main__":
    unittest.main()
