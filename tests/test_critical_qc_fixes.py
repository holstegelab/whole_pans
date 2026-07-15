#!/usr/bin/env python3

import csv
import gzip
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

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


class CompleasmReuseTests(unittest.TestCase):
    def test_post_qc_reuses_persistent_original_compleasm_summaries(self):
        qc_rules = (PROJECT_ROOT / "workflow" / "rules" / "QC.smk").read_text()
        post_rules = (
            PROJECT_ROOT / "workflow" / "rules" / "post_decontamination_QC.smk"
        ).read_text()

        self.assertNotIn("rule post_qc_compleasm:", post_rules)
        self.assertIn("def all_original_compleasm", post_rules)
        self.assertIn(
            'f"{QC_OUTDIR}/compleasm/{{assembly}}/summary.txt"', post_rules
        )
        self.assertIn("--compleasm-results-dir", post_rules)
        self.assertIn(
            'summary=f"{QC_OUTDIR}/compleasm/{{assembly}}/summary.txt"', qc_rules
        )
        self.assertNotIn(
            'summary=temp(f"{QC_OUTDIR}/compleasm/{{assembly}}/summary.txt")',
            qc_rules,
        )


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


if __name__ == "__main__":
    unittest.main()
