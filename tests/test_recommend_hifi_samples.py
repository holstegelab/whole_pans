#!/usr/bin/env python3

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "workflow" / "scripts" / "recommend_hifi_samples.py"


class RecommendHifiSamplesTests(unittest.TestCase):
    def test_large_catalog_field_is_streamed_and_ranked(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            manifest = tmp / "manifest.tsv"
            manifest.write_text(
                "assembly_id\tsample_id\thaplotype\tgraph_member\trescue_tier\t"
                "graph_context\tmate_status\n"
                "sample.hap1\tsample\thap1\tfalse\tfragmented_rescue\t"
                "outside_graph\tPRESENT\n"
            )
            screen = tmp / "screen.tsv"
            screen.write_text(
                "assembly_id\tprimary_callable_fraction\t"
                "sensitivity_callable_fraction\n"
                "sample.hap1\t0.70\t0.80\n"
            )
            catalog = tmp / "catalog.tsv"
            catalog_header = (
                "event_id\tsvtype\tcarrier_samples\tindependent_sample_count\t"
                "validation_status\tconfidence\tdiscovery_methods\tevidence_ids\n"
            )
            catalog_rows = [
                "PSV_00\tDEL\tsample\t1\tPENDING_HIFI\tMEDIUM\tgraph_residual\t"
                + "E" * 200_000
                + "\n",
                "PSV_01\tDEL\tsample\t2\tASSEMBLY_ONLY_REVIEW\tHIGH\t"
                "graph_residual\tEVD_1\n",
            ]
            catalog_rows.extend(
                f"PSV_{index:02d}\tDEL\tsample\t2\tASSEMBLY_ONLY_REVIEW\t"
                f"MEDIUM\tgraph_residual\tEVD_{index}\n"
                for index in range(2, 30)
            )
            catalog.write_text(catalog_header + "".join(catalog_rows))
            evidence = tmp / "evidence.tsv"
            evidence.write_text(
                "sample_id\tdiscovery_method\tconfidence_tier\tsvtype\n"
                "sample\tgraph_residual\tHIGH_CONFIDENCE\tDEL\n"
            )
            output = tmp / "recommended.tsv"

            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--assembly-manifest",
                    str(manifest),
                    "--screen-summary",
                    str(screen),
                    "--catalog",
                    str(catalog),
                    "--evidence",
                    str(evidence),
                    "--output",
                    str(output),
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            with output.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["candidate_count"], "30")
            self.assertEqual(
                rows[0]["important_candidate_ids"],
                ";".join(f"PSV_{index:02d}" for index in range(25)),
            )
            self.assertEqual(rows[0]["priority_tier"], "P1_DISCOVERY_BLIND_SPOT")
            self.assertEqual(
                rows[0]["reason_codes"],
                "P1_DISCOVERY_BLIND_SPOT;P1_VALIDATE_CANDIDATE",
            )
            self.assertEqual(
                rows[0]["candidate_counts_by_method"], "graph_residual=30"
            )

    def test_validation_is_capped_and_graph_member_control_is_reserved(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            manifest = tmp / "manifest.tsv"
            manifest.write_text(
                "assembly_id\tsample_id\thaplotype\tgraph_member\trescue_tier\t"
                "graph_context\tmate_status\n"
                "blind.h1\tblind\thap1\tfalse\tfragmented_rescue\toutside\tFAIL\n"
                "candidate1.h1\tcandidate1\thap1\tfalse\tbest_rescue\toutside\tPASS\n"
                "candidate2.h1\tcandidate2\thap1\tfalse\tbest_rescue\toutside\tPASS\n"
                "control.h1\tcontrol\thap1\ttrue\tgraph_member\tcurrent\tIN_GRAPH\n"
            )
            screen = tmp / "screen.tsv"
            screen.write_text(
                "assembly_id\tprimary_callable_fraction\t"
                "sensitivity_callable_fraction\n"
                "blind.h1\t0.95\t0\n"
                "candidate1.h1\t0.96\t0\n"
                "candidate2.h1\t0.97\t0\n"
                "control.h1\t0.98\t0\n"
            )
            catalog = tmp / "catalog.tsv"
            catalog.write_text(
                "event_id\tsvtype\tcarrier_samples\tvalidation_status\t"
                "confidence\tdiscovery_methods\n"
                "PSV_1\tDEL\tcandidate1\tLINEAR_CALLER_SUPPORTED\tHIGH\t"
                "graph_residual;assembly_dipcall\n"
                "PSV_2\tINV\tcandidate2\tASSEMBLY_ONLY_REVIEW\tMEDIUM\t"
                "graph_residual\n"
            )
            output = tmp / "recommended.tsv"

            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--assembly-manifest",
                    str(manifest),
                    "--screen-summary",
                    str(screen),
                    "--catalog",
                    str(catalog),
                    "--output",
                    str(output),
                    "--validation-count",
                    "1",
                    "--control-count",
                    "1",
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            with output.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 3)
            self.assertEqual(
                [row["priority_tier"] for row in rows],
                [
                    "P1_DISCOVERY_BLIND_SPOT",
                    "P1_VALIDATE_CANDIDATE",
                    "P3_CONTROL",
                ],
            )
            self.assertEqual(rows[1]["sample_id"], "candidate1")
            self.assertEqual(rows[2]["sample_id"], "control")


if __name__ == "__main__":
    unittest.main()
