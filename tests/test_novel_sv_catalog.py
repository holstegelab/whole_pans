#!/usr/bin/env python3

import csv
import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "workflow" / "scripts" / "merge_novel_sv_catalog.py"


class NovelSvCatalogTests(unittest.TestCase):
    def test_graph_and_two_caller_evidence_cluster_with_provenance(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            graph = tmp / "graph.tsv.gz"
            with gzip.open(graph, "wt") as handle:
                handle.write(
                    "event_id\tassembly_id\tsample_id\thaplotype\tstable_source\tstable_position_0\t"
                    "source_rank\tchromosome\tsvtype\tsvlen\tevent_size_bp\tgraph_segment\t"
                    "segment_offset_0\tconfidence_tier\tfilter_reasons\n"
                    "GSV_1\ts.hifi.hifiasm.bp.hap1.p_ctg\ts\thap1\tCHM13.chr1\t100\t0\tchr1\tDEL\t-60\t60\ts1\t100\tHIGH_CONFIDENCE\t\n"
                )
            vcf = tmp / "svim.vcf.gz"
            with gzip.open(vcf, "wt") as handle:
                handle.write(
                    "##fileformat=VCFv4.2\n"
                    "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=type>\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    "chr1\t102\tsvim1\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;SVLEN=-61;END=162\n"
                    "chrX\t500\tsex_artifact\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;SVLEN=-80;END=580\n"
                )
            calls = tmp / "calls.tsv"
            calls.write_text(
                "caller\tcoordinate_system\tsample_id\tassembly_id\thaplotype\tpath\tcallable_bed\n"
                f"svim_asm\tCHM13\ts\ts.hifi.hifiasm.bp.hap1.p_ctg\thap1\t{vcf}\t\n"
            )
            manifest = tmp / "manifest.tsv"
            manifest.write_text(
                "assembly_id\tsample_id\thaplotype\tgraph_member\trescue_tier\n"
                "s.hifi.hifiasm.bp.hap1.p_ctg\ts\thap1\tfalse\tbest_rescue\n"
            )
            catalog = tmp / "catalog.tsv"
            evidence = tmp / "evidence.tsv.gz"

            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--graph-candidates",
                    str(graph),
                    "--call-manifest",
                    str(calls),
                    "--assembly-manifest",
                    str(manifest),
                    "--catalog-output",
                    str(catalog),
                    "--evidence-output",
                    str(evidence),
                    "--temp-dir",
                    str(tmp),
                    "--exclude-contig-regex",
                    "(?:chr)?[XY]",
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            with catalog.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            with gzip.open(evidence, "rt") as handle:
                evidence_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(len(evidence_rows), 2)
            self.assertEqual(rows[0]["confidence"], "HIGH")
            self.assertEqual(rows[0]["discovery_methods"], "assembly_svim_asm;graph_residual")
            self.assertEqual(rows[0]["caller_support"], "minigraph_cigar;svim_asm")
            self.assertEqual(rows[0]["carrier_samples"], "s")
            self.assertEqual(
                rows[0]["sex_chromosome_status"], "NOT_KNOWN_SEX_CHROMOSOME"
            )

    def test_source_specific_filter_keeps_only_par_aware_sex_chromosome_calls(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            graph = tmp / "graph.tsv.gz"
            with gzip.open(graph, "wt") as handle:
                handle.write(
                    "event_id\tassembly_id\tsample_id\thaplotype\tstable_source\t"
                    "stable_position_0\tsource_rank\tchromosome\tsvtype\tsvlen\t"
                    "event_size_bp\tgraph_segment\tsegment_offset_0\tconfidence_tier\t"
                    "filter_reasons\n"
                )
            vcf = tmp / "calls.vcf.gz"
            with gzip.open(vcf, "wt") as handle:
                handle.write(
                    "##fileformat=VCFv4.2\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    "chr1\t100\tauto\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;SVLEN=-60;END=160\n"
                    "chrX\t200\tsex\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;SVLEN=-60;END=260\n"
                )
            calls = tmp / "calls.tsv"
            calls.write_text(
                "caller\tcoordinate_system\tsample_id\tassembly_id\thaplotype\tpath\t"
                "callable_bed\texclude_contig_regex\n"
                f"svim_asm\tCHM13\ts\t\tdiploid\t{vcf}\t\t(?:chr)?[XY]\n"
                f"dipcall\tCHM13\ts\t\tdiploid\t{vcf}\t\t\n"
            )
            manifest = tmp / "manifest.tsv"
            manifest.write_text(
                "assembly_id\tsample_id\thaplotype\tgraph_member\trescue_tier\n"
                "s.hifi.hifiasm.bp.hap1.p_ctg\ts\thap1\tfalse\tbest_rescue\n"
            )
            catalog = tmp / "catalog.tsv"
            evidence = tmp / "evidence.tsv.gz"

            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--graph-candidates",
                    str(graph),
                    "--call-manifest",
                    str(calls),
                    "--assembly-manifest",
                    str(manifest),
                    "--catalog-output",
                    str(catalog),
                    "--evidence-output",
                    str(evidence),
                    "--temp-dir",
                    str(tmp),
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            with catalog.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            with gzip.open(evidence, "rt") as handle:
                evidence_rows = list(csv.DictReader(handle, delimiter="\t"))

            self.assertEqual(len(rows), 2)
            self.assertEqual(len(evidence_rows), 3)
            sex_rows = [row for row in rows if row["chromosome"] == "chrX"]
            self.assertEqual(len(sex_rows), 1)
            self.assertEqual(sex_rows[0]["caller_support"], "dipcall")
            self.assertEqual(
                sex_rows[0]["sex_chromosome_status"], "KNOWN_SEX_CHROMOSOME"
            )

    def test_unprojected_graph_origin_is_filtered_when_known_and_flagged_when_unknown(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            graph = tmp / "graph.tsv.gz"
            with gzip.open(graph, "wt") as handle:
                handle.write(
                    "event_id\tassembly_id\tsample_id\thaplotype\tstable_source\t"
                    "stable_position_0\tsource_rank\tchromosome\tsvtype\tsvlen\t"
                    "event_size_bp\tgraph_segment\tsegment_offset_0\tconfidence_tier\t"
                    "filter_reasons\n"
                    "GSV_X\ta\ts\thap1\tdonor.chrX\t10\t2\tchrX\tINS\t60\t60\tsx\t10\tHIGH_CONFIDENCE\t\n"
                    "GSV_U\ta\ts\thap1\tdonor.contig\t20\t2\tdonor.contig\tINS\t70\t70\tsu\t20\tHIGH_CONFIDENCE\t\n"
                )
            calls = tmp / "calls.tsv"
            calls.write_text(
                "caller\tcoordinate_system\tsample_id\tassembly_id\thaplotype\tpath\t"
                "callable_bed\texclude_contig_regex\n"
            )
            manifest = tmp / "manifest.tsv"
            manifest.write_text(
                "assembly_id\tsample_id\thaplotype\tgraph_member\trescue_tier\n"
                "a\ts\thap1\tfalse\tbest_rescue\n"
            )
            catalog = tmp / "catalog.tsv"
            evidence = tmp / "evidence.tsv.gz"

            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--graph-candidates",
                    str(graph),
                    "--call-manifest",
                    str(calls),
                    "--assembly-manifest",
                    str(manifest),
                    "--catalog-output",
                    str(catalog),
                    "--evidence-output",
                    str(evidence),
                    "--temp-dir",
                    str(tmp),
                    "--exclude-contig-regex",
                    "(?:chr)?[XY]",
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            with catalog.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["graph_segment"], "su")
            self.assertEqual(
                rows[0]["sex_chromosome_status"], "UNRESOLVED_GRAPH_ORIGIN"
            )


if __name__ == "__main__":
    unittest.main()
