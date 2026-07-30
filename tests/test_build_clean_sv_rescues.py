#!/usr/bin/env python3

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "workflow" / "scripts" / "build_clean_sv_rescues.py"
SPEC = importlib.util.spec_from_file_location("build_clean_sv_rescues", SCRIPT)
RESCUE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RESCUE)


def write_tsv(path, fields, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


class RescuePlanningTests(unittest.TestCase):
    def test_default_plan_chooses_eligible_source_and_computes_seed(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            catalog = tmp / "catalog.tsv"
            candidates = tmp / "candidates.tsv"
            assemblies = tmp / "assemblies.tsv"
            reads = tmp / "reads.tsv"
            output = tmp / "plan.tsv"

            write_tsv(
                catalog,
                [
                    "event_id",
                    "coordinate_system",
                    "chromosome",
                    "position_0",
                    "svtype",
                    "svlen",
                    "confidence",
                    "validation_status",
                    "independent_sample_count",
                    "linear_supporting_sample_count",
                    "discovery_methods",
                    "carrier_assemblies",
                ],
                [
                    {
                        "event_id": "PSV_keep",
                        "coordinate_system": "CHM13",
                        "chromosome": "chr1",
                        "position_0": 1000,
                        "svtype": "INS",
                        "svlen": 60,
                        "confidence": "UNCERTAIN",
                        "validation_status": "PENDING_HIFI",
                        "independent_sample_count": 2,
                        "linear_supporting_sample_count": 2,
                        "discovery_methods": (
                            "assembly_dipcall;assembly_svim_asm;graph_residual"
                        ),
                        "carrier_assemblies": "assembly.hap1",
                    },
                    {
                        "event_id": "PSV_excluded",
                        "coordinate_system": "CHM13",
                        "chromosome": "chr2",
                        "position_0": 2000,
                        "svtype": "DEL",
                        "svlen": 80,
                        "confidence": "UNCERTAIN",
                        "validation_status": "PENDING_HIFI",
                        "independent_sample_count": 1,
                        "linear_supporting_sample_count": 1,
                        "discovery_methods": (
                            "assembly_dipcall;assembly_svim_asm;graph_residual"
                        ),
                        "carrier_assemblies": "other.hap1",
                    },
                ],
            )
            write_tsv(
                candidates,
                [
                    "event_id",
                    "assembly_id",
                    "sample_id",
                    "haplotype",
                    "graph_member",
                    "mapping_pass",
                    "query_name",
                    "query_length",
                    "query_position_0",
                    "graph_segment",
                    "segment_offset_0",
                    "stable_source",
                    "stable_position_0",
                    "source_rank",
                    "chromosome",
                    "svtype",
                    "svlen",
                    "event_size_bp",
                    "left_anchor_bp",
                    "right_anchor_bp",
                    "identity",
                    "mapq",
                    "primary",
                    "confidence_tier",
                    "filter_reasons",
                ],
                [
                    {
                        "event_id": "GSV_source",
                        "assembly_id": "assembly.hap1",
                        "sample_id": "sample",
                        "haplotype": "hap1",
                        "graph_member": "false",
                        "mapping_pass": "primary",
                        "query_name": "contig1",
                        "query_length": 300000,
                        "query_position_0": 100000,
                        "graph_segment": "chr1:0-2000",
                        "segment_offset_0": 1000,
                        "stable_source": "chr1",
                        "stable_position_0": 1005,
                        "source_rank": 0,
                        "chromosome": "chr1",
                        "svtype": "INS",
                        "svlen": 60,
                        "event_size_bp": 60,
                        "left_anchor_bp": 50000,
                        "right_anchor_bp": 50000,
                        "identity": 0.999,
                        "mapq": 60,
                        "primary": "true",
                        "confidence_tier": "HIGH_CONFIDENCE",
                        "filter_reasons": "",
                    }
                ],
            )
            write_tsv(
                assemblies,
                [
                    "assembly_id",
                    "sample_id",
                    "haplotype",
                    "cleaned_fasta_path",
                    "rescue_tier",
                ],
                [
                    {
                        "assembly_id": "assembly.hap1",
                        "sample_id": "sample",
                        "haplotype": "hap1",
                        "cleaned_fasta_path": "/assemblies/sample.hap1.fa.gz",
                        "rescue_tier": "fragmented_rescue",
                    }
                ],
            )
            write_tsv(
                reads,
                ["sample_id", "hifi_path"],
                [
                    {
                        "sample_id": "sample",
                        "hifi_path": "/reads/sample.hifi.fastq.gz",
                    }
                ],
            )

            RESCUE.build_plan(
                SimpleNamespace(
                    catalog=catalog,
                    graph_candidates=candidates,
                    assembly_manifest=assemblies,
                    read_manifest=reads,
                    output=output,
                    event_ids=None,
                    confidence="UNCERTAIN",
                    validation_status="PENDING_HIFI",
                    require_methods=(
                        "graph_residual;assembly_dipcall;assembly_svim_asm"
                    ),
                    min_independent_samples=2,
                    max_events=100,
                    source_qc_tiers="fragmented_rescue;not_recommended",
                    breakpoint_distance=500,
                    length_similarity=0.70,
                    flank_bp=50000,
                    min_flank_bp=20000,
                    allow_contig_end=False,
                )
            )

            rows = RESCUE.read_tsv(output)
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertEqual(row["rescue_id"], "PSV_keep")
            self.assertEqual(row["plan_status"], "READY")
            self.assertEqual(row["source_event_id"], "GSV_source")
            self.assertEqual(row["seed_start_0"], "50000")
            self.assertEqual(row["seed_end_0"], "150060")
            self.assertEqual(row["hifi_read_paths"], "/reads/sample.hifi.fastq.gz")

class RescueSequenceTests(unittest.TestCase):
    def test_fastq_read_extraction_and_reverse_complement(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            reads = tmp / "reads.fastq"
            output = tmp / "selected.fa.gz"
            reads.write_text(
                "@read1 description\nACGTN\n+\nIIIII\n" "@read2\nTTAA\n+\nIIII\n"
            )

            count = RESCUE.extract_selected_reads([reads], {"read2"}, output)

            self.assertEqual(count, 1)
            self.assertEqual(RESCUE.parse_fasta(output), {"read2": "TTAA"})
            self.assertEqual(RESCUE.reverse_complement("ACGTRY"), "RYACGT")

    def test_contig_selection_prefers_agreement_at_source_event(self):
        base = {
            "query_length": 100060,
            "query_start": 0,
            "query_end": 100060,
            "strand": "+",
            "target_name": "seed",
            "target_length": 100060,
            "target_start": 0,
            "target_end": 100060,
            "matches": 100000,
            "alignment_length": 100060,
            "mapq": 60,
        }
        discordant = dict(
            base,
            query_name="discordant",
            tags={"cg": "50000M60I50060M"},
        )
        agreeing = dict(
            base,
            query_name="agreeing",
            matches=100060,
            tags={"cg": "100060M"},
        )

        selected = RESCUE.select_local_contig(
            [discordant, agreeing],
            event_start=50000,
            event_end=50060,
            seed_length=100060,
            anchor_bp=10000,
            min_mapq=20,
            min_identity=0.99,
            discordance_window=500,
        )

        self.assertEqual(selected["query_name"], "agreeing")
        self.assertEqual(selected["event_discordance_bp"], 0)
        self.assertEqual(
            RESCUE.cigar_event_discordance("50000M60I50060M", 0, 50000, 50060, 500),
            60,
        )

    def test_read_support_requires_spanning_and_event_agreement(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paf = Path(tmpdir) / "reads.paf"
            paf.write_text(
                "agree\t5000\t0\t5000\t+\tseed\t10000\t2500\t7500\t"
                "5000\t5000\t60\tcg:Z:5000M\n"
                "other_haplotype\t5060\t0\t5060\t+\tseed\t10000\t"
                "2500\t7500\t5000\t5060\t60\tcg:Z:2500M60I2500M\n"
                "short\t2000\t0\t2000\t+\tseed\t10000\t4500\t6500\t"
                "2000\t2000\t60\tcg:Z:2000M\n"
                "secondary\t5000\t0\t5000\t+\tseed\t10000\t2500\t7500\t"
                "5000\t5000\t60\ttp:A:S\tcg:Z:5000M\n"
            )

            selected, spanning, supporting = RESCUE.read_selection_from_paf(
                paf,
                event_start=5000,
                event_end=5001,
                min_mapq=10,
                min_aligned_bp=2000,
                read_anchor_bp=1000,
                max_event_discordance_bp=30,
            )

            self.assertEqual(
                selected, {"agree", "other_haplotype", "short"}
            )
            self.assertEqual(spanning, {"agree", "other_haplotype"})
            self.assertEqual(supporting, {"agree"})


if __name__ == "__main__":
    unittest.main()
