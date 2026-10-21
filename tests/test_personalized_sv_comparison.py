"""Focused synthetic checks for reusing Giraffe alignments in the four-graph comparison."""

import argparse
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "workflow" / "scripts"))
import personalized_sv_statistics as stats
import tmem_whole_graph_target as target


GRAPHS = ("hprc", "local_v1", "local_v2", "hprc_plus_local")


def table(path, rows, fields=None):
    stats.write_tsv(path, rows, fields or list(rows[0]))
    return path


def gaf(name, traversal, mapq=60, secondary=False):
    return (f"{name}\t100\t0\t100\t+\t{traversal}\t1000\t0\t100\t100\t100\t{mapq}"
            f"\ttp:A:{'S' if secondary else 'P'}\tcg:Z:100=\n")


def source_target():
    queries = [{
        "query_name": name, "allele_type": kind, "query_length": "2100",
        "left_flank_bp": "1000", "right_flank_bp": "1000",
        "allele_start": "1000", "allele_end": "1100", "derivation_source": "synthetic",
    } for name, kind in (("REF", "reference"), ("ALT", "alternate"))]

    def alignment(nodes):
        return {"query_length": 2100, "query_start": 0, "query_end": 2100,
                "coverage": 1.0, "identity": 1.0, "mapq": 60, "cigar": "2100=",
                "primary": True, "traversal": nodes}

    return target.derive_target("local_v1", queries, {
        "REF": [alignment([">left", ">ref", ">right"])],
        "ALT": [alignment([">left", "<insert", ">right"])],
    }, 0.95, 0.90, 20, 900)


class TargetTests(unittest.TestCase):
    def test_stable_coordinates_require_explicit_segment_output(self):
        for path in ("CHM13.chr7", ">CHM13#0#chr7:100-200>sample:5-25"):
            with self.assertRaisesRegex(ValueError, "--vc"):
                target.path_tokens(path)
        self.assertEqual(target.canonical_traversal("<b<a", "-"), [">a", ">b"])

    def test_translation_uses_chopped_nodes_and_reverse_orientation(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "translation.tsv"
            path.write_text("T\tleft\t101,102\nT\tref\t900,901\n"
                            "T\tinsert\t300,301,302\nT\tright\t501,502\n")
            row = target.translate_target(source_target(), path)
            self.assertEqual(row["left_junction"], ">102,<302")
            self.assertEqual(row["right_junction"], "<300,>501")
            self.assertEqual(set(row["reference_defining_edges"].split(";")),
                             {">102,>900", ">901,>501"})
            self.assertEqual(row["insertion_nodes"], "300,301,302")
            self.assertEqual(row["edge_namespace"], "vg_node_ids")

    def test_missing_translation_is_an_error(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "translation.tsv"
            path.write_text("T\tleft\t101,102\n")
            with self.assertRaisesRegex(ValueError, "lacks target segments"):
                target.translate_target(source_target(), path)

    def test_inserted_flank_bases_do_not_count_as_aligned(self):
        record = {"query_start": 0, "cigar": "200I1900=", "strand": "+"}
        self.assertEqual(target.flank_aligned_bases(record, 0, 1000), 800)
        record["strand"] = "-"
        self.assertEqual(target.flank_aligned_bases(record, 1100, 2100), 800)


class EvidenceTests(unittest.TestCase):
    def prepare(self, root, callable_target=True):
        row = {
            "graph_id": "local_v1", "target_present": "true",
            "callable_status": ("callable_complete_graph_specific_target" if callable_target else
                                "not_callable_target_ambiguous"),
            "edge_namespace": "vg_node_ids", "left_junction": ">10,>20",
            "right_junction": ">21,>40", "reference_defining_edges": ">10,>30;>30,>40",
            "insertion_edges": ">20,>21", "insertion_nodes": "20,21", "failure_reason": "synthetic",
        }
        target_file = table(root / "target.tsv", [row])
        inventory = table(root / "inventory.tsv", [{
            "target": "local_v1", "source_graph": "source.gfa", "graph_checksum": "checksum",
        }])
        return argparse.Namespace(gaf=root / "reads.gaf", target_table=target_file, inventory=inventory,
                                  sample="S", sample_group="prior", batch="B", target="local_v1",
                                  mapper_version="1.76.1", min_mapq=20, min_junction_reads=3,
                                  min_reference_reads=3, output=root / "evidence.tsv")

    def test_uncallable_counts_are_blank_and_gaf_is_not_opened(self):
        with tempfile.TemporaryDirectory() as directory:
            args = self.prepare(Path(directory), False)
            stats.evidence_stats(args)  # reads.gaf deliberately does not exist
            row = stats.read_tsv(args.output)[0]
            self.assertEqual(row["call_class"], "NO_DECISION")
            self.assertEqual(row["left_alt_junction_reads"], "")
            self.assertEqual(row["reference_path_reads"], "")
            self.assertEqual(row["evidence_available"], "false")

    def test_mates_reverse_edges_and_unknown_mapq(self):
        with tempfile.TemporaryDirectory() as directory:
            args = self.prepare(Path(directory))
            records = []
            for index in range(3):
                records.extend([
                    gaf(f"L{index}/1", ">10>20"), gaf(f"L{index}/2", ">10>20"),
                    gaf(f"R{index}", "<40<21"), gaf(f"F{index}", ">10>30>40"),
                ])
            records.extend([gaf("unrelated", ">1000>1001"), gaf("unknown", ">10>20", 255),
                            gaf("low", ">10>20", 5), gaf("secondary", ">10>20", secondary=True)])
            args.gaf.write_text("".join(records))
            stats.evidence_stats(args)
            row = stats.read_tsv(args.output)[0]
            self.assertEqual(row["call_class"], "HET_LIKE")
            self.assertEqual([row[field] for field in ("left_alt_junction_reads", "right_alt_junction_reads",
                                                       "reference_path_reads")], ["3", "3", "3"])
            self.assertEqual(row["locus_mapq_qualified_fragments"], "9")
            self.assertEqual(row["genome_mapq_qualified_alignments"], "13")

    def test_conflicting_fragments_do_not_create_a_heterozygote(self):
        with tempfile.TemporaryDirectory() as directory:
            args = self.prepare(Path(directory))
            args.gaf.write_text("".join(gaf(f"C{i}/1", ">10>20") + gaf(f"C{i}/2", ">10>30")
                                        for i in range(3)))
            stats.evidence_stats(args)
            row = stats.read_tsv(args.output)[0]
            self.assertEqual(row["call_class"], "CONFLICTING")
            self.assertEqual(row["conflicting_fragments"], "3")
            self.assertEqual(row["left_alt_junction_reads"], "0")
            self.assertEqual(row["reference_path_reads"], "0")


class ComparisonTests(unittest.TestCase):
    def test_read_tsv_accepts_large_graph_catalog_fields(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "catalog.tsv"
            long_segment_list = "segment," * 20000
            path.write_text(f"target\tsegment_ids\nlocal_v1\t{long_segment_list}\n")
            self.assertEqual(stats.read_tsv(path)[0]["segment_ids"], long_segment_list)

    def test_reference_aliases_recover_shared_sites_without_claiming_same_allele(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bubbles = []
            for graph in GRAPHS:
                path = root / f"{graph}.bed"
                reference = "CHM13.chr7" if graph.startswith("local") else "CHM13#0#chr7"
                path.write_text(f"{reference}\t100\t100\t3\t2\t0\t0\t100\t0\t0\t0\ta,b,c\n")
                bubbles.append(f"{graph}={path}")
            output = root / "catalog.tsv"
            stats.build_catalog(argparse.Namespace(bubble=bubbles, reference_sample="CHM13", output=output))
            rows = stats.read_tsv(output)
            self.assertEqual(len(rows), 4)
            self.assertEqual({row["coordinate_key"] for row in rows}, {"chr7:100-100:INS"})
            self.assertEqual({row["shared_graph_count"] for row in rows}, {"4"})
            self.assertEqual({row["allele_equivalence"] for row in rows}, {"not_assessed"})
            self.assertEqual(stats.canonical_reference_path("GRCh38#0#chr7"), "GRCh38#0#chr7")

    def test_histogram_statistics_exclude_missing_mapq(self):
        self.assertEqual(stats.histogram_median(Counter({0: 1, 60: 3})), 60)
        self.assertEqual(stats.histogram_median(Counter({0: 1, 60: 1})), 30)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reads.gaf"
            path.write_text(gaf("a", ">1>2", 0) + gaf("b", ">1>2", 60) + gaf("c", ">1>2", 255))
            result = stats.mapping_from_gaf(path, 3, 300)
            self.assertEqual(result["median_mapq"], 30)
            self.assertEqual(result["mean_mapq"], 30)
            self.assertAlmostEqual(result["mapq_ge_20_fraction"], 1 / 3)

    def test_paired_effect_units_and_missing_pairs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rows = [{"sample_id": "S", "sample_group": "prior", "batch": "B", "target": graph,
                     "input_reads": "1000", "input_bases": "100000",
                     "mapq_ge_20_fraction": "0.91" if graph.startswith("local") else "0.90"}
                    for graph in GRAPHS]
            args = argparse.Namespace(metrics=table(root / "metrics.tsv", rows), output=root / "paired.tsv",
                                      permutations=100, seed=1)
            stats.paired_comparison(args)
            result = next(row for row in stats.read_tsv(args.output) if row["level"] == "sample"
                          and row["contrast"] == "local_v1_minus_hprc" and row["metric"] == "mapq_ge_20_fraction")
            self.assertAlmostEqual(float(result["difference_percentage_points"]), 1.0)
            self.assertEqual(result["read_count_difference"], "10")
            table(args.metrics, rows[:-1])
            with self.assertRaisesRegex(ValueError, "Incomplete four-graph"):
                stats.paired_comparison(args)

    def test_report_does_not_rank_missing_tmem_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sample = {"sample_id": "S", "sample_group": "prior", "batch": "B"}
            inventory = [{"target": graph, "graph_checksum": graph} for graph in GRAPHS]
            metadata = [{**sample, **row} for row in inventory]
            metrics = [{**row, "input_reads": "100", "input_bases": "10000",
                        "mapped_read_fraction": "0.99", "mapq_ge_20_fraction": "0.9"} for row in metadata]
            evidence = [{**row, "call_class": "NO_DECISION", "min_mapq": "20",
                         "min_junction_fragments": "3", "min_reference_fragments": "3"} for row in metadata]
            audits = [{"graph_id": graph, "callable_status": "not_callable_target_ambiguous",
                       "reference_query_coverage": "", "reference_query_identity": "",
                       "alternate_query_coverage": "", "alternate_query_identity": "",
                       "edge_namespace": "vg_node_ids", "failure_reason": "synthetic"} for graph in GRAPHS]
            args = argparse.Namespace(
                samples=table(root / "samples.tsv", [sample]), inventory=table(root / "inventory.tsv", inventory),
                membership=table(root / "members.tsv", [{**row, "is_graph_member": "false", "membership_status": "known"}
                                                       for row in metadata]),
                catalog=table(root / "catalog.tsv", [], ["target", "coordinate_key", "shared_graph_count"]),
                metrics=table(root / "metrics.tsv", metrics), evidence=table(root / "evidence.tsv", evidence),
                comparison=root / "paired.tsv", targets=table(root / "targets.tsv", audits),
                validation=table(root / "validation.tsv", [{**sample, "runtime_validation": "PASS"}]),
                output=root / "report.md",
            )
            stats.paired_comparison(argparse.Namespace(metrics=args.metrics, output=args.comparison, permutations=100, seed=1))
            stats.render_report(args)
            report = args.output.read_text()
            self.assertIn("TMEM comparison unavailable", report)
            self.assertIn("No graph can be ranked", report)
            self.assertIn("NO_DECISION (N/A/N/A/N/A)", report)


if __name__ == "__main__":
    unittest.main()
