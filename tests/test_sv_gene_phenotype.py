"""Tiny server-only regression fixtures; no native bioinformatics tools are called."""

import csv
import io
import json
import runpy
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace as NS
from unittest.mock import patch

SCRIPTS = Path(__file__).absolute().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))
import svgp_inventory as inventory
import svgp_catalogue as catalogue
import svgp_measurement as measurement


class Fixtures(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.evidence = self.root / "review.txt"
        self.evidence.write_text("Reviewed tiny fixture\n")

    def table(self, name, records, fields=None):
        path = self.root / name
        records = list(records)
        fields = fields or list(dict.fromkeys(k for r in records for k in r))
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(records)
        return path

    def read(self, path):
        return list(measurement.rows(path))

    def inventory_args(self):
        person = dict.fromkeys(inventory.SAMPLE_FIELDS, "unknown")
        person.update(person_id="P1", illumina_id="I1", fastq_manifest_id="F1",
                      phenotype_id="PH1", phenotype="case", sex="female", age_definition="at_sampling",
                      ancestry_covariates="PC1-PC10", discovery_donor="no", related_to_discovery_donor="no",
                      analysis_eligible="yes", exclusion_reason="none", covariates_available="yes",
                      snp_data_available="yes", fastq_manifest_complete="yes", unrelated_primary="yes")
        lane = dict.fromkeys(inventory.READ_FIELDS, "unknown")
        lane.update(read_id="R1", person_id="P1", illumina_id="I1", fastq_manifest_id="F1",
                    library_id="LIB", lane="L1", batch="B1")
        for mate in ("r1", "r2"):
            path = self.root / (mate + ".fq")
            path.write_text("@read\nAC\n+\nII\n")
            lane[mate] = str(path)
        return NS(samples=self.table("samples.tsv", [person]), reads=self.table("reads.tsv", [lane]),
                  resources=self.table("resources.tsv", [], inventory.RESOURCE_FIELDS), crosswalk=None,
                  output=self.root / "inventory"), person, lane

    def test_explicit_alias_and_unknown_eligibility(self):
        args, person, lane = self.inventory_args()
        lane.update(person_id="unknown", illumina_id="I1_lane_label")
        args.reads = self.table("reads.tsv", [lane])
        inventory.inventory(args)
        self.assertEqual(self.read(args.output / "metadata/samples.tsv")[0]["analysis_eligible"], "no")
        args.crosswalk = self.table("crosswalk.tsv", [dict(namespace="illumina_id", alias="I1_lane_label", person_id="P1")])
        inventory.inventory(args)
        self.assertEqual(self.read(args.output / "metadata/samples.tsv")[0]["analysis_eligible"], "yes")
        person["related_to_discovery_donor"] = "unknown"
        args.samples = self.table("samples.tsv", [person])
        inventory.inventory(args)
        self.assertEqual(self.read(args.output / "metadata/samples.tsv")[0]["analysis_eligible"], "no")

    def test_ambiguous_or_unknown_person_crosswalk_fails(self):
        args, person, lane = self.inventory_args()
        other = dict(person, person_id="P2", illumina_id="I2", fastq_manifest_id="F2", phenotype_id="PH2")
        args.samples = self.table("samples.tsv", [person, other])
        for alias, target in (("I1", "P2"), ("other", "MISSING")):
            with self.subTest(alias=alias):
                args.crosswalk = self.table("crosswalk.tsv", [dict(namespace="illumina_id", alias=alias, person_id=target)])
                with self.assertRaises(ValueError):
                    inventory.inventory(args)

    def test_fingerprint_confirmation_requires_existing_evidence(self):
        args, person, lane = self.inventory_args()
        person.update(ont_id="ONT1", ont_path=str(self.evidence), ont_fingerprint_status="confirmed",
                      ont_fingerprint_evidence=str(self.root / "missing_fingerprint.txt"))
        args.samples = self.table("samples.tsv", [person])
        inventory.inventory(args)
        pairing = self.read(args.output / "metadata/ont_pairing.tsv")[0]
        self.assertEqual(pairing["status"], "ID_LINK_ONLY")
        person["ont_fingerprint_evidence"] = str(self.evidence)
        args.samples = self.table("samples.tsv", [person])
        inventory.inventory(args)
        pairing = self.read(args.output / "metadata/ont_pairing.tsv")[0]
        self.assertEqual(pairing["status"], "DECLARED_FINGERPRINT_CONFIRMED")

    def test_cohort_release_is_bound_to_consumed_pilot_and_ont(self):
        pilot = self.table("pilot.tsv", [dict(graph_id="main", locus_id="L1", allele_id="A1", callable="yes")])
        ont = self.table("ont.tsv", [dict(locus_id="L1", allele_id="A1", status="RESOLVED_CARRIER")])
        lock = self.table("lock.tsv", [dict(graph_id="main", locus_id="L1", allele_id="A1",
                                             status="LOCKED", technical_status="PASS")])
        review = dict(graph_id="main", status="APPROVED", pilot_evidence=str(pilot), ont_evidence=str(ont),
                      calibration_evidence=str(self.evidence), compute_budget_evidence=str(self.evidence),
                      retention_policy="preserve", reviewer="fixture")
        args = NS(graph="main", pilot=pilot, ont=ont, lock=lock, output=self.root / "release.tsv")
        for field in ("pilot_evidence", "ont_evidence"):
            args.review = self.table("release_review.tsv", [dict(review, **{field: str(self.evidence)})])
            with self.subTest(field=field), self.assertRaises(ValueError):
                measurement.release(args)
        (self.root / "nested").mkdir()
        review["pilot_evidence"] = str(self.root / "nested" / ".." / "pilot.tsv")
        args.review = self.table("release_review.tsv", [review])
        measurement.release(args)
        self.assertEqual(self.read(args.output)[0]["status"], "APPROVED")

    def test_forward_reverse_source_and_zero_length_reference_boundary(self):
        inserted = "AGT"
        reference = self.root / "reference.fa"
        reference.write_text(">chr1\n" + "A" * 100 + "C" * 100 + "\n")
        alleles = self.root / "insertions.fa"
        alleles.write_text(">C1\n" + inserted + "\n")
        for orientation in ("+", "-"):
            with self.subTest(orientation=orientation):
                forward = "A" * 100 + inserted + "C" * 100
                source = self.root / "source.fa"
                source.write_text(">ctg\n" + (forward if orientation == "+" else catalogue.reverse_complement(forward)) + "\n")
                row = dict(candidate_id="C1", sv_type="INS", reference_length_bp=0, provenance_status="PASS",
                           sequence_sha256=catalogue.digest(inserted), source_assembly=str(source), source_name="ctg",
                           source_offset_0based=100, chromosome="chr1", position_0based=100, reference_build="CHM13")
                args = NS(candidates=self.table("candidates.tsv", [row]), alleles=alleles, reference=reference,
                          graph_version="fixture", graph_sha256="a" * 64, flank=100, output=self.root / "catalogue",
                          source_review=None, tmem_fasta=None, tmem_targets=None, controls=None, control_fasta=None)
                catalogue.catalogue(args)
                result = self.read(args.output / "alleles.tsv")[0]
                self.assertEqual((result["source_orientation"], result["status"]), (orientation, "PASS"))
                self.assertEqual((result["reference_start_0"], result["reference_end_0"]), ("100", "100"))
                self.assertEqual((result["inserted_source_start_0"], result["inserted_source_end_0"]), ("100", "103"))
                native = self.read(args.output / "coordinate_crosswalk.tsv")[0]
                self.assertEqual(native["vcf_anchor_pos_1based"], "100")
                self.assertEqual(catalogue.read_fasta(args.output / "alleles.fa")[result["allele_id"] + "|ALT"], forward)

    def test_unmapped_stays_unresolved_and_controls_do_not_inflate_local_counts(self):
        alleles = [dict(allele_id=key, locus_id=key, graph_version="fixture", status=status)
                   for key, status in (("LOCAL", "PASS"), ("CONTROL", "PASS_CONTROL"))]
        self.table("candidate_flow.tsv", [dict(candidate_id="C1", allele_id="LOCAL", status="PASS")])
        empty_paf = self.root / "empty.paf"
        empty_paf.touch()
        searches = [dict(resource=r, alignment_path=str(empty_paf), resource_version="frozen", resource_sha256="a" * 64,
                         available="true", search_complete="true") for r in catalogue.RESOURCES]
        args = NS(catalogue=self.table("alleles.tsv", alleles), queries=self.table("queries.tsv", [], ["query_name"]),
                  search_manifest=self.table("search.tsv", searches), review=self.table("review.tsv", [], ["allele_id", "resource"]),
                  output=self.root / "representation")
        catalogue.representation(args)
        self.assertEqual({r["status"] for r in self.read(args.output / "resource_matrix.tsv")}, {"UNRESOLVED"})
        reviews = [dict(allele_id=a["allele_id"], resource=r, status="NOT_REPRESENTED_AT_RESOLVED_LOCUS",
                        locus_resolved="true", search_adequate="true", structural_configuration_reviewed="true",
                        both_flanks_syntenic="true", comparator_locus_covered="true", evidence=str(self.evidence), reason="reviewed")
                   for a in alleles for r in catalogue.RESOURCES]
        args.review = self.table("review.tsv", reviews)
        catalogue.representation(args)
        flow = {r["metric"]: int(r["count"]) for r in self.read(args.output / "flow.tsv")}
        self.assertEqual(flow["strict_additional_distinct_local_alt_alleles"], 1)
        self.assertEqual(flow["all_attempted_candidate_records"], 1)
        reviews[0]["evidence"] = str(self.root / "missing_resource_review.txt")
        args.review = self.table("review.tsv", reviews)
        with self.assertRaises(ValueError):
            catalogue.representation(args)

    def test_duplicate_molecules_and_overlapping_depth_count_once(self):
        key = dict(locus_id="L1", allele_id="A1", graph_sha256="a" * 64)
        graph = dict(graph_sha256="b" * 64, node_space_id="nodes", reference_build="CHM13")
        base = dict(key, graph_id="main", graph_sha256=graph["graph_sha256"], catalogue_graph_sha256=key["graph_sha256"],
                    node_space_id="nodes", status="VERIFIED", evidence=str(self.evidence))
        targets = [dict(base, feature=f, edge_from=a, edge_to=b) for f, a, b in
                   (("ref", ">1", ">3"), ("left", ">1", ">2"), ("right", ">2", ">3"))]
        targets += [dict(base, feature="flank", node_id="1", node_length=10, start0=a, end0=b)
                    for a, b in ((0, 7), (5, 10))]
        fields = list(dict.fromkeys(k for r in targets for k in r))
        targets = [{k: r.get(k, "") for k in fields} for r in targets]
        mapping = self.root / "mapping"
        mapping.mkdir()
        self.table("mapping/libraries.tsv", [dict(library_id="LIB", gam="unused.gam")])
        alignments = [dict(name=name, mapping_quality=60, path=dict(mapping=[
            dict(position=dict(node_id=node), edit=[dict(from_length=10, to_length=10)]) for node in (1, 2, 3)]))
            for name in ("read/1", "read/2", "duplicate/1")]
        molecules = [dict(person_id="P1", library_id="LIB", read_name=name, molecule_id="M1", status="VERIFIED",
                          method="reviewed_fixture", evidence=str(self.evidence)) for name in ("read", "duplicate")]
        process = NS(stdout=io.StringIO("".join(json.dumps(a) + "\n" for a in alignments)), wait=lambda: 0)
        args = NS(graph_config=json.dumps(graph), samples=self.table("people.tsv", [dict(person_id="P1", analysis_group="independent_primary")]),
                  panel=self.table("panel.tsv", [key]), targets=self.table("targets.tsv", targets), graph="main", person="P1",
                  mapping=mapping, molecules=self.table("molecules.tsv", molecules), anchor_matches=5, min_mapq=20,
                  output=self.root / "extracted.tsv")
        with patch.object(measurement.subprocess, "Popen", return_value=process):
            measurement.extract(args)
        result = self.read(args.output)[0]
        self.assertEqual((result["alt_left_junction_support"], result["alt_right_junction_support"]), ("1", "1"))
        self.assertEqual(float(result["flank_depth"]), 1)
        self.assertEqual(result["dosage"], ".")

    def test_carrier_is_not_diploid_dosage_and_one_sided_support_is_no_call(self):
        base = dict.fromkeys(measurement.EVIDENCE_FIELDS, ".")
        base.update(locus_id="L1", allele_id="A1", reason="UNCALIBRATED", graph_sha256="g", catalogue_graph_sha256="c", node_space_id="nodes",
                    flank_depth=10, mapping_ambiguity=0, ref_support=3, alt_left_junction_support=2,
                    alt_right_junction_support=2, carrier_status="NO_CALL", genotype="./.", dosage=".", callable="no")
        lock = dict(status="LOCKED", graph_sha256="g", catalogue_graph_sha256="c", node_space_id="nodes", min_mapq=20,
                    anchor_matches=5, locked_before_phenotypes="yes", technical_status="PASS", measurement_type="carrier",
                    min_ref_molecules=2, min_alt_left_molecules=2, min_alt_right_molecules=2,
                    min_flank_depth=2, max_mapping_ambiguity=0.1)
        lock.update({k: str(self.evidence) for k in ("calibration_evidence", "threshold_sensitivity_evidence", "depth_batch_review", "technical_review")})
        called = measurement.apply_lock(base, lock, 20, 5)
        self.assertEqual((called["carrier_status"], called["genotype"], called["dosage"]), ("ALT_SUPPORTED", "./.", "."))
        partial = measurement.apply_lock(dict(base, alt_right_junction_support=0), lock, 20, 5)
        self.assertEqual((partial["carrier_status"], partial["dosage"]), ("NO_CALL", "."))
        second_alt = measurement.apply_lock(dict(base, allele_id="A2", alt_left_junction_support=0,
                                                alt_right_junction_support=0, ref_support=0), lock, 20, 5)
        self.assertEqual((second_alt["carrier_status"], second_alt["genotype"], second_alt["dosage"]), ("NO_CALL", "./.", "."))
        unassigned = measurement.apply_lock(dict(base, reason="DUPLICATE_ASSIGNMENTS_INCOMPLETE"), lock, 20, 5)
        self.assertEqual(unassigned["callable"], "no")

    def test_annotation_insertion_boundary_and_tss_distance(self):
        locus = dict(locus_id="L1", allele_id="A1", graph_version="fixture", status="READY",
                     annotation_start_0=100, annotation_end_0=100, inserted_sequence_length="")
        feature = dict(chromosome="chr1", start_0=100, end_0=120, feature_id="EX1", feature_type="exon",
                       gene_id="G1", gene_name="GENE", transcript_id="T1", strand="+", evidence="fixture", build="CHM13", release="frozen")
        overlap = self.root / "overlap.bed"
        overlap.write_text("chr1\t99\t101\tA1|locus\tchr1\t100\t120\tEX1\n")
        empty = self.root / "empty.bed"
        empty.touch()
        sm = NS(params=NS(task="summarize", settings=dict(build="CHM13", release="frozen", assessed_features=["exon"])),
                input=NS(loci=self.table("loci.tsv", [locus]), metadata=self.table("features.tsv", [feature]), overlaps=overlap,
                         nearby=empty, nearest_gene=empty, nearest_tss=empty, links=[], repeatmasker=[]),
                output=NS(context=self.root / "context.tsv", evidence=self.root / "annotation_evidence.tsv", counts=self.root / "counts.tsv"),
                log=[str(self.root / "annotation.log")])
        functions = runpy.run_path(str(SCRIPTS / "svgp_annotation.py"), init_globals={"snakemake": sm})
        self.assertEqual(self.read(sm.output.evidence)[0]["relation"], "INSERTION_AT_FEATURE_BOUNDARY")
        tss = dict(feature, feature_type="tss", start_0=110, end_0=111)
        self.assertEqual(functions["distance"](locus, tss), 10)


if __name__ == "__main__":
    unittest.main()
