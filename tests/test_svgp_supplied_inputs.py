"""Tiny server-only fixtures for supplied clinical and existing-caller imports."""

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace as NS

from openpyxl import Workbook

SCRIPTS = Path(__file__).absolute().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))
import svgp_expansionhunter as eh
import svgp_gcnv as gcnv
import svgp_inventory as inventory
import svgp_phenotypes as phenotype


class SuppliedInputs(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def table(self, name, rows, fields):
        path = self.root / name
        inventory.write_table(path, rows, fields)
        return path

    def read(self, path):
        with Path(path).open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def workbook_args(self, age=72):
        workbook = Workbook()
        sheet = workbook.active
        sheet.title = "Blad1"
        sheet.append(phenotype.SOURCE_FIELDS)
        # Repeated source aliases with conflicting clinical values are not people.
        sheet.append(["SEQ1", "WGS1", "BIO1", "F", age, "svPPA", "detail1", "blood"])
        sheet.append(["SEQ1", "WGS2", "BIO1", "M", "NA", "control", "detail2", "brain"])
        path = self.root / "clinical.xlsx"
        workbook.save(path)
        workbook.close()
        sample = dict.fromkeys(inventory.SAMPLE_FIELDS, "unknown")
        sample.update(person_id="P1", age_definition="unknown")
        return NS(workbook=path, sheet="Blad1", id_column="Sample_ID_WGS_given",
                  phenotype_column="Diagnosis_svPPA",
                  crosswalk=self.table("phenotype_crosswalk.tsv", [{"Sample_ID_WGS_given": "WGS1", "person_id": "P1"}],
                                       ["Sample_ID_WGS_given", "person_id"]),
                  samples=self.table("samples.tsv", [sample], inventory.SAMPLE_FIELDS), output=self.root / "metadata")

    def crosswalk(self, source, label="SAMPLE.markdup"):
        return self.table("caller_crosswalk.tsv", [dict(source_id=source, caller_sample_id=label,
                                                       person_id="P1", review_status="REVIEWED")],
                          ["source_id", "caller_sample_id", "person_id", "review_status"])

    def test_workbook_rows_are_retained_and_sampling_age_is_joined(self):
        args = self.workbook_args()
        phenotype.import_workbook(args)
        rows = self.read(args.output / "phenotypes.tsv")
        self.assertEqual(len(rows), 2)
        self.assertEqual([row["join_status"] for row in rows], ["MAPPED", "UNMAPPED"])
        self.assertEqual((rows[1]["source_age"], rows[1]["age_at_sampling"]), ("NA", "unknown"))
        sample = self.read(args.output / "samples_with_phenotypes.tsv")[0]
        self.assertEqual((sample["phenotype"], sample["sex"], sample["source_age"]), ("svPPA", "F", "72"))
        self.assertEqual((sample["age_at_sampling"], sample["age_definition"]), ("72", "age_at_sampling"))
        self.assertEqual((sample["diagnosis_svppa"], sample["diagnosis_complete"]), ("svPPA", "detail1"))

    def test_missing_sampling_age_stays_missing_and_unmapped_people_are_unchanged(self):
        args = self.workbook_args(age="NA")
        samples = self.read(args.samples)
        unrelated = dict.fromkeys(inventory.SAMPLE_FIELDS, "unknown")
        unrelated.update(person_id="P2", age_definition="age_at_onset", age_at_sampling="55")
        inventory.write_table(args.samples, samples + [unrelated], inventory.SAMPLE_FIELDS + ["age_at_sampling"])
        phenotype.import_workbook(args)
        matched, unmapped = self.read(args.output / "samples_with_phenotypes.tsv")
        self.assertEqual((matched["age_at_sampling"], matched["age_definition"]), ("unknown", "age_at_sampling"))
        self.assertEqual((unmapped["age_at_sampling"], unmapped["age_definition"]), ("55", "age_at_onset"))

    def test_known_curated_age_is_preserved_when_equal_or_source_age_is_missing(self):
        for source_age in (72, "NA", None):
            with self.subTest(source_age=source_age):
                args = self.workbook_args(age=source_age)
                sample = self.read(args.samples)[0]
                sample.update(age_at_sampling="72.0", age_definition="age_at_sampling")
                inventory.write_table(args.samples, [sample], inventory.SAMPLE_FIELDS + ["age_at_sampling"])
                phenotype.import_workbook(args)
                imported = self.read(args.output / "samples_with_phenotypes.tsv")[0]
                self.assertEqual(imported["age_at_sampling"], "72.0")
                self.assertEqual(imported["source_age"], "" if source_age is None else str(source_age))

    def test_conflicting_curated_age_or_definition_requires_review(self):
        for field, value in (("age_at_sampling", "73"), ("age_definition", "age_at_onset")):
            with self.subTest(field=field):
                args = self.workbook_args()
                sample = self.read(args.samples)[0]
                sample[field] = value
                inventory.write_table(args.samples, [sample], inventory.SAMPLE_FIELDS + ["age_at_sampling"])
                with self.assertRaisesRegex(ValueError, f"{field} conflicts"):
                    phenotype.import_workbook(args)

    def test_invalid_source_or_curated_sampling_age_is_rejected(self):
        for value in ("-1", "not an age", "NaN", "Infinity"):
            for source in (True, False):
                with self.subTest(value=value, source=source):
                    args = self.workbook_args(age=value if source else 72)
                    if not source:
                        sample = self.read(args.samples)[0]
                        sample["age_at_sampling"] = value
                        inventory.write_table(args.samples, [sample], inventory.SAMPLE_FIELDS + ["age_at_sampling"])
                    with self.assertRaisesRegex(ValueError, "non-negative finite number"):
                        phenotype.import_workbook(args)

    def test_workbook_does_not_silently_recode_conflicting_sex(self):
        args = self.workbook_args()
        sample = self.read(args.samples)[0]
        sample["sex"] = "female"
        inventory.write_table(args.samples, [sample], inventory.SAMPLE_FIELDS)
        with self.assertRaisesRegex(ValueError, "sex conflicts"):
            phenotype.import_workbook(args)

    def test_workbook_crosswalk_cannot_duplicate_people(self):
        args = self.workbook_args()
        inventory.write_table(args.crosswalk, [{"Sample_ID_WGS_given": key, "person_id": "P1"}
                                               for key in ("WGS1", "WGS2")],
                              ["Sample_ID_WGS_given", "person_id"])
        with self.assertRaises(ValueError):
            phenotype.import_workbook(args)

    def test_helpers_refuse_outputs_that_alias_inputs(self):
        source = self.root / "source.txt"
        source.write_text("preserve\n")
        alias = self.root / "alias.txt"
        alias.symlink_to(source)
        with self.assertRaises(ValueError):
            inventory.protect_inputs([source], [alias])
        with self.assertRaises(ValueError):
            inventory.protect_inputs([source], [self.root / "nested" / ".." / "source.txt"])
        args = self.workbook_args()
        args.samples = args.output / "samples_with_phenotypes.tsv"
        with self.assertRaises(ValueError):
            phenotype.import_workbook(args)
        with self.assertRaises(ValueError):
            eh.prepare_catalog(NS(catalog=source, loci=alias, output=source))
        with self.assertRaises(ValueError):
            eh.extract(NS(vcf=source, catalog=alias, crosswalk=source, run_log="",
                          out_evidence=alias, out_status=self.root / "status.tsv"))
        with self.assertRaises(ValueError):
            gcnv.import_vcf(NS(vcf=source, crosswalk=alias, out_evidence=source,
                               out_status=self.root / "status.tsv"))
        self.assertEqual(source.read_text(), "preserve\n")

    def test_gcnv_keeps_reported_gain_with_missing_gt_and_unassessed_filter(self):
        vcf = self.root / "segments.vcf"
        vcf.write_text("##fileformat=VCFv4.2\n##FORMAT=<ID=NP,Number=1,Type=Integer,Description=Bins>\n"
                       "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
                       "chr1\t101\tgain\tN\t<DUP>\t.\t.\tEND=200\tGT:CN:NP:QS\t./.:3:4:60\n"
                       "chrX\t301\tref\tN\t.\t.\tPASS\tEND=400\tGT:CN:NP\t0:1:2\n"
                       "chr2\t501\tmissing\tN\t<DEL>\t.\t.\tEND=600\tGT:CN:NP\t./.:.:2\n")
        args = NS(vcf=vcf, crosswalk=self.crosswalk("gcnv_standard", "SAMPLE"), source_id="gcnv_standard",
                  expected_sample="SAMPLE", reference_build="custom_GRCh38",
                  out_evidence=self.root / "cnv.tsv", out_status=self.root / "cnv_status.tsv")
        gcnv.import_vcf(args)
        rows = self.read(args.out_evidence)
        self.assertEqual((rows[0]["reported_cn"], rows[0]["raw_gt"], rows[0]["filter"]), ("3", "./.", "."))
        self.assertEqual((rows[0]["start0"], rows[0]["end0"]), ("100", "200"))
        self.assertEqual([row["reported_state"] for row in rows],
                         ["REPORTED_CNV", "REPORTED_REFERENCE_SEGMENT", "MISSING_CN"])
        self.assertEqual(self.read(args.out_status)[0]["records"], "3")

    def repeat_catalog(self):
        catalog = self.root / "catalog.json"
        catalog.write_text(json.dumps([dict(LocusId="LOC", ReferenceRegion=["chr1:100-109", "chr1:110-119"],
                                            VariantType=["Repeat", "Repeat"], LocusStructure="(CAG)*A(CCG)*")]))
        loci = self.table("loci.tsv", [dict(locus_id="LOC", variant_id="LOC_chr1:100-109")], ["locus_id", "variant_id"])
        output = self.root / "selected.tsv"
        eh.prepare_catalog(NS(catalog=catalog, loci=loci, output=output))
        return output

    def test_repeat_catalog_selection_and_empty_selection(self):
        prepared = self.repeat_catalog()
        self.assertEqual(self.read(prepared)[0]["start0"], "100")
        empty = self.table("empty.tsv", [], ["locus_id", "variant_id"])
        with self.assertRaisesRegex(ValueError, "prespecified"):
            eh.prepare_catalog(NS(catalog=self.root / "catalog.json", loci=empty, output=self.root / "empty_selected.tsv"))

    def test_repeat_anchor_partial_uncertainty_and_lowdepth(self):
        metadata = self.read(self.repeat_catalog())[0]
        fields = ["chr1", "100", ".", "C", "<STR3>,<STR4>", ".", "PASS", "END=109;RU=CAG;REF=3;VARID=LOC_chr1:100-109",
                  "GT:REPCN:REPCI:LC", "0/2:3/4:3-3/4-5:20"]
        result = eh.parse_call(fields, metadata)
        self.assertEqual((result["call_status"], result["allele2_repeat_units"]), ("CALLER_PASS", "4"))
        fields[1] = "99"
        self.assertEqual(eh.parse_call(fields, metadata)["call_status"], "CATALOG_COORDINATE_MISMATCH")
        fields[1], fields[9] = "100", "0/2:3/4:.:20"
        self.assertEqual(eh.parse_call(fields, metadata)["call_status"], "PARTIAL_GENOTYPE")
        fields[6] = "LowDepth"
        self.assertEqual(eh.parse_call(fields, metadata)["call_status"], "LOW_DEPTH")

    def test_repeat_header_uses_exact_label_including_suffix_and_crlf(self):
        prepared = self.repeat_catalog()
        vcf = self.root / "sample.vcf"
        vcf.write_bytes(("##fileformat=VCFv4.2\r\n"
                         "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE.markdup\r\n"
                         "chr1\t100\t.\tC\t<STR4>\t.\tPASS\tEND=109;RU=CAG;VARID=LOC_chr1:100-109\tGT:REPCN:REPCI\t0/1:3/4:3-3/4-4\r\n").encode())
        args = NS(vcf=vcf, catalog=prepared, crosswalk=self.crosswalk("expansionhunter"), source_id="expansionhunter",
                  run_log="", reference_build="custom_GRCh38", reference_fasta="unknown", reference_sha256="unknown",
                  out_evidence=self.root / "repeat.tsv", out_status=self.root / "repeat_status.tsv")
        eh.extract(args)
        row = self.read(args.out_evidence)[0]
        self.assertEqual((row["caller_sample_id"], row["person_id"], row["review_status"]),
                         ("SAMPLE.markdup", "P1", "NOT_VALIDATED"))
        args.crosswalk = self.crosswalk("expansionhunter", "SAMPLE")
        eh.extract(args)
        self.assertEqual(self.read(args.out_evidence)[0]["identity_status"], "UNRESOLVED")


if __name__ == "__main__":
    unittest.main()
