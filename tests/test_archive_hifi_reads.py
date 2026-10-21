#!/usr/bin/env python3

import importlib.util
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "workflow" / "scripts" / "archive_hifi_reads.py"
RULES = PROJECT_ROOT / "workflow" / "rules" / "local_sv_graph_v2.smk"
SPEC = importlib.util.spec_from_file_location("archive_hifi_reads", SCRIPT)
ARCHIVE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ARCHIVE)


def copy_args(source, marker, output, log):
    return SimpleNamespace(source=source, marker=marker, output=output, log=log)


class ArchiveBamCopyTests(unittest.TestCase):
    def test_copy_bam_validates_marker_and_publishes_complete_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            source = tmp / "archive" / "sample.bam"
            marker = tmp / "sample.staged"
            output = tmp / "project" / "sample.bam"
            log = tmp / "logs" / "sample.log"
            source.parent.mkdir()
            source.write_bytes(b"BAM contents")
            marker.write_text(
                f"archive_bam\t{source}\n"
                "status\tDUL\n"
                "staged_utc\t2026-07-30T19:16:18+00:00\n"
            )

            ARCHIVE.copy_bam(copy_args(source, marker, output, log))

            self.assertEqual(output.read_bytes(), b"BAM contents")
            self.assertIn(f"SOURCE: {source}", log.read_text())
            self.assertIn("SIZE_BYTES: 12", log.read_text())
            self.assertEqual(list(output.parent.glob(".*.part")), [])

    def test_copy_bam_rejects_marker_for_another_source(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            source = tmp / "sample.bam"
            marker = tmp / "sample.staged"
            output = tmp / "copied.bam"
            log = tmp / "copy.log"
            source.write_bytes(b"BAM")
            marker.write_text(
                "archive_bam\t/archive/another/sample.bam\n"
                "status\tDUL\n"
            )

            with self.assertRaisesRegex(ValueError, "does not match"):
                ARCHIVE.copy_bam(copy_args(source, marker, output, log))

            self.assertFalse(output.exists())

    def test_copy_bam_removes_partial_output_after_copy_failure(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            source = tmp / "sample.bam"
            marker = tmp / "sample.staged"
            output = tmp / "project" / "sample.bam"
            log = tmp / "copy.log"
            source.write_bytes(b"BAM")
            marker.write_text(f"archive_bam\t{source}\nstatus\tREG\n")

            def fail_after_partial_copy(_source, partial):
                Path(partial).write_bytes(b"partial")
                raise OSError("simulated copy failure")

            with mock.patch.object(
                ARCHIVE.shutil, "copyfile", side_effect=fail_after_partial_copy
            ):
                with self.assertRaisesRegex(OSError, "simulated copy failure"):
                    ARCHIVE.copy_bam(copy_args(source, marker, output, log))

            self.assertFalse(output.exists())
            self.assertEqual(list(output.parent.glob(".*.part")), [])


class WorkflowWiringTests(unittest.TestCase):
    def test_archive_bam_is_copied_before_compute_conversion(self):
        rules = RULES.read_text()

        self.assertIn("rule copy_local_v2_hifi_bam:", rules)
        self.assertIn("bam=temp(STAGED_BAM)", rules)
        self.assertIn("partition=HIFI.get(\"archive_partition\", \"archive\")", rules)
        self.assertIn("staged_bam=STAGED_BAM", rules)
        self.assertIn("samtools quickcheck -v {input.staged_bam:q}", rules)
        self.assertIn("samtools fastq", rules)
        self.assertIn(
            'project_work=$(mktemp -d -p "$(dirname {output.reads:q})"',
            rules,
        )
        self.assertIn('> "$project_work/reads.fastq" 2> {log:q}', rules)
        self.assertIn(
            'bgzip -@ {threads} -c "$project_work/reads.fastq"',
            rules,
        )
        self.assertNotIn('| bgzip -@ {threads}', rules)


if __name__ == "__main__":
    unittest.main()
