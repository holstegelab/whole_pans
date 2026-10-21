#!/usr/bin/env python3
"""Prepare and stage archived HiFi BAMs for local SV and benchmark workflows."""

import argparse
import csv
import os
import re
import shlex
import shutil
import subprocess
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path


SAFE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fieldnames, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_ids(path):
    ids = set()
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, 1):
            value = line.strip().split()
            if not value or value[0].startswith("#"):
                continue
            sample = value[0]
            if not SAFE_ID.fullmatch(sample):
                raise ValueError(f"Unsafe sample ID at {path}:{line_number}: {sample}")
            if sample in ids:
                raise ValueError(f"Duplicate sample ID in {path}: {sample}")
            ids.add(sample)
    return ids


def make_manifest(args):
    rows = read_tsv(args.assembly_manifest)
    if not rows or "sample_id" not in rows[0]:
        raise ValueError(f"No sample_id column in {args.assembly_manifest}")

    samples = sorted({row["sample_id"].strip() for row in rows})
    if not samples or "" in samples:
        raise ValueError(f"Missing sample IDs in {args.assembly_manifest}")

    centenarians = read_ids(args.centenarian_samples)
    manifest = []
    for sample in samples:
        if not SAFE_ID.fullmatch(sample):
            raise ValueError(f"Unsafe sample ID: {sample}")
        cohort = "centenarian" if sample in centenarians else "ad"
        archive_bam = os.path.normpath(
            args.archive_pattern.format(
                archive_root=args.archive_root,
                cohort=cohort,
                sample=sample,
            )
        )
        manifest.append(
            {
                "sample_id": sample,
                "cohort": cohort,
                "archive_bam": archive_bam,
                "hifi_path": args.fastq_pattern.format(sample=sample),
            }
        )

    write_tsv(
        args.output,
        ["sample_id", "cohort", "archive_bam", "hifi_path"],
        manifest,
    )


def run_logged(command, log):
    log.write("COMMAND: " + shlex.join(command) + "\n")
    log.flush()
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    log.write((result.stdout or "") + (result.stderr or ""))
    log.flush()
    return result


def chunks(values, size=100):
    for start in range(0, len(values), size):
        yield values[start : start + size]


def stage_bams(args):
    planned_samples = set(args.sample)
    if args.plan:
        planned_samples.update(
            row["sample_id"]
            for row in read_tsv(args.plan)
            if row.get("plan_status") == "READY"
        )
    if not planned_samples:
        raise ValueError("No READY rescues or explicit samples were requested")

    manifest = {row["sample_id"]: row for row in read_tsv(args.read_manifest)}
    missing = sorted(planned_samples - manifest.keys())
    if missing:
        raise ValueError("Samples missing from read manifest: " + ", ".join(missing))

    bam_to_sample = {
        manifest[sample]["archive_bam"]: sample for sample in planned_samples
    }
    archive_bams = sorted(bam_to_sample)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = Path(args.log)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    online = {}
    with log_path.open("w") as log:
        for batch in chunks(archive_bams):
            result = run_logged([args.daget, "-av", *batch], log)
            if result.returncode:
                raise RuntimeError(f"daget failed; see {log_path}")

        deadline = time.monotonic() + args.max_wait_minutes * 60
        while len(online) < len(archive_bams):
            pending = [path for path in archive_bams if path not in online]
            for batch in chunks(pending):
                result = run_logged([args.dals, "-l", *batch], log)
                text = (result.stdout or "") + (result.stderr or "")
                for line in text.splitlines():
                    match = re.search(r"\(([A-Z]{3})\)\s+(.*)$", line)
                    if not match:
                        continue
                    path = "/" + os.path.normpath(match.group(2)).lstrip("/")
                    if path in bam_to_sample and match.group(1) in {"DUL", "REG"}:
                        online[path] = match.group(1)

            log.write(f"ONLINE {len(online)}/{len(archive_bams)}\n")
            log.flush()
            if len(online) == len(archive_bams):
                break
            if time.monotonic() >= deadline:
                missing = [path for path in archive_bams if path not in online]
                raise TimeoutError(
                    f"Timed out waiting for {len(missing)} archived BAM(s); "
                    f"see {log_path}"
                )
            time.sleep(args.poll_seconds)

    staged_rows = []
    for archive_bam in archive_bams:
        if not Path(archive_bam).is_file():
            raise FileNotFoundError(f"Staged BAM is not visible: {archive_bam}")
        sample = bam_to_sample[archive_bam]
        staged_at = datetime.now(timezone.utc).isoformat()
        staged_rows.append(
            {
                "sample_id": sample,
                "archive_bam": archive_bam,
                "status": online[archive_bam],
                "staged_utc": staged_at,
            }
        )
        (output_dir / f"{sample}.staged").write_text(
            f"archive_bam\t{archive_bam}\n"
            f"status\t{online[archive_bam]}\n"
            f"staged_utc\t{staged_at}\n"
        )

    write_tsv(
        output_dir / "staged_bams.tsv",
        ["sample_id", "archive_bam", "status", "staged_utc"],
        staged_rows,
    )


def read_stage_marker(path):
    values = {}
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t", 1)
            if len(fields) != 2 or not fields[0]:
                raise ValueError(f"Malformed stage marker at {path}:{line_number}")
            values[fields[0]] = fields[1]
    return values


def copy_bam(args):
    """Copy one recalled BAM from /archive to compute-visible storage."""
    source = Path(args.source)
    marker = Path(args.marker)
    output = Path(args.output)
    log_path = Path(args.log)
    output.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    partial = None

    try:
        marker_values = read_stage_marker(marker)
        expected_source = marker_values.get("archive_bam")
        if expected_source is None:
            raise ValueError(f"No archive_bam entry in stage marker: {marker}")
        if os.path.normpath(expected_source) != os.path.normpath(str(source)):
            raise ValueError(
                f"Stage marker source does not match requested source: "
                f"{expected_source} != {source}"
            )
        status = marker_values.get("status")
        if status not in {"DUL", "REG"}:
            raise ValueError(f"Archive BAM is not online in {marker}: {status!r}")
        if not source.is_file():
            raise FileNotFoundError(f"Staged BAM is not visible: {source}")
        if not os.access(source, os.R_OK):
            raise PermissionError(f"Staged BAM is not readable: {source}")

        source_size = source.stat().st_size
        if source_size == 0:
            raise ValueError(f"Staged BAM is empty: {source}")

        with log_path.open("w") as log:
            log.write(f"SOURCE: {source}\n")
            log.write(f"DESTINATION: {output}\n")
            log.write(f"ARCHIVE_STATUS: {status}\n")
            log.write(f"SIZE_BYTES: {source_size}\n")
            log.flush()

            descriptor, partial_name = tempfile.mkstemp(
                dir=output.parent,
                prefix=f".{output.name}.",
                suffix=".part",
            )
            os.close(descriptor)
            partial = Path(partial_name)
            shutil.copyfile(source, partial)

            copied_size = partial.stat().st_size
            if copied_size != source_size:
                raise OSError(
                    f"Incomplete BAM copy: expected {source_size} bytes, "
                    f"copied {copied_size} bytes"
                )
            os.replace(partial, output)
            partial = None
            log.write("COPY_COMPLETE\n")
    except Exception as error:
        if partial is not None:
            partial.unlink(missing_ok=True)
        with log_path.open("a") as log:
            log.write(f"ERROR: {error}\n")
        raise


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    manifest = subparsers.add_parser("manifest")
    manifest.add_argument("--assembly-manifest", required=True)
    manifest.add_argument("--centenarian-samples", required=True)
    manifest.add_argument("--archive-root", required=True)
    manifest.add_argument("--archive-pattern", required=True)
    manifest.add_argument("--fastq-pattern", required=True)
    manifest.add_argument("--output", required=True)
    manifest.set_defaults(func=make_manifest)

    stage = subparsers.add_parser("stage")
    stage.add_argument("--plan")
    stage.add_argument("--sample", action="append", default=[])
    stage.add_argument("--read-manifest", required=True)
    stage.add_argument("--output-dir", required=True)
    stage.add_argument("--daget", required=True)
    stage.add_argument("--dals", required=True)
    stage.add_argument("--poll-seconds", type=int, default=30)
    stage.add_argument("--max-wait-minutes", type=int, default=1440)
    stage.add_argument("--log", required=True)
    stage.set_defaults(func=stage_bams)

    copy = subparsers.add_parser("copy")
    copy.add_argument("--source", required=True)
    copy.add_argument("--marker", required=True)
    copy.add_argument("--output", required=True)
    copy.add_argument("--log", required=True)
    copy.set_defaults(func=copy_bam)

    return parser.parse_args()


def main():
    args = parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
