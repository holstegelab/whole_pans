#!/usr/bin/env python3
"""Add validated rescue sequences to a frozen Minigraph rGFA."""

import argparse
import csv
import hashlib
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_metrics(path):
    with Path(path).open(newline="") as handle:
        return {
            row["metric"]: row["value"]
            for row in csv.DictReader(handle, delimiter="\t")
        }


def count_rows(path):
    with Path(path).open() as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def count_fasta(path):
    with Path(path).open() as handle:
        return sum(line.startswith(">") for line in handle)


def graph_metrics(path):
    metrics = {
        "segments": 0,
        "links": 0,
        "paths": 0,
        "segment_bp": 0,
        "rescue_named_segments": 0,
    }
    with Path(path).open("rb") as handle:
        for line in handle:
            if line.startswith(b"S\t"):
                fields = line.rstrip().split(b"\t")
                metrics["segments"] += 1
                metrics["segment_bp"] += len(fields[2])
                metrics["rescue_named_segments"] += any(
                    field.startswith(b"SN:Z:rescue_") for field in fields[3:]
                )
            elif line.startswith(b"L\t"):
                metrics["links"] += 1
            elif line.startswith((b"P\t", b"W\t")):
                metrics["paths"] += 1
    return metrics


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


def version(command):
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    text = (result.stdout or "") + (result.stderr or "")
    return next((line for line in text.splitlines() if line.strip()), "unavailable")


def build(args):
    source = read_metrics(args.frozen_summary)
    inventory = read_metrics(args.frozen_inventory)
    observed_source_hash = sha256(args.frozen_graph)
    if inventory.get("sha256") != observed_source_hash:
        raise ValueError("Frozen graph checksum does not match its inventory")

    rescue_count = count_rows(args.graph_inputs)
    if rescue_count < 1 or rescue_count != count_fasta(args.rescue_fasta):
        raise ValueError("Rescue FASTA and graph-input table are empty or disagree")

    output_graph = Path(args.output_graph)
    output_graph.parent.mkdir(parents=True, exist_ok=True)
    Path(args.log).parent.mkdir(parents=True, exist_ok=True)

    scratch_root = Path(args.scratch_dir or "/tmp")
    with tempfile.TemporaryDirectory(
        prefix="local-sv-graph-v2.",
        dir=scratch_root,
    ) as work:
        temporary_graph = Path(work) / output_graph.name
        command = [
            "minigraph",
            f"-cx{args.preset}",
            "-t",
            str(args.threads),
            args.frozen_graph,
            args.rescue_fasta,
        ]
        with Path(args.log).open("w") as log, temporary_graph.open("w") as output:
            log.write("COMMAND: " + " ".join(command) + "\n")
            result = subprocess.run(command, stdout=output, stderr=log, check=False)
        if result.returncode or not temporary_graph.stat().st_size:
            raise RuntimeError(f"Minigraph failed; see {args.log}")
        log_text = Path(args.log).read_text(errors="replace")
        if "inconsistent rGFA" in log_text or "associated with different ranks" in log_text:
            raise RuntimeError(f"Minigraph reported inconsistent rGFA; see {args.log}")

        final = graph_metrics(temporary_graph)
        if (
            final["segments"] == int(source["segments"])
            and final["links"] == int(source["links"])
        ):
            raise ValueError("Minigraph added no new graph structure")
        shutil.move(temporary_graph, output_graph)

    summary_rows = [
        {"metric": f"source_{key}", "value": source[key]}
        for key in ("segments", "links", "paths", "segment_bp")
    ]
    summary_rows += [
        {"metric": f"final_{key}", "value": final[key]}
        for key in ("segments", "links", "paths", "segment_bp")
    ]
    summary_rows += [
        {
            "metric": f"delta_{key}",
            "value": final[key] - int(source[key]),
        }
        for key in ("segments", "links", "paths", "segment_bp")
    ]
    summary_rows += [
        {"metric": "unique_pass_rescues", "value": rescue_count},
        {
            "metric": "rescue_fasta_bp",
            "value": sum(
                len(line.strip())
                for line in Path(args.rescue_fasta).open()
                if not line.startswith(">")
            ),
        },
        {
            "metric": "rescue_named_segments",
            "value": final["rescue_named_segments"],
        },
    ]
    write_tsv(args.output_summary, ["metric", "value"], summary_rows)

    write_tsv(
        args.output_versions,
        ["tool", "version"],
        [
            {"tool": "minigraph", "version": version(["minigraph", "--version"])},
            {"tool": "gfatools", "version": version(["gfatools"])},
        ],
    )

    provenance_files = {
        "source_graph": args.frozen_graph,
        "source_graph_summary": args.frozen_summary,
        "source_graph_inventory": args.frozen_inventory,
        "rescue_fasta": args.rescue_fasta,
        "graph_inputs": args.graph_inputs,
        "rescue_summary": args.rescue_summary,
        "rescue_plan": args.plan,
        "hifi_read_manifest": args.read_manifest,
        "hifi_header_template": args.header_template,
        "graph_conda_environment": args.graph_environment,
    }
    provenance = [
        {"key": "graph_name", "value": "local_sv_graph_v2"},
        {
            "key": "created_utc",
            "value": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        },
        {"key": "minigraph_preset", "value": args.preset},
    ]
    hashes = {"source_graph": observed_source_hash}
    for key, path in provenance_files.items():
        provenance.append({"key": key, "value": path})
        if key not in hashes:
            hashes[key] = sha256(path)
        provenance.append({"key": f"{key}_sha256", "value": hashes[key]})
    provenance += [
        {"key": "final_graph", "value": str(output_graph)},
        {"key": "final_graph_sha256", "value": sha256(output_graph)},
    ]
    write_tsv(args.output_provenance, ["key", "value"], provenance)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in (
        "frozen_graph",
        "frozen_summary",
        "frozen_inventory",
        "rescue_fasta",
        "graph_inputs",
        "rescue_summary",
        "plan",
        "read_manifest",
        "header_template",
        "graph_environment",
        "output_graph",
        "output_summary",
        "output_versions",
        "output_provenance",
        "log",
    ):
        parser.add_argument(f"--{name.replace('_', '-')}", required=True)
    parser.add_argument("--preset", default="ggs")
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--scratch-dir")
    return parser.parse_args()


def main():
    build(parse_args())


if __name__ == "__main__":
    main()
