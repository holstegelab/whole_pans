#!/usr/bin/env python3
"""Launch or safely resume one immutable Minigraph-Cactus build."""

import argparse
import fcntl
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from minigraph_cactus_common import (atomic_bytes, canonical, check_input_stats,
                                     checked_outputs, file_sha256, identity, read_json)


def now():
    return datetime.now(timezone.utc).isoformat()


def installed_cactus_version(executable):
    binary = shutil.which(executable)
    if not binary:
        raise ValueError(f"Cactus executable unavailable: {executable}")
    # This workflow pins Cactus in its rule-level Conda environment. Check that
    # the selected launcher comes from that environment and its package record.
    prefix = Path(sys.prefix)
    if Path(os.path.abspath(binary)).parent != prefix / "bin":
        raise ValueError(f"Cactus launcher is outside the active Conda environment: {binary}")
    records = list((prefix / "conda-meta").glob("cactus-*.json"))
    versions = set()
    for record in records:
        package = read_json(record)
        if package.get("name") == "cactus":
            versions.add(package.get("version"))
    if len(versions) != 1:
        raise ValueError("cannot determine Cactus version from active Conda package record")
    return binary, versions.pop()


def write_json(path, value):
    atomic_bytes(path, canonical(value) + b"\n")


def require_directory(path):
    if not path.is_absolute() or not path.is_dir() or not os.access(path, os.W_OK | os.X_OK):
        raise ValueError(f"directory must exist and be writable on the server: {path}")


def completion_logged(path, fingerprint):
    if not path.is_file():
        return False
    with open(path, "rb") as handle:
        handle.seek(0, os.SEEK_END)
        handle.seek(max(0, handle.tell() - 4096))
        tail = handle.read()
    return (f"CACTUS_COMPLETED {fingerprint}\n").encode() in tail


def make_command(binary, args, settings):
    return [binary, str(args.jobstore), str(args.seqfile),
            "--outDir", str(args.outdir), "--outName", args.out_name,
            "--reference", "CHM13", "GRCh38", "--vcfReference", "CHM13",
            "--gfa", "--gbz", "--vcf", "--giraffe", "--validate",
            "--batchSystem", str(settings["batch_system"]),
            "--workDir", str(args.work_dir), "--clean=onSuccess",
            "--maxCores", str(settings["max_cores"]),
            "--maxMemory", str(settings["max_memory"]),
            "--maxDisk", str(settings["max_disk"]),
            "--mgCores", str(settings["mg_cores"]),
            "--mgMemory", str(settings["mg_memory"]),
            "--mapCores", str(settings["map_cores"]),
            "--consCores", str(settings["cons_cores"]),
            "--consMemory", str(settings["cons_memory"]),
            "--indexCores", str(settings["index_cores"]),
            "--indexMemory", str(settings["index_memory"])]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("preflight-json", "seqfile", "settings-json", "outdir", "out-name",
                 "jobstore", "work-dir", "marker", "log"):
        parser.add_argument(f"--{name}", required=True)
    parser.add_argument("--mode", choices=("fresh", "restart"), required=True)
    args = parser.parse_args()
    preflight, settings, preflight_hash, settings_hash, fingerprint = identity(
        args.preflight_json, args.settings_json)
    if "run_mode" in settings:
        raise ValueError("run_mode must be excluded from settings fingerprint")
    for key, value in (("outdir", args.outdir), ("jobstore", args.jobstore),
                       ("work_dir", args.work_dir)):
        if not Path(value).is_absolute() or str(value) != str(settings[key]):
            raise ValueError(f"{key} must match an absolute path in settings")
    if args.out_name != settings["out_name"]:
        raise ValueError("output name does not match settings")
    if file_sha256(args.seqfile) != preflight["seqfile_sha256"]:
        raise ValueError("seqfile differs from preflight; choose a new build ID")
    check_input_stats(preflight)
    require_directory(Path(args.work_dir))
    outdir = Path(args.outdir)
    jobstore = Path(args.jobstore)
    marker = Path(args.marker)
    metadata = marker.parent
    metadata.mkdir(parents=True, exist_ok=True)
    with open(metadata / "build.lock", "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if marker.exists():
            raise ValueError("success marker already exists; request the minigraph_cactus verifier target")
        state_path = metadata / "state.json"
        settings_path = metadata / "settings.json"
        locked_preflight = metadata / "preflight.locked.json"
        state = read_json(state_path) if state_path.exists() else None
        if state and state.get("fingerprint") != fingerprint:
            raise ValueError("build ID has different settings or inputs; choose a new build ID")
        if locked_preflight.exists():
            if read_json(locked_preflight) != preflight:
                raise ValueError("locked preflight differs from current inputs; choose a new build ID")
        elif state:
            raise ValueError("build state exists without its locked preflight; choose a new build ID")
        binary, version = installed_cactus_version(settings["cactus_executable"])
        if version != str(settings["cactus_version"]):
            raise ValueError(f"Cactus version {version} does not match configured {settings['cactus_version']}")
        command = make_command(binary, args, settings)
        if state and state.get("status") == "cactus_completed":
            files = checked_outputs(outdir, args.out_name)
            publish(marker, metadata, state, preflight_hash, settings_hash, files)
            return
        log_path = Path(args.log)
        if state and state.get("status") in {"running", "failed"}:
            if completion_logged(log_path, fingerprint):
                state["status"] = "cactus_completed"
                state["completed_at"] = now()
                write_json(state_path, state)
                files = checked_outputs(outdir, args.out_name)
                publish(marker, metadata, state, preflight_hash, settings_hash, files)
                return
        if args.mode == "fresh":
            if state or jobstore.exists() or (outdir.exists() and any(outdir.iterdir())):
                raise ValueError("fresh build requires no prior state, job store, or occupied output directory")
            require_directory(jobstore.parent)
            require_directory(outdir.parent)
        else:
            if not state or state.get("status") not in {"running", "failed"} or not jobstore.exists():
                raise ValueError("restart requires matching running/failed state and an existing Toil job store")
            command.append("--restart")
        if not locked_preflight.exists():
            atomic_bytes(locked_preflight, canonical(preflight) + b"\n")
        write_json(settings_path, settings)
        state = {"fingerprint": fingerprint, "preflight_sha256": preflight_hash,
                 "settings_sha256": settings_hash, "cactus_version": version,
                 "command": command, "status": "running", "started_at": now(),
                 "mode": args.mode}
        write_json(state_path, state)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "a", encoding="utf-8") as log:
            log.write(f"START {now()} fingerprint={fingerprint}\n")
            log.flush()
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                                    check=False)
            if result.returncode == 0:
                log.write(f"CACTUS_COMPLETED {fingerprint}\n")
                log.flush()
                os.fsync(log.fileno())
        state["returncode"] = result.returncode
        state["completed_at"] = now()
        state["status"] = "cactus_completed" if result.returncode == 0 else "failed"
        write_json(state_path, state)
        if result.returncode:
            raise ValueError(f"Cactus exited {result.returncode}; inspect {log_path} and Toil store")
        files = checked_outputs(outdir, args.out_name)
        publish(marker, metadata, state, preflight_hash, settings_hash, files)


def publish(marker, metadata, state, preflight_hash, settings_hash, files):
    metadata_hashes = {name: file_sha256(metadata / name) for name in (
        "seqfile.txt", "input_manifest.tsv", "contig_map.tsv.gz", "preflight.json",
        "preflight.locked.json", "settings.json")}
    write_json(marker, {"fingerprint": state["fingerprint"],
                        "cactus_version": state["cactus_version"],
                        "preflight_sha256": preflight_hash,
                        "settings_sha256": settings_hash,
                        "metadata_sha256": metadata_hashes,
                        "outputs": files, "completed_at": now()})


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError, json.JSONDecodeError) as exc:
        sys.exit(f"Minigraph-Cactus build refused or failed: {exc}")
