"""Small shared helpers for the Minigraph-Cactus build and verifier."""

import hashlib
import json
import os
import tempfile
from pathlib import Path


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_bytes(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def read_json(path):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def identity(preflight_path, settings_text):
    settings = json.loads(settings_text)
    if canonical(settings).decode() != settings_text:
        raise ValueError("settings JSON must be canonical")
    preflight = read_json(preflight_path)
    preflight_hash = sha256(canonical(preflight))
    settings_hash = sha256(settings_text.encode())
    fingerprint = sha256(canonical({"preflight_sha256": preflight_hash,
                                    "settings_sha256": settings_hash}))
    return preflight, settings, preflight_hash, settings_hash, fingerprint


def output_paths(outdir, out_name):
    base = Path(outdir) / out_name
    files = [Path(str(base) + suffix) for suffix in (
        ".gfa.gz", ".gbz", ".vcf.gz",
        ".raw.vcf.gz", ".d2.gbz", ".d2.dist",
        ".d2.shortread.withzip.min", ".d2.shortread.zipcodes",
        ".full.hal", ".stats/graph-stats.tsv",
    )]
    for suffix in (".vcf.gz", ".raw.vcf.gz"):
        tbi = Path(str(base) + suffix + ".tbi")
        csi = Path(str(base) + suffix + ".csi")
        files.append(tbi if tbi.exists() or not csi.exists() else csi)
    return files


def file_record(path, hash_content=False):
    path = Path(path)
    if not path.is_file():
        raise ValueError(f"missing output: {path}")
    stat = path.stat()
    if stat.st_size == 0:
        raise ValueError(f"empty output: {path}")
    row = {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    if hash_content:
        row["sha256"] = file_sha256(path)
    return row


def check_file_record(row):
    path = Path(row["path"])
    if not path.is_file():
        raise ValueError(f"missing output: {path}")
    stat = path.stat()
    if (stat.st_size, stat.st_mtime_ns) != (row["size"], row["mtime_ns"]) or stat.st_size == 0:
        raise ValueError(f"output size or modification time changed: {path}")
    if "sha256" in row and file_sha256(path) != row["sha256"]:
        raise ValueError(f"output content changed: {path}")


def checked_outputs(outdir, out_name, return_paths=False):
    if return_paths:
        return output_paths(outdir, out_name)
    files = []
    for path in output_paths(outdir, out_name):
        files.append(file_record(path))
    return files


def check_input_stats(preflight):
    rows = [preflight["included_list"], *preflight.get("source_manifests", []),
            *preflight["inputs"]]
    if "seed_gfa" in preflight:
        rows.append(preflight["seed_gfa"])
    for row in rows:
        path = Path(row["path"])
        if not path.is_absolute() or not path.is_file():
            raise ValueError(f"input path missing or no longer absolute: {path}; choose a new build ID")
        stat = path.stat()
        if (stat.st_size, stat.st_mtime_ns) != (row["size"], row["mtime_ns"]):
            raise ValueError(f"input changed: {path}; choose a new build ID and rerun preflight")
