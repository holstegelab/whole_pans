#!/usr/bin/env python3
"""Check that a staged Cactus build still uses its prepared inputs and settings."""

import argparse
import json
import sys
from pathlib import Path

from minigraph_cactus_common import atomic_bytes, check_input_stats, file_sha256


def check_lock(path, content, initialize):
    if path.exists():
        if path.read_bytes() != content:
            raise ValueError(f"build settings or inputs changed: {path}; choose a new build ID")
    elif initialize:
        atomic_bytes(path, content)
    else:
        raise ValueError(f"missing build lock: {path}")


def check_graphmap_seqfile(original_path, graphmap_path, fasta_path):
    from cactus.progressive.seqFile import SeqFile
    from cactus.shared.common import makeURL

    original = SeqFile(str(original_path))
    graphmap = SeqFile(str(graphmap_path))
    graph_event = "_MINIGRAPH_"
    expected = set(original.pathMap) | {graph_event}
    if (graph_event in original.pathMap or set(graphmap.pathMap) != expected or
            any(graphmap.pathMap[name] != path for name, path in original.pathMap.items())):
        raise ValueError("graphmap changed an original genome path")
    if graphmap.pathMap[graph_event] not in {str(fasta_path), makeURL(str(fasta_path))}:
        raise ValueError("_MINIGRAPH_ does not point to the graph FASTA")
    leaves = [graphmap.tree.getName(node) for node in graphmap.tree.getLeaves()]
    if len(leaves) != len(expected) or set(leaves) != expected:
        raise ValueError("graphmap seqfile has missing or duplicate tree leaves")
    rows = sum(line.split()[0] == graph_event
               for line in graphmap_path.read_text(encoding="utf-8").splitlines()
               if line.split())
    if rows != 1:
        raise ValueError("graphmap seqfile needs exactly one _MINIGRAPH_ row")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight-json", type=Path, required=True)
    parser.add_argument("--seqfile", type=Path, required=True)
    parser.add_argument("--settings-json", required=True)
    parser.add_argument("--metadata-dir", type=Path, required=True)
    parser.add_argument("--initialize", action="store_true")
    parser.add_argument("--graphmap-copy", type=Path)
    parser.add_argument("--graph-fasta", type=Path)
    args = parser.parse_args()

    preflight_bytes = args.preflight_json.read_bytes()
    preflight = json.loads(preflight_bytes)
    settings = json.loads(args.settings_json)
    versions = set()
    for path in (Path(sys.prefix) / "conda-meta").glob("cactus-*.json"):
        package = json.loads(path.read_text(encoding="utf-8"))
        if package.get("name") == "cactus":
            versions.add(package["version"])
    if versions != {str(settings["minigraph_cactus"]["cactus_version"])}:
        raise ValueError("active Conda Cactus version differs from config")
    check_input_stats(preflight)
    if file_sha256(args.seqfile) != preflight["seqfile_sha256"]:
        raise ValueError("original seqfile changed since preflight")
    check_lock(args.metadata_dir / "preflight.locked.json", preflight_bytes, args.initialize)
    check_lock(args.metadata_dir / "settings.locked.json",
               (args.settings_json + "\n").encode(), args.initialize)
    if args.graphmap_copy:
        if not args.graph_fasta or not args.graph_fasta.is_file() or args.graph_fasta.stat().st_size == 0:
            raise ValueError("graphmap FASTA is missing")
        check_graphmap_seqfile(args.seqfile, args.graphmap_copy, args.graph_fasta)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError, json.JSONDecodeError) as exc:
        sys.exit(f"Minigraph-Cactus input check failed: {exc}")
