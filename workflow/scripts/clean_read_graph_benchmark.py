#!/usr/bin/env python3
"""Safely remove completed intermediates from the read/graph benchmark.

The default is a read-only inventory. Pass --delete only after all Snakemake
processes and cluster jobs using this results directory have stopped.
"""

from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


REPOSITORY_ROOT = Path(os.path.abspath(__file__)).parent.parent.parent


@dataclass(frozen=True)
class Candidate:
    path: Path
    category: str
    downstream: tuple[Path, ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Inventory or delete regenerable benchmark intermediates whose "
            "downstream result already exists and is at least as new."
        )
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=REPOSITORY_ROOT / "config/read_graph_benchmark.yaml",
        help="Benchmark config (default: the config next to this workflow).",
    )
    parser.add_argument(
        "--results",
        type=Path,
        help="Override the results path from the config; relative to the repository root.",
    )
    parser.add_argument(
        "--delete",
        action="store_true",
        help="Delete safe candidates. Without this option the script is read-only.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print every candidate instead of category summaries only.",
    )
    return parser.parse_args()


def config_results(config_path: Path) -> str:
    with config_path.open() as handle:
        for raw_line in handle:
            if raw_line.startswith("results:"):
                value = raw_line.split(":", 1)[1].split("#", 1)[0].strip()
                if len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'":
                    value = value[1:-1]
                if value:
                    return value
    raise SystemExit(f"Could not find a top-level results setting in {config_path}")


def config_section_values(
    config_path: Path,
    section: str,
    defaults: dict[str, str],
) -> dict[str, str]:
    """Read simple scalar keys from one top-level YAML section."""
    values = defaults.copy()
    in_section = False
    with config_path.open() as handle:
        for raw_line in handle:
            content = raw_line.split("#", 1)[0].rstrip()
            if not content:
                continue
            indent = len(content) - len(content.lstrip())
            stripped = content.strip()
            if indent == 0:
                in_section = stripped == f"{section}:"
                continue
            if not in_section or indent != 2 or ":" not in stripped:
                continue
            key, value = stripped.split(":", 1)
            if key in values and value.strip():
                scalar = value.strip()
                if len(scalar) >= 2 and scalar[0] == scalar[-1] and scalar[0] in "\"'":
                    scalar = scalar[1:-1]
                values[key] = scalar
    return values


def nonempty_and_not_older(downstream: Path, intermediate: Path) -> bool:
    try:
        return (
            downstream.is_file()
            and downstream.stat().st_size > 0
            and downstream.stat().st_mtime_ns >= intermediate.stat().st_mtime_ns
        )
    except OSError:
        return False


def human_size(size: int) -> str:
    value = float(size)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024 or unit == "TiB":
            return f"{value:.1f} {unit}"
        value /= 1024
    raise AssertionError("unreachable")


def inventory(root: Path) -> tuple[int, int, int]:
    files = directories = size = 0
    for _, dirnames, filenames in os.walk(root, followlinks=False):
        directories += len(dirnames)
        files += len(filenames)
        for filename in filenames:
            path = Path(_) / filename
            try:
                if path.is_file() and not path.is_symlink():
                    size += path.stat().st_size
            except OSError:
                pass
    return files, directories + 1, size


def discover(results: Path, rescue_graphs: set[str]) -> list[Candidate]:
    candidates: list[Candidate] = []

    for path in results.rglob("*.tmp"):
        if path.is_file() and not path.is_symlink():
            candidates.append(Candidate(path, "incomplete .tmp artifact", ()))

    for path in results.glob("preflight/*.ready"):
        if path.is_file() and not path.is_symlink():
            candidates.append(Candidate(path, "preflight gate marker", ()))

    rescue_table = results / "tables/v2_rescue_capture.tsv"
    for path in results.glob("long_reads/alignments/*/*.gaf"):
        relative = path.relative_to(results / "long_reads/alignments")
        alignment_stem = relative.name[:-4]
        metric = (
            results
            / "long_reads/metrics"
            / relative.parent
            / f"{alignment_stem}.tsv"
        )
        downstream = [metric]
        if relative.name.startswith("minigraph."):
            graph = alignment_stem[len("minigraph."):]
            if graph in rescue_graphs:
                downstream.append(rescue_table)
        candidates.append(Candidate(path, "long-read GAF", tuple(downstream)))

    for path in results.glob("v2_rescue/alignments/*/*.gaf"):
        candidates.append(Candidate(path, "legacy duplicate rescue GAF", (rescue_table,)))

    for path in results.glob("simulation/alignments/*/*.gaf"):
        relative = path.relative_to(results / "simulation/alignments")
        alignment_stem = relative.name[:-4]
        accuracy = (
            results
            / "simulation/accuracy"
            / relative.parent
            / f"{alignment_stem}.tsv"
        )
        candidates.append(Candidate(path, "simulated-read GAF", (accuracy,)))

    for path in results.glob("tmem/alignments/*/*.gaf"):
        relative = path.relative_to(results / "tmem/alignments")
        alignment_stem = relative.name[:-4]
        evidence = (
            results
            / "tmem/evidence"
            / relative.parent
            / f"{alignment_stem}.tsv"
        )
        candidates.append(Candidate(path, "TMEM GAF", (evidence,)))

    merged_metrics = results / "tables/mapping_metrics_per_sample.tsv"
    for path in results.glob("long_reads/metrics/*/*.tsv"):
        candidates.append(Candidate(path, "per-job mapping metric", (merged_metrics,)))

    merged_accuracy = results / "tables/simulated_mapping_accuracy.tsv"
    for path in results.glob("simulation/accuracy/*/*.tsv"):
        candidates.append(Candidate(path, "per-job simulation accuracy", (merged_accuracy,)))

    merged_tmem = results / "tables/tmem_genotype_evidence.tsv"
    for path in results.glob("tmem/evidence/*/*.tsv"):
        candidates.append(Candidate(path, "per-job TMEM evidence", (merged_tmem,)))

    for path in results.glob("long_reads/selection/*.read_names.txt"):
        candidates.append(
            Candidate(
                path,
                "selected-read name list",
                (merged_metrics, rescue_table),
            )
        )

    for path in results.glob("simulation/reads/*.simulated.hifi.fastq.gz"):
        candidates.append(Candidate(path, "simulated FASTQ", (merged_accuracy,)))

    return sorted(candidates, key=lambda item: os.fspath(item.path))


def classify(candidates: list[Candidate]) -> tuple[list[Candidate], list[Candidate]]:
    safe = []
    blocked = []
    for candidate in candidates:
        if candidate.path.is_file() and all(
            nonempty_and_not_older(path, candidate.path)
            for path in candidate.downstream
        ):
            safe.append(candidate)
        else:
            blocked.append(candidate)
    return safe, blocked


def summarize(candidates: list[Candidate]) -> dict[str, tuple[int, int]]:
    summary: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    for candidate in candidates:
        try:
            size = candidate.path.stat().st_size
        except OSError:
            size = 0
        summary[candidate.category][0] += 1
        summary[candidate.category][1] += size
    return {key: (value[0], value[1]) for key, value in summary.items()}


def print_summary(label: str, candidates: list[Candidate]) -> None:
    print(label)
    if not candidates:
        print("  none")
        return
    for category, (count, size) in sorted(summarize(candidates).items()):
        print(f"  {category}: {count} files, {human_size(size)}")


def empty_snakemake_temp_directories(repository_root: Path) -> list[Path]:
    snakemake_root = repository_root / ".snakemake"
    if not snakemake_root.is_dir():
        return []
    empty = []
    for path in snakemake_root.glob("tmp.*"):
        try:
            if path.is_dir() and not path.is_symlink() and not any(path.iterdir()):
                empty.append(path)
        except OSError:
            continue
    return sorted(empty)


def prune_empty_intermediate_directories(results: Path) -> int:
    roots = (
        results / "long_reads/alignments",
        results / "long_reads/metrics",
        results / "simulation/alignments",
        results / "simulation/accuracy",
        results / "tmem/alignments",
        results / "tmem/evidence",
        results / "v2_rescue/alignments",
    )
    removed = 0
    for root in roots:
        if not root.is_dir() or root.is_symlink():
            continue
        directories = sorted(
            (path for path in root.rglob("*") if path.is_dir() and not path.is_symlink()),
            key=lambda path: len(path.parts),
            reverse=True,
        )
        directories.append(root)
        for directory in directories:
            try:
                directory.rmdir()
                removed += 1
            except OSError:
                pass
    return removed


def main() -> int:
    args = parse_args()
    config_path = Path(os.path.abspath(args.config))
    if not config_path.is_file():
        raise SystemExit(f"Config does not exist: {config_path}")

    repository_root = config_path.parent.parent
    configured = args.results if args.results is not None else Path(config_results(config_path))
    results = configured if configured.is_absolute() else repository_root / configured
    results = Path(os.path.abspath(results))

    if results == Path("/") or results.is_symlink() or not results.is_dir():
        raise SystemExit(f"Refusing unsafe or missing results directory: {results}")
    expected_children = {"tables", "logs", "long_reads", "tmem", "simulation", "v2_rescue"}
    if not any((results / name).exists() for name in expected_children):
        raise SystemExit(f"Directory does not look like benchmark results: {results}")

    file_count, directory_count, total_size = inventory(results)
    rescue_config = config_section_values(
        config_path,
        "v2_rescue_validation",
        {"v1_graph_id": "local_v1", "v2_graph_id": "local_v2"},
    )
    rescue_graphs = {rescue_config["v1_graph_id"], rescue_config["v2_graph_id"]}
    candidates = discover(results, rescue_graphs)
    safe, blocked = classify(candidates)
    safe_size = sum(candidate.path.stat().st_size for candidate in safe)
    snakemake_temp_dirs = empty_snakemake_temp_directories(repository_root)

    print(f"Results: {results}")
    print(
        f"Inventory: {file_count} files, {directory_count} directories, "
        f"{human_size(total_size)} logical size"
    )
    print_summary("Safe now (downstream output is complete and not older):", safe)
    print_summary("Kept for now (downstream output is missing, empty, or older):", blocked)
    print(f"Empty .snakemake/tmp.* directories safe to remove: {len(snakemake_temp_dirs)}")
    print("Persistent inputs, final tables, reports, logs, benchmarks, and QC files are never selected.")

    if args.verbose:
        for candidate in safe:
            print(f"SAFE\t{candidate.path}")
        for candidate in blocked:
            downstream = ",".join(os.fspath(path) for path in candidate.downstream)
            print(f"KEEP\t{candidate.path}\twaiting_for={downstream}")

    if not args.delete:
        print(
            f"Dry run: {len(safe)} files ({human_size(safe_size)}) and "
            f"{len(snakemake_temp_dirs)} empty directories would be deleted."
        )
        print("Re-run with --delete only when no Snakemake or cluster job is using these files.")
        return 0

    lock_directories = {
        repository_root / ".snakemake/locks",
        Path.cwd() / ".snakemake/locks",
    }
    active_lock_directories = [
        path for path in lock_directories
        if path.is_dir() and any(path.iterdir())
    ]
    if active_lock_directories:
        lock_list = ", ".join(os.fspath(path) for path in sorted(active_lock_directories))
        print(
            f"Refusing deletion while Snakemake lock entries exist in {lock_list}. "
            "Check jobs and use 'snakemake --unlock' only for a stale lock.",
            file=sys.stderr,
        )
        return 2

    failures = []
    for candidate in safe:
        try:
            candidate.path.unlink()
        except OSError as error:
            failures.append((candidate.path, error))
    deleted_files = len(safe) - len(failures)

    removed_snakemake_directories = 0
    for path in snakemake_temp_dirs:
        try:
            path.rmdir()
            removed_snakemake_directories += 1
        except OSError as error:
            failures.append((path, error))

    removed_directories = prune_empty_intermediate_directories(results)
    print(
        f"Deleted {deleted_files} files, freed up to "
        f"{human_size(safe_size)}, and removed "
        f"{removed_directories + removed_snakemake_directories} empty directories."
    )
    for path, error in failures:
        print(f"Could not delete {path}: {error}", file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
