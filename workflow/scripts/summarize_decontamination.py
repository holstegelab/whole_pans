#!/usr/bin/env python3
"""Create cohort-level contamination reports and graph input list."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--actions", required=True)
    parser.add_argument("--review", required=True)
    parser.add_argument("--graph-list", required=True)
    parser.add_argument("--complete-marker", required=True)
    return parser.parse_args()


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    manifest = read_tsv(args.manifest)
    results = Path(args.results_dir)
    summaries = []
    actions = []
    review = []
    graph_paths = []

    for entry in manifest:
        assembly_id = entry["assembly_id"]
        cleaned_path = results / "cleaned" / f"{assembly_id}.clean.fa.gz"
        rows = read_tsv(results / "decisions" / f"{assembly_id}.contigs.tsv")
        counts = defaultdict(int)
        bases = defaultdict(int)
        for row in rows:
            decision = row["decision"]
            counts[decision] += 1
            bases[decision] += int(row["contig_length_bp"])
            if decision != "KEEP":
                action = dict(row)
                action["source_path"] = entry["source_path"]
                action["cleaned_path"] = str(cleaned_path)
                actions.append(action)
                if decision == "REVIEW":
                    review.append(action)

        summaries.append(
            {
                "assembly_id": assembly_id,
                "source_path": entry["source_path"],
                "cleaned_path": str(cleaned_path),
                "total_contigs": len(rows),
                "keep_contigs": counts["KEEP"],
                "review_contigs": counts["REVIEW"],
                "split_contigs": counts["SPLIT"],
                "removed_contigs": counts["REMOVE"],
                "removed_whole_contig_bp": bases["REMOVE"],
                "split_nonhuman_bp": sum(
                    int(row["nonhuman_covered_bp"])
                    for row in rows
                    if row["decision"] == "SPLIT"
                ),
                "review_bp": bases["REVIEW"],
            }
        )
        graph_paths.append(str(cleaned_path))

    summary_fields = [
        "assembly_id",
        "source_path",
        "cleaned_path",
        "total_contigs",
        "keep_contigs",
        "review_contigs",
        "split_contigs",
        "removed_contigs",
        "removed_whole_contig_bp",
        "split_nonhuman_bp",
        "review_bp",
    ]
    action_fields = list(actions[0]) if actions else [
        "assembly_id",
        "contig",
        "contig_length_bp",
        "decision",
        "reason",
        "source_path",
        "cleaned_path",
    ]
    write_tsv(args.summary, summaries, summary_fields)
    write_tsv(args.actions, actions, action_fields)
    write_tsv(args.review, review, action_fields)
    Path(args.graph_list).parent.mkdir(parents=True, exist_ok=True)
    with open(args.graph_list, "w") as handle:
        for path in graph_paths:
            handle.write(path + "\n")
    Path(args.complete_marker).touch()
    print(
        f"Assemblies={len(summaries)} actions={len(actions)} review={len(review)}"
    )


if __name__ == "__main__":
    main()
