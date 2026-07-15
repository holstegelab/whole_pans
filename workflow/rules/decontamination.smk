DECONTAM_CONFIG = workflow.source_path("../../config/config.yaml")
configfile: DECONTAM_CONFIG

DECONTAM_ENTRY_OUTDIR = config["results"]["decontamination"]


rule decontamination:
    input:
        f"{DECONTAM_ENTRY_OUTDIR}/summary/contamination_summary.tsv",
        f"{DECONTAM_ENTRY_OUTDIR}/summary/contig_actions.tsv",
        f"{DECONTAM_ENTRY_OUTDIR}/summary/review_candidates.tsv",
        f"{DECONTAM_ENTRY_OUTDIR}/summary/graph_cleaned_assemblies.txt",
        complete=lambda wildcards: DECONTAM_COMPLETE_MARKER,


include: "QC.smk"


import csv
import gzip
import hashlib
import json
import os
import sys
from collections import defaultdict
from pathlib import Path


QC_RESULTS_DIR = config["results"]["qc"]
DECONTAM_OUTDIR = config["results"]["decontamination"]
DECONTAM_ENV = f"{WORKFLOW_ROOT}/workflow/envs/tools.yaml"
FILTER_FASTA_SCRIPT = f"{WORKFLOW_ROOT}/workflow/scripts/filter_fasta.py"
SUMMARIZE_DECONTAMINATION_SCRIPT = (
    f"{WORKFLOW_ROOT}/workflow/scripts/summarize_decontamination.py"
)
WORKFLOW_SCRIPTS_DIR = f"{WORKFLOW_ROOT}/workflow/scripts"
if WORKFLOW_SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, WORKFLOW_SCRIPTS_DIR)
from kraken_report import parse_kraken_report

DECONTAM_ASSEMBLY_MANIFEST_SIGNATURE = hashlib.sha256(
    "\n".join(
        f"{assembly_id}\t{ASSEMBLIES[assembly_id]}"
        for assembly_id in ASSEMBLY_IDS
    ).encode()
).hexdigest()
DECONTAM_COMPLETE_MARKER = (
    f"{DECONTAM_OUTDIR}/summary/"
    f"all_assemblies.{DECONTAM_ASSEMBLY_MANIFEST_SIGNATURE[:12]}.complete"
)

DEFAULT_KRAKEN_DB = "/gpfs/work3/0/qtholstg/hg38_res_v2/kraken/pluspf_20230605"
DEFAULT_KRAKEN_TARGET_TAXIDS = {
    "2": "bacteria",
    "10239": "viruses",
}


DECISION_FIELDS = [
    "assembly_id",
    "contig",
    "contig_length_bp",
    "decision",
    "nonhuman_covered_bp",
    "nonhuman_covered_percent",
    "human_covered_bp",
    "human_covered_percent",
    "human_overlap_nonhuman_bp",
    "human_overlap_nonhuman_percent",
    "human_outside_nonhuman_bp",
    "largest_nonhuman_block_bp",
    "best_nonhuman_identity_percent",
    "nonhuman_groups",
    "nonhuman_hit_count",
    "reason",
]


def strip_fasta_suffix(path):
    name = os.path.basename(path)
    for suffix in (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    raise ValueError(f"Unsupported FASTA suffix: {path}")


def decontamination_assembly_input(wildcards):
    return ASSEMBLIES[wildcards.assembly]


def all_decontamination_decisions(wildcards):
    return expand(
        f"{DECONTAM_OUTDIR}/decisions/{{assembly}}.contigs.tsv",
        assembly=ASSEMBLY_IDS,
    )


def all_cleaned_assemblies(wildcards):
    return expand(
        f"{DECONTAM_OUTDIR}/cleaned/{{assembly}}.clean.fa.gz",
        assembly=ASSEMBLY_IDS,
    )


def all_cleaned_stats(wildcards):
    return expand(
        f"{DECONTAM_OUTDIR}/stats/{{assembly}}.clean.seqkit.tsv",
        assembly=ASSEMBLY_IDS,
    )


def resolved_kraken_config():
    values = dict(config.get("kraken", {}))
    target_taxids = values.get("target_taxids", DEFAULT_KRAKEN_TARGET_TAXIDS)
    if not isinstance(target_taxids, dict):
        raise WorkflowError("kraken.target_taxids must be a mapping of taxid to label")
    return {
        "db": values.get("db", DEFAULT_KRAKEN_DB),
        "confidence": float(values.get("confidence", 0.0)),
        "minimum_hit_groups": int(values.get("minimum_hit_groups", 0)),
        "report_minimizer_data": bool(values.get("report_minimizer_data", True)),
        "target_taxids": {
            str(taxid): str(label) for taxid, label in target_taxids.items()
        },
    }


KRAKEN_CONFIG = resolved_kraken_config()


def kraken_extra_args():
    arguments = []
    if KRAKEN_CONFIG["report_minimizer_data"]:
        arguments.append("--report-minimizer-data")
    if KRAKEN_CONFIG["confidence"] > 0:
        arguments.extend(["--confidence", str(KRAKEN_CONFIG["confidence"])])
    if KRAKEN_CONFIG["minimum_hit_groups"] > 0:
        arguments.extend(
            ["--minimum-hit-groups", str(KRAKEN_CONFIG["minimum_hit_groups"])]
        )
    return " ".join(arguments)


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def fasta_lengths(path):
    lengths = {}
    current = None
    length = 0
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith(">"):
                if current is not None:
                    lengths[current] = length
                current = line[1:].strip().split()[0]
                if not current:
                    raise ValueError(f"Empty FASTA identifier at {path}:{line_number}")
                if current in lengths:
                    raise ValueError(f"Duplicate FASTA identifier: {current}")
                length = 0
            else:
                if current is None:
                    raise ValueError(
                        f"Sequence before first FASTA header at {path}:{line_number}"
                    )
                length += len(line.strip())
    if current is not None:
        lengths[current] = length
    if not lengths:
        raise ValueError(f"No sequences found in {path}")
    return lengths


def merge_intervals(intervals):
    if not intervals:
        return []
    merged = []
    start, end = sorted(intervals)[0]
    for next_start, next_end in sorted(intervals)[1:]:
        if next_start <= end:
            end = max(end, next_end)
        else:
            merged.append((start, end))
            start, end = next_start, next_end
    merged.append((start, end))
    return merged


def interval_length(intervals):
    return sum(end - start for start, end in intervals)


def read_human_alignments(paths, min_mapq, min_identity):
    intervals = defaultdict(list)
    for path in paths:
        with open_text(path) as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 12:
                    raise ValueError(f"Malformed PAF line at {path}:{line_number}")
                matches = int(fields[9])
                block_length = int(fields[10])
                mapq = int(fields[11])
                identity = 100.0 * matches / block_length if block_length else 0.0
                if mapq >= min_mapq and identity >= min_identity:
                    intervals[fields[0]].append((int(fields[2]), int(fields[3])))
    return intervals


def read_kraken_classifications(path, taxid_to_group, taxid_to_name, target_taxids):
    classifications = {}
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(
                    f"Malformed Kraken classification line at {path}:{line_number}"
                )
            status, contig, taxid = fields[:3]
            if contig in classifications:
                raise ValueError(f"Duplicate contig in Kraken output: {contig}")
            classifications[contig] = {
                "status": status,
                "taxid": taxid,
                "group": target_taxids.get(taxid, taxid_to_group.get(taxid, "")),
                "taxon_name": taxid_to_name.get(taxid, ""),
            }
    return classifications


def write_kraken_decisions(
    assembly,
    assembly_id,
    human_paf,
    kraken_classification,
    kraken_report,
    resolved_config,
    decisions,
    split_bed,
    remove_list,
    review_list,
):
    with open(resolved_config) as handle:
        resolved = json.load(handle)
    values = resolved["classification"]
    kraken_values = resolved["kraken"]

    lengths = fasta_lengths(assembly)
    human_raw = read_human_alignments(
        human_paf,
        int(values.get("human_min_mapq", 5)),
        float(values.get("human_min_identity", 95.0)),
    )
    target_taxids = {
        str(taxid): str(label)
        for taxid, label in kraken_values["target_taxids"].items()
    }
    taxid_to_group, taxid_to_name = parse_kraken_report(
        kraken_report, target_taxids
    )
    classifications = read_kraken_classifications(
        kraken_classification, taxid_to_group, taxid_to_name, target_taxids
    )

    unknown_queries = (set(human_raw) | set(classifications)) - set(lengths)
    if unknown_queries:
        raise ValueError(
            "Alignment or Kraken query IDs absent from FASTA: "
            + ", ".join(sorted(unknown_queries)[:5])
        )
    missing_classifications = set(lengths) - set(classifications)
    if missing_classifications:
        raise ValueError(
            "Kraken output is missing FASTA records: "
            + ", ".join(sorted(missing_classifications)[:5])
        )

    remove_max_human_percent = float(values.get("remove_max_human_percent", 10.0))
    rows = []
    remove_contigs = []
    review_contigs = []

    for contig, length in lengths.items():
        human = merge_intervals(human_raw.get(contig, []))
        human_bp = interval_length(human)
        human_percent = 100.0 * human_bp / length if length else 0.0
        classification = classifications[contig]
        group = classification["group"]
        taxid = classification["taxid"]
        taxon_name = classification["taxon_name"] or "unknown"

        if group:
            nonhuman_bp = length
            nonhuman_percent = 100.0
            overlap_bp = human_bp
            overlap_percent = human_percent
            human_outside_bp = 0
            largest_block = length
            hit_count = 1
            groups = group
            if human_percent <= remove_max_human_percent:
                decision = "REMOVE"
                reason = (
                    f"Kraken2 classified contig under {group} taxid {taxid} "
                    f"({taxon_name}) with little human-reference support"
                )
                remove_contigs.append(contig)
            else:
                decision = "REVIEW"
                reason = (
                    f"Kraken2 classified contig under {group} taxid {taxid} "
                    f"({taxon_name}) but human-reference support exceeds the "
                    "removal threshold"
                )
                review_contigs.append(contig)
        else:
            nonhuman_bp = 0
            nonhuman_percent = 0.0
            overlap_bp = 0
            overlap_percent = 0.0
            human_outside_bp = human_bp
            largest_block = 0
            hit_count = 0
            groups = ""
            decision = "KEEP"
            if classification["status"] == "C":
                reason = (
                    f"Kraken2 classified contig as taxid {taxid} ({taxon_name}), "
                    "outside bacterial and viral target taxa"
                )
            else:
                reason = "Kraken2 did not classify contig"

        rows.append(
            {
                "assembly_id": assembly_id,
                "contig": contig,
                "contig_length_bp": length,
                "decision": decision,
                "nonhuman_covered_bp": nonhuman_bp,
                "nonhuman_covered_percent": f"{nonhuman_percent:.6f}",
                "human_covered_bp": human_bp,
                "human_covered_percent": f"{human_percent:.6f}",
                "human_overlap_nonhuman_bp": overlap_bp,
                "human_overlap_nonhuman_percent": f"{overlap_percent:.6f}",
                "human_outside_nonhuman_bp": human_outside_bp,
                "largest_nonhuman_block_bp": largest_block,
                "best_nonhuman_identity_percent": "0.000000",
                "nonhuman_groups": groups,
                "nonhuman_hit_count": hit_count,
                "reason": reason,
            }
        )

    for output_path in (decisions, split_bed, remove_list, review_list):
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    with open(decisions, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=DECISION_FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    with open(split_bed, "w"):
        pass
    with open(remove_list, "w") as handle:
        for contig in remove_contigs:
            handle.write(contig + "\n")
    with open(review_list, "w") as handle:
        for contig in review_contigs:
            handle.write(contig + "\n")

    counts = defaultdict(int)
    for row in rows:
        counts[row["decision"]] += 1
    return counts


rule decontamination_assembly_manifest:
    output:
        f"{DECONTAM_OUTDIR}/resources/all_assemblies.tsv"
    params:
        assembly_signature=DECONTAM_ASSEMBLY_MANIFEST_SIGNATURE
    run:
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        with open(output[0], "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["assembly_id", "source_path"])
            for assembly_id in ASSEMBLY_IDS:
                writer.writerow([assembly_id, ASSEMBLIES[assembly_id]])


rule resolved_classification_config:
    output:
        f"{DECONTAM_OUTDIR}/resources/resolved_classification_config.json"
    params:
        content=json.dumps(
            {
                "kraken": KRAKEN_CONFIG,
                "classification": config["classification"],
            },
            sort_keys=True,
        )
    run:
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        with open(output[0], "w") as handle:
            json.dump(json.loads(params.content), handle, indent=2, sort_keys=True)
            handle.write("\n")


rule kraken_contigs:
    input:
        assembly=decontamination_assembly_input
    output:
        report=f"{DECONTAM_OUTDIR}/kraken/{{assembly}}.report.tsv",
        classification=f"{DECONTAM_OUTDIR}/kraken/{{assembly}}.contig_classification.tsv.gz"
    params:
        db=KRAKEN_CONFIG["db"],
        extra_args=kraken_extra_args()
    log:
        f"{DECONTAM_OUTDIR}/logs/kraken/{{assembly}}.log"
    benchmark:
        f"{DECONTAM_OUTDIR}/benchmarks/kraken/{{assembly}}.tsv"
    conda:
        DECONTAM_ENV
    threads: config["resources"].get("kraken_threads", 8)
    resources:
        mem_mb=config["resources"].get("kraken_mem_mb", 65000),
        runtime_min=config["resources"].get("kraken_runtime_min", 720)
    shell:
        r"""
        set -o pipefail
        mkdir -p "$(dirname {output.report:q})" \
                 "$(dirname {output.classification:q})" \
                 "$(dirname {log:q})"

        compression_arg=""
        case {input.assembly:q} in
          *.gz) compression_arg="--gzip-compressed" ;;
          *.bz2) compression_arg="--bzip2-compressed" ;;
        esac

        kraken2 \
          --db {params.db:q} \
          --threads {threads} {params.extra_args} \
          --report {output.report:q} \
          --output >(gzip -c > {output.classification:q}) \
          $compression_arg \
          {input.assembly:q} > {log:q} 2>&1
        """


rule classify_contigs:
    input:
        assembly=decontamination_assembly_input,
        kraken=f"{DECONTAM_OUTDIR}/kraken/{{assembly}}.contig_classification.tsv.gz",
        report=f"{DECONTAM_OUTDIR}/kraken/{{assembly}}.report.tsv",
        chm13=f"{QC_RESULTS_DIR}/alignments/CHM13/{{assembly}}.paf.gz",
        hg38=f"{QC_RESULTS_DIR}/alignments/hg38/{{assembly}}.paf.gz",
        config=f"{DECONTAM_OUTDIR}/resources/resolved_classification_config.json"
    output:
        decisions=f"{DECONTAM_OUTDIR}/decisions/{{assembly}}.contigs.tsv",
        split_bed=f"{DECONTAM_OUTDIR}/decisions/{{assembly}}.split_nonhuman.bed",
        remove_list=f"{DECONTAM_OUTDIR}/decisions/{{assembly}}.remove_contigs.txt",
        review_list=f"{DECONTAM_OUTDIR}/decisions/{{assembly}}.review_contigs.txt"
    log:
        f"{DECONTAM_OUTDIR}/logs/classify/{{assembly}}.log"
    benchmark:
        f"{DECONTAM_OUTDIR}/benchmarks/classify/{{assembly}}.tsv"
    threads: 1
    resources:
        mem_mb=4000,
        runtime_min=60
    run:
        counts = write_kraken_decisions(
            str(input.assembly),
            wildcards.assembly,
            [str(input.chm13), str(input.hg38)],
            str(input.kraken),
            str(input.report),
            str(input.config),
            str(output.decisions),
            str(output.split_bed),
            str(output.remove_list),
            str(output.review_list),
        )
        Path(log[0]).parent.mkdir(parents=True, exist_ok=True)
        with open(log[0], "w") as handle:
            handle.write(
                f"{wildcards.assembly}: KEEP={counts['KEEP']} "
                f"REVIEW={counts['REVIEW']} SPLIT={counts['SPLIT']} "
                f"REMOVE={counts['REMOVE']}\n"
            )


rule clean_assembly:
    input:
        assembly=decontamination_assembly_input,
        decisions=f"{DECONTAM_OUTDIR}/decisions/{{assembly}}.contigs.tsv",
        split_bed=f"{DECONTAM_OUTDIR}/decisions/{{assembly}}.split_nonhuman.bed"
    output:
        cleaned=f"{DECONTAM_OUTDIR}/cleaned/{{assembly}}.clean.fa.gz",
        removed=f"{DECONTAM_OUTDIR}/removed/{{assembly}}.nonhuman.fa.gz",
        review=f"{DECONTAM_OUTDIR}/review/{{assembly}}.review.fa.gz",
        split_map=f"{DECONTAM_OUTDIR}/split_maps/{{assembly}}.split_map.tsv"
    params:
        script=FILTER_FASTA_SCRIPT,
        header_prefix=lambda wildcards: wildcards.assembly
    log:
        f"{DECONTAM_OUTDIR}/logs/clean/{{assembly}}.log"
    benchmark:
        f"{DECONTAM_OUTDIR}/benchmarks/clean/{{assembly}}.tsv"
    conda:
        DECONTAM_ENV
    threads: 1
    resources:
        mem_mb=8000,
        runtime_min=120
    shell:
        r"""
        python {params.script:q} \
          --assembly {input.assembly:q} \
          --decisions {input.decisions:q} \
          --split-bed {input.split_bed:q} \
          --cleaned {output.cleaned:q} \
          --removed {output.removed:q} \
          --review {output.review:q} \
          --split-map {output.split_map:q} \
          --header-prefix {params.header_prefix:q} > {log:q} 2>&1
        """


rule cleaned_stats:
    input:
        f"{DECONTAM_OUTDIR}/cleaned/{{assembly}}.clean.fa.gz"
    output:
        f"{DECONTAM_OUTDIR}/stats/{{assembly}}.clean.seqkit.tsv"
    log:
        f"{DECONTAM_OUTDIR}/logs/seqkit/{{assembly}}.clean.log"
    conda:
        DECONTAM_ENV
    threads: 1
    resources:
        mem_mb=4000,
        runtime_min=60
    shell:
        r"""
        seqkit stats --all --tabular {input:q} > {output:q} 2> {log:q}
        """


rule summarize_decontamination:
    input:
        manifest=f"{DECONTAM_OUTDIR}/resources/all_assemblies.tsv",
        decisions=all_decontamination_decisions,
        cleaned=all_cleaned_assemblies,
        stats=all_cleaned_stats
    output:
        summary=f"{DECONTAM_OUTDIR}/summary/contamination_summary.tsv",
        actions=f"{DECONTAM_OUTDIR}/summary/contig_actions.tsv",
        review=f"{DECONTAM_OUTDIR}/summary/review_candidates.tsv",
        graph_list=f"{DECONTAM_OUTDIR}/summary/graph_cleaned_assemblies.txt",
        complete=DECONTAM_COMPLETE_MARKER
    params:
        results_dir=DECONTAM_OUTDIR,
        script=SUMMARIZE_DECONTAMINATION_SCRIPT
    log:
        f"{DECONTAM_OUTDIR}/logs/summarize.log"
    conda:
        DECONTAM_ENV
    threads: 1
    resources:
        mem_mb=8000,
        runtime_min=60
    shell:
        r"""
        python {params.script:q} \
          --manifest {input.manifest:q} \
          --results-dir {params.results_dir:q} \
          --summary {output.summary:q} \
          --actions {output.actions:q} \
          --review {output.review:q} \
          --graph-list {output.graph_list:q} \
          --complete-marker {output.complete:q} > {log:q} 2>&1
        """
