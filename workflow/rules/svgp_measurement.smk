"""Bounded graph pilot, ONT corroboration and a gated one-graph cohort."""

import csv
import json
import re
import shlex
from pathlib import Path

SGP_M = SGP["measurement"]
SGP_M_RES = SGP_M["resources"]
SGP_M_SCRIPT = str(workflow.source_path("../scripts/svgp_measurement.py"))
SGP_M_ENV = str(workflow.source_path("../envs/svgp_mapping.yaml"))
SGP_PILOT_EVIDENCE = f"{SGP_OUT}/pilot/evidence.tsv.gz"
SGP_PILOT_SUMMARY = f"{SGP_OUT}/pilot/summary.tsv"
SGP_ONT_EVIDENCE = f"{SGP_OUT}/validation/ont_by_sample_locus.tsv"
SGP_ONT_SUMMARY = f"{SGP_OUT}/validation/concordance.tsv"
SGP_COHORT_EVIDENCE = f"{SGP_OUT}/genotypes/evidence.tsv.gz"
SGP_COHORT_QC = f"{SGP_OUT}/genotypes/frequencies.tsv"
SGP_COHORT_GATE = f"{SGP_OUT}/gates/cohort_release.tsv"


def sgp_graph(wildcards):
    if wildcards.graph not in SGP_M["pilot_graphs"]:
        raise ValueError("Only the frozen two pilot graphs are permitted")
    if (wildcards.graph != SGP_M["main_graph"]
            and getattr(wildcards, "person", None) not in (None, *SGP_M["pilot_people"])):
        raise ValueError("Comparator mapping is restricted to the pilot")
    return SGP_M["graphs"][wildcards.graph]


def sgp_extraction(name):
    value = SGP_M["extraction"].get(name)
    if value is None:
        raise ValueError(f"Set measurement.extraction.{name} using development controls first")
    return int(value)


def sgp_pilot_jobs(pattern):
    people = SGP_M["pilot_people"]
    if not 12 <= len(people) <= 20 or len(people) != len(set(people)):
        raise ValueError("Freeze 12-20 unique pilot people; target one mapping directory for a single-person smoke run")
    if set(SGP_M["pilot_graphs"]) != {"hprc_only", "hprc_plus_local"} or len(SGP_M["pilot_graphs"]) != 2:
        raise ValueError("The paired pilot requires exactly hprc_only and hprc_plus_local")
    return expand(pattern, graph=SGP_M["pilot_graphs"], person=people)


def sgp_cohort_jobs(wildcards):
    people = SGP_M["cohort_people"]
    if not people or len(people) != len(set(people)):
        raise ValueError("Freeze a nonempty list of unique technically ready cohort people")
    return expand(f"{SGP_OUT}/measurement/support/{{graph}}/{{person}}.tsv.gz",
                  graph=[SGP_M["main_graph"]], person=people)


def sgp_illumina_reads(wildcards):
    # The checkpoint resolves explicit aliases before Snakemake tracks the raw
    # FASTQs. Changed/missing files are therefore actual mapping dependencies.
    checkpoints.svgp_inventory.get()
    with open(f"{SGP_OUT}/metadata/reads.tsv") as handle:
        lanes = [row for row in csv.DictReader(handle, delimiter="\t")
                 if row["person_id"] == wildcards.person]
    if not lanes or any(row["status"] != "READY" for row in lanes):
        raise ValueError(f"All paired lanes must be READY for {wildcards.person}")
    return [row[mate] for row in lanes for mate in ("r1", "r2")]


rule svgp_pilot_mapping:
    input:
        lambda wc: sgp_pilot_jobs(f"{SGP_OUT}/mapping/{{graph}}/{{person}}")


rule svgp_pilot_supporting_calls:
    input:
        lambda wc: sgp_pilot_jobs(f"{SGP_OUT}/measurement/supporting_calls/{{graph}}/{{person}}")


rule svgp_pilot:
    input:
        SGP_PILOT_EVIDENCE,
        SGP_PILOT_SUMMARY


rule svgp_ont_validation:
    input:
        SGP_ONT_EVIDENCE,
        SGP_ONT_SUMMARY


rule svgp_cohort:
    input:
        SGP_COHORT_EVIDENCE,
        SGP_COHORT_QC,
        f"{SGP_OUT}/genotypes/carrier_matrix.tsv.gz",
        f"{SGP_OUT}/genotypes/dosage_matrix.tsv.gz",
        f"{SGP_OUT}/genotypes/callability.tsv.gz"


rule svgp_audit_graph_bundle:
    input:
        script=SGP_M_SCRIPT,
        manifest=lambda wc: sgp_graph(wc)["manifest"],
        reference_review=lambda wc: ([sgp_graph(wc)["reference_review"]]
                                     if sgp_graph(wc).get("reference_paths") else []),
        components=lambda wc: [sgp_graph(wc)[name] for name in
                                ("gbz", "dist", "minimizer", "zipcodes")]
    output:
        f"{SGP_OUT}/measurement/bundles/{{graph}}.tsv"
    params:
        graph=lambda wc: json.dumps(sgp_graph(wc), sort_keys=True),
        review=lambda wc: sgp_graph(wc).get("reference_review", "")
    threads: 1
    resources:
        mem_mb=2000,
        runtime_min=360
    log: f"{SGP_OUT}/logs/measurement/bundles/{{graph}}.log"
    conda: SGP_M_ENV
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {input.script:q} audit-bundle --manifest {input.manifest:q} \
          --reference-review {params.review:q} --graph {wildcards.graph:q} \
          --graph-config {params.graph:q} --output {output:q} > {log:q} 2>&1
        """


rule svgp_map_illumina:
    input:
        samples=f"{SGP_OUT}/metadata/samples.tsv",
        reads=f"{SGP_OUT}/metadata/reads.tsv",
        fastqs=sgp_illumina_reads,
        bundle=f"{SGP_OUT}/measurement/bundles/{{graph}}.tsv",
        release=lambda wc: ([] if wc.person in SGP_M["pilot_people"] else [SGP_COHORT_GATE]),
        gbz=lambda wc: sgp_graph(wc)["gbz"],
        dist=lambda wc: sgp_graph(wc)["dist"],
        minimizer=lambda wc: sgp_graph(wc)["minimizer"],
        zipcodes=lambda wc: sgp_graph(wc)["zipcodes"]
    output:
        directory(f"{SGP_OUT}/mapping/{{graph}}/{{person}}")
    wildcard_constraints:
        graph="hprc_only|hprc_plus_local",
        person="[A-Za-z0-9_.-]+"
    threads: SGP_M_RES["mapping_threads"]
    resources:
        mem_mb=SGP_M_RES["mapping_mem_mb"],
        runtime_min=SGP_M_RES["mapping_runtime_min"]
    log: f"{SGP_OUT}/logs/measurement/mapping/{{graph}}/{{person}}.log"
    benchmark: f"{SGP_OUT}/benchmarks/measurement/mapping/{{graph}}/{{person}}.tsv"
    conda: SGP_M_ENV
    run:
        allowed = SGP_M["pilot_people"] + SGP_M["cohort_people"]
        if wildcards.person not in allowed or wildcards.person in {".", ".."}:
            raise ValueError("Person was not explicitly selected for pilot/cohort")
        with open(input.samples) as handle:
            samples = {r["person_id"]: r for r in csv.DictReader(handle, delimiter="\t")}
        person = samples[wildcards.person]
        # Giraffe can update a writable distance-index mmap. Never change the
        # shared source index or let parallel jobs write it.
        if Path(input.dist).stat().st_mode & 0o222:
            raise ValueError("Use an existing read-only distance index (no write permission bits); shared input is not modified")
        if person["fastq_manifest_complete"] != "yes":
            raise ValueError("Read manifest completeness must be explicitly confirmed")
        # Technical measurement can include labelled development/sensitivity
        # people. The association rule separately selects the primary cohort.
        with open(input.reads) as handle:
            reads = [r for r in csv.DictReader(handle, delimiter="\t")
                     if r["person_id"] == wildcards.person]
        if not reads or any(r["status"] != "READY" for r in reads):
            raise ValueError("All paired lanes for a selected person must be READY")
        libraries = sorted({r["library_id"] for r in reads})
        Path(output[0]).mkdir(parents=True, exist_ok=True)
        Path(log[0]).parent.mkdir(parents=True, exist_ok=True)
        records = []
        for library in libraries:
            if not re.fullmatch(r"[A-Za-z0-9_.-]+", library) or library in {".", ".."}:
                raise ValueError("Use an explicit safe library ID in the read manifest")
            lanes = sorted((r for r in reads if r["library_id"] == library),
                           key=lambda r: (r["lane"], r["read_id"]))
            r1, r2 = ([r[mate] for r in lanes] for mate in ("r1", "r2"))
            gam = str(Path(output[0]) / f"{library}.gam")
            pipe1, pipe2 = (str(Path(output[0]) / f"{library}.{mate}.fifo")
                            for mate in ("r1", "r2"))
            mapper_threads = max(1, threads - 2)
            shell(r"""
                stream_fastq() {{
                    for fastq in "$@"; do
                        case "$fastq" in
                          *.gz) gzip -cd -- "$fastq" ;;
                          *) cat -- "$fastq" ;;
                        esac || return 1
                    done
                }}
                vg version >> {log:q} 2>&1
                stream_pipe1={pipe1:q}
                stream_pipe2={pipe2:q}
                mkfifo "$stream_pipe1" "$stream_pipe2"
                trap 'rm -f -- "$stream_pipe1" "$stream_pipe2"' EXIT
                stream_fastq {r1:q} > {pipe1:q} 2>> {log:q} & reader1=$!
                stream_fastq {r2:q} > {pipe2:q} 2>> {log:q} & reader2=$!
                if vg giraffe -p -t {mapper_threads} -Z {input.gbz:q} -d {input.dist:q} \
                  -m {input.minimizer:q} -z {input.zipcodes:q} \
                  -N {wildcards.person:q} -R {library:q} \
                  -f {pipe1:q} -f {pipe2:q} > {gam:q} 2>> {log:q}; then
                    reader_status=0
                    wait "$reader1" || reader_status=1
                    wait "$reader2" || reader_status=1
                    test "$reader_status" -eq 0
                else
                    kill "$reader1" "$reader2" 2>/dev/null || true
                    wait "$reader1" "$reader2" 2>/dev/null || true
                    exit 1
                fi
                vg stats -a {gam:q} >> {log:q} 2>&1
            """)
            records.append(dict(library_id=library, gam=gam, lanes=len(lanes),
                                r1=";".join(r1), r2=";".join(r2)))
        with (Path(output[0]) / "libraries.tsv").open("w") as handle:
            writer = csv.DictWriter(handle, list(records[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(records)


rule svgp_extract_molecule_support:
    input:
        script=SGP_M_SCRIPT,
        mapping=f"{SGP_OUT}/mapping/{{graph}}/{{person}}",
        samples=f"{SGP_OUT}/metadata/samples.tsv",
        panel=f"{SGP_OUT}/panel/panel.tsv",
        targets=lambda wc: sgp_graph(wc)["targets"],
        molecules=SGP_M["molecule_assignments"]
    output:
        f"{SGP_OUT}/measurement/support/{{graph}}/{{person}}.tsv.gz"
    params:
        graph=lambda wc: json.dumps(sgp_graph(wc), sort_keys=True),
        mapq=lambda wc: sgp_extraction("min_mapq"),
        anchor=lambda wc: sgp_extraction("anchor_matches")
    threads: 1
    resources:
        mem_mb=SGP_M_RES["summary_mem_mb"],
        runtime_min=SGP_M_RES["summary_runtime_min"]
    log: f"{SGP_OUT}/logs/measurement/support/{{graph}}/{{person}}.log"
    benchmark: f"{SGP_OUT}/benchmarks/measurement/support/{{graph}}/{{person}}.tsv"
    conda: SGP_M_ENV
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {input.script:q} extract --mapping {input.mapping:q} \
          --samples {input.samples:q} \
          --panel {input.panel:q} --targets {input.targets:q} --molecules {input.molecules:q} \
          --person {wildcards.person:q} --graph {wildcards.graph:q} \
          --graph-config {params.graph:q} --min-mapq {params.mapq} \
          --anchor-matches {params.anchor} --output {output:q} > {log:q} 2>&1
        """


rule svgp_supporting_graph_calls:
    input:
        mapping=f"{SGP_OUT}/mapping/{{graph}}/{{person}}",
        bundle=f"{SGP_OUT}/measurement/bundles/{{graph}}.tsv",
        gbz=lambda wc: sgp_graph(wc)["gbz"],
        review=lambda wc: ([sgp_graph(wc)["genotyping_review"]]
                            if sgp_graph(wc).get("genotyping_review") else [])
    output:
        directory(f"{SGP_OUT}/measurement/supporting_calls/{{graph}}/{{person}}")
    threads: SGP_M_RES["mapping_threads"]
    resources:
        mem_mb=SGP_M_RES["mapping_mem_mb"],
        runtime_min=SGP_M_RES["mapping_runtime_min"]
    log: f"{SGP_OUT}/logs/measurement/supporting_calls/{{graph}}/{{person}}.log"
    benchmark: f"{SGP_OUT}/benchmarks/measurement/supporting_calls/{{graph}}/{{person}}.tsv"
    conda: SGP_M_ENV
    run:
        graph = sgp_graph(wildcards)
        Path(output[0]).mkdir(parents=True, exist_ok=True)
        Path(log[0]).parent.mkdir(parents=True, exist_ok=True)
        status = Path(output[0]) / "status.tsv"
        if not input.review or not graph.get("reference_paths"):
            status.write_text("status\treason\nSKIPPED\tverified_reference_paths_and_direct_genotyping_review_required\n")
            Path(log[0]).write_text("No valid direct graph genotyping configuration. Carrier assay remains available.\n")
        else:
            with open(input.review[0]) as handle:
                reviews = list(csv.DictReader(handle, delimiter="\t"))
            if (len(reviews) != 1 or reviews[0]["status"] != "APPROVED"
                    or reviews[0]["mode"] != "direct_graph"
                    or reviews[0]["graph_sha256"] != graph["graph_sha256"]
                    or reviews[0]["node_space_id"] != graph["node_space_id"]
                    or not Path(reviews[0]["allele_call_crosswalk"]).is_file()):
                raise ValueError("Direct graph call configuration/crosswalk not reviewed")
            with (Path(input.mapping) / "libraries.tsv").open() as handle:
                libraries = list(csv.DictReader(handle, delimiter="\t"))
            packs = []
            mapq = sgp_extraction("min_mapq")
            for library in libraries:
                gam = library["gam"]
                pack = str(Path(output[0]) / (library["library_id"] + ".pack"))
                shell("vg pack -x {input.gbz:q} -g {gam:q} -Q {mapq} -t {threads} -o {pack:q} >> {log:q} 2>&1")
                packs.append(pack)
            pack_options = " ".join("-i " + shlex.quote(path) for path in packs)
            reference_options = " ".join("-p " + shlex.quote(path) for path in graph["reference_paths"])
            combined = str(Path(output[0]) / "combined.pack")
            vcf = str(Path(output[0]) / "supporting.vcf.gz")
            shell(r"""
                vg pack -x {input.gbz:q} {pack_options} -t {threads} -o {combined:q} >> {log:q} 2>&1
                vg call {input.gbz:q} -k {combined:q} {reference_options} \
                  -s {wildcards.person:q} -t {threads} 2>> {log:q} \
                  | bcftools view -Oz -o {vcf:q} 2>> {log:q}
                bcftools index --csi {vcf:q} >> {log:q} 2>&1
            """)
            status.write_text("status\treason\nSUPPORTING_ONLY\tpack_is_read_depth_not_duplicate_resolved_molecule_evidence;VCF_not_calibrated_genotypes\n")


rule svgp_summarize_pilot:
    input:
        script=SGP_M_SCRIPT,
        support=lambda wc: sgp_pilot_jobs(f"{SGP_OUT}/measurement/support/{{graph}}/{{person}}.tsv.gz"),
        lock=SGP_M["measurement_lock"]
    output:
        evidence=SGP_PILOT_EVIDENCE,
        summary=SGP_PILOT_SUMMARY,
        paired=f"{SGP_OUT}/pilot/paired_by_person.tsv"
    params:
        mapq=lambda wc: sgp_extraction("min_mapq"),
        anchor=lambda wc: sgp_extraction("anchor_matches")
    threads: 1
    resources:
        mem_mb=SGP_M_RES["summary_mem_mb"],
        runtime_min=SGP_M_RES["summary_runtime_min"]
    log: f"{SGP_OUT}/logs/measurement/pilot_summary.log"
    conda: SGP_M_ENV
    shell:
        r"""
        mkdir -p "$(dirname {output.evidence:q})" "$(dirname {log:q})"
        python {input.script:q} summarize --inputs {input.support:q} --lock {input.lock:q} \
          --min-mapq {params.mapq} --anchor-matches {params.anchor} \
          --output {output.evidence:q} --summary {output.summary:q} \
          --paired-summary {output.paired:q} > {log:q} 2>&1
        """


rule svgp_ont_force_genotype:
    input:
        reference=SGP_M["ont"]["reference"] or [],
        known_sites=SGP_M["ont"]["known_sites_vcf"] or [],
        review=SGP_M["ont"]["known_sites_review"] or [],
        source=lambda wc: ([SGP_M["ont"]["samples"][wc.person]["bam"],
                            SGP_M["ont"]["samples"][wc.person]["bai"]]
                           if SGP_M["ont"]["samples"][wc.person].get("bam")
                           else SGP_M["ont"]["samples"][wc.person]["fastqs"])
    output:
        directory(f"{SGP_OUT}/validation/force_calls/{{person}}")
    threads: SGP_M_RES["ont_threads"]
    resources:
        mem_mb=SGP_M_RES["ont_mem_mb"],
        runtime_min=SGP_M_RES["ont_runtime_min"]
    log: f"{SGP_OUT}/logs/measurement/ont/{{person}}.log"
    benchmark: f"{SGP_OUT}/benchmarks/measurement/ont/{{person}}.tsv"
    conda: SGP_M_ENV
    run:
        if wildcards.person not in SGP_M["ont"]["people"]:
            raise ValueError("ONT person not in the frozen manifest")
        if not input.reference or not input.known_sites or not input.review:
            raise ValueError("Set the ONT reference, known-site VCF and matching-build review")
        sample = SGP_M["ont"]["samples"][wildcards.person]
        with open(input.review) as handle:
            review = list(csv.DictReader(handle, delimiter="\t"))
        if (len(review) != 1 or review[0]["status"] != "VERIFIED"
                or review[0]["reference_build"] != SGP_M["ont"]["reference_build"]
                or review[0]["vcf"] != str(input.known_sites)
                or review[0]["reference"] != str(input.reference)):
            raise ValueError("Known sites and ONT whole-genome reference need matching build review")
        if (sample.get("reference_build") != SGP_M["ont"]["reference_build"]
                or sample.get("read_scope") != "whole_genome"
                or not Path(sample.get("alignment_review", "")).is_file()):
            raise ValueError("ONT requires explicit whole-genome scope and matching-build alignment review")
        Path(output[0]).mkdir(parents=True, exist_ok=True)
        Path(log[0]).parent.mkdir(parents=True, exist_ok=True)
        bam = sample.get("bam", "")
        if bam and sample["bai"] not in {bam + ".bai", str(Path(bam).with_suffix(".bai"))}:
            raise ValueError("Supply a conventional BAM index path discoverable by Sniffles")
        if not bam:
            preset = sample["preset"]
            if preset not in {"map-ont", "lr:hq"}:
                raise ValueError("Choose ONT preset from actual chemistry/basecalling quality")
            fastqs = sample["fastqs"]
            bam = str(Path(output[0]) / "whole_genome.bam")
            # Split the job's threads between mapping and sorting.
            mapper_threads = max(1, threads - 2)
            shell(r"""
                minimap2 --version >> {log:q} 2>&1
                minimap2 -ax {preset:q} -t {mapper_threads} {input.reference:q} {fastqs:q} \
                  2>> {log:q} | samtools sort -@ 1 -m 2G -o {bam:q} - 2>> {log:q}
                samtools index -@ 1 {bam:q} 2>> {log:q}
            """)
        vcf = str(Path(output[0]) / "known_sites.vcf.gz")
        shell(r"""
            sniffles --version >> {log:q} 2>&1
            sniffles --input {bam:q} --reference {input.reference:q} \
              --genotype-vcf {input.known_sites:q} --vcf {vcf:q} \
              --threads {threads} --output-rnames >> {log:q} 2>&1
        """)


rule svgp_summarize_ont:
    input:
        script=SGP_M_SCRIPT,
        calls=expand(f"{SGP_OUT}/validation/force_calls/{{person}}", person=SGP_M["ont"]["people"]),
        panel=f"{SGP_OUT}/panel/panel.tsv",
        samples=f"{SGP_OUT}/metadata/samples.tsv",
        pairing=f"{SGP_OUT}/metadata/ont_pairing.tsv",
        adjudication=SGP_M["ont"]["adjudication"],
        illumina=SGP_PILOT_EVIDENCE
    output:
        evidence=SGP_ONT_EVIDENCE,
        summary=SGP_ONT_SUMMARY
    params:
        people=SGP_M["ont"]["people"],
        graph=SGP_M["main_graph"],
        build=SGP_M["ont"]["reference_build"]
    threads: 1
    resources:
        mem_mb=SGP_M_RES["summary_mem_mb"],
        runtime_min=SGP_M_RES["summary_runtime_min"]
    log: f"{SGP_OUT}/logs/measurement/ont_summary.log"
    conda: SGP_M_ENV
    shell:
        r"""
        mkdir -p "$(dirname {output.evidence:q})" "$(dirname {log:q})"
        python {input.script:q} validate-ont --panel {input.panel:q} \
          --samples {input.samples:q} --pairing {input.pairing:q} \
          --adjudication {input.adjudication:q} --illumina {input.illumina:q} \
          --people {params.people:q} --graph {params.graph:q} --reference-build {params.build:q} \
          --output {output.evidence:q} --summary {output.summary:q} > {log:q} 2>&1
        """


rule svgp_release_cohort:
    input:
        script=SGP_M_SCRIPT,
        review=SGP_M["cohort_release_review"],
        lock=SGP_M["measurement_lock"],
        pilot=SGP_PILOT_EVIDENCE,
        ont=SGP_ONT_EVIDENCE
    output:
        SGP_COHORT_GATE
    params:
        graph=SGP_M["main_graph"]
    threads: 1
    resources:
        mem_mb=1000,
        runtime_min=30
    log: f"{SGP_OUT}/logs/measurement/cohort_release.log"
    conda: SGP_M_ENV
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python {input.script:q} release --review {input.review:q} --graph {params.graph:q} \
          --lock {input.lock:q} --pilot {input.pilot:q} --ont {input.ont:q} \
          --output {output:q} > {log:q} 2>&1
        """


rule svgp_summarize_cohort:
    input:
        script=SGP_M_SCRIPT,
        release=SGP_COHORT_GATE,
        support=sgp_cohort_jobs,
        lock=SGP_M["measurement_lock"]
    output:
        evidence=SGP_COHORT_EVIDENCE,
        summary=SGP_COHORT_QC,
        carriers=f"{SGP_OUT}/genotypes/carrier_matrix.tsv.gz",
        dosage=f"{SGP_OUT}/genotypes/dosage_matrix.tsv.gz",
        callability=f"{SGP_OUT}/genotypes/callability.tsv.gz"
    params:
        directory=f"{SGP_OUT}/genotypes",
        mapq=lambda wc: sgp_extraction("min_mapq"),
        anchor=lambda wc: sgp_extraction("anchor_matches")
    threads: 1
    resources:
        mem_mb=SGP_M_RES["summary_mem_mb"],
        runtime_min=SGP_M_RES["summary_runtime_min"]
    log: f"{SGP_OUT}/logs/measurement/cohort_summary.log"
    conda: SGP_M_ENV
    shell:
        r"""
        mkdir -p {params.directory:q} "$(dirname {log:q})"
        python {input.script:q} summarize --inputs {input.support:q} --lock {input.lock:q} \
          --min-mapq {params.mapq} --anchor-matches {params.anchor} \
          --matrix-dir {params.directory:q} --output {output.evidence:q} \
          --summary {output.summary:q} > {log:q} 2>&1
        """
