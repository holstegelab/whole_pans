"""Small, standard-library reports; study_status never schedules upstream analyses."""

import gzip
import json
from collections import Counter
from datetime import datetime, timezone
from xml.sax.saxutils import escape


def svgp_report_rows(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def svgp_report_section(title, records):
    lines = ["### " + title, ""]
    for record in records:
        lines += ["| Field | Value |", "| --- | --- |"]
        lines += ["| " + str(key) + " | " + str(value).replace("|", "\\|").replace("\n", " ") + " |"
                  for key, value in record.items()]
        lines += [""]
    return lines if records else lines + ["NOT_AVAILABLE / NOT_ASSESSED", ""]


def svgp_report_group(sample):
    flags = [sample.get(key) for key in ("discovery_donor", "related_to_discovery_donor")]
    if flags == ["no", "no"]:
        return "independent_primary" if sample.get("analysis_eligible") == sample.get("unrelated_primary") == "yes" else "non_donor_not_primary"
    return "discovery_overlap" if "yes" in flags else "independence_unknown"


rule study_status:
    input: SGP_OUT + "/report/stages.tsv", SGP_OUT + "/report/report.md", SGP_OUT + "/report/exclusions.tsv"


rule pilot_compute:
    input: SGP_OUT + "/pilot/compute_budget.tsv"


rule svgp_pilot_compute:
    input: lambda wc: sgp_pilot_jobs(f"{SGP_OUT}/mapping/{{graph}}/{{person}}")
    output: SGP_OUT + "/pilot/compute_budget.tsv"
    threads: 1
    resources: mem_mb=1000, runtime_min=10
    log: SGP_OUT + "/logs/pilot_compute.log"
    run:
        records = []
        for graph in SGP["measurement"]["pilot_graphs"]:
            for person in SGP["measurement"]["pilot_people"]:
                benchmark = Path(SGP_OUT) / "benchmarks/measurement/mapping" / graph / (person + ".tsv")
                values = list(svgp_report_rows(benchmark))
                if len(values) != 1:
                    raise ValueError("Use one complete benchmark per pilot person/graph; repeated benchmarks need separate review")
                value = values[0]
                directory_path = Path(SGP_OUT) / "mapping" / graph / person
                records.append(dict(person_id=person, graph_id=graph, benchmark=str(benchmark),
                    wall_seconds=value.get("s", "NOT_MEASURED"), cpu_seconds=value.get("cpu_time", "NOT_MEASURED"),
                    max_rss_mb=value.get("max_rss", "NOT_MEASURED"),
                    retained_output_bytes=sum(path.stat().st_size for path in directory_path.rglob("*") if path.is_file()),
                    temporary_disk_peak="NOT_MEASURED", status="MEASURED_PILOT",
                    reason="Full native read mapping; temporary disk peak requires separate server monitoring"))
        if not records:
            raise ValueError("Configure and run the bounded pilot before estimating compute")
        with open(output[0], "w", newline="") as handle:
            writer = csv.DictWriter(handle, list(records[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(records)
        Path(log[0]).write_text("Recorded actual pilot benchmarks and retained bytes; cohort extrapolation must account for read counts/depth and concurrency.\n")


rule svgp_study_status:
    # Deliberately no analysis inputs: refresh explicitly with --forcerun svgp_study_status.
    output:
        stages=SGP_OUT + "/report/stages.tsv",
        report=SGP_OUT + "/report/report.md",
        exclusions=SGP_OUT + "/report/exclusions.tsv"
    threads: 1
    resources: mem_mb=2000, runtime_min=30
    log: SGP_OUT + "/logs/study_status.log"
    run:
        stages = {
            "phenotypes": ["metadata/phenotypes.tsv", "metadata/phenotype_join_report.tsv"],
            "reported_ont": ["metadata/reported_ont_runs.tsv"],
            "existing_calls_inventory": ["existing_calls/source_files.tsv", "existing_calls/source_summary.tsv"],
            "gcnv_import": ["existing_calls/gcnv/import_status.tsv"],
            "repeat_import": ["existing_calls/expansionhunter/import_status.tsv"],
            "inventory": ["metadata/samples.tsv", "metadata/resources.tsv", "metadata/exclusions.tsv"],
            "catalogue": ["catalogue/alleles.tsv", "catalogue/coordinate_crosswalk.tsv", "catalogue/rejections.tsv"],
            "representation": ["representation/resource_matrix.tsv", "representation/flow.tsv"],
            "panel": ["panel/panel.tsv"], "pilot": ["pilot/summary.tsv"],
            "ont_validation": ["validation/ont_by_sample_locus.tsv", "validation/concordance.tsv"],
            "cohort": ["genotypes/frequencies.tsv"], "annotation": ["annotation/allele_gene_context.tsv"],
            "association": ["association/status.tsv"], "prediction": ["prediction/status.tsv"], "dossiers": ["report/candidate_dossiers.md"]}
        snapshot, lines = [], ["# SV study evidence snapshot", "", datetime.now(timezone.utc).isoformat(), "",
            "OUTPUT_AVAILABLE means an output exists, not that all scientific gates passed. Missing outputs cannot distinguish an unrun job from an unrecorded failure.", "",
            "Refresh after running analyses with `--forcerun svgp_study_status`. This target does not launch mapping.", ""]
        for stage, names in stages.items():
            available = [str(Path(SGP_OUT) / name) for name in names if (Path(SGP_OUT) / name).is_file()]
            missing = [name for name in names if not (Path(SGP_OUT) / name).is_file()]
            states, error = Counter(), ""
            try:
                for path in available:
                    if path.endswith(".tsv"):
                        states.update(row["status"] for row in svgp_report_rows(path) if row.get("status"))
            except (OSError, ValueError, csv.Error) as exc:
                error = str(exc)
            state = "OUTPUT_AVAILABLE" if available else "UNRUN"
            if states and set(states) <= {"COMPLETED", "COMPLETE", "ESTIMATED"} and not missing:
                state = "COMPLETED"
            if states and set(states) <= {"SKIPPED", "NOT_IMPLEMENTED"}:
                state = "SKIPPED"
            if error or set(states) & {"FAILED", "ERROR"}:
                state = "FAILED"
            disabled = stage == "prediction" and not available and not SGP.get("prediction", {}).get("enabled", False)
            if disabled:
                state = "SKIPPED"
            record = dict(stage=stage, status=state, observed_status_counts=json.dumps(dict(states), sort_keys=True),
                          available_outputs=";".join(available), missing_outputs=";".join(missing), reason=error or ("disabled_by_configuration" if disabled else "."))
            snapshot.append(record)
            lines += [f"- **{stage}: {state}**; statuses: {record['observed_status_counts']}; missing: {record['missing_outputs'] or 'none'}"]
        lines += ["", "Scope: insertion-first panel; additional relative to frozen resources. Gene context is a hypothesis, not functional or disease evidence.", "",
                  "Available benchmark files (measured values remain in the source tables):"]
        lines += ["- " + str(path) for path in sorted((Path(SGP_OUT) / "benchmarks").rglob("*.tsv"))]
        Path(output.stages).parent.mkdir(parents=True, exist_ok=True)
        with open(output.stages, "w", newline="") as handle:
            writer = csv.DictWriter(handle, list(snapshot[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(snapshot)
        ledger = ["source_table", "entity", "identifier", "person_id", "locus_id", "allele_id", "graph_id", "status", "reason", "source_record_json"]
        with open(output.exclusions, "w", newline="") as handle:
            writer = csv.DictWriter(handle, ledger, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for name in ("metadata/exclusions.tsv", "catalogue/rejections.tsv", "pilot/evidence.tsv.gz", "genotypes/evidence.tsv.gz"):
                path = Path(SGP_OUT) / name
                if not path.is_file():
                    continue
                for row in svgp_report_rows(path):
                    if name.endswith(".gz") and row.get("callable") == "yes":
                        continue
                    record = {key: row.get(key, "") for key in ledger}
                    record.update(source_table=str(path), source_record_json=json.dumps(row, sort_keys=True))
                    record["entity"] = row.get("entity") or ("sample_allele" if name.endswith(".gz") else "candidate")
                    record["identifier"] = row.get("identifier") or row.get("candidate_id", "")
                    writer.writerow(record)
        Path(output.report).write_text("\n".join(lines) + "\n")
        Path(log[0]).write_text("Recorded existing outputs and content statuses; streamed available exclusions and no-calls.\n")


rule candidate_flow:
    input: SGP_OUT + "/report/candidate_flow.svg"


rule svgp_candidate_flow:
    input: SGP_OUT + "/representation/flow.tsv"
    output: SGP_OUT + "/report/candidate_flow.svg"
    threads: 1
    resources: mem_mb=1000, runtime_min=5
    log: SGP_OUT + "/logs/candidate_flow.log"
    run:
        counts = {row["metric"]: int(row["count"]) for row in svgp_report_rows(input[0])}
        groups = [("Candidate records (all attempted records denominator)", [
            "all_attempted_candidate_records", "rejected_candidate_records", "source_anchor_pending_candidate_records", "strict_additional_attempted_candidate_records"]),
            ("Distinct local ALT alleles (reconstructed alleles denominator)", [
            "reconstructed_distinct_local_alt_alleles", "three_resource_resolved_distinct_local_alt_alleles", "strict_additional_distinct_local_alt_alleles"]),
            ("Distinct local loci (separate unit)", ["strict_additional_distinct_local_loci"])]
        svg, y = ['<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="550" viewBox="0 0 1200 550">',
                  '<rect width="1200" height="550" fill="white"/>', '<g font-family="sans-serif" font-size="15" fill="#182b39">'], 30
        for title, keys in groups:
            svg.append(f'<text x="20" y="{y}" font-weight="bold">{escape(title)}</text>')
            y += 30
            scale = max(1, counts[keys[0]])
            for key in keys:
                count = counts[key]
                svg += [f'<text x="20" y="{y}">{escape(key.replace("_", " "))}</text>',
                        f'<rect x="600" y="{y-15}" width="{440*count/scale:.2f}" height="20" fill="#338d9c"/>',
                        f'<text x="{610+440*count/scale:.2f}" y="{y}">{count}</text>']
                y += 32
            y += 20
        svg += ['<text x="20" y="530">Descriptive counts; rows are not sequential subtraction steps. Units and denominators are separate.</text>', '</g></svg>']
        Path(output[0]).write_text("\n".join(svg) + "\n")
        Path(log[0]).write_text("Plotted existing flow counts with separate record, allele and locus units.\n")


rule dossiers:
    input: SGP_OUT + "/report/candidate_dossiers.md"


rule svgp_dossiers:
    input:
        catalogue=SGP_OUT + "/catalogue/alleles.tsv", coordinates=SGP_OUT + "/catalogue/coordinate_crosswalk.tsv",
        panel=SGP_OUT + "/panel/panel.tsv", lock=SGP["measurement"]["measurement_lock"],
        representation=SGP_OUT + "/representation/resource_matrix.tsv", annotation=SGP_OUT + "/annotation/allele_gene_context.tsv",
        frequencies=SGP_OUT + "/genotypes/frequencies.tsv", evidence=SGP_OUT + "/genotypes/evidence.tsv.gz",
        ont=SGP_OUT + "/validation/ont_by_sample_locus.tsv", samples=SGP_OUT + "/metadata/samples.tsv",
        pairing=SGP_OUT + "/metadata/ont_pairing.tsv"
    output: SGP_OUT + "/report/candidate_dossiers.md"
    params: selected=SGP.get("report", {}).get("dossier_alleles", [])
    threads: 1
    resources: mem_mb=4000, runtime_min=30
    log: SGP_OUT + "/logs/dossiers.log"
    run:
        selected = params.selected
        if len(selected) not in (2, 3) or len(set(selected)) != len(selected):
            raise ValueError("Set report.dossier_alleles to 2–3 distinct, phenotype-blind selected candidate_addition allele IDs")
        tables = {name: list(svgp_report_rows(path)) for name, path in input.items() if name != "evidence"}
        graph = SGP["measurement"]["main_graph"]
        samples = {row["person_id"]: row for row in tables["samples"]}
        pairing = {row["person_id"]: row for row in tables["pairing"]}
        evidence = [row for row in svgp_report_rows(input.evidence) if row["allele_id"] in selected and row["graph_id"] == graph]
        lines = ["# Selected candidate dossiers", "", "Selection is prespecified; represented controls cannot be additional-allele discoveries. Coordinates are 0-based, half-open; insertion boundaries have zero reference span.", ""]
        for allele in selected:
            matched = {name: [row for row in records if row.get("allele_id") == allele] for name, records in tables.items()}
            locks = [row for row in matched["lock"] if row["graph_id"] == graph]
            if (len(matched["catalogue"]) != 1 or len(matched["panel"]) != 1 or len(locks) != 1
                    or matched["panel"][0]["panel_role"] != "candidate_addition"
                    or matched["catalogue"][0]["status"] != "PASS" or locks[0]["status"] != "LOCKED"
                    or locks[0].get("technical_status") != "PASS" or locks[0].get("locked_before_phenotypes") != "yes"
                    or len(matched["representation"]) != 3 or {r["resource"] for r in matched["representation"]} != {"HPRC", "GRCh38", "CHM13"}
                    or any(r["status"] != "NOT_REPRESENTED_AT_RESOLVED_LOCUS" for r in matched["representation"])):
                raise ValueError(f"{allele}: requires passing source, frozen candidate_addition panel, phenotype-blind technical lock and strict three-resource status")
            lines += ["## " + allele, ""]
            lines += svgp_report_section("Frozen measurement lock — " + str(input.lock), locks)
            for name in ("catalogue", "coordinates", "representation", "annotation", "frequencies"):
                lines += svgp_report_section(name + " — " + str(input[name]), matched[name])
            cells = [row for row in evidence if row["allele_id"] == allele]
            counts = []
            for label in ("independent_primary", "non_donor_not_primary", "discovery_overlap", "independence_unknown"):
                group = [r for r in cells if svgp_report_group(samples[r["person_id"]]) == label]
                called = [r for r in group if r["callable"] == "yes"]
                counts.append(dict(group=label, attempted=len(group), callable=len(called), carriers=sum(r["carrier_status"] == "ALT_SUPPORTED" for r in called), no_call=len(group)-len(called)))
            lines += svgp_report_section("Illumina callable carrier counts — " + str(input.evidence), counts)
            ont = [dict(row, analysis_group=svgp_report_group(samples[row["person_id"]]), discovery_donor=samples[row["person_id"]].get("discovery_donor", "unknown"), related_to_discovery_donor=samples[row["person_id"]].get("related_to_discovery_donor", "unknown"),
                        fingerprint_status=pairing.get(row["person_id"], {}).get("fingerprint_status", "unknown"),
                        fingerprint_evidence=pairing.get(row["person_id"], {}).get("fingerprint_evidence", "")) for row in matched["ont"]]
            lines += svgp_report_section("ONT resolved and unresolved evidence — " + str(input.ont), ont)
            for name in ("counts", "effects"):
                path = Path(SGP_OUT) / "association" / (name + ".tsv")
                records = [r for r in svgp_report_rows(path) if r["allele_id"] == allele and r.get("graph_id") == graph] if path.is_file() else []
                lines += svgp_report_section("Phenotype " + name + " (optional existing output)", records)
            lines += ["### Interpretation and future experiment", "", "Gene/regulatory overlap and proximity are context only; no functional effect is established. ONT cells without a resolved ALT carrier do not validate ALT carriage. Discovery-donor support is not independent-carrier validation. Competing explanations include linked SNPs, paralog/repeat mapping, ascertainment/batch effects and no functional effect. A decisive future experiment is native-locus allele removal/restoration in a relevant cell model, with target-gene expression and the proposed regulatory effect measured; this has not been performed.", ""]
        Path(output[0]).write_text("\n".join(lines) + "\n")
        Path(log[0]).write_text("Assembled selected passing candidate evidence; optional phenotype outputs were read only if already available.\n")
