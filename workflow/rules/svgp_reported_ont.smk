"""Retain reported ONT sample/run assignments before read and identity review."""


rule reported_ont:
    input: SGP_OUT + "/metadata/reported_ont_runs.tsv"


rule svgp_reported_ont_inventory:
    input:
        manifest=SGP["inventory"].get("ont_runs", "config/svgp_ont_runs.tsv")
    output:
        SGP_OUT + "/metadata/reported_ont_runs.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime_min=5
    log: SGP_OUT + "/logs/reported_ont_runs.log"
    run:
        import hashlib
        import os

        source, target = Path(input.manifest), Path(output[0])
        if (os.path.abspath(source) == os.path.abspath(target)
                or (target.exists() and source.samefile(target))):
            raise ValueError("Reported ONT inventory output would overwrite its source manifest")
        required = ["sample_id", "ont_run", "assignment_source", "reported_date"]
        with source.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if (not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames))
                    or not set(required).issubset(reader.fieldnames)):
                raise ValueError(f"Missing or duplicate reported ONT manifest columns: {source}")
            rows, seen = [], set()
            for number, row in enumerate(reader, 2):
                if None in row or any(value is None for value in row.values()):
                    raise ValueError(f"{source}:{number}: malformed TSV row")
                if any(not row[field].strip() for field in required):
                    raise ValueError(f"{source}:{number}: empty sample, run or provenance field")
                key = (row["sample_id"], row["ont_run"])
                if key in seen:
                    raise ValueError(f"{source}:{number}: duplicate sample/run assignment: {key}")
                seen.add(key)
                rows.append({field: row[field] for field in required})
        manifest_sha256 = hashlib.sha256(source.read_bytes()).hexdigest()
        for row in rows:
            row.update(availability_status="USER_REPORTED", person_id="unknown",
                       read_path="unknown", identity_status="UNRESOLVED",
                       read_status="NOT_PROVIDED", fingerprint_status="NOT_ASSESSED",
                       measurement_status="NOT_READY", source_manifest=str(source),
                       source_manifest_sha256=manifest_sha256)
        fields = required + ["availability_status", "person_id", "read_path", "identity_status", "read_status",
                             "fingerprint_status", "measurement_status", "source_manifest",
                             "source_manifest_sha256"]
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        Path(log[0]).parent.mkdir(parents=True, exist_ok=True)
        Path(log[0]).write_text(
            f"Reported sample labels: {len({row['sample_id'] for row in rows})}\n"
            f"Reported sample/run assignments: {len(rows)}\n"
            f"Distinct run labels: {len({row['ont_run'] for row in rows})}\n"
            "People, read files and cross-platform fingerprints remain unverified.\n")
