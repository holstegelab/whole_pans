#!/usr/bin/env python3
"""Freeze the assemblies and source graph for an existing-GFA Cactus build."""

import argparse
import csv
import gzip
import io
import os
import sys
import tempfile
from pathlib import Path

from minigraph_cactus_common import atomic_bytes, canonical, file_sha256, read_json, sha256
from prepare_minigraph_cactus_seqfile import ASSEMBLY, required_path, scan_fasta


def file_info(path, content_hash=True):
    path = required_path(str(path))
    before = path.stat()
    digest = file_sha256(path) if content_hash else None
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"file changed while reading: {path}")
    return {"path": str(path), "size": after.st_size, "mtime_ns": after.st_mtime_ns,
            "sha256": digest}


def local_rows(manifest):
    rows = []
    prefixes = {}
    with open(manifest, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"order", "role", "assembly_id", "path", "header_prefix"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"local V1 manifest lacks {sorted(required)}")
        for line_no, row in enumerate(reader, 2):
            if row["role"] != "assembly":
                continue
            path = required_path(row["path"])
            match = ASSEMBLY.fullmatch(path.name)
            if not match:
                raise ValueError(f"unexpected local assembly path at {manifest}:{line_no}: {path}")
            sample, hap = match.group(1), match.group(2)
            if "." in sample or any(char.isspace() for char in sample):
                raise ValueError(f"unsafe local sample name at {manifest}:{line_no}")
            label = f"{sample}.{hap}"
            prefix = row["header_prefix"]
            if not prefix or "#" in prefix or any(char.isspace() for char in prefix):
                raise ValueError(f"unsafe local GFA prefix at {manifest}:{line_no}")
            if prefix in prefixes or any(existing[0] == label for existing in rows):
                raise ValueError(f"duplicate local sample or GFA prefix: {label}, {prefix}")
            prefixes[prefix] = label
            rows.append((label, sample, hap, path, line_no))
    if not rows:
        raise ValueError("local V1 manifest contains no assemblies")
    return rows, prefixes


def hprc_rows(manifest, downloads):
    rows = []
    with open(manifest, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Sample ID", "Haplotype", "Filename"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"HPRC assembly manifest lacks {sorted(required)}")
        for line_no, row in enumerate(reader, 2):
            sample, hap = row["Sample ID"].strip(), row["Haplotype"].strip()
            filename = row["Filename"].strip()
            if not sample or "." in sample or "#" in sample or any(c.isspace() for c in sample):
                raise ValueError(f"unsafe HPRC sample at {manifest}:{line_no}")
            if hap not in {"0", "1", "2"} or Path(filename).name != filename:
                raise ValueError(f"unexpected HPRC haplotype or filename at {manifest}:{line_no}")
            label = sample if hap == "0" else f"{sample}.{hap}"
            rows.append((label, sample, hap, required_path(str(downloads / filename)), line_no))
    labels = [row[0] for row in rows]
    if len(labels) != len(set(labels)) or labels.count("CHM13") != 1 or labels.count("GRCh38") != 1:
        raise ValueError("HPRC manifest needs unique haplotypes, CHM13 and GRCh38")
    references = [next(row for row in rows if row[0] == label) for label in ("CHM13", "GRCh38")]
    return references + [row for row in rows if row[0] not in {"CHM13", "GRCh38"}]


def write_preflight(args):
    source = file_info(args.seed_gfa)
    local_manifest = file_info(args.local_manifest)
    local, prefixes = local_rows(args.local_manifest)
    local_prefixes = {label: prefix for prefix, label in prefixes.items()}
    manifests = [local_manifest]
    if args.hprc_manifest:
        if not args.hprc_downloads:
            raise ValueError("HPRC downloads directory is required with the HPRC manifest")
        hprc_manifest = file_info(args.hprc_manifest)
        manifests.append(hprc_manifest)
        rows = hprc_rows(args.hprc_manifest, Path(args.hprc_downloads)) + local
        cohort = "hprc_plus_local_v1"
    else:
        if not args.chm13 or not args.grch38:
            raise ValueError("local V1 requires CHM13 and GRCh38 FASTAs")
        rows = [("CHM13", "CHM13", "0", required_path(args.chm13), 0),
                ("GRCh38", "GRCh38", "0", required_path(args.grch38), 0)] + local
        cohort = "local_v1"
    labels = [row[0] for row in rows]
    if len(labels) != len(set(labels)):
        raise ValueError("HPRC and local seqfile sample names collide")
    paths = [str(row[3]) for row in rows]
    if len(paths) != len(set(paths)):
        raise ValueError("the same FASTA is assigned to more than one sample")
    seqfile = "".join(f"{label} {path}\n" for label, _, _, path, _ in rows).encode()
    metadata = Path(args.preflight_json).parent
    metadata.mkdir(parents=True, exist_ok=True)
    fd, staged_map = tempfile.mkstemp(prefix=".contig-map.", dir=metadata)
    os.close(fd)
    inputs = []
    refs = {}
    try:
        with open(staged_map, "wb") as raw_map:
            with gzip.GzipFile(fileobj=raw_map, mode="wb", filename="", mtime=0) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text_map:
                    writer = csv.writer(text_map, delimiter="\t")
                    writer.writerow(["seqfile_name", "original_header", "cactus_contig_id", "length_bp",
                                     "gfa_source_name", "pansn_source_name"])
                    source_names = set()
                    for label, sample, hap, path, line_no in rows:
                        class SourceWriter:
                            def writerow(self, values):
                                _, original, effective, length = values
                                raw_name = original.split()[0]
                                if label in local_prefixes and not raw_name.startswith(local_prefixes[label] + "."):
                                    raise ValueError(f"local FASTA contig lacks frozen V1 prefix: {raw_name}")
                                if cohort == "hprc_plus_local_v1" and label not in local_prefixes:
                                    if not raw_name.startswith(f"{sample}#{hap}#"):
                                        raise ValueError(f"HPRC FASTA contig has wrong PanSN sample: {raw_name}")
                                source_name = (f"{label}.{raw_name}" if cohort == "local_v1" and
                                               label in {"CHM13", "GRCh38"} else raw_name)
                                if source_name in source_names:
                                    raise ValueError(f"duplicate GFA source name: {source_name}")
                                source_names.add(source_name)
                                pansn = f"{sample}#{hap}#{effective}"
                                writer.writerow([label, original, effective, length, source_name, pansn])

                        stats, contigs = scan_fasta(path, label, SourceWriter())
                        inputs.append({"name": label, "sample": sample, "haplotype": hap,
                                       "path": str(path), "source_line": line_no, **stats})
                        if label in {"CHM13", "GRCh38"}:
                            refs[label] = contigs
        preflight = {"preflight_version": 3, "cohort": cohort,
                     "included_list": local_manifest, "source_manifests": manifests[1:],
                     "seed_gfa": source, "seed_prefixes": prefixes,
                     "inputs": inputs, "seqfile_sha256": sha256(seqfile),
                     "haplotype_count": len(rows) - 2, "local_haplotype_count": len(local),
                     "reference_contigs": refs,
                     "total_input_bases": sum(row["total_bases"] for row in inputs)}
        preflight["identity_sha256"] = sha256(canonical(preflight))
        preflight_hash = sha256(canonical(preflight))
        locked = metadata / "preflight.locked.json"
        states = [*metadata.glob("*.state.json"), *metadata.glob("*.success.json")]
        if locked.exists() and read_json(locked) != preflight:
            raise ValueError("seeded build ID is occupied by different inputs")
        if states and not locked.exists():
            raise ValueError("stage state exists without locked preflight")
        for state in states:
            if read_json(state).get("preflight_sha256") != preflight_hash:
                raise ValueError("existing stage state has different inputs")
        manifest_text = io.StringIO()
        writer = csv.writer(manifest_text, delimiter="\t", lineterminator="\n")
        writer.writerow(["seqfile_name", "sample", "haplotype", "source_path", "source_list_line",
                         "size", "mtime_ns", "sha256_uncompressed", "contig_count", "total_bases"])
        for row in inputs:
            writer.writerow([row["name"], row["sample"], row["haplotype"], row["path"],
                             row["source_line"], row["size"], row["mtime_ns"], row["sha256"],
                             row["contig_count"], row["total_bases"]])
        atomic_bytes(args.seqfile, seqfile)
        atomic_bytes(args.input_manifest, manifest_text.getvalue().encode())
        os.replace(staged_map, args.contig_map)
        atomic_bytes(args.preflight_json, canonical(preflight) + b"\n")
    finally:
        if os.path.exists(staged_map):
            os.unlink(staged_map)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("seed-gfa", "local-manifest", "seqfile", "input-manifest", "contig-map", "preflight-json"):
        parser.add_argument(f"--{key}", required=True)
    for key in ("hprc-manifest", "hprc-downloads", "chm13", "grch38"):
        parser.add_argument(f"--{key}")
    write_preflight(parser.parse_args())


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, EOFError, KeyError, UnicodeError, gzip.BadGzipFile) as exc:
        sys.exit(f"Minigraph-Cactus seed preflight failed: {exc}")
