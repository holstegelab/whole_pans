#!/usr/bin/env python3
"""Parse Kraken2 taxonomy reports with or without minimizer columns."""


def report_taxon_fields(fields, path, line_number):
    """Return (taxid, indented name) from a Kraken2 report row."""
    if len(fields) == 6:
        return fields[4], fields[5]
    if len(fields) == 8:
        return fields[6], fields[7]
    raise ValueError(
        f"Expected 6 or 8 columns in Kraken report at {path}:{line_number}, "
        f"observed {len(fields)}"
    )


def parse_kraken_report(path, target_taxids):
    """Map report taxids to configured ancestor groups and taxon names."""
    taxid_to_group = {}
    taxid_to_name = {}
    stack = []

    with open(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            taxid, raw_name = report_taxon_fields(fields, path, line_number)
            depth = (len(raw_name) - len(raw_name.lstrip(" "))) // 2
            if depth > len(stack):
                stack.extend([""] * (depth - len(stack)))
            else:
                stack = stack[:depth]

            inherited_group = stack[-1] if stack else ""
            group = target_taxids.get(taxid, inherited_group)
            stack.append(group)

            taxid_to_name[taxid] = raw_name.strip()
            if group:
                taxid_to_group[taxid] = group

    return taxid_to_group, taxid_to_name
