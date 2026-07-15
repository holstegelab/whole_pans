SV_PANGENOME_CONFIG = workflow.source_path("../../config/config.yaml")
configfile: SV_PANGENOME_CONFIG


import os
from pathlib import Path
from pathlib import PurePosixPath


SV_PANGENOME = config.get("pangenome", {})
SV_PANGENOME_ASSEMBLIES = SV_PANGENOME.get("assemblies", {})
SV_PANGENOME_REFERENCES = SV_PANGENOME.get("references", {})
SV_PANGENOME_PROJECT_ROOT = Path(workflow.current_basedir).resolve().parents[1]
workdir: str(SV_PANGENOME_PROJECT_ROOT)

SV_PANGENOME_PATHS = {
    "results": SV_PANGENOME.get("results", "results/sv_pangenome"),
    "assemblies.included_list": SV_PANGENOME_ASSEMBLIES.get(
        "included_list",
        "results/post_decontamination_qc/summary/graph_included_assemblies.txt",
    ),
    "references.CHM13": SV_PANGENOME_REFERENCES.get("CHM13", "../refs/CHM13.fa.gz"),
    "references.hg38": SV_PANGENOME_REFERENCES.get(
        "hg38", "../refs/GRCh38_primary_chr1-22_XY.fa"
    ),
    "vg": SV_PANGENOME.get("vg", "../vg"),
}
for key, value in SV_PANGENOME_PATHS.items():
    if os.path.isabs(str(value)):
        raise WorkflowError(
            f"pangenome.{key} must be relative to the Snakefile, got: {value}"
        )

SV_PANGENOME_OUTDIR = str(PurePosixPath(str(SV_PANGENOME_PATHS["results"])))
SV_PANGENOME_INCLUDED_LIST = str(
    PurePosixPath(str(SV_PANGENOME_PATHS["assemblies.included_list"]))
)
SV_PANGENOME_CHM13 = str(PurePosixPath(str(SV_PANGENOME_PATHS["references.CHM13"])))
SV_PANGENOME_HG38 = str(PurePosixPath(str(SV_PANGENOME_PATHS["references.hg38"])))
SV_PANGENOME_VG = str(PurePosixPath(str(SV_PANGENOME_PATHS["vg"])))
SV_PANGENOME_ENV = str(workflow.source_path("../envs/pangenome.yaml"))

SV_PANGENOME_ORDER_BY_MASH = str(
    SV_PANGENOME.get("order_by_mash", "ascending")
).lower()
if SV_PANGENOME_ORDER_BY_MASH not in {"ascending", "descending"}:
    raise WorkflowError("pangenome.order_by_mash must be 'ascending' or 'descending'")

# The 1KCP construction code sorts Mash distance to CHM13 ascending before
# Minigraph; hg38 is still forced into position 2 for GRCh38 compatibility.
SV_PANGENOME_SORT_REVERSE = (
    "r" if SV_PANGENOME_ORDER_BY_MASH == "descending" else ""
)


rule build_sv_pangenome_graph:
    input:
        chm13=SV_PANGENOME_CHM13,
        hg38=SV_PANGENOME_HG38,
        assemblies=SV_PANGENOME_INCLUDED_LIST
    output:
        gfa=f"{SV_PANGENOME_OUTDIR}/graphs/sv_pangenome.minigraph.gfa",
        ordered=f"{SV_PANGENOME_OUTDIR}/metadata/sv_pangenome.ordered_assemblies.tsv",
        mash_distances=f"{SV_PANGENOME_OUTDIR}/metadata/chm13.mash_distances.tsv",
        versions=f"{SV_PANGENOME_OUTDIR}/metadata/tool_versions.tsv",
        summary=f"{SV_PANGENOME_OUTDIR}/metadata/sv_pangenome.graph_summary.tsv"
    params:
        vg=SV_PANGENOME_VG,
        minigraph_preset=SV_PANGENOME.get("minigraph_preset", "ggs"),
        mash_kmer_size=int(SV_PANGENOME.get("mash_kmer_size", 21)),
        mash_sketch_size=int(SV_PANGENOME.get("mash_sketch_size", 100000)),
        sort_reverse=SV_PANGENOME_SORT_REVERSE
    log:
        f"{SV_PANGENOME_OUTDIR}/logs/build_sv_pangenome_graph.log"
    benchmark:
        f"{SV_PANGENOME_OUTDIR}/benchmarks/build_sv_pangenome_graph.tsv"
    conda:
        SV_PANGENOME_ENV
    threads: int(SV_PANGENOME.get("threads", 16))
    resources:
        mem_mb=int(SV_PANGENOME.get("mem_mb", 240000)),
        runtime_min=int(SV_PANGENOME.get("runtime_min", 2880))
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {output.gfa:q})" \
                 "$(dirname {output.ordered:q})" \
                 "$(dirname {output.mash_distances:q})" \
                 "$(dirname {output.versions:q})" \
                 "$(dirname {output.summary:q})" \
                 "$(dirname {log:q})"
        : > {log:q}

        if [ ! -s {input.chm13:q} ]; then
            echo "CHM13 FASTA is missing or empty: {input.chm13}" >&2
            exit 1
        fi
        if [ ! -s {input.hg38:q} ]; then
            echo "hg38 FASTA is missing or empty: {input.hg38}" >&2
            exit 1
        fi
        if [ ! -s {input.assemblies:q} ]; then
            echo "Assembly inclusion list is missing or empty: {input.assemblies}" >&2
            exit 1
        fi

        scratch_root="${{TMPDIR:-$(dirname {output.gfa:q})}}"
        if [ ! -d "$scratch_root" ] || [ ! -w "$scratch_root" ]; then
            echo "TMPDIR is not a writable directory: $scratch_root" >&2
            exit 1
        fi
        workdir=$(mktemp -d -p "$scratch_root" "sv-pangenome.XXXXXX")
        trap 'rm -rf "$workdir"' EXIT
        echo "SV pangenome scratch directory: $workdir" >> {log:q}

        {{
            printf 'tool\tversion\n'
            printf 'minigraph\t'
            {{ minigraph --version 2>&1 || minigraph -V 2>&1 || minigraph 2>&1 || true; }} \
              | awk 'NF {{print; found=1; exit}} END {{if (!found) print "unavailable"}}'
            printf 'mash\t'
            {{ mash --version 2>&1 || mash 2>&1 || true; }} \
              | awk 'NF {{print; found=1; exit}} END {{if (!found) print "unavailable"}}'
            printf 'vg\t'
            {{ {params.vg:q} version 2>&1 || {params.vg:q} --version 2>&1 || true; }} \
              | awk 'NF {{print; found=1; exit}} END {{if (!found) print "unavailable"}}'
        }} > {output.versions:q}

        fasta_id() {{
            local name
            name=$(basename "$1")
            case "$name" in
                *.fasta.gz) printf '%s\n' "${{name%.fasta.gz}}" ;;
                *.fna.gz) printf '%s\n' "${{name%.fna.gz}}" ;;
                *.fa.gz) printf '%s\n' "${{name%.fa.gz}}" ;;
                *.fasta) printf '%s\n' "${{name%.fasta}}" ;;
                *.fna) printf '%s\n' "${{name%.fna}}" ;;
                *.fa) printf '%s\n' "${{name%.fa}}" ;;
                *) printf '%s\n' "$name" ;;
            esac
        }}

        assembly_paths="$workdir/assembly_paths.txt"
        : > "$assembly_paths"
        declare -A seen_paths
        declare -A seen_ids
        assembly_count=0
        while IFS= read -r assembly || [ -n "$assembly" ]; do
            assembly="${{assembly%$'\r'}}"
            [ -z "$assembly" ] && continue
            case "$assembly" in \#*) continue ;; esac
            if [ "$assembly" = {input.chm13:q} ] || [ "$assembly" = {input.hg38:q} ]; then
                continue
            fi
            if [ -n "${{seen_paths[$assembly]+x}}" ]; then
                echo "Duplicate assembly path in {input.assemblies}: $assembly" >&2
                exit 1
            fi
            assembly_id=$(fasta_id "$assembly")
            if [ -n "${{seen_ids[$assembly_id]+x}}" ]; then
                echo "Duplicate assembly ID after FASTA suffix stripping: $assembly_id" >&2
                exit 1
            fi
            if [ ! -s "$assembly" ]; then
                echo "Assembly FASTA is missing or empty: $assembly" >&2
                exit 1
            fi
            seen_paths[$assembly]=1
            seen_ids[$assembly_id]=1
            printf '%s\n' "$assembly" >> "$assembly_paths"
            assembly_count=$((assembly_count + 1))
        done < {input.assemblies:q}

        if [ "$assembly_count" -eq 0 ]; then
            echo "No non-reference assemblies found in {input.assemblies}" >&2
            exit 1
        fi
        echo "Assemblies to add after CHM13 and hg38: $assembly_count" >> {log:q}

        mash_prefix="$workdir/chm13"
        mash sketch \
          -p {threads} \
          -k {params.mash_kmer_size} \
          -s {params.mash_sketch_size} \
          -o "$mash_prefix" \
          {input.chm13:q} >> {log:q} 2>&1

        printf 'assembly_id\tassembly_path\tmash_distance\tp_value\tmatching_hashes\n' \
          > {output.mash_distances:q}
        while IFS= read -r assembly; do
            assembly_id=$(fasta_id "$assembly")
            mash dist "$mash_prefix.msh" "$assembly" 2>> {log:q} \
              | awk -v id="$assembly_id" -v path="$assembly" 'BEGIN {{OFS="\t"}} NF >= 5 {{print id, path, $3, $4, $5; found=1; exit}} END {{if (!found) exit 1}}' \
              >> {output.mash_distances:q}
        done < "$assembly_paths"

        printf 'order\trole\tassembly_id\tpath\tmash_distance\theader_prefix\n' > {output.ordered:q}
        printf '1\treference\tCHM13\t%s\tNA\tCHM13\n' {input.chm13:q} >> {output.ordered:q}
        printf '2\tcompatibility_path\thg38\t%s\tNA\tGRCh38\n' {input.hg38:q} >> {output.ordered:q}
        order=3
        tail -n +2 {output.mash_distances:q} \
          | sort -t $'\t' -k3,3g{params.sort_reverse} \
          | while IFS=$'\t' read -r assembly_id assembly_path mash_distance p_value matching_hashes; do
                header_prefix="${{assembly_id%.clean}}"
                printf '%s\tassembly\t%s\t%s\t%s\t%s\n' \
                  "$order" "$assembly_id" "$assembly_path" "$mash_distance" "$header_prefix" \
                  >> {output.ordered:q}
                order=$((order + 1))
            done

        mapfile -t fasta_order < <(awk -F '\t' 'NR > 1 {{print $4}}' {output.ordered:q})
        if [ "${{#fasta_order[@]}}" -lt 3 ]; then
            echo "Ordered FASTA list must contain CHM13, hg38, and at least one assembly" >&2
            exit 1
        fi

        prefix_reference_fasta() {{
            local source="$1"
            local prefix="$2"
            local destination="$3"
            case "$source" in
                *.gz) gzip -cd "$source" ;;
                *) cat "$source" ;;
            esac | awk -v prefix="$prefix" '
                /^>/ {{
                    header = substr($0, 2)
                    split(header, fields, /[[:space:]]/)
                    suffix = substr(header, length(fields[1]) + 1)
                    print ">" prefix "." fields[1] suffix
                    next
                }}
                {{print}}
            ' > "$destination"
            test -s "$destination"
        }}

        graph_chm13="$workdir/CHM13.graph.fa"
        graph_hg38="$workdir/GRCh38.graph.fa"
        prefix_reference_fasta {input.chm13:q} CHM13 "$graph_chm13"
        prefix_reference_fasta {input.hg38:q} GRCh38 "$graph_hg38"
        fasta_order[0]="$graph_chm13"
        fasta_order[1]="$graph_hg38"

        tmp_gfa="$workdir/sv_pangenome.minigraph.gfa"
        minigraph -cx{params.minigraph_preset} -t {threads} "${{fasta_order[@]}}" \
          > "$tmp_gfa" 2>> {log:q}
        test -s "$tmp_gfa"
        if grep -Eq 'inconsistent rGFA|associated with different ranks' {log:q}; then
            echo "Minigraph reported inconsistent rGFA source names; refusing graph output" \
              >> {log:q}
            exit 1
        fi
        mv "$tmp_gfa" {output.gfa:q}

        awk -F '\t' '
            BEGIN {{OFS="\t"; print "metric", "value"}}
            $1 == "S" {{segments++; segment_bp += length($3)}}
            $1 == "L" {{links++}}
            $1 == "P" || $1 == "W" {{paths++}}
            END {{
                print "segments", segments + 0
                print "links", links + 0
                print "paths", paths + 0
                print "segment_bp", segment_bp + 0
            }}
        ' {output.gfa:q} > {output.summary:q}
        test -s {output.summary:q}
        """
