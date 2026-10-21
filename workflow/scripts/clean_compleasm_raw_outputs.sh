#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: clean_compleasm_raw_outputs.sh [OPTIONS]

Remove obsolete per-assembly Compleasm primates_odb12 directories while
preserving each summary.txt and the aggregate assembly_qc.tsv.

The default is a read-only inventory. Nothing is removed unless --delete is
given.

Options:
  --compleasm-dir PATH  Override the Compleasm results directory.
  --delete              Delete validated primates_odb12 directories.
  --verbose             Print every validated directory.
  -h, --help            Show this help message.
EOF
}

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
repository_root=$(cd -- "$script_dir/../.." && pwd -P)
compleasm_dir="$repository_root/../whole_pangenome/assembly_qc/results/compleasm"
delete=false
verbose=false

while (($# > 0)); do
    case "$1" in
        --compleasm-dir)
            if (($# < 2)); then
                echo "ERROR: --compleasm-dir requires a path" >&2
                exit 2
            fi
            compleasm_dir=$2
            shift 2
            ;;
        --delete)
            delete=true
            shift
            ;;
        --verbose)
            verbose=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [[ ! -d "$compleasm_dir" ]]; then
    echo "ERROR: Compleasm results directory does not exist: $compleasm_dir" >&2
    exit 1
fi

compleasm_dir=$(cd -- "$compleasm_dir" && pwd -P)
assembly_qc_root=$(cd -- "$compleasm_dir/../.." && pwd -P)
aggregate_qc="$assembly_qc_root/results/summary/assembly_qc.tsv"

if [[ "$compleasm_dir" == "/" || "$compleasm_dir" == "$repository_root" ]]; then
    echo "ERROR: refusing unsafe Compleasm directory: $compleasm_dir" >&2
    exit 1
fi

if [[ ! -s "$aggregate_qc" ]]; then
    echo "ERROR: aggregate QC table is missing or empty: $aggregate_qc" >&2
    exit 1
fi

aggregate_header=$(head -n 1 -- "$aggregate_qc")
if [[ "${aggregate_header%%$'\t'*}" != "assembly_id" ]]; then
    echo "ERROR: first column of $aggregate_qc is not assembly_id" >&2
    exit 1
fi

lock_dirs=(
    "$repository_root/.snakemake/locks"
    "$assembly_qc_root/.snakemake/locks"
)

assert_no_snakemake_locks() {
    local lock_dir lock_entry
    for lock_dir in "${lock_dirs[@]}"; do
        [[ -d "$lock_dir" ]] || continue
        lock_entry=$(find "$lock_dir" -mindepth 1 -maxdepth 1 -print -quit)
        if [[ -n "$lock_entry" ]]; then
            echo "ERROR: Snakemake lock found: $lock_entry" >&2
            echo "Stop the workflow and resolve the lock before deleting files." >&2
            return 1
        fi
    done
}

assembly_is_in_aggregate() {
    local assembly_id=$1
    awk -F '\t' -v assembly_id="$assembly_id" '
        NR > 1 && $1 == assembly_id { found = 1; exit }
        END { exit(found ? 0 : 1) }
    ' "$aggregate_qc"
}

validate_candidate() {
    local candidate=$1
    local sample_dir sample_id summary

    if [[ ! -d "$candidate" || -L "$candidate" ]]; then
        echo "ERROR: candidate is not a normal directory: $candidate" >&2
        return 1
    fi
    if [[ "${candidate##*/}" != "primates_odb12" ]]; then
        echo "ERROR: unexpected candidate name: $candidate" >&2
        return 1
    fi

    sample_dir=${candidate%/*}
    if [[ "${sample_dir%/*}" != "$compleasm_dir" ]]; then
        echo "ERROR: candidate is outside the expected sample level: $candidate" >&2
        return 1
    fi

    sample_id=${sample_dir##*/}
    summary="$sample_dir/summary.txt"
    if [[ ! -s "$summary" ]]; then
        echo "ERROR: preserving candidate because summary.txt is missing or empty: $candidate" >&2
        return 1
    fi
    if ! assembly_is_in_aggregate "$sample_id"; then
        echo "ERROR: preserving candidate because $sample_id is absent from $aggregate_qc" >&2
        return 1
    fi
}

assert_no_snakemake_locks

mapfile -d '' candidates < <(
    find "$compleasm_dir" \
        -mindepth 2 -maxdepth 2 \
        -type d -name primates_odb12 \
        -print0 | LC_ALL=C sort -z
)

if ((${#candidates[@]} == 0)); then
    echo "No per-assembly primates_odb12 directories found under: $compleasm_dir"
    exit 0
fi

validated=()
total_entries=0
for candidate in "${candidates[@]}"; do
    validate_candidate "$candidate"
    validated+=("$candidate")

    read -r entry_count < <(find "$candidate" -printf '.' | wc -c)
    total_entries=$((total_entries + entry_count))
    if [[ "$verbose" == true ]]; then
        printf 'SAFE\t%s entries\t%s\n' "$entry_count" "$candidate"
    fi
done

printf 'Validated %d raw Compleasm directories containing %d filesystem entries.\n' \
    "${#validated[@]}" "$total_entries"
printf 'Preserved per-sample summary.txt files and aggregate table: %s\n' "$aggregate_qc"

if [[ "$delete" != true ]]; then
    echo "DRY RUN: nothing deleted. Re-run with --delete after all related jobs have stopped."
    exit 0
fi

assert_no_snakemake_locks

removed=0
for candidate in "${validated[@]}"; do
    # Revalidate immediately before each destructive operation.
    validate_candidate "$candidate"
    assert_no_snakemake_locks
    printf 'DELETE\t%s\n' "$candidate"
    rm -rf -- "$candidate"
    removed=$((removed + 1))
done

printf 'Deleted %d raw Compleasm directories (%d inventoried filesystem entries).\n' \
    "$removed" "$total_entries"
