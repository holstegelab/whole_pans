# Sourced by explicit Snakemake shell rules; never launches bioinformatics tools.
# Each job owns exactly one mktemp directory, and cleanup removes only that one.
one_kcp_scratch() {
    local scratch_root="$1" required_mb="$2" usage_file="$3"
    mkdir -p "$scratch_root" "$(dirname "$usage_file")"
    local free_kb
    free_kb=$(df -Pk "$scratch_root" | awk 'NR==2 {print $4}')
    if (( free_kb < required_mb * 1000 )); then
        echo "Insufficient scratch: need ${required_mb} MB at $scratch_root" >&2
        return 1
    fi
    ONE_KCP_TMP=$(mktemp -d "$scratch_root/one-kcp.XXXXXXXXXX")
    export TMPDIR="$ONE_KCP_TMP"
    ONE_KCP_USAGE="$usage_file"
    : > "$ONE_KCP_USAGE.samples"
    (
        while :; do
            du -sk "$ONE_KCP_TMP" | cut -f1 >> "$ONE_KCP_USAGE.samples"
            sleep 30
        done
    ) &
    ONE_KCP_MONITOR=$!
    trap 'one_kcp_cleanup' EXIT
    trap 'exit 130' INT
    trap 'exit 143' TERM
}
one_kcp_cleanup() {
    local status=$?
    trap - EXIT
    kill "$ONE_KCP_MONITOR" 2>/dev/null || true
    wait "$ONE_KCP_MONITOR" 2>/dev/null || true
    du -sk "$ONE_KCP_TMP" | cut -f1 >> "$ONE_KCP_USAGE.samples" || true
    awk 'BEGIN {print "sampled_peak_kib"} {if ($1>m) m=$1} END {print m+0}' "$ONE_KCP_USAGE.samples" > "$ONE_KCP_USAGE"
    rm -f -- "$ONE_KCP_USAGE.samples"
    rm -rf -- "$ONE_KCP_TMP"
    exit "$status"
}
