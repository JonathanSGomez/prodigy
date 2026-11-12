#!/usr/bin/env bash

# ============================================================
# Parallel runner for make_truth_areas.py
# Robust version: no associative PID array and no "set -u".
# ============================================================

set -o pipefail

resolution="3k"
smass_cut="1.0e9"
mass_label="Mstar1e9"
step="2"

# Number of regions running simultaneously.
# You can override it when launching:
#   MAX_JOBS=50 ./run_make_truth_areas_all_fixed.sh
max_jobs="${MAX_JOBS:-50}"

# The Python code is assumed to be in the same directory as this script.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
code="${script_dir}/make_truth_areas.py"

logdir="${script_dir}/logfiles/truth_areas/${resolution}/${mass_label}/step${step}"
mkdir -p "$logdir"

# Set to true to recalculate existing outputs.
overwrite=true

# Prevent NumPy/BLAS from using several threads inside each Python process.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Valid 7k regions: 1...324, excluding region 107.
regions=()
for region in $(seq 1 324); do
    if [[ "$region" -ne 107 ]]; then
        regions+=("$region")
    fi
done

status_dir="${logdir}/run_status"
mkdir -p "$status_dir"

success_file="${status_dir}/success.txt"
failed_file="${status_dir}/failed.txt"

: > "$success_file"
: > "$failed_file"

echo ""
echo "=============================================="
echo "make_truth_areas parallel run"
echo "resolution : $resolution"
echo "smass_cut  : $smass_cut"
echo "step       : $step"
echo "N regions  : ${#regions[@]}"
echo "max_jobs   : $max_jobs"
echo "code       : $code"
echo "logdir     : $logdir"
echo "=============================================="
echo ""

if [[ ! -f "$code" ]]; then
    echo "ERROR: no existe el código:"
    echo "  $code"
    exit 1
fi

start_time=$(date +%s)
n_started=0

running_jobs() {
    jobs -pr | wc -l
}

launch_region() {
    local region="$1"
    local region4
    local logfile

    region4=$(printf "%04d" "$region")
    logfile="${logdir}/truth_areas.${resolution}.${mass_label}.step${step}.region${region4}.log"

    (
        cmd=(
            python "$code" "$resolution"
            --regions "$region"
            --smass-cut "$smass_cut"
            --step "$step"
        )

        if [[ "$overwrite" == "true" ]]; then
            cmd+=(--overwrite)
        fi

        "${cmd[@]}" > "$logfile" 2>&1
        status=$?

        if [[ "$status" -eq 0 ]]; then
            echo "$region" >> "$success_file"
        else
            echo "$region" >> "$failed_file"
        fi

        exit "$status"
    ) &

    pid=$!
    n_started=$((n_started + 1))

    printf "started region %4d | pid=%d | running=%d/%d\n" \
        "$region" "$pid" "$(running_jobs)" "$max_jobs"
}

for region in "${regions[@]}"; do

    while [[ "$(running_jobs)" -ge "$max_jobs" ]]; do
        sleep 1
    done

    launch_region "$region"
done

# Wait for all remaining jobs.
wait

end_time=$(date +%s)
elapsed=$((end_time - start_time))

n_success=$(wc -l < "$success_file")
n_failed=$(wc -l < "$failed_file")
n_finished=$((n_success + n_failed))

echo ""
echo "=============================================="
echo "END"
echo "started   : $n_started"
echo "finished  : $n_finished"
echo "successful: $n_success"
echo "failed    : $n_failed"
echo "elapsed   : ${elapsed} s"
echo "logs      : $logdir"

if [[ "$n_failed" -gt 0 ]]; then
    echo -n "failed regions: "
    tr '\n' ' ' < "$failed_file"
    echo ""
else
    echo "failed regions: none"
fi

echo "=============================================="
