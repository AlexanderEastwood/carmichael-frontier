#!/usr/bin/env bash
# run_excl.sh <exponent> [KMAX] : complete partitioned k=64 exclusion below 10^exponent
set -u
cd "$(dirname "$0")"
E="${1:?exponent}"; KMAX="${2:-20}"
B=1$(printf "0%.0s" $(seq 1 "$E"))
JD="jobs_1e${E}"; rm -rf "$JD"; mkdir -p "$JD"
python3 gen_jobs.py > "$JD/JOBS"
NJ=$(wc -l < "$JD/JOBS")
echo "exponent=$E kmax=$KMAX jobs=$NJ start=$(date -u +%FT%TZ)" > "$JD/PROGRESS"
run_one(){ f="$JD/$(echo "$1"|tr , _).out"; DYNT_KMAX="$KMAX" ./cm_dynt 64 "$B" --prefix "$1" > "$f" 2>&1; }
export -f run_one; export B JD KMAX
xargs -P 30 -I{} bash -c 'run_one "$@"' _ {} < "$JD/JOBS"
hits=$(grep -l "^n=[0-9]" "$JD"/*.out 2>/dev/null)
{ if [ -n "$hits" ]; then echo "RESULT=HIT"; echo "$hits"; else echo "RESULT=ALL_NONE"; fi; } > "$JD/DONE"
echo "done=$(date -u +%FT%TZ) $(head -1 "$JD/DONE")" >> "$JD/PROGRESS"
