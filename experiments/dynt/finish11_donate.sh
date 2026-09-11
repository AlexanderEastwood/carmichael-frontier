#!/usr/bin/env bash
# Donate idle cores to the last outstanding depth-18 chain child WITHOUT killing its running workers:
# run its feasibility-pruned +4 children (SUBJOBS11) in parallel under DYNT_KMAX=0 (measured faster than 20 on
# these dense slices). The verifier accepts either coverage (own n=NONE or all children n=NONE). No kills.
set -u; cd ~/carmichael-frontier/experiments/dynt; B=1$(printf "0%.0s" $(seq 1 147)); JD=jobs_1e147; KMAX=${KMAX:-0}; NP=${NP:-22}
[ -s "$JD/SUBJOBS11" ] || { echo "no SUBJOBS11" >&2; exit 1; }
echo "subjobs11=$(wc -l < "$JD/SUBJOBS11") kmax=$KMAX np=$NP start=$(date -u +%FT%TZ)" >> "$JD/PROGRESS"
run_one(){ f="$JD/$(echo "$1"|tr , _).out"; [ -s "$f" ] && grep -q "^n=" "$f" && return 0; DYNT_KMAX="$KMAX" ./cm_dynt 64 "$B" --prefix "$1" > "$f" 2>&1; }
export -f run_one; export B JD KMAX
xargs -P "$NP" -I{} bash -c "run_one \"\$@\"" _ {} < "$JD/SUBJOBS11"
echo "frontier11 done=$(date -u +%FT%TZ)" >> "$JD/PROGRESS"
