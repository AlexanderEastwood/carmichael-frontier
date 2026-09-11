#!/usr/bin/env bash
# Donate idle cores to the 5 hard depth-14 chain slices WITHOUT killing their running workers:
# run their feasibility-pruned +4 children (SUBJOBS9, generated once by sub_gen_feasible.py) in
# parallel. The verifier accepts either coverage (own n=NONE or all children n=NONE). No kills.
set -u; cd ~/carmichael-frontier/experiments/dynt; B=1$(printf "0%.0s" $(seq 1 147)); JD=jobs_1e147; KMAX=${KMAX:-20}; NP=${NP:-22}
[ -s "$JD/SUBJOBS9" ] || { echo "no SUBJOBS9" >&2; exit 1; }
echo "subjobs9=$(wc -l < "$JD/SUBJOBS9") kmax=$KMAX np=$NP start=$(date -u +%FT%TZ)" >> "$JD/PROGRESS"
run_one(){ f="$JD/$(echo "$1"|tr , _).out"; [ -s "$f" ] && grep -q "^n=" "$f" && return 0; DYNT_KMAX="$KMAX" ./cm_dynt 64 "$B" --prefix "$1" > "$f" 2>&1; }
export -f run_one; export B JD KMAX
xargs -P "$NP" -I{} bash -c "run_one \"\$@\"" _ {} < "$JD/SUBJOBS9"
echo "frontier9 done=$(date -u +%FT%TZ)" >> "$JD/PROGRESS"
