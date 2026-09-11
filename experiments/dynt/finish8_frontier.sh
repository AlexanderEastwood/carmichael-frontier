#!/usr/bin/env bash
# Persistent feasibility-pruned frontier for the exclusion chain: generated ONCE, run to completion, no restarts.
set -u; cd ~/carmichael-frontier/experiments/dynt; B=1$(printf "0%.0s" $(seq 1 147)); JD=jobs_1e147; KMAX=${KMAX:-20}
for p in $(pgrep -f "[f]inish7_1e147"); do kill -9 "$p"; done; for p in $(pgrep -f "[x]args -P 28"); do kill -9 "$p"; done; sleep 1
P10="19,23,31,37,41,43,53,59,61,67 19,29,31,37,41,43,47,53,61,67 19,31,37,41,43,47,53,59,61,67 19,31,37,41,43,47,59,61,67,71 19,31,41,43,47,53,59,61,67,71"
for p in $P10; do for q in $(pgrep -f "[p]refix $p,"); do kill -9 "$q"; done; done   # depth-11 chain workers (children of the depth-10 prefixes)
sleep 2; A=""; for p in $P10; do A="$A $p:4"; done
python3 sub_gen_feasible.py $A > "$JD/SUBJOBS8"
echo "subjobs8=$(wc -l < "$JD/SUBJOBS8") kmax=$KMAX start=$(date -u +%FT%TZ)" >> "$JD/PROGRESS"
run_one(){ f="$JD/$(echo "$1"|tr , _).out"; DYNT_KMAX="$KMAX" ./cm_dynt 64 "$B" --prefix "$1" > "$f" 2>&1; }
export -f run_one; export B JD KMAX
xargs -P 28 -I{} bash -c "run_one \"\$@\"" _ {} < "$JD/SUBJOBS8"
echo "frontier8 done=$(date -u +%FT%TZ)" >> "$JD/PROGRESS"
