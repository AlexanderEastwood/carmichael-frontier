#!/usr/bin/env bash
set -u; cd ~/carmichael-frontier/experiments/dynt; B=1$(printf "0%.0s" $(seq 1 147)); JD=jobs_1e147; KMAX=20
for p in $(pgrep -f "[f]inish2_1e147"); do kill -9 "$p"; done; sleep 1        # old runner first: no false tally
LAG="19,23,31,37,41,43:2 19,29,31,37,41,43:2 19,29,31,41,43,47:2 19,31,37,41,43,47:2 19,31,41,43,47,53:2"
for p in 19,23,31,37,41,43 19,29,31,37,41,43 19,29,31,41,43,47 19,31,37,41,43,47 19,31,41,43,47,53; do
  for q in $(pgrep -f "[p]refix $p\$"); do kill -9 "$q"; done; rm -f "$JD/$(echo $p|tr , _).out"; done
rm -f "$JD/DONE"; sleep 2
python3 sub_gen.py $LAG > "$JD/SUBJOBS4"
echo "subjobs4=$(wc -l < "$JD/SUBJOBS4") start=$(date -u +%FT%TZ)" >> "$JD/PROGRESS"
run_one(){ f="$JD/$(echo "$1"|tr , _).out"; DYNT_KMAX="$KMAX" ./cm_dynt 64 "$B" --prefix "$1" > "$f" 2>&1; }
export -f run_one; export B JD KMAX
xargs -P 28 -I{} bash -c "run_one \"\$@\"" _ {} < "$JD/SUBJOBS4"
missing=0; hits=0; total=0
for f in "$JD"/*.out; do total=$((total+1)); if grep -q "^n=NONE" "$f"; then :; elif grep -q "^n=[0-9]" "$f"; then hits=$((hits+1)); echo "HIT $f" >> "$JD/DONE.tmp"; else missing=$((missing+1)); echo "NOVERDICT $f" >> "$JD/DONE.tmp"; fi; done
{ echo "total_out=$total hits=$hits noverdict=$missing"; [ -f "$JD/DONE.tmp" ] && cat "$JD/DONE.tmp"; if [ "$hits" = 0 ] && [ "$missing" = 0 ]; then echo "RESULT=ALL_NONE"; else echo "RESULT=INCOMPLETE_OR_HIT"; fi; } > "$JD/DONE"
rm -f "$JD/DONE.tmp"; echo "done=$(date -u +%FT%TZ) $(tail -1 "$JD/DONE")" >> "$JD/PROGRESS"
