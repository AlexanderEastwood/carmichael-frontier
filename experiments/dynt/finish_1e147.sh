#!/usr/bin/env bash
set -u
cd "$(dirname "$0")"
B=1$(printf "0%.0s" $(seq 1 147)); JD=jobs_1e147; KMAX=20
# 1) stop the original runner + the laggard workers so no false tally is written
pkill -9 -f run_excl.sh 2>/dev/null; pkill -9 -f "cm_dynt 64" 2>/dev/null; sleep 2
# 2) laggard prefixes (single classes get +1 level -> depth 2; triples get +1 -> depth 4)
LAG="11:1 17:1 31:1 19,23,29:1 19,23,31:1 19,29,31:1 19,31,37:1 19,31,41:1"
for p in 11 17 31 19_23_29 19_23_31 19_29_31 19_31_37 19_31_41; do rm -f "$JD/$p.out"; done
python3 sub_gen.py $LAG > "$JD/SUBJOBS"
echo "subjobs=$(wc -l < "$JD/SUBJOBS") start=$(date -u +%FT%TZ)" >> "$JD/PROGRESS"
run_one(){ f="$JD/$(echo "$1"|tr , _).out"; DYNT_KMAX="$KMAX" ./cm_dynt 64 "$B" --prefix "$1" > "$f" 2>&1; }
export -f run_one; export B JD KMAX
xargs -P 30 -I{} bash -c 'run_one "$@"' _ {} < "$JD/SUBJOBS"
# 3) STRICT tally: every .out must contain exactly a NONE verdict; flag anything else
missing=0; hits=0; total=0
for f in "$JD"/*.out; do
  total=$((total+1))
  if grep -q "^n=NONE" "$f"; then :; elif grep -q "^n=[0-9]" "$f"; then hits=$((hits+1)); echo "HIT $f" >> "$JD/DONE.tmp"; else missing=$((missing+1)); echo "NOVERDICT $f" >> "$JD/DONE.tmp"; fi
done
{ echo "total_out=$total hits=$hits noverdict=$missing"; [ -f "$JD/DONE.tmp" ] && cat "$JD/DONE.tmp"
  if [ "$hits" = 0 ] && [ "$missing" = 0 ]; then echo "RESULT=ALL_NONE"; else echo "RESULT=INCOMPLETE_OR_HIT"; fi; } > "$JD/DONE"
rm -f "$JD/DONE.tmp"; echo "done=$(date -u +%FT%TZ) $(tail -1 "$JD/DONE")" >> "$JD/PROGRESS"
