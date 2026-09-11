#!/usr/bin/env bash
set -u; cd "$(dirname "$0")"; B=1$(printf "0%.0s" $(seq 1 147)); JD=jobs_1e147; KMAX=20
pkill -9 -f finish_1e147.sh 2>/dev/null; sleep 1
LAG="19,23,31,37:1 19,29,31,37:1 19,29,31,41:1 19,31,37,41:1 19,31,41,43:1"
for p in 19,23,31,37 19,29,31,37 19,29,31,41 19,31,37,41 19,31,41,43; do pkill -9 -f "prefix $p\$" 2>/dev/null; rm -f "$JD/$(echo $p|tr , _).out"; done
rm -f "$JD/DONE"; sleep 2
python3 sub_gen.py $LAG > "$JD/SUBJOBS2"
echo "subjobs2=$(wc -l < "$JD/SUBJOBS2") start=$(date -u +%FT%TZ)" >> "$JD/PROGRESS"
run_one(){ f="$JD/$(echo "$1"|tr , _).out"; DYNT_KMAX="$KMAX" ./cm_dynt 64 "$B" --prefix "$1" > "$f" 2>&1; }
export -f run_one; export B JD KMAX
xargs -P 28 -I{} bash -c "run_one \"\$@\"" _ {} < "$JD/SUBJOBS2"
missing=0; hits=0; total=0
for f in "$JD"/*.out; do total=$((total+1)); if grep -q "^n=NONE" "$f"; then :; elif grep -q "^n=[0-9]" "$f"; then hits=$((hits+1)); echo "HIT $f" >> "$JD/DONE.tmp"; else missing=$((missing+1)); echo "NOVERDICT $f" >> "$JD/DONE.tmp"; fi; done
{ echo "total_out=$total hits=$hits noverdict=$missing"; [ -f "$JD/DONE.tmp" ] && cat "$JD/DONE.tmp"; if [ "$hits" = 0 ] && [ "$missing" = 0 ]; then echo "RESULT=ALL_NONE"; else echo "RESULT=INCOMPLETE_OR_HIT"; fi; } > "$JD/DONE"
rm -f "$JD/DONE.tmp"; echo "done=$(date -u +%FT%TZ) $(tail -1 "$JD/DONE")" >> "$JD/PROGRESS"
