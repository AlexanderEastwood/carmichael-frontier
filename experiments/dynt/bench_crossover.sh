#!/usr/bin/env bash
# Crossover experiment (GPT-6's recommendation, 2026-09-11): at the ENDGAME bound N (Webster's
# incumbent) with the interval (Y, N] -- everything <= Y = 10^147 already excluded -- compare
#   A  baseline descent, no interval           (cm_dynt_base, production switch)
#   B  interval only                           (cm_dynt_iv, DYNT_LO=Y, production switch)
#   C  interval + early AP completion KMAX=1e2 (cm_dynt_iv, DYNT_LO=Y, DYNT_KMAX=100)
#   D  interval + early AP completion KMAX=1e4 (cm_dynt_iv, DYNT_LO=Y, DYNT_KMAX=10000)
#   E  interval + early AP completion KMAX=1e6 (cm_dynt_iv, DYNT_LO=Y, DYNT_KMAX=1000000)
# on the same sample of hard all-smallest-primes chain prefixes. Every arm runs the SAME prefixes
# sequentially (arm by arm) with the same parallelism, so timings are comparable.
# Correctness gate: every arm must report the same completions in (Y, N].
set -u
cd ~/carmichael-frontier/experiments/dynt || exit 1
N=${N:-5243945587436457869157283735696333365758212589223889025728267613522676393978172071040517944862906263750800337318732066878122067505358422938650976001}
Y=${Y:-1$(printf "0%.0s" $(seq 1 147))}
SAMPLE=${SAMPLE:-jobs_1e147/SUBJOBS11}
NJOBS=${NJOBS:-24}
NP=${NP:-8}
OUT=${OUT:-$HOME/bench_crossover}
mkdir -p "$OUT"
head -n "$NJOBS" "$SAMPLE" > "$OUT/prefixes.txt"
echo "bench start $(date -u +%FT%TZ) jobs=$NJOBS np=$NP bound=N($(echo -n $N | wc -c) digits) lo=Y(10^147)" >> "$OUT/log"

run_arm() {   # $1 = arm name, $2 = binary, $3.. = env assignments
  local arm="$1" bin="$2"; shift 2
  local d="$OUT/$arm"; mkdir -p "$d"
  run_one() {
    local pfx="$1" f="$D_ARM/$(echo "$pfx" | tr , _).out"
    [ -s "$f" ] && return 0
    /usr/bin/time -f "TIME %e" env $ENV_ARM DYNT_DUMP=1 DYNT_FREEZE=1 ./"$BIN_ARM" 64 "$N_ARM" --prefix "$pfx" > "$f" 2>&1
  }
  export -f run_one; export D_ARM="$d" BIN_ARM="$bin" N_ARM="$N" ENV_ARM="$*"
  local t0=$(date +%s.%N)
  xargs -P "$NP" -I{} bash -c 'run_one "$@"' _ {} < "$OUT/prefixes.txt"
  local t1=$(date +%s.%N)
  echo "$arm wall_total=$(echo "$t1 - $t0" | bc) env=[$*]" >> "$OUT/log"
}

run_arm A cm_dynt_base X=0
run_arm B cm_dynt_iv   DYNT_LO="$Y"
run_arm C cm_dynt_iv   DYNT_LO="$Y" DYNT_KMAX=100
run_arm D cm_dynt_iv   DYNT_LO="$Y" DYNT_KMAX=10000
run_arm E cm_dynt_iv   DYNT_LO="$Y" DYNT_KMAX=1000000
echo "bench done $(date -u +%FT%TZ)" >> "$OUT/log"
