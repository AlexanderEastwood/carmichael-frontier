#!/bin/bash
# KMAX=0 duplicates of the remaining chain slices (no kills; the KMAX=20 frontier workers keep running).
# Each duplicate writes to ~/kmax0_dup/<prefix>.out and, on a clean n=NONE verdict, installs that file as
# jobs_1e147/<prefix>.out if no verdict exists there yet -- the coverage verifier then sees it.
set -u; cd ~/carmichael-frontier/experiments/dynt || exit 1
B=1$(printf "0%.0s" $(seq 1 147)); JD=$PWD/jobs_1e147; OUT=~/kmax0_dup; mkdir -p $OUT
# 1) the finished KMAX=0 test verdict
P0=19,23,31,37,41,43,53,59,61,67,71,73,79,89,97,101,103,109; F0=$JD/$(echo $P0 | tr , _).out
if grep -q "^n=NONE" ~/kmax_test2/out_k0.txt && ! grep -q "^n=" "$F0" 2>/dev/null; then cp ~/kmax_test2/out_k0.txt "$F0"; echo "installed kmax0 verdict for $P0" >> $OUT/log; fi
# 2) duplicates
for P in 19,29,31,37,41,43,47,53,61,67,71,73,79,89,97,101,103,109 19,31,37,41,43,47,53,59,61,67,71,73,79,89,97,101,103,109 19,31,41,43,47,53,59,61,67,71,73,79,89,97,101,103,109,113 19,29,31,37,41,43,47,53,61,67,71,73,79,89 19,31,37,41,43,47,53,59,61,67,71,73,79,89 19,31,41,43,47,53,59,61,67,71,73,79,89,97; do
  F=$OUT/$(echo $P | tr , _).out; T=$JD/$(echo $P | tr , _).out
  nohup bash -c "s=\$(date +%s); DYNT_KMAX=0 ./cm_dynt 64 $B --prefix $P > $F 2>&1; echo \"KMAX=0 $P secs=\$(( \$(date +%s) - s )) \$(grep -m1 '^n=' $F)\" >> $OUT/log; if grep -q '^n=NONE' $F && ! grep -q '^n=' $T 2>/dev/null; then cp $F $T; echo \"installed $P\" >> $OUT/log; fi" > /dev/null 2>&1 < /dev/null &
done
echo "launched 6 KMAX=0 duplicates at $(date -u +%FT%TZ)" >> $OUT/log
