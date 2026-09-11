#!/bin/bash
# Bare KMAX=0 duplicates of the 6 remaining chain slices: plain nohup cm_dynt, output only to ~/kmax0_dup.
# No wrappers, no copying, no kills. verify_1e147.py reads ~/kmax0_dup/<prefix>.out as a secondary verdict.
set -u; cd ~/carmichael-frontier/experiments/dynt || exit 1
B=1$(printf "0%.0s" $(seq 1 147)); OUT=~/kmax0_dup
for P in 19,29,31,37,41,43,47,53,61,67,71,73,79,89,97,101,103,109 19,31,37,41,43,47,53,59,61,67,71,73,79,89,97,101,103,109 19,31,41,43,47,53,59,61,67,71,73,79,89,97,101,103,109,113 19,29,31,37,41,43,47,53,61,67,71,73,79,89 19,31,37,41,43,47,53,59,61,67,71,73,79,89 19,31,41,43,47,53,59,61,67,71,73,79,89,97; do
  F=$OUT/$(echo $P | tr , _).out
  DYNT_KMAX=0 nohup ./cm_dynt 64 $B --prefix $P > $F 2>&1 < /dev/null &
done
sleep 2
echo "relaunched 6 bare KMAX=0 duplicates at $(date -u +%FT%TZ); cm_dynt=$(pgrep -fc '[c]m_dynt 64')" >> $OUT/log
