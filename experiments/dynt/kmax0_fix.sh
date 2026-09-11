#!/bin/bash
# Remove the auto-copy: kill only the bash wrappers of the KMAX=0 duplicates (their cm_dynt children
# keep running and write to ~/kmax0_dup/<prefix>.out). Never touches jobs_1e147.
OUT=~/kmax0_dup
for p in $(pgrep -f "[k]max0_dup/"); do kill "$p"; done
sleep 1
cp ~/kmax_test2/out_k0.txt $OUT/19_23_31_37_41_43_53_59_61_67_71_73_79_89_97_101_103_109.out
echo "wrappers killed; kmax_test2 KMAX=0 verdict placed in $OUT at $(date -u +%FT%TZ); cm_dynt=$(pgrep -fc '[c]m_dynt 64')" >> $OUT/log
