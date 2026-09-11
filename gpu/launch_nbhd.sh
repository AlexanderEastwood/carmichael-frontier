#!/bin/bash
# Experiment A: lambda-neighbourhood discovery pass (ksmall base, radii 1..min(rmax,6), richest moduli first),
# two runners sharing the GPU. Bound = current incumbent (Webster's N). Records: gpu/records/nbhd_s{0,1}_U5243.../
cd ~/carmichael-frontier || exit 1
PY=~/erdos/.venv/bin/python
DL=$(date -u -d '+6 hours' +%FT%TZ)
for i in 0 1; do
  FAM_TAG=nbhd_s$i FAMILY_JSON=gpu/nbhd_p300.json FAM_ORDER=file FAM_SHARD=$i/2 FAM_R_MAX=6 IDUMP_T=6 DEADLINE=$DL \
    nohup $PY -u gpu/run_mstar_family.py > gpu/nbhd_s$i.stdout 2>&1 < /dev/null &
done
sleep 2
echo "launched 2 shards at $(date -u +%FT%TZ), deadline $DL, pids $(pgrep -f '[r]un_mstar_family.py' | tr '\n' ' ')" >> gpu/nbhd_launch.log
