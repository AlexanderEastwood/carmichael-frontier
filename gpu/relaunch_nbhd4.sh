#!/bin/bash
# Re-split the neighbourhood discovery pass into 4 shards (GPU now free; pass is CPU-generation-bound).
cd ~/carmichael-frontier || exit 1
PY=~/erdos/.venv/bin/python
for p in $(pgrep -f "[r]un_mstar_family.py"); do kill "$p"; done
sleep 1
for p in $(pgrep -f "[m]itm_gpu_mstar_driver.py"); do kill "$p"; done
sleep 3
R=gpu/records; DD=$R/nbhd_s0_U524394558743,$R/nbhd_s1_U524394558743,$R/nbhd_s2_U524394558743,$R/nbhd_s3_U524394558743
DL=$(date -u -d '+6 hours' +%FT%TZ)
for i in 0 1 2 3; do
  FAM_TAG=nbhd_s$i FAMILY_JSON=gpu/nbhd_p300.json FAM_ORDER=file FAM_SHARD=$i/4 FAM_R_MAX=6 IDUMP_T=8 DEADLINE=$DL \
  FAM_SKIP_DONE=1 FAM_DONE_DIRS=$DD nohup $PY -u gpu/run_mstar_family.py >> gpu/nbhd_s$i.stdout 2>&1 < /dev/null &
done
sleep 2
echo "relaunched 4 shards at $(date -u +%FT%TZ), deadline $DL, pids $(pgrep -f '[r]un_mstar_family.py' | tr '\n' ' ')" >> gpu/nbhd_launch.log
