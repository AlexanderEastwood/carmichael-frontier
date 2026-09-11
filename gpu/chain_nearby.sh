#!/bin/bash
# Nearby-modulus expanded-pool searches (Astra 11:40Z), then r=8 block 361-400.
# Sequential on the single GPU; detached; never touches ~/erdos.
cd ~/carmichael-frontier || exit 1
LOG=gpu/chain_nearby.log
PY=~/erdos/.venv/bin/python
echo "[$(date -u +%FT%TZ)] chain_nearby start" >> $LOG
# keep the completed M* checkpoint (the driver keys its checkpoint file by modulus and would overwrite it)
cp -n gpu/checkpoint_mstar.json gpu/checkpoint_mstar_Mstar_done.json 2>/dev/null
for MOD in 589416059232000 252606882528000 21304194912000; do
  [ "$(pgrep -f '[p]ython -u gpu/mitm_gpu_mstar' | wc -l)" = "0" ] || { echo "driver still running?!" >> $LOG; exit 2; }
  rm -f gpu/checkpoint_mstar.json
  echo "[$(date -u +%FT%TZ)] MODULUS=$MOD radii 1..8 start" >> $LOG
  MODULUS=$MOD R_MIN=1 R_MAX=8 INS_PRIME_CAP=10000 INS_MAX=256 SPLIT=4000 WALL_CAP=7200 \
    $PY -u gpu/mitm_gpu_mstar_driver.py > gpu/run_mstar_M$MOD.log 2>&1 < /dev/null
  cp gpu/checkpoint_mstar.json gpu/checkpoint_mstar_M$MOD.json 2>/dev/null
  echo "[$(date -u +%FT%TZ)] MODULUS=$MOD exit=$? | $(tail -n 1 gpu/run_mstar_M$MOD.log | cut -c1-120)" >> $LOG
  grep -q "NEW BEST" gpu/run_mstar_M$MOD.log && echo "[$(date -u +%FT%TZ)] *** NEW BEST at M=$MOD ***" >> $LOG
done
# r=8 portfolio block 361-400 (fresh range: keep checkpoint incumbent, just set status running)
python3 - <<'PY'
import json; c=json.load(open("gpu/checkpoint_r8.json")); c["status"]="running"
json.dump(c,open("gpu/checkpoint_r8.json","w"),indent=1); print("checkpoint: last_completed_idx", c["last_completed_idx"])
PY
echo "[$(date -u +%FT%TZ)] launching r8 block 361-400" >> $LOG
STREAM=1 R=8 IDX_MIN=361 IDX_MAX=400 BASES=S64,ksmall NCHUNKS_D=16 IDUMP_T=8 IDUMP_TIMEOUT=1800 WALL_CAP=28800 \
  nohup $PY -u gpu/mitm_gpu_r8_driver.py >> gpu/run_r8.log 2>&1 < /dev/null &
echo "[$(date -u +%FT%TZ)] r8 driver pid $!" >> $LOG
