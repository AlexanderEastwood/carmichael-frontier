#!/bin/bash
# Wait for the r=8 block 241-360 sweep to finish, then launch the M* expanded-pool driver.
# Detached on titan; never touches ~/erdos.
cd ~/carmichael-frontier || exit 1
LOG=gpu/chain_mstar.log
echo "[$(date -u +%FT%TZ)] chain armed: waiting for r8 driver to exit" >> $LOG
until [ "$(pgrep -f '[m]itm_gpu_r8_driver' | wc -l)" = "0" ]; do sleep 60; done
echo "[$(date -u +%FT%TZ)] r8 driver exited; tail:" >> $LOG
tail -n 3 gpu/run_r8.log >> $LOG
if ! grep -q "sweep_complete" gpu/run_r8.log; then
  echo "[$(date -u +%FT%TZ)] WARNING: no sweep_complete in run_r8.log (driver died?) - launching M* anyway, r8 needs relaunch" >> $LOG
fi
sleep 20
R_MIN=5 R_MAX=8 INS_PRIME_CAP=10000 INS_MAX=256 SPLIT=4000 \
  nohup ~/erdos/.venv/bin/python -u gpu/mitm_gpu_mstar_driver.py > gpu/run_mstar.log 2>&1 < /dev/null &
echo "[$(date -u +%FT%TZ)] M* driver launched pid $!" >> $LOG
