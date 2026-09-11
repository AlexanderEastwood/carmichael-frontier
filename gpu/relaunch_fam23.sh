#!/bin/bash
# Kill any running family runner / M* driver (by PID) and relaunch the 23-modulus exact run with the fixed driver.
cd ~/carmichael-frontier || exit 1
for p in $(pgrep -f "[r]un_mstar_family.py"); do kill "$p"; done
sleep 1
for p in $(pgrep -f "[m]itm_gpu_mstar_driver.py"); do kill "$p"; done
sleep 2
[ -f gpu/run_mstar_family.log ] && mv gpu/run_mstar_family.log "gpu/run_mstar_family23_attempt_$(date -u +%H%M%S).log"
S_ONLY_LE=3 FAM_WALL_CAP=3000 FAM_IDUMP_TIMEOUT=2400 DEADLINE=2026-09-11T16:08:00+00:00 \
  nohup ~/erdos/.venv/bin/python -u gpu/run_mstar_family.py > gpu/run_mstar_family23.stdout 2>&1 < /dev/null &
echo "relaunched pid $! at $(date -u +%FT%TZ)" >> gpu/relaunch_fam23.log
