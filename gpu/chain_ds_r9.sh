#!/bin/bash
# D_s family (s<=3, 22 moduli besides M*) at radii 9-10 with the filtered deletion side. Detached.
cd ~/carmichael-frontier || exit 1
LOG=gpu/chain_ds_r9.log; PY=~/erdos/.venv/bin/python
echo "[$(date -u +%FT%TZ)] chain_ds_r9 start" >> $LOG
for MOD in 589416059232000 252606882528000 21304194912000 1014485472000 2367132768000 3043456416000 4260838982400 5326048728000 7101398304000 10652097456000 28067431392000 50521376505600 63151720632000 84202294176000 117883211846400 126303441264000 147354014808000 196472019744000 294708029616000 353649635539200 442062044424000 884124088848000; do
  rm -f gpu/checkpoint_mstar_r9.json
  MODULUS=$MOD DFILTER=1 CKPT_TAG=r9 R_MIN=9 R_MAX=10 BASES=S64,ksmall INS_PRIME_CAP=10000 INS_MAX=256 SPLIT=4000 IDUMP_T=8 IDUMP_TIMEOUT=3600 WALL_CAP=3000 \
    $PY -u gpu/mitm_gpu_mstar_driver.py > gpu/run_mstar_r9_M$MOD.log 2>&1 < /dev/null
  echo "[$(date -u +%FT%TZ)] M=$MOD | $(tail -n 1 gpu/run_mstar_r9_M$MOD.log | cut -c1-100)" >> $LOG
  grep -q "NEW BEST" gpu/run_mstar_r9_M$MOD.log && echo "[$(date -u +%FT%TZ)] *** NEW BEST at M=$MOD ***" >> $LOG
done
echo "[$(date -u +%FT%TZ)] chain_ds_r9 DONE" >> $LOG
