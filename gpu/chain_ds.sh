#!/bin/bash
# Astra's small-change modulus family D_s = { D | M* : #{p | N' : (p-1) ∤ D} <= 3 } — the 19 moduli
# not yet run (M*, M*/3, M*/7, M*/83 done). Expanded-pool driver, radii 1-8, both bases, sequential
# on the GPU after the r=8 driver exits. Detached; never touches ~/erdos.
cd ~/carmichael-frontier || exit 1
LOG=gpu/chain_ds.log; PY=~/erdos/.venv/bin/python
echo "[$(date -u +%FT%TZ)] chain_ds armed: waiting for r8 driver to exit" >> $LOG
until [ "$(pgrep -f '[m]itm_gpu_r8_driver' | wc -l)" = "0" ]; do sleep 30; done
echo "[$(date -u +%FT%TZ)] r8 exited: $(tail -n 1 gpu/run_r8.log | cut -c1-100)" >> $LOG
for MOD in 1014485472000 2367132768000 3043456416000 4260838982400 5326048728000 7101398304000 10652097456000 28067431392000 50521376505600 63151720632000 84202294176000 117883211846400 126303441264000 147354014808000 196472019744000 294708029616000 353649635539200 442062044424000 884124088848000; do
  rm -f gpu/checkpoint_mstar.json
  MODULUS=$MOD R_MIN=1 R_MAX=8 INS_PRIME_CAP=10000 INS_MAX=256 SPLIT=4000 WALL_CAP=3600 \
    $PY -u gpu/mitm_gpu_mstar_driver.py > gpu/run_mstar_M$MOD.log 2>&1 < /dev/null
  cp gpu/checkpoint_mstar.json gpu/checkpoint_mstar_M$MOD.json 2>/dev/null
  echo "[$(date -u +%FT%TZ)] M=$MOD | $(tail -n 1 gpu/run_mstar_M$MOD.log | cut -c1-100)" >> $LOG
  grep -q "NEW BEST" gpu/run_mstar_M$MOD.log && echo "[$(date -u +%FT%TZ)] *** NEW BEST at M=$MOD ***" >> $LOG
done
echo "[$(date -u +%FT%TZ)] chain_ds DONE" >> $LOG
