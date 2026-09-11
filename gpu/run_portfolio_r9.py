#!/usr/bin/env python3
"""Filtered r=9-10 exchange-MITM over the Stage-B portfolio moduli (the r=8 sweep's list), one
mitm_gpu_mstar_driver run per modulus (DFILTER=1, both bases, pool<10000). Skips moduli already
covered at r=9-10 (M* and the 22 D_s moduli). Stops at DEADLINE (UTC ISO) so it never runs past
the authorised window. Log: gpu/run_portfolio_r9.log"""
import os, sys, time, subprocess, importlib.util, datetime as dt
HERE = os.path.dirname(os.path.abspath(__file__)); FR = os.path.dirname(HERE)
spec = importlib.util.spec_from_file_location("stageB_driver", os.path.join(FR, "stageB_driver.py"))
assert spec is not None and spec.loader is not None
SB = importlib.util.module_from_spec(spec); spec.loader.exec_module(SB)
DEADLINE = dt.datetime.fromisoformat(os.environ.get("DEADLINE", "2026-09-11T16:00:00+00:00"))
IDX_MIN = int(os.environ.get("IDX_MIN", "1")); IDX_MAX = int(os.environ.get("IDX_MAX", "400"))
DONE_MODULI = {1768248177696000, 589416059232000, 252606882528000, 21304194912000, 1014485472000, 2367132768000,
    3043456416000, 4260838982400, 5326048728000, 7101398304000, 10652097456000, 28067431392000, 50521376505600,
    63151720632000, 84202294176000, 117883211846400, 126303441264000, 147354014808000, 196472019744000,
    294708029616000, 353649635539200, 442062044424000, 884124088848000}
PY = sys.executable; DRIVER = os.path.join(HERE, "mitm_gpu_mstar_driver.py")
def log(*a):
    with open(os.path.join(HERE, "run_portfolio_r9.log"), "a") as f: print(dt.datetime.now(dt.timezone.utc).strftime("%H:%M:%SZ"), *a, file=f, flush=True)
def main() -> None:
    portfolio = SB.build_portfolio(SB.load_facts()); seen = set(DONE_MODULI); log(f"portfolio {len(portfolio)} moduli; idx {IDX_MIN}..{IDX_MAX}; deadline {DEADLINE.isoformat()}")
    for idx, (M, _dmod, prov) in enumerate(portfolio):
        if idx < IDX_MIN or idx > IDX_MAX: continue
        M = int(M)
        if M in seen: log(f"idx {idx} M={M} already covered; skip"); continue
        seen.add(M)
        if dt.datetime.now(dt.timezone.utc) > DEADLINE: log("DEADLINE reached; stopping"); break
        ck = os.path.join(HERE, "checkpoint_mstar_pf.json")
        if os.path.exists(ck): os.remove(ck)
        env = dict(os.environ, MODULUS=str(M), DFILTER="1", CKPT_TAG="pf", R_MIN="9", R_MAX="10", BASES="S64,ksmall",
                   INS_PRIME_CAP="10000", INS_MAX="256", SPLIT="4000", IDUMP_T="8", IDUMP_TIMEOUT="1800", WALL_CAP="1500")
        t = time.time(); lf = os.path.join(HERE, f"run_pf_r9_idx{idx}.log")
        with open(lf, "w") as out: rc = subprocess.call([PY, "-u", DRIVER], stdout=out, stderr=subprocess.STDOUT, env=env, stdin=subprocess.DEVNULL)
        txt = open(lf).read(); last = txt.strip().splitlines()[-1][:110] if txt.strip() else ""
        log(f"idx {idx} M={M} prov={prov} rc={rc} {time.time()-t:.0f}s | {last}")
        if "NEW BEST" in txt: log(f"*** NEW BEST at idx {idx} M={M} -> gpu/results_k64_gpu_mstar_pf.json ***")
    log("portfolio r9-10 pass finished")
if __name__ == "__main__": main()
