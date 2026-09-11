#!/usr/bin/env python3
"""EXACT search of the M*-divisor family (Astra, 2026-09-11): every Carmichael n with exactly 64
prime factors, n < N', and lambda(n) | M* has lambda(n) = D for one of the 111 divisors D of M*
whose 64 smallest eligible primes multiply to less than N' (gpu/mstar_divisor_family.json, recomputed
independently: 111 survivors). For each D the affordable pool is {q eligible : q <= cap(D)} with
cap = (N'-1)/T63, and any 64-subset with product < N' differs from the 64 smallest (ksmall) in at
most rmax(D) places (size-forced). So ksmall + radii 1..rmax over the FULL capped pool is exhaustive
at D. This runner does the 88 moduli outside the s<=3 neighbourhood (rmax <= 10, |INS| <= 206);
radius 0 (ksmall itself) is checked directly with the frozen oracle. Log: gpu/run_mstar_family.log"""
import os, sys, json, time, subprocess, datetime as dt
HERE = os.path.dirname(os.path.abspath(__file__)); FR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(FR, "ref"))
import ref_carmichael as ref
PY = sys.executable; DRIVER = os.path.join(HERE, "mitm_gpu_mstar_driver.py")
DEADLINE = dt.datetime.fromisoformat(os.environ.get("DEADLINE", "2026-09-11T16:05:00+00:00"))
S_MAX_SKIP = int(os.environ.get("S_MAX_SKIP", "3"))   # skip the neighbourhood already run (s <= 3)
def log(*a):
    with open(os.path.join(HERE, "run_mstar_family.log"), "a") as f: print(dt.datetime.now(dt.timezone.utc).strftime("%H:%M:%SZ"), *a, file=f, flush=True)
def pool(D, cap):
    out = []; q = 3
    while q <= cap:
        if D % (q - 1) == 0 and D % q != 0 and ref.is_prime(q): out.append(q)
        q += 2
    return out
def main() -> None:
    fam = json.load(open(os.environ.get("FAMILY_JSON", os.path.join(HERE, "mstar_divisor_family.json"))))
    if isinstance(fam, dict): fam = fam["survivors"]          # divisor_family_census.py format
    inc = json.load(open(os.path.join(FR, "results_k64_best_global.json"))); U = int(inc["n"])
    todo = [r for r in fam if r["s"] > S_MAX_SKIP]
    if os.environ.get("S_ONLY_LE"):                      # e.g. S_ONLY_LE=3: the 23 neighbourhood moduli at ALL feasible radii, full capped pools
        todo = [r for r in fam if r["s"] <= int(os.environ["S_ONLY_LE"])]
    log(f"family {len(fam)} survivors; running {len(todo)}; U={len(str(U))} digits")
    n_r0 = 0; n_hit = 0
    for row in sorted(todo, key=lambda r: r["rmax"]):
        D, cap, rmax = int(row["D"]), int(row["cap"]), int(row["rmax"])
        if dt.datetime.now(dt.timezone.utc) > DEADLINE: log("DEADLINE reached; stopping"); break
        P = pool(D, cap); assert len(P) == row["pool_capped"], (D, len(P), row["pool_capped"])
        base = P[:64]; ok, n, _ = ref.verify_certificate(base); n_r0 += 1
        if ok and n < U: log(f"*** r=0 HIT: ksmall at D={D} is Carmichael below U ***"); n_hit += 1
        if rmax == 0: log(f"D={D} s={row['s']} pool={len(P)} rmax=0: base checked, nothing else affordable"); continue
        ck = os.path.join(HERE, "checkpoint_mstar_fam.json")
        if os.path.exists(ck): os.remove(ck)
        env = dict(os.environ, MODULUS=str(D), DFILTER="1", CKPT_TAG="fam", R_MIN="1", R_MAX=str(rmax), BASES="ksmall",
                   INS_PRIME_CAP=str(cap + 1), INS_MAX=os.environ.get("FAM_INS_MAX", "1024"), SPLIT="4000", IDUMP_T="8",
                   IDUMP_TIMEOUT=os.environ.get("FAM_IDUMP_TIMEOUT", "1800"), WALL_CAP=os.environ.get("FAM_WALL_CAP", "1500"))
        t = time.time(); lf = os.path.join(HERE, f"run_fam_D{D}.log")
        with open(lf, "w") as out: rc = subprocess.call([PY, "-u", DRIVER], stdout=out, stderr=subprocess.STDOUT, env=env, stdin=subprocess.DEVNULL)
        txt = open(lf).read(); inst = txt.count("dfilter"); skipped = txt.count("no stratum"); inc_ = txt.count("INCOMPLETE")
        log(f"D={D} s={row['s']} pool={len(P)} |INS|={len(P)-64} rmax={rmax} rc={rc} {time.time()-t:.0f}s instances={inst} bound-skipped={skipped} incomplete={inc_} complete={'DONE complete' in txt}")
        if "NEW BEST" in txt: n_hit += 1; log(f"*** NEW BEST at D={D} -> gpu/results_k64_gpu_mstar_fam.json ***")
    log(f"family pass finished: r0_checked={n_r0} hits={n_hit}")
if __name__ == "__main__": main()
