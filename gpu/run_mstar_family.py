#!/usr/bin/env python3
"""EXACT search of a divisor family (Astra, 2026-09-11): every Carmichael n with exactly 64 prime
factors, n < U, and lambda(n) | M has lambda(n) = D for one of the survivor divisors D of M listed in
the census file (gpu/divisor_family_census.py). For each D the affordable pool is
{q eligible : q <= cap(D)} with cap = (U-1)/T63(D), and any 64-subset with product < U differs from
the 64 smallest (ksmall) in at most rmax(D) places (size-forced). So ksmall + radii 1..rmax over the
FULL capped pool is exhaustive at D; radius 0 (ksmall itself) is checked directly with the frozen
oracle.

Records: one file per modulus under gpu/records/<FAM_TAG>_U<first 12 digits of U>/D<D>.log, holding
the driver's output for radii 1..rmax plus a '[runner] r=0 ...' line for the base check. The same
bound U (U_OVERRIDE if set) is used for the base checks and for the exchange instances. The pass
ends with 'family pass finished' only if every requested modulus completed every requested radius;
otherwise it prints 'FAMILY PASS INCOMPLETE' and exits 1. Verify a finished pass with
    python3 gpu/verify_family_coverage.py <census.json> <records dir> [<bound>]

env: FAMILY_JSON (census file; default gpu/mstar_divisor_family.json), FAM_TAG (default: census
file stem), U_OVERRIDE (frozen bound), S_MAX_SKIP / S_ONLY_LE (row filters, default: all rows),
FAM_R_MIN / FAM_R_MAX (radius window, default 1..rmax), FAM_WALL_CAP, FAM_IDUMP_TIMEOUT, FAM_INS_MAX,
DEADLINE (UTC ISO; default none)."""
import os, sys, json, time, subprocess, datetime as dt
HERE = os.path.dirname(os.path.abspath(__file__)); FR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(FR, "ref"))
import ref_carmichael as ref
PY = sys.executable; DRIVER = os.path.join(HERE, "mitm_gpu_mstar_driver.py")
DEADLINE = dt.datetime.fromisoformat(os.environ["DEADLINE"]) if os.environ.get("DEADLINE") else None
S_MAX_SKIP = int(os.environ.get("S_MAX_SKIP", "-1"))
FAMILY_JSON = os.environ.get("FAMILY_JSON", os.path.join(HERE, "mstar_divisor_family.json"))
FAM_TAG = os.environ.get("FAM_TAG", os.path.splitext(os.path.basename(FAMILY_JSON))[0])
LOG = os.path.join(HERE, f"run_family_{FAM_TAG}.log")

def log(*a):
    with open(LOG, "a") as f: print(dt.datetime.now(dt.timezone.utc).strftime("%H:%M:%SZ"), *a, file=f, flush=True)

def pool(D, cap):
    out = []; q = 3
    while q <= cap:
        if D % (q - 1) == 0 and D % q != 0 and ref.is_prime(q): out.append(q)
        q += 2
    return out

def main() -> None:
    fam = json.load(open(FAMILY_JSON))
    if isinstance(fam, dict): fam = fam["survivors"]          # divisor_family_census.py format
    inc = json.load(open(os.path.join(FR, "results_k64_best_global.json")))
    U = int(os.environ["U_OVERRIDE"]) if os.environ.get("U_OVERRIDE") else int(inc["n"])
    rec_dir = os.path.join(HERE, "records", f"{FAM_TAG}_U{str(U)[:12]}"); os.makedirs(rec_dir, exist_ok=True)
    todo = [r for r in fam if r["s"] > S_MAX_SKIP]
    if os.environ.get("S_ONLY_LE"): todo = [r for r in fam if r["s"] <= int(os.environ["S_ONLY_LE"])]
    if os.environ.get("FAM_ORDER", "rmax") == "file": order = lambda rows: rows            # keep the census file's order (e.g. richest first)
    else: order = lambda rows: sorted(rows, key=lambda r: r["rmax"])
    todo = order(todo)
    if os.environ.get("FAM_SHARD"):                      # FAM_SHARD=i/n : this runner takes rows with index % n == i
        i, n = (int(x) for x in os.environ["FAM_SHARD"].split("/")); todo = [r for k, r in enumerate(todo) if k % n == i]
    log(f"family {FAM_TAG}: {len(fam)} survivors; running {len(todo)}; U={len(str(U))} digits (U_OVERRIDE={'yes' if os.environ.get('U_OVERRIDE') else 'no'}); records {rec_dir}")
    n_r0 = 0; n_hit = 0; failures = []; skipped = 0
    done_dirs = [d for d in os.environ.get("FAM_DONE_DIRS", "").split(",") if d]   # FAM_SKIP_DONE=1: skip moduli already complete there
    def already_done(D):
        for d in done_dirs:
            f = os.path.join(d, f"D{D}.log")
            if os.path.exists(f):
                t = open(f, errors="replace").read()
                if "DONE complete" in t and "INCOMPLETE" not in t: return True
                if "[runner] no exchange radii requested" in t: return True
        return False
    for row in todo:
        D, cap, rmax = int(row["D"]), int(row["cap"]), int(row["rmax"])
        if os.environ.get("FAM_SKIP_DONE") == "1" and already_done(D): skipped += 1; continue
        if DEADLINE and dt.datetime.now(dt.timezone.utc) > DEADLINE: log("DEADLINE reached; stopping"); failures.append(("deadline", D)); break
        P = pool(D, cap); assert len(P) == row["pool_capped"], (D, len(P), row["pool_capped"])
        base = P[:64]; ok, n, _ = ref.verify_certificate(base); n_r0 += 1
        lf = os.path.join(rec_dir, f"D{D}.log")
        with open(lf, "w") as out:
            out.write(f"[runner] family={FAM_TAG} D={D} U={U} pool={len(P)} cap={cap} rmax={rmax}\n")
            out.write(f"[runner] r=0 base=ksmall({len(base)} primes, largest {base[-1]}) carmichael={ok} below_U={ok and n < U}\n")
        if ok and n < U: log(f"*** r=0 HIT: ksmall at D={D} is Carmichael below U ***"); n_hit += 1
        rlo = int(os.environ.get("FAM_R_MIN", "1")); rhi = min(rmax, int(os.environ.get("FAM_R_MAX", "99")))
        if rhi < rlo:
            with open(lf, "a") as out: out.write(f"[runner] no exchange radii requested (rmax={rmax})\n")
            log(f"D={D} s={row['s']} pool={len(P)} rmax={rmax}: base checked, no exchange radii"); continue
        env = dict(os.environ, MODULUS=str(D), DFILTER="1", CKPT_TAG=f"fam_{FAM_TAG}", R_MIN=str(rlo), R_MAX=str(rhi), BASES="ksmall",
                   INS_PRIME_CAP=str(cap + 1), INS_MAX=os.environ.get("FAM_INS_MAX", "1024"), SPLIT="4000", IDUMP_T="8",
                   IDUMP_TIMEOUT=os.environ.get("FAM_IDUMP_TIMEOUT", "1800"), WALL_CAP=os.environ.get("FAM_WALL_CAP", "3000"), U_OVERRIDE=str(U))
        ck = os.path.join(HERE, f"checkpoint_mstar_fam_{FAM_TAG}.json")
        if os.path.exists(ck): os.remove(ck)
        t = time.time()
        with open(lf, "a") as out: rc = subprocess.call([PY, "-u", DRIVER], stdout=out, stderr=subprocess.STDOUT, env=env, stdin=subprocess.DEVNULL)
        txt = open(lf).read(); inc_ = txt.count("INCOMPLETE"); done = "DONE complete" in txt
        radii_seen = {int(l.split("r=")[1].split()[0]) for l in txt.splitlines() if " r=" in l and ("dfilter" in l or "host-exact" in l or "nothing to search" in l)}
        missing = [r for r in range(rlo, rhi + 1) if r not in radii_seen]
        ok_run = rc == 0 and inc_ == 0 and done and not missing
        if not ok_run: failures.append((D, f"rc={rc} incomplete={inc_} done={done} missing_radii={missing}"))
        log(f"D={D} s={row['s']} pool={len(P)} |INS|={len(P)-64} radii {rlo}..{rhi} rc={rc} {time.time()-t:.0f}s incomplete={inc_} missing_radii={missing} complete={ok_run}")
        if "NEW BEST" in txt: n_hit += 1; log(f"*** NEW BEST at D={D} -> gpu/results_k64_gpu_mstar_fam_{FAM_TAG}.json ***")
    if failures:
        log(f"FAMILY PASS INCOMPLETE: {len(failures)} failures: {failures[:10]}"); sys.exit(1)
    log(f"family pass finished: moduli={len(todo)} skipped_done={skipped} r0_checked={n_r0} hits={n_hit} records={rec_dir}")

if __name__ == "__main__": main()
