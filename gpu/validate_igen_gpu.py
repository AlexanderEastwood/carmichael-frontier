#!/usr/bin/env python3
"""Validate igen_gpu against the exact-integer C++ dumper (src/mitm64_idump_stream2, IDUMP_RANK=1).
For each test instance: the GPU record set must be a SUPERSET of the CPU set, and every extra record
must decode to a subset whose exact product is >= Icap (i.e. only the float-margin boundary cases).
Instances: the ksmall base at M0 with pool<4000 (r=3..8), and at M0 with the full 461-insertion pool
(r=3..8, per-radius cap as in the driver), all with U frozen at N0 or at the previous incumbent."""
import os, sys, json, subprocess, time
import numpy as np, cupy as cp
HERE = os.path.dirname(os.path.abspath(__file__)); FR = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(FR, "ref"))
import ref_carmichael as ref, igen_gpu
from math import comb
IDUMP = os.path.join(FR, "src", "mitm64_idump_stream2")

def pool(M, cap):
    out = []; q = 3
    while q < cap:
        if M % (q - 1) == 0 and M % q != 0 and ref.is_prime(q): out.append(q)
        q += 2
    return out

def unrank(v, r, n):
    v = int(v); idx = []; x = n - 1
    for j in range(r - 1, -1, -1):
        while comb(x, j + 1) > v: x -= 1
        idx.append(x); v -= comb(x, j + 1); x -= 1
    assert v == 0; return tuple(reversed(idx))

def cpu_records(M, r, P0, Icap, B0, INS, split, hmax):
    inp = ("\n".join([f"{M} {r} {P0 % M} {Icap} 1 8 0", "64 " + " ".join(map(str, B0)), f"{len(INS)} " + " ".join(map(str, INS))]) + "\n").encode()
    env = dict(os.environ, IDUMP_SPLIT=str(split), IDUMP_HMAX=str(hmax), IDUMP_RANK="1")
    t = time.time(); pr = subprocess.run([IDUMP], input=inp, capture_output=True, env=env, check=True); dt = time.time() - t
    a = np.frombuffer(pr.stdout, dtype=np.uint64).reshape(-1, 2)
    return a, dt

def run_case(M, U, cap_pool, r, split=4000):
    P = pool(M, cap_pool); B0 = P[:64]; P0 = 1
    for p in B0: P0 *= p
    INS = P[64:]; ptb = 1
    for p in sorted(B0)[-r:]: ptb *= p
    topI = 1
    for p in INS[-r:]: topI *= p
    Icap = max(1, min(topI + 1, (U * ptb) // P0 + 1))
    Imin1 = 1
    for p in INS[:r - 1]: Imin1 *= p
    cap_r = (Icap - 1) // Imin1; INS = [q for q in INS if q <= cap_r]
    if comb(len(INS), r) >= 2 ** 64: print(f"r={r}: |INS|={len(INS)} rank ids overflow; skipped"); return True
    large = [p for p in INS if p > split]; small = [p for p in INS if p <= split]; Dmax = ptb
    hmax = -1
    for h in range(0, r + 1):
        if h > len(large) or r - h > len(small): continue
        v = P0
        for p in small[:r - h]: v *= p
        for p in large[:h]: v *= p
        if v < U * Dmax: hmax = h
    if hmax < 0: print(f"r={r}: no stratum can improve; skipped"); return True
    est = comb(len(INS), r)
    if est > 3_000_000_000: print(f"r={r}: |INS|={len(INS)} C(n,r)={est:.2e} too large for an in-memory comparison; skipped"); return True
    cpu, tc = cpu_records(M, r, P0, Icap, B0, INS, split, hmax)
    if len(cpu) > 300_000_000: print(f"r={r}: {len(cpu):,} CPU records too many for an in-memory comparison; skipped"); return True
    t = time.time(); k, i, cnt = igen_gpu.generate(INS, r, M, P0, Icap, split, hmax); cp.cuda.Device().synchronize(); tg = time.time() - t
    gk = k.get() if cnt else np.zeros(0, np.uint64); gi = i.get() if cnt else np.zeros(0, np.uint64); del k, i
    ck = np.ascontiguousarray(cpu[:, 0]); ci = np.ascontiguousarray(cpu[:, 1]); del cpu
    # rank ids are unique per subset: compare id sets, then keys on the common ids
    oc = np.argsort(ci, kind="stable"); ci_s = ci[oc]; ck_s = ck[oc]
    og = np.argsort(gi, kind="stable"); gi_s = gi[og]; gk_s = gk[og]
    assert len(np.unique(ci_s)) == len(ci_s) and len(np.unique(gi_s)) == len(gi_s), "duplicate ids"
    missing_ids = np.setdiff1d(ci_s, gi_s, assume_unique=True); extra_ids = np.setdiff1d(gi_s, ci_s, assume_unique=True)
    common = np.intersect1d(ci_s, gi_s, assume_unique=True)
    kc = ck_s[np.searchsorted(ci_s, common)]; kg = gk_s[np.searchsorted(gi_s, common)]; key_mismatch = int((kc != kg).sum())
    bad_extra = 0
    for rid in extra_ids[:100000]:
        idx = unrank(rid, r, len(INS)); pr_ = 1
        for x in idx: pr_ *= INS[x]
        if pr_ < Icap: bad_extra += 1
    ok = len(missing_ids) == 0 and key_mismatch == 0 and bad_extra == 0 and len(gi) == cnt
    print(f"r={r} |INS|={len(INS)} hmax={hmax}: cpu {len(ci):,} ({tc:.2f}s)  gpu {cnt:,} ({tg:.2f}s)  missing={len(missing_ids)} extra={len(extra_ids)} key_mismatch={key_mismatch} bad_extra={bad_extra}  {'OK' if ok else 'FAIL'}")
    return ok

def main():
    incJ = json.load(open(os.path.join(FR, "results_k64_best_global.json")))
    N0 = int(json.load(open(os.path.join(FR, "results_k64_best_global_prev_148_539839443357.json")))["n"])
    prevN = int(json.load(open(os.path.join(FR, "results_k64_best_global_prev_148_N.json")))["n"])
    M0 = 1768248177696000; allok = True
    print("== M0, pool<4000, U = previous incumbent (rediscovery instance)")
    for r in range(3, 9): allok &= run_case(M0, prevN, 4000, r)
    print("== M0, full pool (cap 82178), U = N0")
    for r in range(3, 9): allok &= run_case(M0, N0, 82178, r)
    print("== lambda(N), pool<65287, U = N0")
    for r in range(3, 7): allok &= run_case(int(incJ["modulus"]), N0, 65287, r)
    print("ALL OK" if allok else "FAILURES PRESENT"); sys.exit(0 if allok else 1)

if __name__ == "__main__": main()
