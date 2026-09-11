#!/usr/bin/env python3
"""lambda-neighbourhood census (experiment A, 2026-09-11).

Our exchange search is exhaustive inside a family lambda(n) | M but cannot propose a modulus; Webster's
N lies one step from our M0: lambda(N) = (M0 / 3) * 263. This script enumerates that kind of step
systematically and finds every NEW candidate modulus that could carry a 64-factor Carmichael number
below the bound U.

Neighbours of a seed modulus M:  M' = (M / q^j) * p  for every prime power q^e || M, 0 <= j <= e, and
every new prime p ∤ M with p <= PMAX. New candidate lambdas are the divisors D | M' with p | D (divisors
without p are divisors of M and already covered by M's own family). D can carry such an n only if the
64 smallest primes of Q(D) = {r prime : (r-1) | D, r ∤ D} multiply to < U ("survivor"); this test is
monotone (D' | D and D' survivor => D survivor), so survivors are found top-down with pruning.

Output: JSON rows compatible with run_mstar_family.py (D, s, pool, cap, pool_capped, rmax, digits_T64)
plus provenance (seed, q, j, p), sorted by digits_T64 then rmax -- richest pools first.

usage: python3 lambda_neighbourhood.py --seeds M1,M2 --pmax 3000 --out nbhd.json [--procs N]
"""
import argparse, json, os, sys, math
from math import isqrt
from multiprocessing import Pool as MPPool
HERE = os.path.dirname(os.path.abspath(__file__)); FR = os.path.dirname(HERE)
SIEVE_BOUND = 2_000_000
_S = None  # sieve, built lazily per process
def sieve(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, isqrt(n) + 1):
        if s[i]: s[i * i::i] = bytearray(len(s[i * i::i]))
    return s
def factor(m):
    f = {}; d = 2
    while d * d <= m:
        while m % d == 0: f[d] = f.get(d, 0) + 1; m //= d
        d += 1
    if m > 1: f[m] = f.get(m, 0) + 1
    return f
def divisors_of(fm):
    ds = [1]
    for p, e in fm.items(): ds = [d * p ** k for d in ds for k in range(e + 1)]
    return ds
def pool(D, fm):
    """Q(D) ∩ [3, SIEVE_BOUND] from the divisors of D (q = d+1)."""
    out = []
    for d in divisors_of(fm):
        q = d + 1
        if 3 <= q <= SIEVE_BOUND and D % q != 0 and _S[q]: out.append(q)
    return sorted(out)
def pool_relaxed(fm):
    """{q prime : (q-1) | D} WITHOUT the coprimality condition -- monotone under divisibility of D, so it is a
    valid (weaker) pool for pruning descendants (Astra 2026-09-11: Q(D) itself is not monotone, e.g. 3 in Q(80)
    but not in Q(240))."""
    out = []
    for d in divisors_of(fm):
        q = d + 1
        if 3 <= q <= SIEVE_BOUND and _S[q]: out.append(q)
    return sorted(out)
M0_ = 1768248177696000; LN_ = 155016423578016000
def relaxed_survives(fm, U):
    """False only if even the relaxed pool cannot produce a qualifying 64-product below U (then no divisor can
    either): both the cheapest 64-product and the cheapest two-witness 64-product (Theorems 2-3) are monotone
    in D when computed over the relaxed pool, so pruning on them is sound."""
    P = pool_relaxed(fm)
    if len(P) < 64:
        lb = 1
        for p in P: lb *= p
        return lb * (SIEVE_BOUND + 1) ** (64 - len(P)) < U
    T64 = 1
    for p in P[:64]: T64 *= p
    if T64 >= U: return False
    w = two_witness_min(P, U, M0_, LN_)
    return w is not None and w < U
def two_witness_min(P, U, M0, LN):
    """Exact minimum 64-product from pool P containing a prime q with (q-1) ∤ M0 AND a prime q' with (q'-1) ∤ LN
    (Theorems 2-3: every 64-factor Carmichael below N has both). Returns (digits, product) or None if impossible."""
    E0 = [q for q in P if M0 % (q - 1) != 0]; EN = [q for q in P if LN % (q - 1) != 0]
    both = [q for q in P if M0 % (q - 1) != 0 and LN % (q - 1) != 0]
    best = None
    def fill(req):
        rest = [q for q in P if q not in req]
        if len(rest) < 64 - len(req): return None
        v = 1
        for q in req + rest[:64 - len(req)]: v *= q
        return v
    if both: v = fill([both[0]]); best = v if best is None or (v is not None and v < best) else best
    e0 = [q for q in E0 if LN % (q - 1) == 0]; en = [q for q in EN if M0 % (q - 1) == 0]
    if e0 and en:
        v = fill(sorted([e0[0], en[0]])); best = v if best is None or (v is not None and v < best) else best
    return best
def census_row(D, fm, U, Ufactors):
    P = pool(D, fm)
    if len(P) < 64:
        lb = 1
        for p in P: lb *= p
        lb *= (SIEVE_BOUND + 1) ** (64 - len(P))
        return None if lb >= U else ("SIEVE", D)
    T64 = 1
    for p in P[:64]: T64 *= p
    if T64 >= U: return None
    T63 = T64 // P[63]; cap = (U - 1) // T63
    if cap >= SIEVE_BOUND: return ("SIEVE", D)
    Pc = [q for q in P if q <= cap]; rmax = 0
    for r in range(1, 65):
        if 64 + r > len(Pc): break
        num = T64
        for q in Pc[64:64 + r]: num *= q
        den = 1
        for q in P[64 - r:64]: den *= q
        if num // den < U: rmax = r
        else: break
    s = sum(1 for p in Ufactors if D % (p - 1) != 0)
    return {"D": D, "s": s, "pool": len(P), "cap": cap, "pool_capped": len(Pc), "rmax": rmax, "digits_T64": len(str(T64))}

def survivors_with_p(Mp_fm, p, U, Ufactors, seen_global=None):
    """Top-down DFS over exponent vectors of M' (p's exponent fixed at 1): survivors are downward-closed
    in the divisor lattice's complement, i.e. if D is NOT a survivor no divisor of D is."""
    primes = [q for q in Mp_fm if q != p]; emax = [Mp_fm[q] for q in primes]
    seen = {} if seen_global is None else seen_global; out = []; problems = []
    def D_of(vec):
        D = p
        for q, e in zip(primes, vec): D *= q ** e
        return D
    def rec(vec):
        D = D_of(vec); key = D                       # key by the modulus itself so variants share visits
        if key in seen: return seen[key]
        fm = {q: e for q, e in zip(primes, vec) if e} ; fm[p] = 1
        if not relaxed_survives(fm, U): seen[key] = False; return False      # sound prune: monotone relaxed pool
        row = census_row(D, fm, U, Ufactors)
        seen[key] = True
        if isinstance(row, tuple): problems.append(row)
        elif row is not None: out.append(row)
        for i in range(len(vec)):
            if vec[i] > 0:
                v2 = list(vec); v2[i] -= 1; rec(v2)
        return True
    rec(list(emax)); return out, problems

def work(args):
    global _S
    if _S is None: _S = sieve(SIEVE_BOUND)
    seed, seed_fm, p, U, Ufactors = args
    rows = {}; problems = []
    if seed % p == 0: return rows, problems
    # M' = M*p covers every (M/q^j)*p (their divisors are divisors of M*p and the top-down census reaches all
    # relaxed survivors), so only the pure addition and the exponent-increase neighbours M*q*p are distinct roots.
    variants = [(None, 0, dict(seed_fm))]
    for q, e in seed_fm.items():
        fm = dict(seed_fm); fm[q] = e + 1; variants.append((q, -1, fm))
    shared = {}
    for q, j, fm in variants:
        fm2 = dict(fm); fm2[p] = 1
        out, prob = survivors_with_p(fm2, p, U, Ufactors, shared); problems += prob
        for r in out:
            if r["D"] not in rows: rows[r["D"]] = dict(r, seed=str(seed), q=q, j=j, p=p)
    return rows, problems

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--seeds", required=True); ap.add_argument("--pmax", type=int, default=3000)
    ap.add_argument("--out", required=True); ap.add_argument("--procs", type=int, default=os.cpu_count()); ap.add_argument("--incumbent", default=os.path.join(FR, "results_k64_best_global.json"))
    a = ap.parse_args()
    inc = json.load(open(a.incumbent)); U = int(inc["n"]); Uf = sorted(int(x) for x in inc["factors"])
    S = sieve(a.pmax); primes = [p for p in range(3, a.pmax + 1) if S[p]]
    seeds = [int(x) for x in a.seeds.split(",")]
    jobs = [(M, factor(M), p, U, Uf) for M in seeds for p in primes]
    rows = {}; problems = []; covered = set()
    for M in seeds: covered |= set(divisors_of(factor(M)))   # already-exhausted families (the seeds' own divisors)
    with MPPool(a.procs) as mp:
        for r, pr in mp.imap_unordered(work, jobs, chunksize=4):
            problems += pr
            for D, row in r.items():
                if D in covered: continue
                if D not in rows or row["digits_T64"] < rows[D]["digits_T64"]: rows[D] = row
    # exact two-witness filter (Theorems 2-3): drop moduli whose cheapest qualifying 64-product is >= U
    global _S
    if _S is None: _S = sieve(SIEVE_BOUND)
    M0 = 1768248177696000; LN = 155016423578016000; kept = []
    for row in rows.values():
        D = row["D"]; fm = factor(D); P = [q for q in pool(D, fm) if q <= row["cap"]]
        w = two_witness_min(P, U, M0, LN)
        if w is None or w >= U: continue
        row["w2_digits"] = len(str(w))
        kept.append(row)
    out = sorted(kept, key=lambda r: (r["w2_digits"], r["digits_T64"], -r["rmax"], r["D"]))
    json.dump({"U": str(U), "seeds": [str(m) for m in seeds], "pmax": a.pmax, "sieve_problems": problems[:50], "survivors": out}, open(a.out, "w"), indent=1)
    print(f"seeds={seeds} pmax={a.pmax}: new survivor moduli={len(rows)}; after two-witness filter={len(out)}; sieve-bound problems={len(problems)}")
    for r in out[:25]: print(r["D"], "p=%d q=%s j=%d" % (r["p"], r["q"], r["j"]), "T64 digits", r["digits_T64"], "pool", r["pool_capped"], "rmax", r["rmax"], "s", r["s"])
    print("rmax distribution:", {k: sum(1 for r in out if r["rmax"] == k) for k in sorted(set(r["rmax"] for r in out))})

if __name__ == "__main__": main()
