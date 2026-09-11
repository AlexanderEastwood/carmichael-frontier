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
import argparse, json, os, sys
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

def survivors_with_p(Mp_fm, p, U, Ufactors):
    """Top-down DFS over exponent vectors of M' (p's exponent fixed at 1): survivors are downward-closed
    in the divisor lattice's complement, i.e. if D is NOT a survivor no divisor of D is."""
    primes = [q for q in Mp_fm if q != p]; emax = [Mp_fm[q] for q in primes]
    seen = {}; out = []; problems = []
    def D_of(vec):
        D = p
        for q, e in zip(primes, vec): D *= q ** e
        return D
    def rec(vec):
        key = tuple(vec)
        if key in seen: return seen[key]
        D = D_of(vec); fm = {q: e for q, e in zip(primes, vec) if e} ; fm[p] = 1
        row = census_row(D, fm, U, Ufactors)
        if row is None: seen[key] = False; return False
        if isinstance(row, tuple): problems.append(row); seen[key] = True; return True
        seen[key] = True; out.append(row)
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
    variants = [(None, 0, dict(seed_fm))]
    for q, e in seed_fm.items():
        for j in range(1, e + 1):
            fm = dict(seed_fm); fm[q] = e - j
            if fm[q] == 0: del fm[q]
            variants.append((q, j, fm))
    for q, j, fm in variants:
        fm2 = dict(fm); fm2[p] = 1
        out, prob = survivors_with_p(fm2, p, U, Ufactors); problems += prob
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
    out = sorted(rows.values(), key=lambda r: (r["digits_T64"], -r["rmax"], r["D"]))
    json.dump({"U": str(U), "seeds": [str(m) for m in seeds], "pmax": a.pmax, "sieve_problems": problems[:50], "survivors": out}, open(a.out, "w"), indent=1)
    print(f"seeds={seeds} pmax={a.pmax}: new survivor moduli={len(out)}; sieve-bound problems={len(problems)}")
    for r in out[:25]: print(r["D"], "p=%d q=%s j=%d" % (r["p"], r["q"], r["j"]), "T64 digits", r["digits_T64"], "pool", r["pool_capped"], "rmax", r["rmax"], "s", r["s"])
    print("rmax distribution:", {k: sum(1 for r in out if r["rmax"] == k) for k in sorted(set(r["rmax"] for r in out))})

if __name__ == "__main__": main()
