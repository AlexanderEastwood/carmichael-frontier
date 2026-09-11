#!/usr/bin/env python3
"""Divisor-family census for a k=64 incumbent (stdlib only; independent of the GPU code).

For an incumbent U with lambda(U) = M: every Carmichael n with omega(n)=64, n < U and lambda(n) | M
has lambda(n) = D for a divisor D of M, all primes in Q(D) = {q prime : (q-1) | D, q ∤ D}, and
n ≡ 1 (mod D). D can carry such an n only if the 64 smallest primes of Q(D) multiply to < U
("survivor"). For a survivor, any 64-subset of Q(D) with product < U uses only primes
q <= cap(D) = (U-1)/T63(D) and differs from the 64 smallest primes of Q(D) in at most rmax(D) places
(rmax = the largest r for which swapping the r largest of the 64 smallest for the r smallest primes
outside stays below U; the ratio grows with r, so this is size-forced).

usage: python3 divisor_family_census.py [results.json] [out.json]
"""
import json, sys, os
from math import isqrt

def factor(m):
    f = {}; d = 2
    while d * d <= m:
        while m % d == 0: f[d] = f.get(d, 0) + 1; m //= d
        d += 1
    if m > 1: f[m] = f.get(m, 0) + 1
    return f

def divisors(fm):
    ds = [1]
    for p, e in fm.items(): ds = [d * p ** k for d in ds for k in range(e + 1)]
    return sorted(ds)

SIEVE_BOUND = 2_000_000          # every survivor's cap is checked to lie below this (else the census aborts)
def sieve(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, isqrt(n) + 1):
        if s[i]: s[i * i::i] = bytearray(len(s[i * i::i]))
    return s
_IS_PRIME = sieve(SIEVE_BOUND)

def pool(D):
    """Q(D) ∩ [3, SIEVE_BOUND]: primes q with (q-1) | D and q ∤ D, generated from the divisors of D (q = d+1)."""
    out = []
    for d in divisors(factor(D)):
        q = d + 1
        if 3 <= q <= SIEVE_BOUND and D % q != 0 and _IS_PRIME[q]: out.append(q)
    return sorted(out)

def census(U, M, U_factors):
    fm = factor(M); rows = []
    for D in divisors(fm):
        P = pool(D)
        if len(P) < 64:
            # finite-sieve completeness guard (Astra 2026-09-11): the class may only be discarded if no completion
            # using primes beyond the sieve could be < U; every omitted prime is > SIEVE_BOUND, so the product of the
            # s known primes times (SIEVE_BOUND+1)^(64-s) is a lower bound on any 64-product from Q(D).
            lb = 1
            for p in P: lb *= p
            lb *= (SIEVE_BOUND + 1) ** (64 - len(P))
            assert lb >= U, f"D={D}: only {len(P)} eligible primes below the sieve bound and the completion bound {len(str(lb))} digits < U; raise SIEVE_BOUND"
            continue
        T64 = 1
        for p in P[:64]: T64 *= p
        if T64 >= U: continue
        T63 = T64 // P[63]; cap = (U - 1) // T63
        assert cap < SIEVE_BOUND, f"D={D}: cap {cap} exceeds the sieve bound; raise SIEVE_BOUND"
        Pc = [q for q in P if q <= cap]
        rmax = 0
        for r in range(1, 65):
            if 64 + r > len(Pc): break
            num = T64
            for q in Pc[64:64 + r]: num *= q
            den = 1
            for q in P[64 - r:64]: den *= q
            if num // den < U: rmax = r
            else: break
        s = sum(1 for p in U_factors if D % (p - 1) != 0)   # incumbent factors forced out by D
        rows.append({"D": D, "s": s, "pool": len(P), "cap": cap, "pool_capped": len(Pc), "rmax": rmax, "digits_T64": len(str(T64))})
    return fm, rows

def main():
    here = os.path.dirname(os.path.abspath(__file__)); fr = os.path.dirname(here)
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(fr, "results_k64_best_global.json")
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "divisor_family.json")
    inc = json.load(open(src)); U = int(inc["n"]); F = sorted(int(p) for p in inc["factors"])
    from math import lcm
    M = lcm(*[p - 1 for p in F]); assert M == int(inc["modulus"]), "incumbent modulus is not lambda(U)"
    fm, rows = census(U, M, F)
    print(f"U = {len(str(U))} digits, lambda(U) = {M} = {fm}; divisors {len(divisors(fm))}; survivors {len(rows)}")
    print("D s |pool| cap |pool<=cap| rmax digits(T64)")
    for r in sorted(rows, key=lambda r: (r["s"], -r["D"])): print(r["D"], r["s"], r["pool"], r["cap"], r["pool_capped"], r["rmax"], r["digits_T64"])
    json.dump({"U": str(U), "lambda": str(M), "lambda_factorization": {str(k): v for k, v in fm.items()}, "survivors": rows}, open(out, "w"), indent=1)
    print("wrote", out)

if __name__ == "__main__": main()
