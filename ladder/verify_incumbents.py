#!/usr/bin/env python3
"""Standalone, dependency-free checker for the upper bounds U_k in ladder_manifest.json.

Trusts NOTHING stored: for every row it re-derives the upper endpoint from its
factor list using only exact integer arithmetic --
  * every listed factor is prime by exact trial division (they are all small),
  * the factors are distinct and there are exactly k of them,
  * their product equals the recorded decimal U_k,
  * Korselt: (p-1) | (U_k - 1) for every factor  (so U_k is a Carmichael number,
    being squarefree with >= 3 prime factors),
  * U_k has the recorded digit count.
Python 3 standard library only.  Usage:  python3 ladder/verify_incumbents.py
"""
import json, os, sys
from math import isqrt

def is_prime(n: int) -> bool:
    if n < 2: return False
    if n % 2 == 0: return n == 2
    for d in range(3, isqrt(n) + 1, 2):
        if n % d == 0: return False
    return True

def check_row(r):
    k = r["k"]; U = int(r["U_k"]); fs = [int(x) for x in r["U_k_factors"]]
    problems = []
    if len(fs) != k: problems.append(f"{len(fs)} factors, expected {k}")
    if len(set(fs)) != len(fs): problems.append("repeated factor")
    bad = [p for p in fs if not is_prime(p)]
    if bad: problems.append(f"non-prime factors {bad[:5]}")
    prod = 1
    for p in fs: prod *= p
    if prod != U: problems.append("product of factors != U_k")
    if any((U - 1) % (p - 1) for p in fs): problems.append("Korselt fails")
    if len(str(U)) != r["U_k_digits"]: problems.append("digit count mismatch")
    return problems

def main():
    here = os.path.dirname(os.path.abspath(__file__))
    m = json.load(open(os.path.join(here, "ladder_manifest.json")))
    rows = [r for r in m["rows"] if r.get("U_k") and r.get("U_k_factors")]
    failures = 0
    for r in rows:
        probs = check_row(r)
        if probs: failures += 1; print(f"k={r['k']}: FAIL -- " + "; ".join(probs))
    print(f"checked {len(rows)} upper bounds (k={rows[0]['k']}..{rows[-1]['k']}): "
          + ("ALL VERIFIED -- each U_k is a Carmichael number with exactly k distinct prime factors" if not failures else f"{failures} FAILURES"))
    sys.exit(1 if failures else 0)

if __name__ == "__main__": main()
