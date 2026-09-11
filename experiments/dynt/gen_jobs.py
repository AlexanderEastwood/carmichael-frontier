#!/usr/bin/env python3
"""Emit prefix partition for a complete k=64 exclusion below a bound.
Union of all emitted --prefix jobs = the full search (every smallest prime, the
19-class split by 2nd and (for heavy 2nd primes) 3rd prime for load balance)."""
from math import isqrt
def sieve(n):
    s=[True]*(n+1); s[0]=s[1]=False
    for i in range(2,isqrt(n)+1):
        if s[i]:
            for j in range(i*i,n+1,i): s[j]=False
    return [i for i in range(3,n+1) if s[i]]
P=sieve(6000)
jobs=[]
for p in P:                                   # smallest-prime classes p != 19
    if p!=19: jobs.append(str(p))
for q in [q for q in P if q>19 and (q-1)%19]: # 19-class by admissible 2nd prime
    if q in (23,29,31,37,41,43):              # heavy -> also split by 3rd prime
        for r in [r for r in P if r>q and (r-1)%19 and (r-1)%q]:
            jobs.append(f"19,{q},{r}")
    else:
        jobs.append(f"19,{q}")
print("\n".join(jobs))
