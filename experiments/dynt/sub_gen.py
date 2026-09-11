#!/usr/bin/env python3
"""Emit one-level-deeper admissible extensions of given prefixes.
args: "p1,p2,...:levels" ... -> for each, all admissible prefixes adding `levels` primes."""
import sys
from math import isqrt
def sieve(n):
    s=[True]*(n+1); s[0]=s[1]=False
    for i in range(2,isqrt(n)+1):
        if s[i]:
            for j in range(i*i,n+1,i): s[j]=False
    return [p for p in range(3,n+1) if s[p]]
PR=sieve(6000)
def admissible(p, chosen): return all((p-1)%c for c in chosen)
def extend(prefix, levels, out):
    if levels==0: out.append(",".join(map(str,prefix))); return
    for p in PR:
        if p>prefix[-1] and admissible(p, prefix): extend(prefix+[p], levels-1, out)
out=[]
for arg in sys.argv[1:]:
    ps, lv = arg.split(":"); extend([int(x) for x in ps.split(",")], int(lv), out)
print("\n".join(out))
