#!/usr/bin/env python3
"""Feasibility-pruned deeper frontier for the k=64 exclusion below 10^147.

Like sub_gen.py, extends given prefixes by admissible primes -- but only along
children that can still complete: a child prefix (p_1..p_d) is emitted/expanded
only if  prod(prefix) * prod(the 64-d smallest primes > p_d)  <  10^147.
That product uses the plain next primes (ignoring admissibility among them),
which is <= the true minimal completion, so the prune is SAFE: any prefix it
drops has NO 64-factor completion below the bound at all, and therefore needs
no coverage. This keeps a +L-level frontier to the branches that carry real
work instead of ~30^L mostly-empty leaves.

args: "p1,p2,...:levels" ...   (each prefix extended by exactly `levels` primes)
"""
import sys
from math import isqrt
B = 10**147; K = 64
def sieve(n):
    s = [True]*(n+1); s[0] = s[1] = False
    for i in range(2, isqrt(n)+1):
        if s[i]:
            for j in range(i*i, n+1, i): s[j] = False
    return [p for p in range(3, n+1) if s[p]]
PR = sieve(6000); IDX = {p: i for i, p in enumerate(PR)}
PREF = [1]
for p in PR: PREF.append(PREF[-1]*p)             # PREF[i] = product of PR[:i]
def min_completion(prod, last_idx, t):            # prod * next t primes after PR[last_idx]
    a = last_idx + 1
    return prod * (PREF[a+t] // PREF[a]) if a + t <= len(PR) else None
def admissible(p, chosen): return all((p-1) % c for c in chosen)
def extend(prefix, prod, levels, out):
    d = len(prefix); t = K - d
    if levels == 0:
        out.append(",".join(map(str, prefix))); return
    for i in range(IDX[prefix[-1]] + 1, len(PR)):
        p = PR[i]
        if not admissible(p, prefix): continue
        np_ = prod * p; mc = min_completion(np_, i, t - 1)
        if mc is None or mc >= B: break            # PR ascending: later children only larger
        extend(prefix + [p], np_, levels - 1, out)
out = []
for arg in sys.argv[1:]:
    ps, lv = arg.split(":"); pref = [int(x) for x in ps.split(",")]
    pr = 1
    for p in pref: pr *= p
    extend(pref, pr, int(lv), out)
print("\n".join(out))
