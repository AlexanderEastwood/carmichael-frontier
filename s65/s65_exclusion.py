#!/usr/bin/env python3
"""Certify S65 >= 10**148; optionally calibrate 10**149 with a time budget.

Run alongside the bundled s65_transfer_certificate.py. The setup recomputes
the global minimum A64 of admissible 64-prime products. It does not assume
that an admissible cofactor is Carmichael. For n<Y, every 65-factor Carmichael
has all its prime factors <=floor((Y-1)/A64).

The DFS keeps the complete remaining pool in an increasing prime bitmask.
Its cardinality and smallest-product bounds are necessary conditions only.
When P*lcm(p-1)>Y, the cofactor progression has at most one term below Y/P.
The first term at or above the safe remaining-product bound is checked by
exact factorization over the entire allowed prime pool and full Korselt.
Completions must have precisely the required number of distinct factors.
Timeouts report INCOMPLETE, never an exclusion.

The native verifier in this bundle does not use the progression switch: it
visits all 22,287,138 admissible full products below 10**148 and tests Korselt.
Both searches share the cofactor reduction, but no implementation code.
"""

import argparse
import json
import sys
from math import lcm
from time import perf_counter
from s65_transfer_certificate import (
    cofactor63, enumerate64_and_extend, sieve_odd_primes, require, Y64,
)


def exclude(k, Y, limit, seconds):
    setup_start = perf_counter()
    ps = sieve_odd_primes(limit)
    index = {p: i for i, p in enumerate(ps)}
    future = [sum(1 << j for j in range(i+1, len(ps))
                  if (ps[j]-1) % p) for i, p in enumerate(ps)]
    sys.setrecursionlimit(max(sys.getrecursionlimit(), len(ps)+300))
    table_seconds = perf_counter()-setup_start
    nodes = inversions = ap_count = max_tail = 0
    hits = []
    started = perf_counter()

    def check_completion(P, L, t, R, mask):
        factors = []
        rem, bits = R, mask
        while bits and rem > 1:
            bit = bits & -bits
            bits -= bit
            p = ps[bit.bit_length()-1]
            if p*p > rem:
                j = index.get(rem)
                if j is None or not (mask >> j & 1):
                    return
                factors.append(rem)
                rem = 1
                break
            if rem % p == 0:
                rem //= p
                factors.append(p)
                if rem % p == 0 or len(factors) > t:
                    return
        if rem != 1 or len(factors) != t:
            return
        n = P*R
        if (n-1) % L == 0 and all((n-1) % (p-1) == 0 for p in factors):
            hits.append(n)

    def dfs(mask, t, P, L):
        nonlocal nodes, inversions, ap_count, max_tail
        nodes += 1
        if seconds is not None and nodes % 16384 == 0:
            if perf_counter()-started > seconds:
                raise TimeoutError
        if t == 0:
            if P < Y and (P-1) % L == 0:
                hits.append(P)
            return
        if mask.bit_count() < t:
            return
        bits, Rmin = mask, 1
        for _ in range(t):
            bit = bits & -bits
            bits -= bit
            Rmin *= ps[bit.bit_length()-1]
        if P*Rmin >= Y:
            return
        if P*L > Y:
            inversions += 1
            max_tail = max(max_tail, t)
            R = Rmin + (pow(P, -1, L)-Rmin) % L
            if R <= (Y-1)//P:
                ap_count += 1
                check_completion(P, L, t, R, mask)
            return
        bit = mask & -mask
        i = bit.bit_length()-1
        dfs(mask & future[i], t-1, P*ps[i], lcm(L, ps[i]-1))
        dfs(mask ^ bit, t, P, L)

    complete = False
    try:
        dfs((1 << len(ps))-1, k, 1, 1)
        complete = True
    except TimeoutError:
        pass
    return {
        "complete": complete, "k": k, "strict_bound": str(Y),
        "prime_cap": limit, "odd_primes": len(ps), "nodes": nodes,
        "inversions": inversions, "AP_candidates": ap_count,
        "maximum_tail_count": max_tail,
        "Carmichael": [str(n) for n in hits],
        "table_seconds": table_seconds, "search_seconds": perf_counter()-started,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exponent', type=int, choices=(148,149), default=148)
    parser.add_argument('--seconds', type=float,
                        help='Search-only budget; timeout is INCOMPLETE')
    args = parser.parse_args()
    if args.seconds is not None and args.seconds <= 0:
        parser.error('--seconds must be positive')
    started = perf_counter()
    C63, classes = cofactor63()
    base, extended, metadata = enumerate64_and_extend(C63, 317*Y64)
    A64 = base.minimum
    Y = 10**args.exponent
    setup_seconds = perf_counter()-started
    result = exclude(65, Y, (Y-1)//A64, args.seconds)
    if args.exponent == 148 and result['complete']:
        require(result['nodes'] == 64355 and result['inversions'] == 5012
                and result['AP_candidates'] == 66, '148 regression mismatch')
    if result['Carmichael']:
        conclusion = 'COUNTEREXAMPLE FOUND: claimed lower bound fails'
    elif result['complete']:
        conclusion = 'CERTIFIED: S65 >= 10**'+str(args.exponent)
    else:
        conclusion = 'INCOMPLETE: no exclusion established by this run'
    print(json.dumps({
        'conclusion': conclusion, 'A64': str(A64),
        'cofactor_classes': classes, 'base64_products': base.count,
        'base64_leaf_sha256': base.digest.hexdigest(),
        'cofactor_setup_seconds': setup_seconds,
        'search': result, 'total_seconds': perf_counter()-started,
    }, indent=2))


if __name__ == '__main__':
    main()
