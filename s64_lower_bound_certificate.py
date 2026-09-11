#!/usr/bin/env python3
"""Reproduce the computational theorem S_64 >= 10**145 (146 digits).

Python 3.10+, standard library only. Default invocation:
    python s64_lower_bound_certificate.py

This file is self-contained: it first recomputes an unconditional cofactor
lower bound, then independently enumerates all possible competitors twice.
No Carmichael oracle, candidate modulus, input prime table, or engine code
is imported. Integer arithmetic is exact; no logarithmic pruning is used.

PROOF OF COVERAGE
1. Every Carmichael number with 64 distinct prime factors is odd. If it were
   even, an odd factor p would make the even integer p-1 divide the odd n-1.
2. Its primes are pairwise admissible: for p<q dividing n, p cannot divide
   q-1, since q-1 divides n-1 whereas p divides n.
3. Remove any factor. For each possible intersection S of the remaining 63
   factors with the odd primes through 97, those remaining factors have
   product at least product(S) times the product of the cheapest 63-|S|
   primes above 97 individually compatible with S. Ignoring conflicts among
   these tail primes can only LOWER that product. Enumerate all internally
   admissible S. For every S, check that the finite prime table contains
   enough eligible tail primes, so its truncation cannot omit a smaller one.
   Taking the minimum gives the unconditional bound C printed below.
4. Thus a competing n<Y has every factor q <= floor((Y-1)/C). At Y=10**145
   this is 5518, containing 727 odd primes (largest 5507).
5. Enumerate every increasing admissible 64-subset of this universe with
   product <Y. Both enumerators prune only when there are too few remaining
   primes, or when the product of the smallest required remaining primes
   already reaches Y. Both are necessary conditions.
6. For every complete product, test Korselt. Enumerator A uses the exact
   incrementally maintained lcm. Enumerator B separately generates primes
   by trial division, uses an ordered-choice traversal and direct divisibility
   by each individual p-1, without using the lcm or the bitmask graph.

DEFAULT RESULT
    cofactor classes: 407760
    universe: 727 odd primes, q <= 5518
    binary enumeration nodes: 327049
    ordered enumeration nodes: 163525
    admissible complete products: 127092 (both enumerators)
    Carmichael products: 0 (both enumerators)
    conclusion: S_64 >= 10**145

The matching leaf digests are an audit aid, not a replacement for the
coverage argument. This is a reproducible computational proof, not a
formally verified proof-assistant development or a proof of minimality.

Optional --exponent changes Y. --seconds bounds EACH enumeration and a
timeout explicitly produces INCOMPLETE, never an exclusion. Setup and the
cofactor-bound computation are outside this limit. With --one-enumerator,
only A runs; its mathematical coverage is unchanged but the independent
cross-check is omitted. These options support bounded follow-up experiments.
"""

import argparse
import hashlib
import json
import sys
from math import isqrt, lcm, prod
from time import perf_counter


EXPECTED_C = int(
    "1812140595092971804464959059612576188399713525008997620073689836850955613086867158799858934921557347504370550240044086773665823716701403213549"
)


def trial_primes(limit):
    return [p for p in range(3, limit + 1, 2)
            if all(p % d for d in range(2, isqrt(p) + 1))]


def cofactor_lower_bound():
    ps = trial_primes(10000)
    head = [p for p in ps if p <= 97]
    tail = [p for p in ps if p > 97]
    allow = [sum(1 << i for i, q in enumerate(tail) if (q - 1) % p)
             for p in head]
    best = None
    classes = 0

    def visit(start, chosen, P, mask):
        nonlocal best, classes
        classes += 1
        need = 63 - len(chosen)
        if mask.bit_count() < need:
            raise RuntimeError("Insufficient tail table; no bound certified")
        value, bits = P, mask
        for _ in range(need):
            bit = bits & -bits
            bits -= bit
            value *= tail[bit.bit_length() - 1]
        best = value if best is None else min(best, value)
        for i in range(start, len(head)):
            q = head[i]
            if all((q - 1) % p for p in chosen):
                visit(i + 1, chosen + (q,), P * q, mask & allow[i])

    visit(0, (), 1, (1 << len(tail)) - 1)
    if best != EXPECTED_C or classes != 407760:
        raise RuntimeError("Cofactor-bound regression mismatch")
    return best, classes


def sieve_primes(limit):
    s = bytearray([1]) * (limit + 1)
    s[0:2] = b"\0\0"
    for p in range(2, isqrt(limit) + 1):
        if s[p]:
            s[p * p:limit + 1:p] = b"\0" * ((limit - p * p) // p + 1)
    return [p for p in range(3, limit + 1) if s[p]]


class Results:
    def __init__(self, seconds):
        self.nodes = 0
        self.leaves = 0
        self.minimum = None
        self.matches = []
        self.digest = hashlib.sha256()
        self.started = perf_counter()
        self.seconds = seconds
        self.complete = False

    def node(self):
        self.nodes += 1
        if (self.seconds is not None and self.nodes % 4096 == 0
                and perf_counter() - self.started > self.seconds):
            raise TimeoutError

    def leaf(self, n, passes):
        self.leaves += 1
        self.minimum = n if self.minimum is None else min(self.minimum, n)
        self.digest.update((str(n) + "\n").encode("ascii"))
        if passes:
            self.matches.append(n)

    def export(self):
        return {
            "complete": self.complete,
            "nodes": self.nodes,
            "admissible_products": self.leaves,
            "minimum_admissible_product": str(self.minimum),
            "carmichael_products": [str(n) for n in self.matches],
            "leaf_sha256": self.digest.hexdigest(),
            "enumeration_seconds": perf_counter() - self.started,
        }


def binary_enumerator(Y, limit, seconds):
    ps = sieve_primes(limit)
    future = [sum(1 << j for j in range(i + 1, len(ps))
                  if (ps[j] - 1) % p) for i, p in enumerate(ps)]
    sys.setrecursionlimit(max(sys.getrecursionlimit(), len(ps) + 200))
    out = Results(seconds)

    def walk(mask, t, P, L):
        out.node()
        if t == 0:
            if P >= Y:
                raise RuntimeError("Leaf outside strict product bound")
            out.leaf(P, (P - 1) % L == 0)
            return
        if mask.bit_count() < t:
            return
        value, bits = P, mask
        for _ in range(t):
            bit = bits & -bits
            bits -= bit
            value *= ps[bit.bit_length() - 1]
        if value >= Y:
            return
        bit = mask & -mask
        i = bit.bit_length() - 1
        walk(mask & future[i], t - 1, P * ps[i], lcm(L, ps[i] - 1))
        walk(mask ^ bit, t, P, L)

    try:
        walk((1 << len(ps)) - 1, 64, 1, 1)
        out.complete = True
    except TimeoutError:
        pass
    return out.export(), len(ps)


def ordered_enumerator(Y, limit, seconds):
    ps = tuple(trial_primes(limit))
    out = Results(seconds)

    def walk(chosen, P, pool):
        out.node()
        t = 64 - len(chosen)
        if t == 0:
            if P >= Y:
                raise RuntimeError("Leaf outside strict product bound")
            out.leaf(P, all((P - 1) % (p - 1) == 0 for p in chosen))
            return
        for i in range(len(pool) - t + 1):
            # Sliding products increase, so failure also excludes later i.
            if P * prod(pool[i:i + t]) >= Y:
                break
            p = pool[i]
            following = tuple(q for q in pool[i + 1:] if (q - 1) % p)
            walk(chosen + (p,), P * p, following)

    try:
        walk((), 1, ps)
        out.complete = True
    except TimeoutError:
        pass
    return out.export()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exponent", type=int, default=145)
    parser.add_argument("--seconds", type=float, default=None)
    parser.add_argument("--one-enumerator", action="store_true")
    args = parser.parse_args()
    if not 144 <= args.exponent <= 147:
        parser.error("Supported experimental exponents: 144 through 147")
    if args.seconds is not None and args.seconds <= 0:
        parser.error("--seconds must be positive")
    started = perf_counter()
    C, classes = cofactor_lower_bound()
    Y = 10 ** args.exponent
    limit = (Y - 1) // C
    first, count = binary_enumerator(Y, limit, args.seconds)
    result = {
        "strict_upper_bound": str(Y),
        "cofactor_lower_bound": str(C),
        "cofactor_classes": classes,
        "largest_prime_bound_inclusive": limit,
        "odd_primes_in_universe": count,
        "binary": first,
    }
    complete = first["complete"]
    if not args.one_enumerator and complete:
        second = ordered_enumerator(Y, limit, args.seconds)
        result["ordered_independent"] = second
        complete = second["complete"]
        if complete:
            for key in ("admissible_products", "minimum_admissible_product",
                        "carmichael_products", "leaf_sha256"):
                if first[key] != second[key]:
                    raise RuntimeError("Independent enumeration mismatch: " + key)
    if complete and not first["carmichael_products"]:
        result["conclusion"] = "CERTIFIED: S_64 >= 10**" + str(args.exponent)
    elif first["carmichael_products"]:
        result["conclusion"] = "COUNTEREXAMPLE(S) FOUND; proposed lower bound fails"
    else:
        result["conclusion"] = "INCOMPLETE: no exclusion certified by this run"
    result["total_seconds"] = perf_counter() - started
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
