#!/usr/bin/env python3
"""Transfer the S64 admissible-product exhaustion to a rigorous S65 bound.

Python 3.10+, standard library only:
    python s65_transfer_certificate.py

Expected conclusion: S65 >= 317 * 10**145 = 3.17 * 10**147.
This does not determine S65 or establish novelty of the bound.

PROOF OF COVERAGE
1. Every divisor of a Carmichael number is a squarefree product of odd
   primes with p not dividing q-1 for p<q. Call such products admissible.
2. Recompute the original lower bound C63 for ANY admissible 63-prime
   product, by partitioning its primes through 97 and ignoring conflicts
   among its cheapest eligible tail primes. The prime table is complete
   through 10000; insufficient tails cause failure, never exclusion.
3. ANY admissible 64-prime product D < Y=10**145 has every prime factor
   at most floor((Y-1)/C63)=5518. This statement does not require D itself
   to be Carmichael. Enumerate all these D with ordered choices and safe
   cardinality/product pruning. The count and digest reproduce the S64
   certificate. Record their minimum A64.
4. A64 is the GLOBAL minimum admissible 64-prime product: A64<Y, and
   every smaller admissible product was included in step 3.
5. A 65-factor Carmichael n has largest factor q >=317 (the 65th odd
   prime). If n<317Y, its 64-prime cofactor D=n/q is <Y, so it occurs in
   step 3. Test every prime extension maxprime(D)<q<317Y/D. Each 65-prime
   product appears once, at its unique largest prime. Direct Korselt and
   a separate residue-based extension check must give identical answers.
6. Cross-check with a direct 65-prime enumeration, without extending the
   recorded 64-prime products. Step 4 bounds its universe by
   floor((317Y-1)/A64). Its product count, minimum, digest, and Carmichael
   list must agree with step 5.
7. Both finish with 285 admissible 65-prime products and no Carmichael.

The minimum A64 also bounds EVERY 64-prime cofactor of any 65-factor
Carmichael. Consequently n<U implies P+(n)<=floor((U-1)/A64).
The 65-product minimum recorded here is global for the same reason as
step 4, but is a minimum among admissible products, not Carmichael numbers.

No floating point, probable-prime tests, candidate modulus, or external
oracle is used. The cross-checks share the mathematical reductions and
some utilities; they are not three independent end-to-end proofs.
"""

import hashlib
import json
from math import isqrt, lcm, prod
from time import perf_counter

Y64 = 10**145
EXPECTED_C63 = int(
    "1812140595092971804464959059612576188399713525008997620073689836850955613086867158799858934921557347504370550240044086773665823716701403213549"
)
EXPECTED_DIGEST64 = (
    "e05f62598eea784b3f60f025c6ccbf30700115ecbef42702e24154dcece32321"
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def trial_odd_primes(limit):
    return [n for n in range(3, limit + 1, 2) if prime(n)]


def sieve_odd_primes(limit):
    flags = bytearray([1]) * (limit + 1)
    flags[:2] = b"\0\0"
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            flags[p*p:limit+1:p] = b"\0" * ((limit-p*p)//p + 1)
    return [p for p in range(3, limit + 1, 2) if flags[p]]


def cofactor63():
    ps = trial_odd_primes(10000)
    head = [p for p in ps if p <= 97]
    tail = [p for p in ps if p > 97]
    allowed = [sum(1 << i for i, q in enumerate(tail) if (q-1) % p)
               for p in head]
    best = None
    classes = 0

    def visit(start, chosen, P, mask):
        nonlocal best, classes
        classes += 1
        need = 63 - len(chosen)
        require(mask.bit_count() >= need, "Insufficient prime table")
        value, bits = P, mask
        for _ in range(need):
            bit = bits & -bits
            bits -= bit
            value *= tail[bit.bit_length() - 1]
        best = value if best is None else min(best, value)
        for i in range(start, len(head)):
            p = head[i]
            if all((p-1) % q for q in chosen):
                visit(i+1, chosen+(p,), P*p, mask & allowed[i])

    visit(0, (), 1, (1 << len(tail)) - 1)
    require(best == EXPECTED_C63 and classes == 407760,
            "Cofactor computation mismatch")
    return best, classes


class Record:
    def __init__(self):
        self.count = 0
        self.minimum = None
        self.hits = []
        self.digest = hashlib.sha256()

    def add(self, n, passes):
        self.count += 1
        self.minimum = n if self.minimum is None else min(self.minimum, n)
        self.digest.update((str(n) + "\n").encode("ascii"))
        if passes:
            self.hits.append(n)

    def export(self):
        return {"admissible_products": self.count,
                "minimum": str(self.minimum),
                "carmichael_products": [str(n) for n in self.hits],
                "leaf_sha256": self.digest.hexdigest()}


def enumerate64_and_extend(C63, Y65):
    base = Record()
    extensions = Record()
    nodes = extendable = ap_candidates = 0
    ap_hits = []
    ps = tuple(sieve_odd_primes((Y64-1)//C63))

    def walk(chosen, P, L, pool):
        nonlocal nodes, extendable, ap_candidates
        nodes += 1
        t = 64 - len(chosen)
        if t == 0:
            require(P < Y64, "64-product out of range")
            base.add(P, (P-1) % L == 0)
            lo, hi = chosen[-1]+1, (Y65-1)//P
            if lo > hi:
                return
            extendable += 1
            r = pow(P, -1, L)
            first = r + ((lo-r+L-1)//L)*L
            for q in range(first, hi+1, L):
                ap_candidates += 1
                if prime(q) and (P-1) % (q-1) == 0:
                    ap_hits.append(P*q)
            for q in range(lo, hi+1):
                if prime(q) and all((q-1) % p for p in chosen):
                    n = P*q
                    extensions.add(n, all((n-1) % (p-1) == 0
                                          for p in chosen+(q,)))
            return
        for i in range(len(pool)-t+1):
            if P*prod(pool[i:i+t]) >= Y64:
                break
            p = pool[i]
            following = tuple(q for q in pool[i+1:] if (q-1) % p)
            walk(chosen+(p,), P*p, lcm(L, p-1), following)

    walk((), 1, 1, ps)
    require(base.count == 127092, "64-product count mismatch")
    require(base.digest.hexdigest() == EXPECTED_DIGEST64,
            "64-product digest mismatch")
    require(base.minimum < Y64, "Missing witness for global A64")
    require(sorted(ap_hits) == sorted(extensions.hits),
            "Residue/direct extension disagreement")
    return base, extensions, {"nodes": nodes, "extendable_prefixes": extendable,
                              "progression_candidates": ap_candidates}


def direct65(Y65, A64):
    # Different traversal: binary include/exclude over a bitmask graph.
    ps = trial_odd_primes((Y65-1)//A64)
    following = [sum(1 << j for j in range(i+1, len(ps))
                     if (ps[j]-1) % p) for i, p in enumerate(ps)]
    out = Record()
    nodes = 0

    def walk(mask, t, P, L):
        nonlocal nodes
        nodes += 1
        if t == 0:
            require(P < Y65, "65-product out of range")
            out.add(P, (P-1) % L == 0)
            return
        if mask.bit_count() < t:
            return
        value, bits = P, mask
        for _ in range(t):
            bit = bits & -bits
            bits -= bit
            value *= ps[bit.bit_length()-1]
        if value >= Y65:
            return
        bit = mask & -mask
        i = bit.bit_length()-1
        walk(mask & following[i], t-1, P*ps[i], lcm(L, ps[i]-1))
        walk(mask ^ bit, t, P, L)

    walk((1 << len(ps))-1, 65, 1, 1)
    return out, {"nodes": nodes, "prime_cap": (Y65-1)//A64,
                 "odd_primes": len(ps)}


def main():
    started = perf_counter()
    odd_primes = trial_odd_primes(1000)
    q_min = odd_primes[64]
    require(q_min == 317, "65th odd prime mismatch")
    Y65 = q_min*Y64
    C63, classes = cofactor63()
    base, extended, extension_meta = enumerate64_and_extend(C63, Y65)
    direct, direct_meta = direct65(Y65, base.minimum)
    require(extended.export() == direct.export(),
            "Direct 65-enumeration/extension mismatch")
    require(extended.count == 285, "65-product regression mismatch")
    require(not base.hits and not extended.hits, "Carmichael found")
    print(json.dumps({
        "conclusion": "CERTIFIED: S65 >= 317 * 10**145",
        "strict_65_exclusion_bound": str(Y65),
        "cofactor63_bound": str(C63), "cofactor_classes": classes,
        "base64": base.export(),
        "extended65": extended.export(), "extension_metadata": extension_meta,
        "direct65": direct.export(), "direct_metadata": direct_meta,
        "largest_prime_caps_for_n_below": {
            "10**"+str(e): (10**e-1)//base.minimum for e in (148,149,150,151)},
        "seconds": perf_counter()-started,
        "scope": "Global lower bound and admissible-product minima; not S65 itself"
    }, indent=2))


if __name__ == "__main__":
    main()
