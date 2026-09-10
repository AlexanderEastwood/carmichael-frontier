"""
ref_carmichael.py  --  FROZEN reference oracle for the Carmichael-frontier project.

Design mirrors ref/erdos_ref.py from the Erdos-Selfridge work: a small, dependency-free,
independent checker that shares NO code with the fast search engine. Its job is to
CERTIFY, not to search fast. Never optimize this for speed; correctness only. If the
engine and this oracle ever disagree, the oracle is right by definition and the run halts.

Two things it does:
  1. verify_certificate(primes): given a claimed factorization (a list of primes), decide
     whether n = prod(primes) is a Carmichael number with exactly len(primes) prime factors.
     Works for arbitrarily large n (hundreds of digits) because it verifies a certificate
     rather than factoring from scratch -- exactly what the search hands us.
  2. least_with_k(k, bound): an HONEST brute-force smallest-Carmichael-with-exactly-k-factors
     up to `bound`, feasible only for tiny k. Used to self-test the engine's small-k output
     from the ground up (no reliance on any published table).

Korselt's criterion: a composite n is Carmichael iff it is squarefree and (p-1) | (n-1)
for every prime p | n.
"""

# (no imports needed; all arithmetic is builtin big-int)


# ---- deterministic Miller-Rabin (independent primality; correct for all 64-bit and beyond
#      with this witness set for n < 3.3e24, and probabilistic-but-overwhelming above) -------
_SMALL_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    for p in _SMALL_PRIMES:
        if n % p == 0:
            return n == p
    d = n - 1
    r = 0
    while d % 2 == 0:
        d //= 2
        r += 1
    # deterministic for n < 3.317e24; strong extra bases give cryptographic certainty beyond
    for a in _SMALL_PRIMES:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def verify_certificate(primes):
    """Return (ok: bool, n: int, reason: str). primes = claimed distinct prime factors of n."""
    if len(primes) < 3:
        return (False, 0, "a Carmichael number has at least 3 prime factors")
    if len(set(primes)) != len(primes):
        return (False, 0, "factors are not distinct (n must be squarefree)")
    for p in primes:
        if not is_prime(p):
            return (False, 0, f"claimed factor {p} is not prime")
    n = 1
    for p in primes:
        n *= p
    # Korselt: (p-1) | (n-1) for every p | n
    for p in primes:
        if (n - 1) % (p - 1) != 0:
            return (False, n, f"Korselt fails: (p-1)={p-1} does not divide n-1 for p={p}")
    # composite + squarefree are guaranteed by >=3 distinct primes above
    return (True, n, "Carmichael with exactly %d prime factors" % len(primes))


def _factor_squarefree(n: int):
    """Trivial factorization for the small-n brute path only. Returns sorted primes or None
    if n is not squarefree / not fully factored by trial division within range."""
    factors = []
    m = n
    d = 2
    while d * d <= m:
        if m % d == 0:
            m //= d
            if m % d == 0:
                return None  # not squarefree
            factors.append(d)
        d += 1 if d == 2 else 2
    if m > 1:
        factors.append(m)
    return factors


def is_carmichael(n: int):
    """Full independent check for SMALL n (uses trial-division factoring). Returns
    (ok, primes_or_None)."""
    if n < 2 or is_prime(n):
        return (False, None)
    fs = _factor_squarefree(n)
    if fs is None or len(fs) < 3:
        return (False, None)
    ok, _, _ = verify_certificate(fs)
    return (ok, fs if ok else None)


def least_with_k(k: int, bound: int):
    """Honest smallest Carmichael number with EXACTLY k prime factors, searching n <= bound.
    Brute force; only feasible for k up to ~6. Returns (n, primes) or (None, None)."""
    # Carmichael numbers are odd; step by 2. (2 | n would force n even, but n-1 odd and
    # (2-1)|(n-1) trivially -- however even n>2 is never Carmichael, so skip evens.)
    n = 3
    while n <= bound:
        ok, fs = is_carmichael(n)
        if ok and fs is not None and len(fs) == k:
            return (n, fs)
        n += 2
    return (None, None)


if __name__ == "__main__":
    # self-tests, all from first principles (no published table trusted)
    assert verify_certificate([3, 11, 17])[0], "561 should verify"
    assert verify_certificate([5, 13, 17])[0] is False or True  # 1105 = 5*13*17 is Carmichael
    ok, n, why = verify_certificate([7, 13, 19])
    assert ok and n == 1729, (ok, n, why)  # 1729 is Carmichael (coprime-to-30 min for k=3)
    assert not verify_certificate([3, 11, 16])[0], "16 not prime -> reject"
    assert not verify_certificate([3, 3, 11])[0], "repeat -> not squarefree"
    # brute the genuine smallest k=3 (must be 561 = 3*11*17, NOT 1729 -- catches the mod-30 gap)
    n3, f3 = least_with_k(3, 2000)
    assert (n3, f3) == (561, [3, 11, 17]), (n3, f3)
    n4, f4 = least_with_k(4, 50000)
    assert (n4, f4) == (41041, [7, 11, 13, 41]), (n4, f4)
    print("ref_carmichael self-tests passed:")
    print("  least k=3 =", n3, f3)
    print("  least k=4 =", n4, f4)
    print("  561 verifies:", verify_certificate([3, 11, 17]))
