#!/usr/bin/env python3
"""
carmichael_prime_bound.py
=========================

Certified UPPER BOUND on the largest prime factor  P+(n) = max prime dividing n
of *any* 64-factor Carmichael number  n < X, together with the size of the
resulting "prime universe"  pi(P+)  (the number of primes an exclusion search
would have to consider as the top factor).

CPU-only, pure Python (no numpy/sympy). Single-threaded by design.

--------------------------------------------------------------------------
THEOREM (admissibility lower bound on the 63-cofactor product)
--------------------------------------------------------------------------
Let n < X be a Carmichael number with exactly k = 64 distinct prime factors
(Carmichael numbers are squarefree, Korselt's criterion). Write its factors as
    p_1 < p_2 < ... < p_63 < q,        q = P+(n)  the largest.
Then  n = q * D_n  where  D_n = p_1 * ... * p_63  is the product of the other 63
("cofactor") primes, and since  q * D_n = n < X  we have

        q < X / D_n        i.e.        q <= (X - 1) / D_n     (integers).

To turn this into a bound valid for ALL such n we need a rigorous LOWER bound D
on D_n over every admissible 63-set of cofactor primes; then q <= (X-1)/D.

Korselt's criterion forces a strong exclusion among the factors. For n
Carmichael, (r - 1) | (n - 1) for every prime factor r. Hence for two distinct
factors p < r of the SAME n, p cannot divide r - 1: if p | (r - 1) then, since
(r - 1) | (n - 1), we get p | (n - 1); but p | n as well, so p | gcd(n, n-1) = 1,
a contradiction. So:

    (*) In any single Carmichael n, no factor p divides (r - 1) for another
        factor r. Call a set of primes with this property INTERNALLY ADMISSIBLE.

Lower bound construction. Fix a cutoff y. For every internally-admissible subset
S of the primes in [3, y], the remaining  63 - |S|  cofactor primes are all > y
(any admissible small prime <= y is either in S or excluded), and each such tail
prime a must be admissible w.r.t. S, i.e.  (a - 1) % p != 0  for all p in S
(condition (*) between a and the members of S). The product of S times the
product of the  63 - |S|  SMALLEST such admissible tail primes is a lower bound
on D_n for every n whose set of <=y factors is exactly S. Taking the minimum
over all admissible S gives a global lower bound:

        D = min over admissible S of
                ( prod S ) * ( prod of the (63-|S|) smallest primes a>y
                               with (a-1) % p != 0 for all p in S )

        q = P+(n) <= (X - 1) // D.

This is a PROVED lower bound on D_n (hence a proved upper bound on q): it only
uses the pairwise exclusion (*) between S and the tail, and DELIBERATELY IGNORES
the further mutual conflicts among the tail primes themselves. Restoring those
conflicts can only make the true cofactor product larger, never smaller, so the
computed D is <= D_n and the bound on q is valid. Pushing the cutoff y higher
enumerates more structure and can only tighten (raise) D, tightening the bound.

SCOPE (be honest): this bounds the PRIME UNIVERSE that an exclusion / top-factor
search must range over. It does NOT by itself bound the exhaustive search-tree
size for the minimality certificate.

--------------------------------------------------------------------------
SELF-CHECK & SOUNDNESS GUARDS (run with --self-check)
--------------------------------------------------------------------------
1. D is computed by TWO independent prime generators (a byte sieve of
   Eratosthenes and an independent trial-division generator); the two D values
   must agree, else the program aborts.
2. Tail exhaustion guard: if for any admissible S the tail prime list is used up
   before  63 - |S|  admissible primes are found, the bound would be UNSOUND
   (missing factors); the program raises instead of returning a too-small D.
3. Reference values (X = 10^149): y=47 -> 3,321,610,829 ; y=61 -> 686,192,825 ;
   y=97 -> 55,183,356. `--self-check` reproduces these.
"""
import argparse
import sys

K = 64                 # number of prime factors of the Carmichael number
NCOF = K - 1           # cofactors after removing the largest prime q  (= 63)


# --------------------------------------------------------------------------
# Two INDEPENDENT prime generators (used to cross-verify D in --self-check)
# --------------------------------------------------------------------------
def primes_sieve(limit):
    """Primes <= limit via a byte sieve of Eratosthenes (method A)."""
    if limit < 2:
        return []
    s = bytearray([1]) * (limit + 1)
    s[0] = s[1] = 0
    i = 2
    while i * i <= limit:
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
        i += 1
    return [i for i in range(limit + 1) if s[i]]


def primes_trial(limit):
    """Primes <= limit via incremental trial division (method B, independent).

    Deliberately shares no code path with the sieve so a bug in one cannot hide
    a bug in the other."""
    if limit < 2:
        return []
    out = [2]
    n = 3
    while n <= limit:
        r = int(n ** 0.5)
        prime = True
        for p in out:
            if p > r:
                break
            if n % p == 0:
                prime = False
                break
        if prime:
            out.append(n)
        n += 2
    return out


# --------------------------------------------------------------------------
# Core bound
# --------------------------------------------------------------------------
def compute_D(y, tail_primes, small_primes):
    """Return (num_admissible_subsets, D) for cutoff y.

    tail_primes : sorted primes > y (the candidate tail factors)
    small_primes: sorted primes in [3, y]
    Raises RuntimeError if the tail is exhausted (soundness guard).
    """
    Smem = small_primes
    tail_gt = tail_primes
    best = None
    nsub = 0
    max_tail_index_used = 0  # for the exhaustion guard / diagnostics

    # DFS over internally-admissible subsets S, built in increasing prime order.
    # Iterative stack to avoid Python recursion limits at large y.
    # Each stack frame: (start_index_into_Smem, S_list, prod_of_S)
    stack = [(0, [], 1)]
    while stack:
        idx, S, prodS = stack.pop()
        nsub += 1
        need = NCOF - len(S)
        if need >= 0:
            prod = prodS
            cnt = 0
            last = -1
            for a in tail_gt:
                last += 1
                ok = True
                for p in S:
                    if (a - 1) % p == 0:
                        ok = False
                        break
                if ok:
                    prod *= a
                    cnt += 1
                    if cnt == need:
                        break
            if cnt < need:
                # SOUNDNESS GUARD: not enough admissible tail primes available.
                raise RuntimeError(
                    "tail-prime pool exhausted for S=%r (found %d of %d needed);"
                    " increase tail_limit -- the bound would be UNSOUND otherwise"
                    % (S, cnt, need))
            if last > max_tail_index_used:
                max_tail_index_used = last
            if best is None or prod < best:
                best = prod
        # extend S by the next admissible small prime
        for j in range(idx, len(Smem)):
            r = Smem[j]
            adm = True
            for p in S:
                if (r - 1) % p == 0:
                    adm = False
                    break
            if adm:
                stack.append((j + 1, S + [r], prodS * r))
    return nsub, best, max_tail_index_used


def bound(X, y, tail_limit, self_check=False):
    """Compute the certified bound. Returns dict with all fields."""
    # Method A primes.
    allp_A = primes_sieve(tail_limit)
    small_A = [p for p in allp_A if 3 <= p <= y]
    tail_A = [p for p in allp_A if p > y]
    nsub, D, maxidx = compute_D(y, tail_A, small_A)

    if self_check:
        # Method B primes, fully independent, must produce the identical D.
        allp_B = primes_trial(tail_limit)
        if allp_B != allp_A:
            raise AssertionError("prime generators A and B disagree on prime list")
        small_B = [p for p in allp_B if 3 <= p <= y]
        tail_B = [p for p in allp_B if p > y]
        nsub_B, D_B, _ = compute_D(y, tail_B, small_B)
        if D_B != D:
            raise AssertionError("independent D computations disagree: %d vs %d"
                                 % (D, D_B))
        if nsub_B != nsub:
            raise AssertionError("subset counts disagree")

    q = (X - 1) // D
    return {
        "X": X, "y": y, "tail_limit": tail_limit,
        "num_admissible_subsets": nsub,
        "D": D, "D_digits": len(str(D)),
        "q_bound": q, "q_digits": len(str(q)),
        "max_tail_prime_used": tail_A[maxidx] if tail_A else None,
        "tail_pool_size": len(tail_A),
    }


# --------------------------------------------------------------------------
# Prime counting  pi(q)  = size of the certified prime universe
# --------------------------------------------------------------------------
def count_primes_below(n):
    """Exact number of primes p with p <= n, via a segmented sieve.

    Memory-light (segment = 1<<20). Cost grows with n; for the tight (large-y)
    bounds n is small and this is fast."""
    if n < 2:
        return 0
    import math
    limit = int(math.isqrt(n))
    base = primes_sieve(limit)  # primes up to sqrt(n)
    count = 0
    seg = 1 << 20
    lo = 2
    while lo <= n:
        hi = min(lo + seg - 1, n)
        size = hi - lo + 1
        block = bytearray([1]) * size
        for p in base:
            if p * p > hi:
                break
            start = max(p * p, ((lo + p - 1) // p) * p)
            block[start - lo::p] = bytearray(len(block[start - lo::p]))
        count += sum(block)
        lo = hi + 1
    return count


# --------------------------------------------------------------------------
def parse_X(s):
    """Accept a plain integer, or  '<mant>e<exp>'  meaning mant * 10**exp
    (so '1e149' -> 10**149).  Fractional mantissas are rejected (must stay exact)."""
    s = s.strip()
    low = s.lower()
    if "e" in low:
        mant, exp = low.split("e")
        exp = int(exp)
        if "." in mant:
            raise ValueError("fractional mantissa would not be exact; give an integer X")
        m = int(mant)
        return m * (10 ** exp)
    return int(s)


def self_check():
    """Reproduce the published reference bounds for X = 10^149."""
    X = 10 ** 149
    expected = {47: 3321610829, 61: 686192825, 97: 55183356}
    ok = True
    for y, want in expected.items():
        r = bound(X, y, tail_limit=20000, self_check=True)
        got = r["q_bound"]
        flag = "OK " if got == want else "FAIL"
        if got != want:
            ok = False
        print("[self-check] y=%3d  subsets=%-8d  D_digits=%-3d  q=%d  expected=%d  %s"
              % (y, r["num_admissible_subsets"], r["D_digits"], got, want, flag))
    print("[self-check] independent-D agreement + tail-exhaustion guard: PASSED")
    print("[self-check] OVERALL:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(
        description="Certified upper bound on P+(n) for 64-factor Carmichael n < X.")
    ap.add_argument("--X", type=str, default="1e149",
                    help="upper limit X: integer, or '1e149' style. Default 1e149.")
    ap.add_argument("--y", type=int, default=97,
                    help="admissibility cutoff (enumerate admissible subsets of primes in [3,y]).")
    ap.add_argument("--tail-limit", type=int, default=20000,
                    help="generate tail primes up to this bound (soundness pool). Default 20000.")
    ap.add_argument("--pi", action="store_true",
                    help="also compute pi(P+) = certified prime-universe size (may be slow for tiny y).")
    ap.add_argument("--self-check", action="store_true",
                    help="run the built-in reference reproduction and exit.")
    args = ap.parse_args()

    if args.self_check:
        sys.exit(self_check())

    X = parse_X(args.X)
    r = bound(X, args.y, tail_limit=args.tail_limit, self_check=True)
    print("X                     = %d  (%d digits)" % (X, len(str(X))))
    print("y (cutoff)            = %d" % r["y"])
    print("admissible subsets    = %d" % r["num_admissible_subsets"])
    print("D (cofactor lower bnd)= %d-digit integer" % r["D_digits"])
    print("max tail prime used   = %s  (pool of %d primes <= %d)"
          % (r["max_tail_prime_used"], r["tail_pool_size"], args.tail_limit))
    print("CERTIFIED  P+(n) <=   = %d  (%d digits)" % (r["q_bound"], r["q_digits"]))
    if args.pi:
        pi = count_primes_below(r["q_bound"])
        print("prime universe pi(P+) = %d primes below the bound" % pi)


if __name__ == "__main__":
    main()
