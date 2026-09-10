# src/ — M1 baseline + M2 R_min prune

Our own preproduct-tabulation search for the least Carmichael number with
**exactly k prime factors** (OEIS A006931). Correctness before speed.

`carmichael.cpp` is a single source. As of 2026-09-10 the **default build (no
flags) is M2.5F** — the fully-tightened production engine: R_min prune +
base-case R\*_min + the **rejection-first Fermat-on-R leaf filter** (see the
M2.5F section at the bottom). Every earlier engine stays reachable behind an
**opt-out** flag, so the A/B harnesses still run on byte-identical source:
`-DNO_FERMAT` = M2.5 (no leaf filter), `-DNO_PRUNE_STAR` = M2, `-DNO_PRUNE` = M1,
`-DPRUNE_STAR_REC` = M2.5rF. `carmichael_par.cpp` (M3) also carries base-case
R\*_min + the Fermat filter by default.

## Files
- `carmichael.cpp` — the engine (C++17 + GMP). Compile flag `-DPRUNE` = M2.
- `validate.py` — thin driver over the FROZEN oracle `../ref/ref_carmichael.py`
  (`verify` a certificate; `brute` = oracle's independent `least_with_k`). Shares
  no code with the engine.
- `check.sh` — regression harness: for each k, run engine, oracle-verify the
  certificate, compare to the OEIS b-file (`../data/b006931.txt`), and for k<=5
  compare to the oracle's independent brute.
- `ab.sh` — M1-vs-M2 A/B on an identical explicit bound `v+1` per k (single-pass,
  clean counters): prints N_E/N_I/N_F/N_H, wall, oracle verdict, and the N_E/N_I
  reduction factors. `bash ab.sh 12 16 20 23 25 26`.

## Build
```
g++ -O2 -std=c++17                  carmichael.cpp -o carmichael          -lgmp -lgmpxx  # M2.5F (DEFAULT)
g++ -O2 -std=c++17 -DNO_FERMAT      carmichael.cpp -o carmichael_star     -lgmp -lgmpxx  # M2.5 (no leaf filter)
g++ -O2 -std=c++17 -DNO_PRUNE_STAR  carmichael.cpp -o carmichael_m2       -lgmp -lgmpxx  # M2
g++ -O2 -std=c++17 -DNO_PRUNE       carmichael.cpp -o carmichael_m1       -lgmp -lgmpxx  # M1
g++ -O2 -std=c++17 -DPRUNE_STAR_REC carmichael.cpp -o carmichael_star_rec -lgmp -lgmpxx  # M2.5rF
```
NOTE: the default binary is now M2.5F, not M1. `ab.sh` still compares the
pre-built `carmichael` (now M2.5F) vs `carmichael_m2` (M2); rebuild those names
per the lines above if you want a specific A/B pairing (e.g. `carmichael_star`
= M2.5 vs `carmichael` = M2.5F to isolate the Fermat filter).

## M2 — the R_min exact-k prune (correctness argument)
At a node with preproduct `P`, largest prime `pmax`, and `t` factors still to
place, let `a_1<…<a_t` be the `t` smallest primes each **individually admissible**
to `P` (prime, `> pmax`, and no prime dividing `P` divides `a_i − 1`). Set
`R_min = ∏ a_i` and **skip the subtree when `P·R_min > bound`.**

`R_min` is a genuine LOWER BOUND on the cofactor `R` of any completion, because
in ANY Carmichael number no prime factor `p` divides `(q−1)` of another prime
factor `q`: if it did, `p | (q−1) | λ(n) | (n−1)` while `p | n`, forcing
`p | gcd(n,n−1)=1`. So the `t` true remaining factors are `t` distinct primes
that are each `> pmax` AND admissible to `P`; their product is therefore `≥` the
product of the `t` smallest such primes `= R_min`. (Those `t` smallest admissible
primes may be mutually incompatible — one dividing another's minus-one — so
`R_min` itself need not be realizable; that only makes the true product larger,
so the bound stays valid.) The prune only skips subtrees whose every completion
already exceeds the current bound/incumbent, so **no answer can change**; it is
applied in the recursive step (per child, before the descent is counted) and at
the top of each base case (before the modular inverse and the AP scan).
Deps: GMP dev (`libgmp`, `libgmpxx`).
Uses `__rdtsc` (x86-64) for the phase timers.

## Run
```
./carmichael k [bound]      # bound omitted or 0 => bound-doubling (no table used)
```
Prints `n`, `FACTORS: p1 p2 ...` (pipe to the oracle), the four counters, wall
time, and the enumeration/first-rejection/deep time split.

Regression sweep + oracle checks:
```
bash check.sh 3 22 5        # klo khi brute_max
```

## Algorithm (what it does, and the one deliberate difference)
Recursive tree over **preproducts** P = a squarefree product of distinct primes
with `gcd(P, lambda(P)) = 1` (lambda = Carmichael function = lcm of (p-1)).

- **Recursive case** (`t` = factors still to place >= 2 and `P*lambda(P) <= bound`):
  for each admissible prime `q` in `(pmax, (bound/P)^(1/t)]`, recurse on `Pq`.
  `q` is *admissible* iff no prime already in P divides `q-1` (this preserves
  `gcd(P, lambda(P)) = 1`, so P is invertible mod lambda(P)).
- **Base case** (`P*lambda(P) > bound`, or `t <= 1`): every completion has
  `R = P^{-1} (mod lambda(P)) + j*lambda(P)` with `R <= bound/P` (because
  `lambda(P) | (n-1)` is necessary). For each such `R`: cheap screen (does R have
  a prime factor <= pmax? — reject fast), then full factorization of R into
  exactly `t` distinct primes all `> pmax`, then a definitive Korselt check on
  n = P*R. Factors are required `> pmax` so each Carmichael number is generated
  exactly once (at the unique prefix where the base case fires).

**The deliberate difference from Jon's tool:** a FULL wheel. We admit every odd
prime (3, 5, 7, ...), so we return the TRUE minima. (Jon's mod-30 wheel skips
2/3/5 and returns 1729 for k=3; ours returns **561**.) 2 is excluded by proof,
not by a wheel shortcut: a Carmichael number is odd — if 2|n then n is even with
an odd prime factor p, and Korselt needs (p-1)|(n-1) with p-1 even, n-1 odd:
impossible.

## Design choices
- **All big arithmetic is GMP `mpz_class`.** No uint64/uint128 fast paths, no
  pair-step, no mod-30 wheel. Those are speed work (M2+); M1 is about an
  obviously-correct reference we trust.
- **Bound-doubling** when no bound is given (start 1000, x4) so the true global
  least is found with NO reliance on any published table; the incumbent shrinks
  the working bound within each pass. (Cost: the final pass may overshoot up to
  4x, and earlier passes re-walk — visible as extra `enum` time.)
- **`factor_R` uses trial division** by the prime sieve (up to sqrt R) plus a
  Miller-Rabin cofactor test. Simple and obviously correct. NOTE: this is a
  different cost model from Jon's polynomial `cn_query` (repeated Fermat-witness
  modexp+gcd); the `deep` time fraction reflects OUR factorizer, not his.
- **Sieve** to 5e6 (odd primes). A hard guard aborts if a needed prime interval
  ever exceeds the sieve (never triggered for k<=25).

## Instrumentation (the four counters)
- `N_E` — extensions examined (recursive descents into a child preproduct).
- `N_I` — terminal inversions (base-case entries; one modular inverse each).
- `N_F` — first rejection tests (AP candidates R that get the cheap screen).
- `N_H` — candidates needing full factorization + Korselt certification.

Time split via `__rdtsc` accumulators (calibrated against `steady_clock`):
`t_deep` = factor+Korselt, `t_first` = cheap screen, `t_enum` = wall - the two.

---

# M2.5 — the R*_min residue-class tightening (DONE 2026-09-10)

Astra's tightening of the base-case R_min bound. At a base-case node the cofactor
`R` must satisfy **both** `R >= R_min` **and** `R ≡ a (mod lambda)` where
`a = P^{-1} mod lambda`. The least feasible cofactor is therefore
`R*_min = R_min + ((a - R_min) mod lambda)` — the least value `>= R_min` in the
forced residue class — and every completion has `R >= R*_min`, so it is a strictly
stronger (still sound) lower bound than `R_min`. It cannot drop the true minimum:
any real completion's `R` is already `>= R*_min`.

Compile flags (all imply the ones below; `PRUNE_STAR` implies `PRUNE`):
```
g++ -O2 -std=c++17 -DPRUNE_STAR     carmichael.cpp -o carmichael_star     -lgmp -lgmpxx  # M2.5
g++ -O2 -std=c++17 -DPRUNE_STAR_REC carmichael.cpp -o carmichael_star_rec -lgmp -lgmpxx  # M2.5r
```
- **M2.5 (`-DPRUNE_STAR`) — base case only.** Reuses the modular inverse the base
  case already computes, so it adds **zero** extra inverses over M2. It (i) starts
  the AP scan at `R*_min`, skipping every candidate below `R_min` that M2 still
  screened, and (ii) tightens the whole-node prune. `N_E`/`N_I` are **unchanged**
  vs M2 (same tree, same base-case entries); `N_F`/`N_H` drop ~3x.
- **M2.5r (`-DPRUNE_STAR_REC`) — also at recursive children.** The residue
  constraint holds at every node, so `R*_min` is a valid stronger child prune too,
  but a recursive node has no inverse in hand: this pays **one extra GMP modular
  inverse per surviving child**. Measured, not assumed.

## A/B (explicit bound V+1, single thread; correctness: all == b006931.txt, oracle OK)
```
k   variant  N_E        N_I       N_F      N_H     wall(s)
20  M2       89511      58121     5385     2784    0.081
20  M2.5     89511      58121     1820      935    0.046
20  M2.5r    30437       1820     1820      935    0.044
25  M2       2458080    1669081   146458   73222   6.45
25  M2.5     2458080    1669081   48977    24456   3.43
25  M2.5r    758680     48977     48977    24456   3.33
27  M2       8441721    5748405   488817   242940  41.5
27  M2.5     8441721    5748405   162655   80786   19.0
27  M2.5r    2579977    162655    162655   80786   18.6
29  M2       15063713   9966462   868596   433627  142.3
29  M2.5     15063713   9966462   284933   142270  60.9
29  M2.5r    4871585    284933    284933   142270  59.3
```
**Verdict: adopt M2.5 (base-case R*_min) always.** ~1.75x–2.34x wall speedup over
M2 (grows with k), at zero extra inverse cost, `N_F`/`N_H` cut ~3x. **Do NOT adopt
M2.5r:** despite huge counter reductions (`N_I` down ~35x, `N_E` ~3x) its wall is
only ~2-3% better than M2.5 — the wall is dominated by base-case deep factor/Korselt
work, which base-case `R*_min` already trims 3x, and the extra recursive inverses
almost exactly offset the residual savings. Correctness re-verified: k=3..29 identical
values, all oracle-verified and == b006931.txt, for both M2.5 and M2.5r. Every M1/M2
code path is byte-for-byte unchanged (all new logic under `#ifdef PRUNE_STAR*`).

---

# M3 — single-box parallel with a SHARED GLOBAL INCUMBENT (DONE)

`carmichael_par.cpp` builds `carmichael_par`: the M2 search math verbatim (sieve,
admissibility, R_min prune, factor_R, Korselt) plus (1) coarse prefix-subtree job
partitioning and (2) a shared, cheaply-read incumbent bound. It achieves POSITIVE
scaling — the failure mode of Jon's `dws`.

Build / run:
    g++ -O2 -std=c++17 carmichael_par.cpp -o carmichael_par -lgmp -lgmpxx -pthread
    ./carmichael_par k [bound] [workers] [jobs_per_worker]
    # bound 0/omitted => bound-doubling (no table). workers default = nproc(=32).
    # jobs_per_worker (default 64) tunes partition granularity; 1024 used below.

## Why dws scales negatively, and what M3 does instead
dws workers deep-copy a GMP big-int preproduct + a bitvector under a mutex on every
work-steal AND do not share a tightening incumbent, so they redo speculative work a
global bound would prune. More workers => more redundant deep work + mutex traffic =>
slower.

M3:
1. PARTITION — enumerate preproduct prefixes down to a split depth d: a frontier cut
   of the search tree = many independent subtree jobs. Running rec() on every frontier
   node covers the single-proc tree exactly once, so the parallel minimum is provably
   the M2 minimum. Granularity is adaptive (deepen d until >= jobs_per_worker*workers
   jobs); coarse jobs + a dynamic queue (one atomic job index), never fine-grained
   stealing. No GMP big-int or bitvector is ever deep-copied across workers.
2. SHARED INCUMBENT — essentially ONE integer (+ its winning factorization) behind a
   mutex. The hot path NEVER locks: each worker holds a thread-local copy of the bound
   and a version stamp, and per node does ONE relaxed atomic load of a global version
   counter; it re-copies the mpz (under the mutex) only on the rare event that some
   worker lowered the bound. A worker writes only when it finds a strictly-smaller n.
   A stale local bound is always >= the true shared bound, so it only ever prunes LESS
   — it can NEVER remove the true minimum. Correct for any worker count / any timing;
   the reported witness is deterministic (the true minimum is unique).

## Correctness (parallel result == M2 result for every k tested)
- k=3..22 bound-doubling (NO table): every value oracle-verified
  (ref_carmichael.verify_certificate) and == b006931.txt, identical at 4 and 32 workers.
- k=3,4,5,10,15,20,23,24,26 bound-doubling on the production binary: oracle-verified +
  b-file MATCH (k<=5 also == oracle brute).
- k=25,27..37 explicit bound V+1 (ab.sh methodology): all == b006931.txt.
- k=25,28,29 identical result at nw=1,2,4,8,16,32 (deterministic across worker counts).

## Scaling (AMD Ryzen 9 9950X, 16 physical cores / 32 SMT threads; explicit
## bound V+1; jobs_per_worker=1024). eff = speedup / workers.
    k=25:  1w 6.53s(1.00x)  2w 1.95x  4w 3.72x  8w 6.84x  16w 9.97x  32w 8.26x
    k=27:  1w 42.5s(1.00x)  2w 1.96x  4w 3.84x  8w 7.30x  16w 11.77x 32w 12.42x
    k=28:  1w 46.5s(1.00x)  2w 1.94x  4w 3.76x  8w 7.40x  16w 12.08x 32w 12.28x
    k=29:  1w 140.1s(1.00x) 2w 1.98x  4w 3.82x  8w 7.30x  16w 12.27x 32w 12.97x
Near-linear to 8 cores (>90% eff), ~12x on 16 physical cores (74-77% eff). Threads
17-32 are SMT siblings and add little on compute-bound GMP (the 16->32 plateau), plus
all-core turbo throttles vs the 1-thread baseline. This is unambiguous POSITIVE scaling
— the exact opposite of dws. k=25 is too small (6.5s total) to feed 32 threads and
regresses past 16; the larger k do not.

## New practical-k ceiling (nw=32, explicit bound V+1, all == b006931.txt)
M2 single-proc ceiling was k=29 in ~142s. M3 at 32 workers, all UNDER 3 min:
    k30 11.2s  k31 0.6s  k32 33.5s  k33 56.1s  k34 35.8s  k35 130.3s  k36 8.3s  k37 79.8s
Wall is highly non-monotonic in k (each minimum's structure + how fast the tight bound
prunes it). New sub-3-min ceiling >= k=37 (from k=29).

## Notes for later milestones
- Load-balance floor: with per-nw granularity the single heaviest job (maxjob_s, printed)
  no longer gates wall; the residual gap at 16 cores is the tail + turbo throttle, at 32
  it is SMT. enum_s (serial job-build) stayed <0.4s — not a bottleneck. jemalloc (LD_PRELOAD)
  made no measurable difference: allocator contention is NOT the limiter.
- M4 (GPU leaf kernel): the wall is ~90%+ base-case work (factor_R + Korselt over the AP
  scan). That is the batched, fixed-width leaf the CGBN-style kernel should target.
- Multi-node: the shared incumbent is one integer + factorization — trivially
  shippable between nodes (broadcast on improvement); jobs are already independent coarse
  subtrees, so a node-level dynamic queue over the same Job list extends this directly.

---

# M2.5F — the rejection-first Fermat-on-R leaf filter (DONE 2026-09-10) + default flip

Profiling M2.5 showed the base-case leaf spent ~63–84% of wall (rising with k) in
`factor_R`'s exhaustive trial division (dividing each surviving cofactor R by the
~5M-prime sieve up to √R) + the Miller–Rabin cofactor test — NOT in modexp. The
fix (expert consult) is a **rejection-first leaf**: filter cheaply *before*
factoring.

## The filter (correctness)
At a base-case candidate `n = P·R` (R odd, R>1) that has passed the cheap
small-factor screen, compute
```
    2^(n-1) mod R      (exponent n-1 = P·R − 1, modulus R; one mpz_powm)
```
and **REJECT if it is ≠ 1**. This is sound because Korselt ⇒ if `n` is Carmichael
then `2^(n-1) ≡ 1 (mod n)`, hence `≡ 1 (mod R)` for every `R | n`. So the test is
a **necessary condition only**: it may REJECT but must **never ACCEPT**. Every
survivor still runs the *full* exact `factor_R` + `korselt` check — a Fermat pass
never short-circuits acceptance (verified in code: the accept path is unchanged).
Valid for **all** `t = k − ω(P)`: when `t > 1` the cofactor `R` is *expected
composite*, so we deliberately do NOT use a compositeness/BPSW test (that would be
the wrong test); Fermat-on-R with the full `n−1` exponent is correct for every t.
Exponent must be `n−1`, NOT `R−1` (the latter is a primality test on R — wrong).

## Default flip
The default `carmichael` binary is now **M2.5F** (all prior engines behind opt-out
flags; see Build). `carmichael_par.cpp` (M3) carries base-case R\*_min **and** the
Fermat filter by default too. All correctness re-verified: **k=3..29 identical
values, every one oracle-verified (`verify_certificate`) and == `b006931.txt`**,
single-thread and parallel (par also deterministic at nw=1,4,8,16,32). New counter
`N_G` = candidates passing the Fermat filter (get the full factor+Korselt);
`survivor_frac = N_G / N_H`.

## Op-mix + speedup (single thread, explicit bound V+1; M2.5 = `carmichael_star`)
```
k   M2.5 wall  M2.5F wall  speedup | N_H       N_G  survivor s   M2.5 t_deep   M2.5F t_fermat  M2.5F t_deep
25  3.42s      1.25s       2.73x   | 24456     22   9.0e-4       2.16s (63%)   0.012s         5.6e-5s
27  19.28s     4.71s       4.09x   | 80786     35   4.3e-4       14.51s (75%)  0.044s         1.2e-4s
29  60.66s     9.53s       6.37x   | 142270    41   2.9e-4       50.96s (84%)  0.083s         3.9e-4s
```
The filter **eliminates the trial-division leaf**: at k=29 it replaces ~51s of
`factor_R`+Korselt (over 142270 candidates) with 0.083s of Fermat modexp (over the
same 142270) + 0.0004s of full factoring on the **41 survivors**. Survivor fraction
`s ≈ 0.03%` and falls with k. Post-filter the wall is **~99% enumeration**
(`frac_enum` 0.987→0.989; the leaf is now `t_fermat+t_deep < 1%`) — i.e. the tree
walk + the base-case modular inverses (`N_I`, ~10M at k=29), NOT the leaf.

## New single-thread ceiling (bound V+1, all == b006931.txt, oracle OK)
M2.5 single-thread ceiling was k=29 (~61s). M2.5F single-thread, all under 3 min:
```
    k30 7.8s  k31 0.17s  k32 15.9s  k33 23.5s  k34 14.1s  k35 44.8s  k36 2.6s  k37 25.0s
```
Highly non-monotonic (minimum structure + how fast the tight bound prunes). k=38
exceeds 185s single-thread. **New reliable single-thread sub-3-min ceiling: k=37**
(was k=29) — an ~8-k jump on ONE core, before any parallelism.

## Implication for the k=64 plan (IMPORTANT — changes M4)
The leaf that motivated the **M4 GPU Fermat kernel is now cheap on the CPU**: the
Fermat filter is <1% of wall and `factor_R` on survivors is negligible. A batched
fixed-width Fermat GPU kernel would accelerate the already-tiny part → the same
mistake the assessment warned against for a naïve tree-port. The new single-thread
bottleneck is **enumeration + the ~10M base-case modular inverses (`N_I`)**, i.e.
`t_enum`. Next real wins are therefore (1) cutting `N_I`/enumeration cost
(cheaper/fewer inverses, tighter node pruning — note M2.5r already cuts `N_I` ~35×
but its extra recursive inverses eat the gain; worth revisiting now that the leaf
is free), and (2) the M3 CPU fleet (par already POSITIVE-scaling) + M5
meet-in-the-middle. Re-profile GPU value against `t_enum`, not the leaf.
