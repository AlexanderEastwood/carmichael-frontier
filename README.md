# carmichael-frontier

An independent, correctness-first implementation of the search for **the smallest
Carmichael number with exactly _k_ prime factors** ([OEIS A006931](https://oeis.org/A006931)),
built to help push the frontier past the current record — and used here to find a **verified
64-factor Carmichael number**.

> **Result (updated 2026-09-11):** an explicit interval for S₆₄,
> $$10^{145} \le S_{64} \le N < 10^{148},$$
> where `N` is an oracle-verified **148-digit** Carmichael number with exactly 64 prime factors.
> So **S₆₄ has 146, 147, or 148 digits.** The lower bound is a finite exhaustion (127,092
> admissible products below 10¹⁴⁵, zero Carmichael) by three separate enumeration implementations,
> with the universe bound independently recomputed. This is a global bound on the minimum —
> **not** a determination of it. Full write-up:
> [`INTERVAL_THEOREM.md`](INTERVAL_THEOREM.md) and [`paper/interval_theorem.tex`](paper/interval_theorem.tex).

Everything here is our own code, developed beside — not committed to — Jonathan Webster's
[`small-carmichael-numbers`](https://github.com/jewebste/small-carmichael-numbers), whose
published table (k = 3–63, and a k = 65 candidate; k = 36–63 from Butler University) is the
authoritative reference for the known values and the neighbours cited below.

## The k = 64 interval

**Upper bound — a 148-digit incumbent** (improves the earlier 149-digit candidate), the product of
64 distinct primes (largest **3697**):

```
19·29·31·37·41·43·47·53·61·67·71·73·79·89·97·101·103·109·113·127·131·137·139·151·157·163·167·181·193·197·199·211·239·241·251·257·271·277·281·307·313·331·337·353·379·397·401·421·433·449·461·463·491·541·547·577·599·631·673·811·829·883·1951·3697
```

- Value + metadata in [`results_k64_best_global.json`](results_k64_best_global.json)
  (modulus `1768248177696000`). **Re-derive it by multiplying the factor list — don't trust a
  pasted decimal.** `ref/ref_carmichael.py::verify_certificate` confirms squarefree, 64 distinct
  primes (independent Miller–Rabin), Korselt `(p−1)|(N−1)` for every `p`; `10^147 ≤ N < 10^148`.
- Sits between Webster's neighbours `N₆₃` (145 digits) and the `N₆₅` candidate (151), where a true
  `S₆₄` must lie. Still an **upper bound**, not a certified minimum — the exchange search explored
  only radius `r ≤ 5` over part of the modulus portfolio.
- The original 149-digit candidate and its full write-up remain in
  [`results_k64.json`](results_k64.json) / [`RESULTS_k64.md`](RESULTS_k64.md).

**Lower bound — `S₆₄ ≥ 10¹⁴⁵`.** Korselt's pairwise condition forces the 63 smallest factors of any
64-factor Carmichael to multiply to at least a 142-digit constant `C`, so any such `n < 10¹⁴⁵` has
all prime factors `≤ 5518` (a universe of 727 odd primes). Exhausting every admissible 64-subset of
that universe below 10¹⁴⁵ gives 127,092 complete products and **zero** Carmichael numbers. See
[`INTERVAL_THEOREM.md`](INTERVAL_THEOREM.md) for the three independent reproductions and
[`paper/interval_theorem.tex`](paper/interval_theorem.tex) for the proofs.

## Method

**Small k (exact tabulation).** A full-wheel Korselt search over preproducts, with:
- an exact-_k_ lower-bound prune `R_min = ∏(t smallest admissible primes)` (reject a subtree once
  `P·R_min > B`), and its residue-class lift `R*_min` at the base case (no extra modular inverse);
- a **Fermat-on-R leaf filter** — for odd cofactor `R`, a completion can be Carmichael only if
  `2^(PR−1) ≡ 1 (mod R)` (small modulus `R`, exponent `PR−1`, **not** `R−1`); factor only survivors;
- a **shared-incumbent** parallelization (one integer + its factorization behind a mutex; the hot
  path never locks) that scales positively to all physical cores.

These reproduce A006931 for k = 3…29, each value oracle-verified.

**k = 64 (meet-in-the-middle).** Exact tabulation is infeasible at k = 64, so we use an
**exchange MITM**: represent a candidate as `S = (B₀ \ D) ∪ I` with `|D| = |I| = r`, seed a
neighbour-derived **RICH modulus portfolio** (perturbed exponents of the neighbours' λ plus a
Sophie–Germain add-pool), and for each modulus `M` match `∏` over the eligible pool `Q(M)` so that
`λ(S) = lcm(q−1) | M | (∏S − 1)` — Carmichael by Korselt, then re-verified by the oracle. The
gate prototype (`src/mitm.cpp`) recovers the true minimum at k = 29, 39, 61 given the modulus;
`src/mitm64.cpp` is the parallel engine used for the k = 64 run, validated by re-recovering N₆₃
exactly before the search.

## Layout

```
ref/ref_carmichael.py    FROZEN independent oracle (Korselt certificate checker + tiny-k brute force)
src/carmichael.cpp       full-wheel tabulation engine (R_min/R*_min prunes, Fermat-on-R leaf)
src/carmichael_par.cpp   shared-incumbent parallel engine
src/mitm.cpp             exchange-MITM gate prototype (recovers k=29/39/61 given the modulus)
src/mitm64.cpp           parallel exchange-MITM engine (the k=64 search)
src/mitm64_driver.py     modulus-portfolio driver
src/modsel.py            modulus-selection portfolio (RICH rule + Sophie-Germain add-pool)
src/README.md            engine build/usage notes, flags, and scaling
results_k64.json         the k=64 candidate: n, factors, modulus, provenance
RESULTS_k64.md           full k=64 write-up (method, stats, correctness caveats)
```

## Build & verify

The C++ engines need GMP (`libgmp`, `libgmpxx`); the Python oracle and drivers need only a
standard Python 3. See [`src/README.md`](src/README.md) for the exact build lines and flags.

Verify the k = 64 candidate independently:

```bash
python3 - <<'PY'
import sys; sys.path.insert(0, "ref")
import ref_carmichael as R
primes = [19,29,31,37,41,43,47,53,61,67,71,73,79,89,97,101,103,109,113,127,131,137,139,151,
157,163,167,181,193,197,199,211,239,241,257,271,277,281,307,313,331,337,379,397,409,421,443,
461,463,487,491,499,521,599,617,641,673,691,701,769,859,1151,1321,1933]
ok, n, why = R.verify_certificate(primes)
print(ok, why, f"({len(str(n))} digits)")
PY
```

## Correctness discipline

1. **A frozen reference oracle** (`ref/ref_carmichael.py`) that shares no code with the fast
   engines certifies every claimed factorization via Korselt and brute-forces the true
   least-with-k for tiny k from first principles. If an engine and the oracle ever disagree, the
   oracle wins and the run halts.
2. **Correctness before speed:** the engines reproduce the known A006931 table (including entries
   with factors 2/3/5) before any optimization is trusted.
3. Every number reported here was re-derived from an oracle-verified factor list, never copied as
   a decimal.

## Provenance & credit

Known values, the k = 65 candidate, and the neighbour factorizations are Jonathan Webster's /
Butler University's work, published at
[jewebste/small-carmichael-numbers](https://github.com/jewebste/small-carmichael-numbers) — cited,
not reproduced here. The implementation and analysis in this repository were carried out by an AI
system (Claude, Anthropic) under my direction, with a second model consulted for review; all
numerical claims are gated by the independent frozen oracle above and re-checked.

— Alexander Eastwood

## License

MIT — see [`LICENSE`](LICENSE).
