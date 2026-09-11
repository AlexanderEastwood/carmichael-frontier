# An interval for S₆₄ (update, 2026-09-11)

**Result.** We now bound the smallest 64-factor Carmichael number *from both sides*:

$$10^{145} \;\le\; S_{64} \;\le\; N \;<\; 10^{148},$$

so **S₆₄ has 146, 147, or 148 decimal digits** (triple-verified core), and the remaining exclusion
interval is `[10^145, N)`. Here `N` is an explicit, oracle-verified 148-digit Carmichael number with
exactly 64 prime factors. A search-engine strengthening (below) excludes everything `≤ 10^146`, so
in fact **S₆₄ has 147 or 148 digits** — see the verification tiering in
[Strengthening](#strengthening--no-64-factor-carmichael--10¹⁴⁶-2026-09-11). This is a global statement about `S₆₄` — *not* a determination of it (that would
require excluding every 64-factor Carmichael in `[10^145, N)`).

Full write-up with proofs: [`paper/interval_theorem.tex`](paper/interval_theorem.tex) (5-page PDF).
Everything below is at repository revision `f872e32`.

---

## Upper bound — a 148-digit incumbent (improves the earlier 149-digit candidate)

`N` is the product of these 64 distinct primes (largest is **3697**):

```
19·29·31·37·41·43·47·53·61·67·71·73·79·89·97·101·103·109·113·127·131·137·139·151·157·163·167·181·193·197·199·211·239·241·251·257·271·277·281·307·313·331·337·353·379·397·401·421·433·449·461·463·491·541·547·577·599·631·673·811·829·883·1951·3697
```

- **Re-derive `N` by multiplying the factor list** — don't trust a pasted decimal. Value + metadata
  in [`results_k64_best_global.json`](results_k64_best_global.json).
- **Verification.** Every factor is `≤ 3697`, so its primality is certified by exact trial division
  (deterministic at this size). Squarefree, 64 distinct primes, and Korselt `(p−1)|(N−1)` for every
  `p` — equivalently `λ(N) = lcm(p−1) = 589416059232000` and `N ≡ 1 (mod λ(N))`. `10^147 ≤ N < 10^148`.
  Checked by `ref/ref_carmichael.py::verify_certificate` (shares no code with the search).
- **Found** by exchange meet-in-the-middle at modulus `M = 3·λ(N) = 1768248177696000` (so `N ≡ 1
  mod M`; `M` is the *search modulus*, not `λ(N)`). An **upper bound**, not a certified minimum.
- Consistent with the neighbours *by digit length + the proved interval*, not by any general
  monotonicity of `Sₖ`: Webster's `N₆₃` has 145 digits (so `N₆₃ < 10^145 ≤ S₆₄`) and the `N₆₅`
  candidate has 151 digits (so it exceeds `N`). The original 149-digit candidate is in
  [`results_k64.json`](results_k64.json) / [`RESULTS_k64.md`](RESULTS_k64.md).

## Lower bound — no 64-factor Carmichael below 10¹⁴⁵ (finite exhaustion)

**1. Certified prime universe.** For every 64-factor Carmichael `n` and every prime factor `q | n`,
the remaining 63 factors satisfy `n/q ≥ C`, a **142-digit constant** `C = 1812…3549`. The bound is
`C = min_A (∏_{p∈A} p)·∏ᵢ bᵢ(A)`, where `A` ranges over internally admissible subsets of the odd
primes ≤ 97 and the `bᵢ(A)` are the increasing primes > 97 individually compatible with `A`; the
actual cofactor tail `rᵢ ≥ bᵢ(A)`, so this is a genuine lower bound. Three points: the `bᵢ(A)` may
conflict with *each other* (ignoring that only lowers the product); the finite prime table (odd
primes ≤ 10⁴) is checked to hold enough eligible primes for every class; the cofactor need not
itself be Carmichael. `C` is a *certified* lower bound, not necessarily the smallest possible
admissible 63-factor product. Hence any `n < 10^145` has largest prime factor
`P⁺(n) ≤ ⌊(10^145−1)/C⌋ = 5518` — a universe of **727 odd primes** (largest 5507).
Tool: `carmichael_prime_bound.py`.

**2. Exhaustion.** Enumerate every increasing pairwise-admissible 64-subset of those 727 primes with
product `< 10^145`, pruning only on necessary conditions (too few primes remain; or current product
× the `64−h` smallest still-available primes individually compatible with the `h` already chosen
reaches `10^145`). Pairwise admissibility (if `p | q−1`, drop `q` after choosing `p`) is enforced
during descent and discards no genuine Carmichael. There are exactly **127,092** such products
(distinct, by unique factorization), and **none is Carmichael.** Therefore `S₆₄ ≥ 10^145`.

> The count means precisely: increasing, pairwise-admissible 64-subsets of the certified 727-prime
> universe with product below 10¹⁴⁵. Keep that qualification — it is *not* "all products of 64
> distinct primes below 10¹⁴⁵."

### Reproduction — three enumeration implementations agree; the cofactor bound independently recomputed

| Implementation | Nodes | Korselt verdict from | Result |
|---|---|---|---|
| `s64_lower_bound_certificate.py`, binary enumerator | 327,049 | incremental `lcm(pᵢ−1)` | 127,092 products, 0 Carmichael |
| same file, ordered enumerator | 163,525 | per-prime `(n−1)%(pᵢ−1)` | 127,092 products, 0 Carmichael |
| `s64_lower_bound_oracle_check.py`, clean-room, **oracle-gated** | 163,525 | **frozen oracle** `verify_certificate` | 127,092 products, 0 Carmichael |

- The two enumerators in the certificate use **different prime generators, traversal logic, and
  Korselt tests**, while **sharing** the cofactor-bound computation and the results/hash reporting
  code. Their identical order-sensitive leaf digest
  `sha256 = e05f62598eea784b3f60f025c6ccbf30700115ecbef42702e24154dcece32321` is a **reproducibility
  aid**, not an extra independent correctness check (the hashing is shared). Correctness comes from
  the algorithmic coverage argument.
- The third enumerator is separate clean-room code that defers every leaf verdict to the frozen
  oracle; it was validated first against the oracle's brute-force `least_with_k`: `k=3 → 561`,
  `k=4 → 41041`, `k=5 → 825265`.
- Separately, **`carmichael_prime_bound.py`** independently recomputes `C` and `P⁺ ≤ 5518`. It does
  *not* enumerate products or establish the absence of Carmichael numbers — it certifies only the
  universe.

Commands: `python3 s64_lower_bound_certificate.py` · `python3 carmichael_prime_bound.py --X 1e145 --pi`
· `python3 s64_lower_bound_oracle_check.py`.

## Strengthening — no 64-factor Carmichael ≤ 10¹⁴⁶ (2026-09-11)

A complete exact-k preproduct search (`experiments/dynt/`, the M2.5F production engine) over all
64-factor products `≤ 10¹⁴⁶` returns **`n = NONE`**: no such Carmichael number exists. Hence

$$10^{146} \;<\; S_{64} \;\le\; N, \qquad\text{so } S_{64}\text{ has } 147 \text{ or } 148 \text{ digits.}$$

This is *far* cheaper than the ≥38.9-billion admissible-product **count** below 10¹⁴⁶, because the
engine rejects prefixes (via the exact `R_min` bound and the `P·λ(P)>B` / AP switch) without
materializing each product — the full search ran in ~71 s single-threaded.

**Verification status (honest tiering).** The `S₆₄ ≥ 10¹⁴⁵` bound above is the rigorously
triple-verified core (two certificate enumerators + a frozen-oracle-gated enumerator, plus the
independently recomputed universe). The `> 10¹⁴⁶` strengthening currently rests on the **search
engine**, corroborated by: (i) two independent switch strategies (production `P·λ>B` and the
dynamic-t early switch) both returning `NONE`; and (ii) the same engine reproducing `NONE` at 10¹⁴⁵,
where the result is independently certified — i.e. the engine is a *fourth* independent confirmation
of the 10¹⁴⁵ bound before it is trusted at 10¹⁴⁶. A fully independent 10¹⁴⁶ certificate (an
admissibility enumerator at that bound, oracle-gated) is the outstanding follow-up. The engine is
also validated end-to-end against OEIS A006931: it reproduces `S_10 … S_29` exactly, oracle-confirmed.

## What's next

- **Independent 10¹⁴⁶ certificate** to raise `> 10¹⁴⁶` to the triple-verified tier.
- Pushing toward 10¹⁴⁷ (each extra digit is ≈10⁵× more work) and, ultimately, the full `[…,N]`
  certification that would prove `S₆₄ = N` — a machine-scale endgame (≈ thread-years), not a rerun.

## Companion result: S₆₅ ≥ 10¹⁴⁸ (≥ 149 digits)

The same cofactor minimum transfers one factor count up, with **no new exhaustion**. Let
`A₆₄ = 4411…943377` (145 digits) be the least product of 64 pairwise-admissible odd primes — it is
exactly the **minimum of the 127,092 admissible 64-products** recorded above (same `e05f62…` digest).
For any 65-factor Carmichael `n` and any prime `q | n`, the other 64 factors are a pairwise-admissible
64-product, so `n/q ≥ A₆₄`, giving `P⁺(n) ≤ ⌊(n−1)/A₆₄⌋`. Below `10¹⁴⁸` that caps the universe at
**334 odd primes** (P⁺ ≤ 2266). An exhaustive search finds **no** 65-factor Carmichael there:

$$S_{65} \ge 10^{148}, \qquad\text{so the smallest 65-factor Carmichael number has} \ge 149 \text{ digits.}$$

**Reproduced four ways** (all 0 Carmichael): the certificate [`s65/s65_exclusion.py`](s65/s65_exclusion.py)
(64,355 nodes, 66 progression candidates); a separate native binary enumeration
[`s65/verify65.cpp`](s65/verify65.cpp) checking all **22,287,138** admissible 65-products below 10¹⁴⁸;
and our preproduct engine `cm_dynt 65` in both switch modes — which reports `n=NONE` with **matching**
counts (5012 inversions, 66 progression candidates). A companion
[`s65/s65_transfer_certificate.py`](s65/s65_transfer_certificate.py) proves the weaker
`S₆₅ ≥ 3.17×10¹⁴⁷` two more independent ways (285 admissible 65-products, 0 Carmichael). All rest on
the same triple-verified `A₆₄`. The universe cap ladder: P⁺ ≤ 2266 / 22666 / 226668 / **2,266,687**
for a 65-factor Carmichael below 10¹⁴⁸ / 10¹⁴⁹ / 10¹⁵⁰ / 10¹⁵¹.

**Upper bound.** Webster reports a 151-digit k=65 candidate; it is **not independently verified here**.
Confirming it would pin `S₆₅` to **149–151 digits**. (`S₆₅` does not depend on determining `S₆₄`.)

## Provenance / credit

Known values `3 ≤ k ≤ 63` (as of Sept 2026): R. G. E. Pinch (earlier `k`) and the Butler University
computation (`k = 36–63`), collected in Webster's
[`small-carmichael-numbers`](https://github.com/jewebste/small-carmichael-numbers) — the repository
cited for the known values and the `N₆₃`/`N₆₅` neighbours. Two existence results of different scope,
not to be merged: **Alford–Grantham–Hayman–Shallue** constructed examples for every
`3 ≤ k ≤ 19,565,220`; **Larsen–Wright** proved existence for every sufficiently large factor count.
Neither gives a minimum; our contribution is the explicit *interval* for `S₆₄`.

AI assistance (OpenAI **GPT-6**, the model designation recorded for the session) contributed the
lower-bound argument, generated the initial certificate code (`s64_lower_bound_certificate.py`), and
ran the first computation. Independently of that script, the author recomputed the cofactor bound
(`carmichael_prime_bound.py`) and wrote the separate oracle-gated enumerator
(`s64_lower_bound_oracle_check.py`); rerunning the supplied certificate is not counted as a
separately implemented proof. The author retains responsibility for the mathematics and code, and
treats the model as an assistant, not a mathematical authority.
