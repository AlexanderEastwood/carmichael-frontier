# An interval for S₆₄ (update, 2026-09-11)

**Result.** We now bound the smallest 64-factor Carmichael number *from both sides*:

$$10^{145} \;\le\; S_{64} \;\le\; N \;<\; 10^{148},$$

so **S₆₄ has 146, 147, or 148 decimal digits.** Here `N` is an explicit, oracle-verified
148-digit Carmichael number with exactly 64 prime factors. This is a global statement about
`S₆₄` — *not* a determination of it (that would require excluding every 64-factor Carmichael in
the open interval `(10^145, N)`).

The full write-up with proofs is [`paper/interval_theorem.tex`](paper/interval_theorem.tex)
(compiles to a 4-page PDF).

---

## Upper bound — a 148-digit incumbent (improves the earlier 149-digit candidate)

`N` is the product of these 64 distinct primes (largest is **3697**):

```
19·29·31·37·41·43·47·53·61·67·71·73·79·89·97·101·103·109·113·127·131·137·139·151·157·163·167·181·193·197·199·211·239·241·251·257·271·277·281·307·313·331·337·353·379·397·401·421·433·449·461·463·491·541·547·577·599·631·673·811·829·883·1951·3697
```

- **Re-derive `N` by multiplying the factor list** — don't trust a pasted decimal. The value and
  metadata are in `results_k64_best_global.json` (modulus `1768248177696000`).
- **Independently verified** by `ref/ref_carmichael.py::verify_certificate`: squarefree, 64
  distinct primes (independent Miller–Rabin), Korselt `(p−1)|(N−1)` for every `p`. `10^147 ≤ N < 10^148`.
- Sits between Webster's known neighbours `N₆₃` (145 digits) and the `N₆₅` candidate (151 digits),
  where a true `S₆₄` must lie. Still an **upper bound**, not a certified minimum.

## Lower bound — no 64-factor Carmichael below 10¹⁴⁵ (finite exhaustion)

Two steps, both exact-integer:

1. **Certified prime universe.** Korselt's pairwise condition (if `p | q−1` then `p`, `q` can't both
   divide a Carmichael `n`) forces the 63 smallest ("cofactor") primes of any 64-factor Carmichael
   to multiply to at least a **142-digit constant** `C = 1812…3549`, independent of any modulus.
   Hence any such `n < 10^145` has largest prime factor `P⁺(n) ≤ ⌊(10^145−1)/C⌋ = 5518` — a universe
   of just **727 odd primes** (largest 5507). Tool: `carmichael_prime_bound.py`.
2. **Exhaustion.** Enumerate every admissible increasing 64-subset of those 727 primes with product
   `< 10^145` (pruning only on necessary conditions; pairwise admissibility enforced during descent,
   which discards no genuine Carmichael). There are exactly **127,092** such complete products, and
   **none is Carmichael.** Therefore `S₆₄ ≥ 10^145`.

### Reproduced three independent ways (all agree: 127,092 products, 0 Carmichael)

| # | Implementation | Nodes | Korselt verdict from | Result |
|---|---|---|---|---|
| 1 | `s64_lower_bound_certificate.py`, binary enumerator | 327,049 | incremental `lcm(pᵢ−1)` | 0 hits |
| 1 | same file, ordered independent enumerator | 163,525 | per-prime `(n−1)%(pᵢ−1)` | 0 hits |
| 2 | `carmichael_prime_bound.py` | — | (recomputes `C`, universe) | `P⁺ ≤ 5518` |
| 3 | our enumerator, oracle-gated | 163,525 | **frozen oracle** `verify_certificate` | 0 hits |

- Enumerators (1) share no code path — different prime generators, traversals, and Korselt tests —
  yet agree on the count, the minimum product, and an **identical order-sensitive leaf digest**
  `sha256 = e05f62598eea784b3f60f025c6ccbf30700115ecbef42702e24154dcece32321`.
- Enumerator (3) takes every one of the 127,092 leaf verdicts from the frozen reference oracle, not
  from its own arithmetic, and was validated first against the oracle's brute-force `least_with_k`:
  `k=3 → 561`, `k=4 → 41041`, `k=5 → 825265` (all exact).

## What's next

Running the same certified machinery at `Y = 10^146` would either surface a smaller Carmichael
number or prove `S₆₄ ≥ 10^146`, narrowing the interval to **147–148 digits**. Whether the direct
enumerator suffices at `10^146` or the exchange decomposition is needed is open — that's the current
work.

## Provenance / credit

Known values and neighbours (`k = 3–63`, `k = 65` candidate) are Webster's
[`small-carmichael-numbers`](https://github.com/jewebste/small-carmichael-numbers) (Butler
University for `k = 36–63`). The lower-bound certificate's design benefited from consultation with a
large-language reasoning model (OpenAI GPT-6); the result here was independently reproduced by the
three implementations above, one of them gated by the frozen oracle. Existence for all large factor
counts is Alford–Grantham–Hayman–Shallue and Larsen–Wright; our contribution is the explicit
*interval* bounding the minimum.
