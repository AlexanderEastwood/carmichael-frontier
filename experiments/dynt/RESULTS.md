# Dynamic-t crossover experiment

A controlled A/B of **where** the exact-k search switches from prefix *descent* to the
arithmetic-progression (AP) *cofactor scan*. Motivated by Shallue–Webster §4.4/§5.3 (preproduct
construction overtaking completion work at high k) and by GPT-6 (Astra) review guidance to **measure**
the crossover rather than adopt Algorithm 3's `P·λ(P)² > B` rule unchanged.

## What the engine already does, and the one change

`src/carmichael.cpp` (production, M2.5F) switches to the AP scan (`base_case`) when
`t ≤ 1 || P·λ(P) > B` — the *simpler* §7.4 rule. `base_case` already handles arbitrary `t`
(the AP scan `R ≡ P⁻¹ mod λ(P)`, `R ≤ B/P`, factor `R` into `t` primes, Korselt).

`carmichael_dynt.cpp` is that engine with **one flag-guarded change** (`-DDYNAMIC_T`): also switch
early when a **cheap** AP-length upper bound
`⌊(B/P)/λ(P)⌋ ≤ DYNT_KMAX` — no modular inverse needed to decide (Astra's "screen switching before
paying for an inverse"). `DYNT_KMAX=0` recovers production **exactly**. Added knobs:
`DYNT_FREEZE` (pin the incumbent so both arms traverse identical space and collect the same
completion set) and `--prefix p1,p2,…` (run one frozen subtree). Everything is inert without the
flag/env, so the default build is behaviorally identical to `src/carmichael.cpp`.

By Shallue–Webster Theorem 11 the switch point affects **efficiency, not completeness** — the
harness asserts this by reporting the full completion set (`completions=…`) in every arm.

## Correctness

- `DYNT_KMAX=0` reproduces the known minima exactly: S₂₀, S₂₅, S₂₉ (= OEIS A006931 / frozen oracle).
- Every early-switch arm finds the **identical** minimum (full search) / identical completion set
  (frozen prefix). Completeness preserved across all KMAX.

## k = 29 (full search, bound-doubling, titan Ryzen 9 9950X single core)

| DYNT_KMAX | descent N_E | leaf N_H | enum-frac | wall (s) | finds S₂₉ |
|---:|---:|---:|---:|---:|:--:|
| 0 (production) | 16,024,440 | 150,659 | 0.99 | 10.58 | ✓ |
| 5 | 11,927,662 | 678,923 | 0.94 | **9.12** | ✓ |
| 20 | 9,834,835 | 1,950,116 | 0.83 | **9.03** | ✓ |
| 100 | 7,757,070 | 7,439,322 | 0.53 | 12.18 | ✓ |

≈ **15 %** wall improvement at KMAX ≈ 5–20; past it, leaf work dominates and it regresses.

## k = 64 (frozen representative subtree; prefix = first 52 primes of the incumbent N, t = 12, bound = N)

| DYNT_KMAX | descent N_E | AP N_F | leaf N_H | completions | max R bits | wall (s) |
|---:|---:|---:|---:|:--:|---:|---:|
| 0 (production) | 34,240,783 | 290,427 | 140,884 | 0 | 79 | 40.61 |
| 5 | 13,028,196 | 2,471,333 | 1,194,377 | 0 | 79 | 19.47 |
| **20** | 8,818,214 | 4,410,964 | 2,146,730 | **0** | 80 | **16.28** |
| 100 | 6,571,991 | 8,347,696 | 4,092,922 | 0 | 80 | 18.46 |
| 1000 | 4,450,439 | 43,487,888 | 21,212,911 | 0 | 80 | 59.44 |

≈ **2.5× speedup** at KMAX ≈ 20 (40.6 s → 16.3 s), completion set unchanged (0 throughout).

## Findings

1. **The crossover is real and its value grows with k** — ~15 % at k=29, ~2.5× at k=64 — consistent
   with the §4.4 diagnosis that preproduct construction increasingly dominates. This is the paper's
   actionable opportunity, using only the existing completion machinery (no new GPU layer).
2. **There is a clear optimum (~KMAX = 20 here) and a hard regression past it** (KMAX = 1000: 59 s > 40 s
   production). So a *fixed* aggressive rule is wrong; the switch must be tuned/measured — exactly
   Astra's warning against adopting the squared rule blind.
3. **Astra's "early switch inflates R / factoring explodes" caution does not bite at this depth**:
   the widest cofactor R actually factored stays **79–80 bits** across all KMAX (the leaf cost per
   candidate is roughly flat; the AP simply processes more candidates). The regression is candidate
   *count* (N_F/N_H), not operand width.
4. `DYNT_KMAX` uses a cheap **upper bound** on the AP length; a per-node cost model (candidate count ×
   measured leaf cost vs. subtree size) would place the switch better, and is the natural next step.

## Reproduce

```
g++ -O3 -std=c++17 -DDYNAMIC_T carmichael_dynt.cpp -o cm_dynt -lgmp -lgmpxx

# k=29 full-search crossover
for KM in 0 5 20 100; do DYNT_KMAX=$KM ./cm_dynt 29; done

# k=64 frozen-subtree A/B (PFX = first 52 primes of N; NN = the 148-digit incumbent)
for KM in 0 5 20 100 1000; do DYNT_FREEZE=1 DYNT_KMAX=$KM ./cm_dynt 64 <NN> --prefix <PFX>; done
```

`<NN>` and `<PFX>` (comma-separated) come from `results_k64_best_global.json`.

## Caveats / next

- One representative subtree (t=12) at k=64; sweep more prefixes and depths for a switch *profile*.
- The interval floor `L = 10¹⁴⁵` (Astra's `a = max(R_min, ⌈L/P⌉)`) is not yet imposed in the AP start;
  it only removes work equally from both arms here, but the real `[10¹⁴⁵, N)` certification wants it.
- Replace the fixed KMAX with an adaptive per-node cost comparison; then parallelize across prefixes
  (the M3 shared-incumbent design) and re-measure.
