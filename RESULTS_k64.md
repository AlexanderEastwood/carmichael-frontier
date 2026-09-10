# k=64 push — RESULT (exchange-MITM route, M5)

Status: **POSITIVE — a verified 64-factor Carmichael number was found** (a strong
S_64 *candidate*, not a certified minimum). Run 2026-09-10 on a single 32-thread workstation (AMD Ryzen 9 9950X),
wall-capped at 900s.

## The candidate

- **decimal (149 digits):**
23197018055291196952520993997769543710196257486734152608091882983852073570129396414304979908337590233782784227789170068607803571949267286133678726401
- **factorization (64 distinct primes, 19 .. 1933):**
19 * 29 * 31 * 37 * 41 * 43 * 47 * 53 * 61 * 67 * 71 * 73 * 79 * 89 * 97 * 101 * 103 * 109 * 113 * 127 * 131 * 137 * 139 * 151 * 157 * 163 * 167 * 181 * 193 * 197 * 199 * 211 * 239 * 241 * 257 * 271 * 277 * 281 * 307 * 313 * 331 * 337 * 379 * 397 * 409 * 421 * 443 * 461 * 463 * 487 * 491 * 499 * 521 * 599 * 617 * 641 * 673 * 691 * 701 * 769 * 859 * 1151 * 1321 * 1933
- Found at modulus **M = 353649635539200** (portfolio provenance ((62, (173,), (83,)), 'N62', 5)),
  B0 = N62-derived 64-set, exchange **r = 5**.
- **Independently re-verified by the FROZEN oracle** ref_carmichael.verify_certificate:
  squarefree, 64 distinct primes, Korselt (p-1)|(n-1) holds for every p. (Every one of
  the 903 matches the run produced passed the oracle; 0 rejects.)

## Comparison to neighbors (all oracle-verified)

| k  | digits | largest prime |
|----|--------|---------------|
| 63 | 145    | 1201          |
| **64 (this candidate)** | **149** | **1933** |
| 65 | 151    | 1993          |

- **N63 (145d) < S64cand (149d) < N65 (151d)** — sits between the known neighbors, as a
  true S_64 should. S64/N63 ~ 3.2e3 ; N65/S64 ~ 2.0e2.
- Shares 53/64 primes with N63 and 51/64 with N65.
- This is almost certainly **NOT the true minimum S_64**: we explored only exchanges
  r<=5 from a handful of base sets over 54 of 120 portfolio moduli. It is a valid
  upper-bound incumbent U.

## How it was found (method)

Meet-in-the-middle over the exchange representation S = (B0\D) U I, |D|=|I|=r
(src/mitm64.cpp, the parallel successor to the gate prototype src/mitm.cpp):
1. Neighbor-seeded RICH modulus portfolio (src/mitm64_driver.py) from lambda(N62/63/65):
   core-exponent +-1 perturbations, drop<=2 / add<=3 mids from the Sophie-Germain +
   neighbor-mid pool. Deduped by integer value, ordered best-first by L1 distance to a
   neighbor lambda. 120-modulus cap this run (54 reached before the wall).
2. For each modulus M: pool Q(M) = {q prime : (q-1)|M, gcd(q,M)=1} (~150-185 primes).
   B0 variants: k-smallest, and N62/N63/N65-derived 64-sets. Insertion pool = the rest,
   bounded < 4000 and to 130 primes.
3. Parallel exchange-MITM mod the FULL M: recD builds a residue->max-prodD map
   (single thread); recI (the bottleneck) is split across 32 threads probing it. Every
   exact match S has lambda(S)=lcm(q-1) | M | (prod(S)-1) => Carmichael by Korselt.
4. Each match re-verified by the oracle; global incumbent U = smallest verified.

## Search statistics

- portfolio: 120 moduli (capped; best-first), **54 tried** before the 900s wall.
- pool sizes: ~150-185 primes per modulus.
- **guaranteed modulus G collapses to ~1 digit** (pool ~170 >> k=64 => r=|pool|-64 ~110
  deletions allowed => the (r+1)-th largest valuation of nearly every prime is 0). So the
  mod-G completeness route gives essentially NO pruning here; we matched mod the FULL M
  instead (sound: sufficient for Korselt).
- **r reached: 5** (RAM-guarded: Dres = C(64,r); C(64,5)=7.6M ~0.6GB; C(64,6)=7.5e7 ~6GB
  was skipped). recD ~2.0s to build 7.6M records; recI ~0.3-1.5s over 30-134M leaves at
  ~100M leaves/s across 32 threads.
- **903 exact matches, all verified Carmichael, 0 oracle rejects.** Incumbent 149 digits.

## Correctness caveats (and the Astra/GPT-6 flag)

- **We match mod the full candidate modulus M, NOT the coarse guaranteed-G modulus.** So
  every match is Carmichael *by construction* and the specific failure mode Astra flagged
  (a per-G-residue dedup silently dropping a candidate that would pass the true-lambda test)
  **does not arise in this pipeline** — there is no separate true-lambda test; all matches
  already satisfy lambda(S)|M.
- Dedup: recD keeps **max-prodD per EXACT residue key (prodD mod M, a u64 — an exact
  residue, not a hash fingerprint)**. For the min-product objective this is provably
  loss-free: min prod(S) = min over I of P0*prodI / (max prodD at its residue), so the
  global minimum is never dropped. **Empirically confirmed:** re-running the winning
  modulus with the truncation-free binary (mitm64_fixed: no per-thread cap, exact u128
  cross-product sort, all matches emitted) returns the *identical* 149-digit minimum.
- Output ordering is exact u128 (no float in acceptance or in the min); acceptance is
  exact residue equality + oracle. 64-bit counters throughout (r<=5 => Dres < 2^31 anyway).
- **Completeness gap (why this is a candidate, not a certified minimum):** matching mod
  full M MISSES any solution whose lambda(S) is a proper divisor of M (prod(S)=1 mod
  lambda(S) but !=1 mod M); and we only searched r<=5 from a few B0 over 54 moduli. A
  certified minimum is out of reach here (as the task notes).

## Recommended next widening (to push toward the true minimum)

1. **Better/searched B0** is the single biggest lever: the exchange distance from
   k-smallest to a real solution is ~7-8, but r is RAM-capped at 5. Seed B0 from a lattice
   of near-optimal 64-sets (interpolate N63<->N65) to bring the optimum within r<=4-5.
2. **Raise r past 5 by shrinking the D side** (swap roles so the smaller side is deletions),
   or stream Dres to disk / use a compact open-addressed table to reach r=6-7 within RAM.
3. **Finish the portfolio** (54/120 this run) and widen it (more core perturbations,
   the drop3/add4 tier) — cheap, best-first ordering means diminishing returns but may
   beat 149 digits.
4. **The mod-G / true-lambda route** would need a modulus with pool much closer to 64
   (fewer, larger primes) so G does not collapse; not the case for these smooth neighbor
   moduli.

## Where a GPU join would help a next, larger pass

The CPU bottleneck here is the **recI insertion scan** (30-134M leaves/call at r=5, and it
would be C(130,6)~1e10 at r=6). That is the batched, fixed-width, embarrassingly-parallel
join the assessment earmarked for CGBN-style GPU: generate insertion r-subset products +
residues on-GPU, sort, and match against the (host-built) residue->max-prodD table. Unlike
the earlier Fermat-leaf idea (which the M2.5F CPU filter already made cheap), THIS join is a
genuinely large uniform-width kernel and is the right place to spend a GPU. It raises the
feasible r (hence reachable exchange distance) rather than accelerating an already-cheap leaf.

Artifacts: src/mitm64.cpp, src/mitm64_driver.py, src/mitm64 (canonical) / mitm64_fixed,
results_k64.json, checkpoint_k64.json, run_k64.log.
