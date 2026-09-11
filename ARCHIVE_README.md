# carmichael-v1 — frozen archive (2026-09-11 23:20Z)

All Carmichael compute was stopped on this date and everything was archived here. `carmichael-v2/`
is the fresh start; nothing in this tree is running.

## Where things are

| Path | Contents |
|---|---|
| (this directory) | the git repo (`origin` = github.com/AlexanderEastwood/carmichael-frontier), HEAD `65cdca9` |
| `archive/logs/` | 504 run logs: GPU family/neighbourhood sweeps, r=8 portfolio blocks, exclusion drivers |
| `archive/checkpoints/` | 47 driver checkpoints and result JSONs from titan |
| `archive/records/gpu_records.tar.gz` | 108,966 per-modulus coverage records (12 families), + `INDEX.txt` |
| `archive/jobs/jobs_1e147.tar.gz` | the whole 10^147 exclusion job tree: 1,032,034 verdict files |
| `archive/jobs/control/` | JOBS / SUBJOBS1–11 / PROGRESS / DONE* kept loose for inspection |
| `archive/scratch/` | KMAX side experiments, the KMAX=0 duplicate runs, the crossover benchmark, the kill log |
| `archive/titan_source_tree.tar.gz` | titan's working mirror of the repo (binaries, scripts, un-pushed files) |
| `archive/ARCHIVE.log` | what the archive run copied, with sizes |

The identical archive exists on titan at `~/carmichael-v1/` (121 MB).

## State at freeze

**Results that stand (all pushed to GitHub):**
- `10^145 <= S_64 <= N < 10^148` with N = Webster's 148-digit number (communicated 2026-09-11,
  independently factored and oracle-verified here); engine exclusion gives `10^146 < S_64`.
- **Theorem 2:** N_0 (our GPU incumbent, 148 digits, largest prime 1249) is the least 64-factor
  Carmichael number with `lambda(n) | lambda(N_0)` — 111 moduli, 684 exchange instances, covered.
- **Theorem 3:** N is the least 64-factor Carmichael number with `lambda(n) | lambda(N)` —
  118 moduli, 728 instances, covered.
- The certified ladder `C_k(97) <= S_k <= U_k` for k = 64…144 (81 rows, all U_k oracle-verified).
- `S_65 in [10^148, U_65]`.
- Paper: `paper/interval_theorem.tex` / `.pdf` (12 pp, v5, reviewed by GPT-6 across four rounds).

**Unfinished when stopped:**
- *10^147 exclusion*: 1,021,546 of 1,021,547 leaves verdict `n=NONE`, **0 hits anywhere**. One
  all-smallest-primes slice outstanding (`19,31,37,41,43,47,53,59,61,67,71,73,79,89`), whose
  donated frontier `SUBJOBS11` was at 2,673/2,680 leaves. Finishing it would give
  `10^147 < S_64 <= N < 10^148`, i.e. **S_64 has exactly 148 digits**. To resume: untar
  `archive/jobs/jobs_1e147.tar.gz` into `experiments/dynt/`, run `finish11_donate.sh`, then
  `python3 verify_1e147.py` (verdict is `DONE.verified` only).
- *Experiment A (lambda-neighbourhood discovery)*: finished and negative — 50,142 moduli,
  215,544 exchange instances at radii <= 6, **0 hits**; 6,100 moduli >= 2^64 were not searched
  (u64 residue keys). Bounded negative, not a global exclusion.
- *Crossover benchmark* (interval + early AP completion vs baseline descent at bound N): built and
  gated (completion sets identical to a pre-patch baseline on 8 (k, bound, lo) cases) but killed
  before any arm completed. Harness: `experiments/dynt/bench_crossover.sh`.

**Open correspondence:** `drafts/reply_to_jon_2026-09-11.html` — draft reply to Jonathan Webster
retracting a mischaracterisation of his preproduct exhaustion. Not sent; Alex is handling it.

## Honest scope

The divisor-family method certifies minima only within a chosen `lambda` family; it does not cover
all Carmichael functions below N (every competitor has `lambda(n) | lcm{q-1 : q <= 2,893,785}`, which
has ~10^10030 divisors) and therefore certifies no global minimality. Determining S_64 rests on an
exhaustive preproduct search — Webster's, which was running when this work stopped. §2.1 of the
paper says exactly this.
