#!/usr/bin/env python3
"""
mitm_gpu_mstar_driver.py -- expanded-insertion-pool exchange-MITM at ONE modulus
(the incumbent's own, M* = 1768248177696000), radii 5..8, GPU D side + streamed I side.

Why (review, 2026-09-11): the r=8 portfolio sweep improved the incumbent only at
M*, and at M* the 4000-prime insertion cutoff is BINDING (eligible/insertions:
<4000 -> 171/107, <10000 -> 254/190, whole product-bounded pool -> 525/461) while
the 160-insertion cap is not. So: keep M*, raise the prime cutoff, and stratify
by the number h of insertion primes above SPLIT (4000): with D_max_r = product of
the r largest base factors, a_j the small insertion primes ascending, c_j the
large ones ascending, every completion with h large insertions satisfies
    n >= (P_B / D_max_r) * prod_{j<=r-h} a_j * prod_{j<=h} c_j ,
so h is skipped whenever that bound is >= the incumbent (at r=8 this gives
h <= 2 for both bases). Exact radii 1..4 over the full pool were already
searched exhaustively (no improvement), so radii start at 5.

Fixes from the same review vs the block driver: (i) completion identity is
(M, base tag, sha256(INS), r), never just (r, idx); (ii) the I side uses
src/mitm64_idump_stream2 (cheapest-completion bound + loop bound; validated
record-identical to v1 and ~10x faster), with IDUMP_SPLIT/IDUMP_HMAX per instance.

Every join pair is re-checked on the host in exact big-int and by the frozen
oracle before it can touch the incumbent. Own checkpoint/results files.
Run:  ~/erdos/.venv/bin/python gpu/mitm_gpu_mstar_driver.py
Env:  R_MIN (5) R_MAX (8) INS_PRIME_CAP (10000) INS_MAX (256) SPLIT (4000)
      BASES (S64,ksmall) NCHUNKS_D (16) I_CHUNK (200000000) IDUMP_T (8)
      IDUMP_TIMEOUT (3600) WALL_CAP (14400) MODULUS (from results_k64_best_global.json)
"""
import os, sys, json, time, hashlib, subprocess
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__)); FR = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(FR, "ref"))
import ref_carmichael as ref
import mitm_gpu2 as G
import numpy as np, cupy as cp

R_MIN = int(os.environ.get("R_MIN", "5")); R_MAX = int(os.environ.get("R_MAX", "8"))
INS_PRIME_CAP = int(os.environ.get("INS_PRIME_CAP", "10000")); INS_MAX = int(os.environ.get("INS_MAX", "256"))
SPLIT = int(os.environ.get("SPLIT", "4000")); BASES = [b for b in os.environ.get("BASES", "S64,ksmall").split(",") if b]
NCHUNKS_D = int(os.environ.get("NCHUNKS_D", "16")); I_CHUNK = int(os.environ.get("I_CHUNK", "200000000"))
IDUMP_T = int(os.environ.get("IDUMP_T", "8")); IDUMP_TIMEOUT = float(os.environ.get("IDUMP_TIMEOUT", "3600"))
WALL_CAP = float(os.environ.get("WALL_CAP", "14400"))
IDUMP = os.path.join(FR, "src", "mitm64_idump_stream2")
CKPT_TAG = os.environ.get("CKPT_TAG", "")        # e.g. "r9" -> checkpoint_mstar_r9.json / results_k64_gpu_mstar_r9.json
CKPT = os.path.join(HERE, f"checkpoint_mstar{('_' + CKPT_TAG) if CKPT_TAG else ''}.json")
RESULTS = os.path.join(HERE, f"results_k64_gpu_mstar{('_' + CKPT_TAG) if CKPT_TAG else ''}.json")
DFILTER = os.environ.get("DFILTER", "0") == "1"   # product-filtered deletion side (src/mitm64_ddump), resident on GPU
RANK = os.environ.get("RANK", "0") == "1"         # force colex-rank insertion ids (automatic for r > 8)
DDUMP = os.path.join(FR, "src", "mitm64_ddump")
K = 64; T0 = time.time()
def log(*a): print(f"[{time.time()-T0:8.1f}s]", *a, flush=True)

inc = json.load(open(os.path.join(FR, "results_k64_best_global.json")))
U0 = int(inc["n"]); U0_factors = sorted(int(p) for p in inc["factors"])
if os.environ.get("U_OVERRIDE"):                  # validation only: search below a frozen (older) bound
    U0 = int(os.environ["U_OVERRIDE"]); assert U0 > 10**140
M = int(os.environ.get("MODULUS") or inc["modulus"])
state = {"modulus": str(M), "orig_incumbent": str(U0), "incumbent": None, "incumbent_factors": None,
         "done": [], "instances": 0, "pairs_total": 0, "completions_below_U": 0, "oracle_rejects": 0,
         "best_digits_seen": len(str(U0)), "status": "running", "started": time.strftime("%Y-%m-%d %H:%M:%S")}
if os.path.exists(CKPT):
    old = json.load(open(CKPT))
    if old.get("modulus") == str(M):
        for k in ("incumbent", "incumbent_factors", "done", "instances", "pairs_total", "completions_below_U", "oracle_rejects", "best_digits_seen"):
            if k in old: state[k] = old[k]
        log("RESUMING:", len(state["done"]), "instances already completed")
def cur_U(): return int(state["incumbent"]) if state["incumbent"] else U0
def save():
    state["elapsed_s"] = round(time.time() - T0, 1); tmp = CKPT + ".tmp"
    json.dump(state, open(tmp, "w"), indent=1); os.replace(tmp, CKPT)
def gpu_used_gb(): f, t = cp.cuda.Device().mem_info; return (t - f) / 2**30

def pool(M, cap):
    out = []; q = 3
    while q < cap:
        if M % (q - 1) == 0 and M % q != 0 and ref.is_prime(q): out.append(q)
        q += 2
    return out

def hmax_for(B0, INS, r, U):
    """Largest h (number of insertion primes > SPLIT) that could still improve on U; -1 if none."""
    small = [p for p in INS if p <= SPLIT]; large = [p for p in INS if p > SPLIT]
    PB = 1
    for p in B0: PB *= p
    Dmax = 1
    for p in sorted(B0)[-r:]: Dmax *= p
    best = -1
    for h in range(0, r + 1):
        if h > len(large) or r - h > len(small): continue
        v = PB
        for p in small[:r - h]: v *= p
        for p in large[:h]: v *= p
        if v < U * Dmax: best = h            # bound below the incumbent: this stratum can improve
    return best

def build_D(B0, M, r, Dmin=None):
    """Sorted D-side chunks. DFILTER: only deletion subsets with prodD >= Dmin (src/mitm64_ddump,
    CPU, resident on GPU); else the full C(64,r) via the GPU kernel in NCHUNKS_D pieces."""
    if DFILTER:
        assert Dmin is not None and os.path.exists(DDUMP), DDUMP
        inp = (f"{M} {r} {Dmin}\n64 " + " ".join(map(str, sorted(B0))) + "\n").encode()
        pr = subprocess.run([DDUMP], input=inp, capture_output=True, check=True)
        a = np.frombuffer(pr.stdout, dtype=np.uint64).reshape(-1, 2)
        nD = int([l for l in pr.stderr.decode().splitlines() if l.startswith("TOTAL")][0].split()[1]); assert a.shape[0] == nD
        chunks = []
        for part in np.array_split(np.arange(nD), max(1, min(NCHUNKS_D, nD // 50_000_000 + 1))):
            if part.size == 0: continue
            k = cp.asarray(np.ascontiguousarray(a[part, 0])); i = cp.asarray(np.ascontiguousarray(a[part, 1]))
            chunks.append(G.sort_side(k, i)); del k, i
        del a; cp.get_default_memory_pool().free_all_blocks(); return chunks, nD
    tot = comb(K, r); b = [tot * i // NCHUNKS_D for i in range(NCHUNKS_D + 1)]; chunks = []
    for a, c in zip(b, b[1:]):
        k, i = G.generate(B0, r, M, 1, bits=6, lo=a, hi=c); chunks.append(G.sort_side(k, i)); del k, i
    cp.get_default_memory_pool().free_all_blocks(); return chunks, tot

def unrank_colex(v, r, n):
    """Inverse of the dumper's colex rank id = sum_j C(i_j, j+1), i_0 < ... < i_{r-1} < n."""
    v = int(v); idx = []; x = n - 1
    for j in range(r - 1, -1, -1):
        while comb(x, j + 1) > v: x -= 1
        idx.append(x); v -= comb(x, j + 1); x -= 1
    assert v == 0; return tuple(reversed(idx))

def _read_block(f, nbytes):
    buf = bytearray(nbytes); mv = memoryview(buf); got = 0
    while got < nbytes:
        n = f.readinto(mv[got:])
        if not n: break
        got += n
    return buf[:got]

def consider(prov, factors, lever):
    if len(factors) != K or len(set(factors)) != K: return False
    ok, n, _why = ref.verify_certificate(sorted(factors))
    if not ok: state["oracle_rejects"] += 1; return False
    dg = len(str(n)); state["best_digits_seen"] = min(state["best_digits_seen"], dg)
    if n < cur_U():
        state.update(incumbent=str(n), incumbent_factors=sorted(factors))
        log(f"  *** NEW BEST k=64 CARMICHAEL {dg} digits (< {len(str(U0))}) prov {prov} lever {lever} ***")
        json.dump({"n": str(n), "factors": sorted(factors), "modulus": str(M), "provenance": str(prov), "digits": dg, "lever": lever,
                   "improves_over": str(U0), "orig_digits": len(str(U0)), "oracle": "Carmichael with exactly 64 prime factors",
                   "when": time.strftime("%Y-%m-%d %H:%M:%S")}, open(RESULTS, "w"), indent=1)
        save(); return True
    return False

def run_instance_host(tag, B0, INS, r, hmax, P0, U, Icap, Dmin, n_full):
    """Exact host-side join for tiny instances or Icap beyond 126 bits (r=14 at M* needs 134-bit
    products): D side from ddump (filtered), I side enumerated in Python with the same product,
    cheapest-completion and h-stratum bounds, joined by residue mod M."""
    from itertools import combinations
    t0 = time.time(); B0set = set(B0)
    inp = (f"{M} {r} {Dmin}\n64 " + " ".join(map(str, sorted(B0))) + "\n").encode()
    pr = subprocess.run([DDUMP], input=inp, capture_output=True, check=True)
    a = np.frombuffer(pr.stdout, dtype=np.uint64).reshape(-1, 2); nD = a.shape[0]
    Dmap = {}
    for kkey, v in a:
        d = tuple(b for b in range(K) if (int(v) >> b) & 1) if r > 10 else tuple((int(v) >> (6 * q)) & 63 for q in range(r))
        Dmap.setdefault(int(kkey), []).append(d)
    nI = 0; pairs = 0; below = 0; nb = 0; rmin = None; P0m = P0 % M
    for c in combinations(range(len(INS)), r):
        if sum(1 for i in c if INS[i] > SPLIT) > hmax: continue
        prodI = 1
        for i in c: prodI *= INS[i]
        if prodI >= Icap: continue
        nI += 1; key = (P0m * (prodI % M)) % M
        for d in Dmap.get(key, ()):
            pairs += 1; D = [B0[x] for x in d]; prodD = 1
            for p in D: prodD *= p
            if (P0 * prodI) % prodD: continue
            n = P0 * prodI // prodD
            if n >= cur_U(): continue
            below += 1; fac = sorted((B0set - set(D)) | set(INS[i] for i in c))
            if consider((M, tag, r, hmax), fac, f"mstar_{tag}_r{r}_h<={hmax}_host"): nb += 1
            if rmin is None or n < rmin: rmin = n
    state["instances"] += 1; state["pairs_total"] += pairs; state["completions_below_U"] += below
    log(f"  {tag} r={r} h<={hmax} |INS|={len(INS)}/{n_full} host-exact: D {nD:,} | I {nI:,} rec {time.time()-t0:.0f}s pairs={pairs} below_U={below} "
        f"min_digits={len(str(rmin)) if rmin else None} improved={nb}")
    return True

def run_instance(tag, B0, INS, r, hmax):
    B0set = set(B0); B0s = sorted(B0); P0 = 1
    for p in B0: P0 *= p
    U = cur_U(); topI = 1
    for p in INS[-r:]: topI *= p
    ptb = 1
    for p in B0s[-r:]: ptb *= p
    Icap = max(1, min(topI, (U * ptb) // P0 + 1)); pdm = 1
    for p in B0s[:r]: pdm *= p
    Imin = 1
    for p in INS[:r]: Imin *= p
    Dmin = (P0 * Imin) // U + 1                    # improving completion needs prodD >= Dmin (P_B*P_I/P_D < U)
    # per-radius insertion cap: in any r-subset with product < Icap the largest prime q satisfies
    # q * (product of the r-1 smallest insertion primes) < Icap  =>  the affordable pool shrinks with r
    Imin1 = 1
    for p in INS[:r - 1]: Imin1 *= p
    cap_r = (Icap - 1) // Imin1; INS_full = INS; INS = [q for q in INS if q <= cap_r]
    if len(INS) < r:
        log(f"  {tag} r={r}: affordable pool has {len(INS)} < r primes (cap_r={cap_r}); nothing to search"); state["instances"] += 1; return True
    rank = RANK or r > 8 or len(INS) > 256         # 8-bit packed I ids only hold r <= 8 and pools <= 256 primes
    if Icap.bit_length() > 126 or comb(len(INS), r) < 200_000:      # tiny or beyond u128: exact host join
        return run_instance_host(tag, B0, INS, r, hmax, P0, U, Icap, Dmin, len(INS_full))
    t = time.time(); Dch, nD = build_D(B0, M, r, Dmin); tD = time.time() - t
    inp = ("\n".join([f"{M} {r} {P0 % M} {Icap} {pdm} {IDUMP_T} 0", "64 " + " ".join(map(str, B0)), f"{len(INS)} " + " ".join(map(str, INS))]) + "\n").encode()
    env = dict(os.environ, IDUMP_SPLIT=str(SPLIT), IDUMP_HMAX=str(hmax), IDUMP_RANK="1" if rank else "0")
    pr = subprocess.Popen([IDUMP], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)
    assert pr.stdin and pr.stdout and pr.stderr
    pr.stdin.write(inp); pr.stdin.close()
    t0 = time.time(); pairs = 0; below = 0; nb = 0; rmin = None; read_rec = 0; timed_out = False
    while True:
        blk = _read_block(pr.stdout, I_CHUNK * 16)
        if not blk: break
        a = np.frombuffer(blk, dtype=np.uint64).reshape(-1, 2); read_rec += a.shape[0]
        Ik = cp.asarray(np.ascontiguousarray(a[:, 0])); Iid = cp.asarray(np.ascontiguousarray(a[:, 1])); del a, blk
        Ik_s, Iid_s = G.sort_side(Ik, Iid); del Ik, Iid
        for Dk_s, Did_s in Dch:
            pd, pi = G.join_sorted(Dk_s, Did_s, Ik_s, Iid_s)
            if pd.size == 0: continue
            pairs += int(pd.size)
            Idec = [unrank_colex(v, r, len(INS)) for v in cp.asnumpy(pi)] if rank else G.decode_ids(pi, r, 8)
            Ddec = [tuple(b for b in range(K) if (int(v) >> b) & 1) for v in cp.asnumpy(pd)] if r > 10 else G.decode_ids(pd, r, 6)
            for d, i in zip(Ddec, Idec):
                D = [B0[x] for x in d]; I = [INS[x] for x in i]
                prodD = 1
                for p in D: prodD *= p
                prodI = 1
                for p in I: prodI *= p
                if (P0 * prodI) % prodD: continue
                n = P0 * prodI // prodD
                if n >= cur_U(): continue
                below += 1; fac = sorted((B0set - set(D)) | set(I))
                if consider((M, tag, r, hmax), fac, f"mstar_{tag}_r{r}_h<={hmax}"): nb += 1
                if rmin is None or n < rmin: rmin = n
        del Ik_s, Iid_s; cp.get_default_memory_pool().free_all_blocks()
        if time.time() - t0 > IDUMP_TIMEOUT: timed_out = True; pr.kill(); break
    pr.wait(); err = pr.stderr.read().decode(errors="replace"); nI = 0
    for line in err.splitlines():
        if line.startswith("TOTAL"): nI = int(line.split()[1])
    del Dch; cp.get_default_memory_pool().free_all_blocks(); dt = time.time() - t0; state["instances"] += 1
    if timed_out or read_rec != nI or pr.returncode != 0:
        log(f"  {tag} r={r} h<={hmax}: INCOMPLETE ({'timeout' if timed_out else f'read {read_rec:,} != TOTAL {nI:,} rc={pr.returncode}'}) after {dt:.0f}s; NOT counted")
        return False
    state["pairs_total"] += pairs; state["completions_below_U"] += below
    log(f"  {tag} r={r} h<={hmax} |INS|={len(INS)}{' dfilter' if DFILTER else ''}{' rank' if rank else ''}: D {nD:,} {tD:.1f}s | I {nI:,} rec {dt:.0f}s pairs={pairs} below_U={below} "
        f"min_digits={len(str(rmin)) if rmin else None} improved={nb} | GPU {gpu_used_gb():.1f}GB")
    return True

def main():
    assert os.path.exists(IDUMP), IDUMP
    P = pool(M, INS_PRIME_CAP); Pset = set(P)
    bases = []
    for tag in BASES:
        if tag == "ksmall": bases.append(("ksmall", sorted(P)[:K]))
        elif tag == "S64":
            b = [p for p in U0_factors if p in Pset]; extra = [p for p in sorted(P) if p not in set(b)]
            bases.append(("S64", sorted(b + extra[:K - len(b)])))
    log(f"M={M} pool<{INS_PRIME_CAP}: {len(P)} eligible primes; bases={[t for t,_ in bases]}; radii {R_MIN}..{R_MAX}; SPLIT={SPLIT}; U={len(str(cur_U()))} digits")
    save()
    for r in range(R_MIN, R_MAX + 1):
        for tag, B0 in bases:
            if time.time() - T0 > WALL_CAP: state["status"] = "wall_cap_reached"; save(); log("WALL CAP"); return
            INS = sorted([p for p in P if p not in set(B0)])[:INS_MAX]
            key = f"{M}|{tag}|{hashlib.sha256(','.join(map(str,INS)).encode()).hexdigest()[:16]}|r{r}"
            if key in state["done"]: log(f"  skip {tag} r={r} (identity {key} already completed)"); continue
            hmax = hmax_for(B0, INS, r, cur_U())
            if hmax < 0: log(f"  {tag} r={r}: no stratum can improve on U (bound); skipped"); state["done"].append(key); save(); continue
            if run_instance(tag, B0, INS, r, hmax): state["done"].append(key)
            else: state["failed"] = state.get("failed", 0) + 1
            save()
    state["status"] = "complete" if not state.get("failed") else f"INCOMPLETE ({state['failed']} instances failed)"; save()
    log("DONE", state["status"], "best_digits_seen:", state["best_digits_seen"], "incumbent improved:", state["incumbent"] is not None)

if __name__ == "__main__": main()
