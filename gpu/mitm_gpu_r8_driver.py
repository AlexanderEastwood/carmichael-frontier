#!/usr/bin/env python3
"""
mitm_gpu_r8_driver.py -- exchange-MITM incumbent search at r=8 (and beyond the
CPU's r<=7 cap) with the D side on the GPU.

Why the GPU lifts the cap: at r=8 the D side has C(64,8) ~ 4.4e9 residue records;
the CPU engine keeps its D side in an in-memory hash, which is infeasible there.
Here the GPU generates the D side in rank-chunks (gpu/mitm_gpu2.generate,
bits=6), sorts each chunk (~0.3 s per 1e9 keys), and keeps them resident. The I
side stays product-pruned on the CPU (src/mitm64_idump, validated record-exact
against the engine's I_leaves) and is streamed through the GPU sort-merge join
against every resident D chunk. Every equal-residue (D,I) pair is emitted; the
completion S=(B0\\D)uI is then re-checked on the host in exact big-int (exactly 64
distinct primes, exact product, and Korselt via the FROZEN oracle) before it can
touch the incumbent. Nothing is trusted without the oracle.

Search definition (portfolio, bases, insertion pools, Icap) is imported UNCHANGED
from stageB_driver.py so results are directly comparable with the CPU sweep.
Own checkpoint/results files; never touches results_k64.json / _best_global.json.
Bounded: an I-side dump larger than MAX_IDUMP_GB is skipped and logged (r=8 is
tried where it is affordable, never promised everywhere).

Run:  ~/erdos/.venv/bin/python gpu/mitm_gpu_r8_driver.py
Env:  R (8) IDX_MIN (0) IDX_MAX (120) BASES (S64,ksmall) NCHUNKS_D (8)
      I_CHUNK (200000000 records) IDUMP_T (8) MAX_IDUMP_GB (40) IDUMP_TIMEOUT (1800)
      SCRATCH (/tmp/mitm_r8) WALL_CAP (36000)
"""
import os, sys, json, time, glob, subprocess, importlib.util
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__)); FR = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(FR, "ref"))
import ref_carmichael as ref
import mitm_gpu2 as G
import numpy as np, cupy as cp

R = int(os.environ.get("R", "8")); IDX_MIN = int(os.environ.get("IDX_MIN", "0")); IDX_MAX = int(os.environ.get("IDX_MAX", "120"))
BASES = [b for b in os.environ.get("BASES", "S64,ksmall").split(",") if b]
NCHUNKS_D = int(os.environ.get("NCHUNKS_D", "8")); I_CHUNK = int(os.environ.get("I_CHUNK", "200000000"))
IDUMP_T = int(os.environ.get("IDUMP_T", "8")); MAX_IDUMP_GB = float(os.environ.get("MAX_IDUMP_GB", "40"))
IDUMP_TIMEOUT = float(os.environ.get("IDUMP_TIMEOUT", "1800")); SCRATCH = os.environ.get("SCRATCH", "/tmp/mitm_r8")
WALL_CAP = float(os.environ.get("WALL_CAP", "36000"))
IDUMP = os.path.join(FR, "src", "mitm64_idump")
STREAM = os.environ.get("STREAM", "1") != "0"                      # stream the I side straight into the join (no disk)
IDUMP_STREAM = os.path.join(FR, "src", "mitm64_idump_stream")
CKPT = os.path.join(HERE, f"checkpoint_r{R}.json"); RESULTS = os.path.join(HERE, f"results_k64_gpu_r{R}.json")
T0 = time.time()
def log(*a): print(f"[{time.time()-T0:8.1f}s]", *a, flush=True)

# ---- import the Stage-B search definition without running its sweep ----
spec = importlib.util.spec_from_file_location("stageB_driver", os.path.join(FR, "stageB_driver.py"))
assert spec is not None and spec.loader is not None, "stageB_driver.py not found"
SB = importlib.util.module_from_spec(spec); spec.loader.exec_module(SB)   # __name__ != "__main__": main() not run
K = 64; assert R * 6 <= 64 and R <= 8, "D side packs 6 bits/idx; I side 8 bits/idx (r<=8)"

inc = json.load(open(os.path.join(FR, "results_k64_best_global.json")))
U0 = int(inc["n"]); U0_factors = sorted(int(p) for p in inc["factors"])
state = {"incumbent": None, "incumbent_factors": None, "incumbent_modulus": None, "orig_incumbent": str(U0), "orig_digits": len(str(U0)),
         "r": R, "instances": 0, "skipped_big_I": 0, "pairs_total": 0, "completions_below_U": 0, "oracle_rejects": 0,
         "best_digits_seen": len(str(U0)), "last_completed_idx": IDX_MIN - 1, "status": "running", "started": time.strftime("%Y-%m-%d %H:%M:%S")}
if os.path.exists(CKPT):
    old = json.load(open(CKPT))
    if old.get("r") == R:
        for k in ("incumbent", "incumbent_factors", "incumbent_modulus", "instances", "skipped_big_I", "pairs_total",
                  "completions_below_U", "oracle_rejects", "best_digits_seen", "last_completed_idx"):
            if k in old: state[k] = old[k]
        log("RESUMING from idx", state["last_completed_idx"])
def cur_U(): return int(state["incumbent"]) if state["incumbent"] else U0
def save():
    state["elapsed_s"] = round(time.time() - T0, 1); tmp = CKPT + ".tmp"
    json.dump(state, open(tmp, "w"), indent=1); os.replace(tmp, CKPT)

def gpu_used_gb(): f, t = cp.cuda.Device().mem_info; return (t - f) / 2**30

def build_D(B0, M):
    """C(64,R) D-side records in NCHUNKS_D sorted rank-chunks, resident on GPU."""
    tot = comb(K, R); b = [tot * i // NCHUNKS_D for i in range(NCHUNKS_D + 1)]; chunks = []
    for a, c in zip(b, b[1:]):
        k, i = G.generate(B0, R, M, 1, bits=6, lo=a, hi=c); chunks.append(G.sort_side(k, i)); del k, i
    cp.get_default_memory_pool().free_all_blocks(); return chunks

def run_idump(M, B0, INS, P0modM, Icap, pdm, prefix):
    """Stream the product-pruned I side to per-thread .bin files; enforce size cap. Returns (ok, nrecords, reason)."""
    for f in glob.glob(prefix + ".t*.bin"): os.remove(f)
    inp = "\n".join([f"{M} {R} {P0modM} {Icap} {pdm} {IDUMP_T} 0", "64 " + " ".join(map(str, B0)), f"{len(INS)} " + " ".join(map(str, INS))]) + "\n"
    pr = subprocess.Popen([IDUMP, prefix], stdin=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert pr.stdin is not None and pr.stderr is not None
    pr.stdin.write(inp); pr.stdin.close(); t = time.time()
    while pr.poll() is None:
        time.sleep(2)
        sz = sum(os.path.getsize(f) for f in glob.glob(prefix + ".t*.bin")) / 2**30
        if sz > MAX_IDUMP_GB: pr.kill(); pr.wait(); return False, 0, f"I-side dump exceeded {MAX_IDUMP_GB} GB"
        if time.time() - t > IDUMP_TIMEOUT: pr.kill(); pr.wait(); return False, 0, f"I-side dump exceeded {IDUMP_TIMEOUT}s"
    err = pr.stderr.read(); n = 0
    for line in err.splitlines():
        if line.startswith("TOTAL"): n = int(line.split()[1])
    sz = sum(os.path.getsize(f) for f in glob.glob(prefix + ".t*.bin")) / 2**30   # re-check AFTER exit (poll gap)
    if sz > MAX_IDUMP_GB: return False, n, f"I-side dump {sz:.1f} GB exceeds cap {MAX_IDUMP_GB} GB"
    return True, n, (err.strip().splitlines()[-1] if err.strip() else "") + f" [{sz:.1f} GB]"

def consider(prov, factors, M, lever):
    if len(factors) != K or len(set(factors)) != K: return False
    ok, n, _why = ref.verify_certificate(sorted(factors))
    if not ok: state["oracle_rejects"] += 1; return False
    dg = len(str(n))
    if dg < state["best_digits_seen"]: state["best_digits_seen"] = dg
    if n < cur_U():
        state.update(incumbent=str(n), incumbent_factors=sorted(factors), incumbent_modulus=str(M), incumbent_prov=str(prov), incumbent_digits=dg)
        log(f"  *** NEW BEST k=64 CARMICHAEL {dg} digits (< {len(str(U0))}) prov {prov} lever {lever} ***")
        json.dump({"n": str(n), "factors": sorted(factors), "modulus": str(M), "provenance": str(prov), "digits": dg, "lever": lever,
                   "improves_over": str(U0), "orig_digits": len(str(U0)), "oracle": "Carmichael with exactly 64 prime factors",
                   "when": time.strftime("%Y-%m-%d %H:%M:%S")}, open(RESULTS, "w"), indent=1)
        save(); return True
    return False

def _read_block(f, nbytes):
    """Read up to nbytes from a pipe, filling as much as possible (short only at EOF)."""
    buf = bytearray(nbytes); mv = memoryview(buf); got = 0
    while got < nbytes:
        n = f.readinto(mv[got:])
        if not n: break
        got += n
    return buf[:got]

def instance_stream(idx, tag, prov, M, B0, B0set, INS, P0, Icap, pdm, Dch, tD):
    """STREAM mode: pipe the dumper's records straight into the GPU join, one I_CHUNK at a time.
    The I side never touches disk; an instance is bounded only by IDUMP_TIMEOUT. Integrity: the
    number of records consumed must equal the dumper's TOTAL, else the instance is NOT counted."""
    inp = ("\n".join([f"{M} {R} {P0 % M} {Icap} {pdm} {IDUMP_T} 0", "64 " + " ".join(map(str, B0)),
                      f"{len(INS)} " + " ".join(map(str, INS))]) + "\n").encode()
    pr = subprocess.Popen([IDUMP_STREAM], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert pr.stdin is not None and pr.stdout is not None and pr.stderr is not None
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
            for d, i in zip(G.decode_ids(pd, R, 6), G.decode_ids(pi, R, 8)):
                D = [B0[x] for x in d]; I = [INS[x] for x in i]
                prodD = 1
                for p in D: prodD *= p
                prodI = 1
                for p in I: prodI *= p
                if (P0 * prodI) % prodD: continue
                n = P0 * prodI // prodD
                if n >= cur_U(): continue
                below += 1; fac = sorted((B0set - set(D)) | set(I))
                if consider((prov, tag, R), fac, M, f"{tag}_r{R}"): nb += 1
                if rmin is None or n < rmin: rmin = n
        del Ik_s, Iid_s; cp.get_default_memory_pool().free_all_blocks()
        if time.time() - t0 > IDUMP_TIMEOUT: timed_out = True; pr.kill(); break
    pr.wait(); err = pr.stderr.read().decode(errors="replace"); nI = 0
    for line in err.splitlines():
        if line.startswith("TOTAL"): nI = int(line.split()[1])
    state["instances"] += 1; dt = time.time() - t0
    if timed_out or read_rec != nI or pr.returncode != 0:
        state["skipped_big_I"] += 1
        why = "timeout" if timed_out else (f"rc={pr.returncode}" if pr.returncode != 0 else f"read {read_rec:,} != TOTAL {nI:,}")
        log(f"  M#{idx} {tag} r={R}: INCOMPLETE ({why}) after {dt:.0f}s, {read_rec:,} records seen; instance NOT counted | D {comb(K,R):,} {tD:.1f}s")
        return
    state["pairs_total"] += pairs; state["completions_below_U"] += below
    log(f"  M#{idx} {tag} r={R}: D {comb(K,R):,} {tD:.1f}s | I {nI:,} rec streamed+joined {dt:.0f}s pairs={pairs} below_U={below} "
        f"min_digits={len(str(rmin)) if rmin else None} improved={nb} | GPU {gpu_used_gb():.1f}GB")

def main():
    os.makedirs(SCRATCH, exist_ok=True); assert os.path.exists(IDUMP), IDUMP
    if STREAM: assert os.path.exists(IDUMP_STREAM), IDUMP_STREAM
    F = SB.load_facts()
    for j in (62, 63, 65): assert ref.verify_certificate(F[j])[0]
    portfolio = SB.build_portfolio(F)
    log(f"portfolio {len(portfolio)}; r={R} idx {IDX_MIN}..{IDX_MAX} bases={BASES} D-chunks={NCHUNKS_D} I-chunk={I_CHUNK:,} idump_T={IDUMP_T} cap={MAX_IDUMP_GB}GB")
    log(f"incumbent U = {len(str(cur_U()))} digits; GPU used {gpu_used_gb():.1f} GB"); save()
    for idx, (M, dmod, prov) in enumerate(portfolio):
        if idx < IDX_MIN or idx <= state["last_completed_idx"]: continue
        if idx > IDX_MAX: break
        if time.time() - T0 > WALL_CAP: state["status"] = "wall_cap_reached"; save(); log("WALL CAP"); break
        pool = SB.build_pool(M, max(3000, SB.INS_PRIME_CAP))
        if len(pool) < K: state["last_completed_idx"] = idx; save(); continue
        poolset = set(pool); bases = []
        for tag in BASES:
            if tag == "ksmall": bases.append(("ksmall", sorted(pool)[:K]))
            elif tag == "S64":
                base = [p for p in U0_factors if p in poolset]; extra = [p for p in sorted(pool) if p not in set(base)]
                b0 = sorted(base + extra[:K - len(base)]);
                if len(b0) == K: bases.append(("S64", b0))
        seen = set()
        for tag, B0 in bases:
            if tuple(B0) in seen: continue
            seen.add(tuple(B0)); B0set = set(B0); B0s = sorted(B0)
            INS = sorted([p for p in pool if p not in B0set and p < SB.INS_PRIME_CAP][:SB.INS_MAX])
            if len(INS) < R or len(INS) > 256: continue
            P0 = 1
            for p in B0: P0 *= p
            U = cur_U(); topI = 1
            for p in INS[-R:]: topI *= p
            ptb = 1
            for p in B0s[-R:]: ptb *= p
            Icap = max(1, min(topI, (U * ptb) // P0 + 1)); pdm = 1
            for p in B0s[:R]: pdm *= p
            # --- D side on GPU ---
            t = time.time(); Dch = build_D(B0, M); tD = time.time() - t
            # --- I side: STREAM mode pipes the dumper straight into the join (no disk cap) ---
            if STREAM:
                instance_stream(idx, tag, prov, M, B0, B0set, INS, P0, Icap, pdm, Dch, tD)
                del Dch; cp.get_default_memory_pool().free_all_blocks(); continue
            # --- FILE mode (legacy): dump to scratch, then join ---
            prefix = os.path.join(SCRATCH, f"i_{idx}_{tag}"); t = time.time()
            ok, nI, reason = run_idump(M, B0, INS, P0 % M, Icap, pdm, prefix); tI = time.time() - t
            state["instances"] += 1
            if not ok:
                state["skipped_big_I"] += 1; log(f"  M#{idx} {tag} r={R}: SKIP ({reason}) | D {comb(K,R):,} in {tD:.1f}s")
                for f in glob.glob(prefix + ".t*.bin"): os.remove(f)
                del Dch; cp.get_default_memory_pool().free_all_blocks(); continue
            pairs = 0; below = 0; nb = 0; rmin = None; t = time.time(); read_rec = 0
            for f in sorted(glob.glob(prefix + ".t*.bin")):
                mm = np.memmap(f, dtype=np.uint64, mode="r"); nrec = mm.size // 2; read_rec += nrec
                for s in range(0, nrec, I_CHUNK):
                    blk = np.asarray(mm[2 * s: 2 * min(nrec, s + I_CHUNK)]).reshape(-1, 2)
                    Ik = cp.asarray(np.ascontiguousarray(blk[:, 0])); Iid = cp.asarray(np.ascontiguousarray(blk[:, 1])); del blk
                    Ik_s, Iid_s = G.sort_side(Ik, Iid); del Ik, Iid
                    for Dk_s, Did_s in Dch:
                        pd, pi = G.join_sorted(Dk_s, Did_s, Ik_s, Iid_s)
                        if pd.size == 0: continue
                        pairs += int(pd.size)
                        for d, i in zip(G.decode_ids(pd, R, 6), G.decode_ids(pi, R, 8)):
                            D = [B0[x] for x in d]; I = [INS[x] for x in i]
                            prodD = 1
                            for p in D: prodD *= p
                            prodI = 1
                            for p in I: prodI *= p
                            if (P0 * prodI) % prodD: continue
                            n = P0 * prodI // prodD
                            if n >= cur_U(): continue
                            below += 1; fac = sorted((B0set - set(D)) | set(I))
                            if consider((prov, tag, R), fac, M, f"{tag}_r{R}"): nb += 1
                            if rmin is None or n < rmin: rmin = n
                    del Ik_s, Iid_s
                del mm
                cp.get_default_memory_pool().free_all_blocks()
            for f in glob.glob(prefix + ".t*.bin"): os.remove(f)
            del Dch; cp.get_default_memory_pool().free_all_blocks()
            if read_rec != nI:   # integrity: a truncated dump (e.g. disk full) must never pass as a completed instance
                state["skipped_big_I"] += 1
                log(f"  M#{idx} {tag} r={R}: INCOMPLETE -- read {read_rec:,} records but dumper reported {nI:,}; instance NOT counted")
                continue
            state["pairs_total"] += pairs; state["completions_below_U"] += below
            log(f"  M#{idx} {tag} r={R}: D {comb(K,R):,} {tD:.1f}s | I {nI:,} rec {tI:.0f}s | join+verify {time.time()-t:.0f}s pairs={pairs} below_U={below} "
                f"min_digits={len(str(rmin)) if rmin else None} improved={nb} | GPU {gpu_used_gb():.1f}GB")
        state["last_completed_idx"] = idx; save()
        log(f"progress idx {idx}/{IDX_MAX} best_digits={state['best_digits_seen']} improved={state['incumbent'] is not None} skipped={state['skipped_big_I']}")
    if state["status"] == "running": state["status"] = "sweep_complete"
    save(); log("DONE", state["status"], "best_digits_seen:", state["best_digits_seen"], "incumbent improved:", state["incumbent"] is not None)

if __name__ == "__main__": main()
