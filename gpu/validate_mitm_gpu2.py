#!/usr/bin/env python3
"""Bit-exact validation of gpu/mitm_gpu2.py against the CPU engine src/mitm64_r7.

Picks one portfolio modulus M with its base B0 (64 primes) and a SMALL insertion
pool INS (so the unpruned I side is enumerable), runs BOTH:
  CPU: mitm64_r7 with the Stage-B stdin format -> MATCH (D,I) pairs
       (engine dedups D by residue keeping MAX prodD -> its matches are the
        dominant representative per residue).
  GPU: mitm_gpu2 generate(D, bits=6, chunked) + generate(I, bits=8) + join
       (emits EVERY equal-residue pair; a superset).
Checks: (1) every CPU MATCH pair is present in the GPU pair set; (2) the minimum
completion product P0*prodI/prodD over each side's pairs is IDENTICAL; (3) the
GPU pair count is >= the CPU count and every GPU pair has exactly equal residue
(recomputed on host in Python ints). Also exercises rank-chunking by generating
the D side in several chunks and asserting the union equals the single-shot set.
Run with the cupy interpreter:  ~/erdos/.venv/bin/python gpu/validate_mitm_gpu2.py
"""
import os, sys, json, subprocess
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__)); FR = os.path.join(HERE, "..")
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(FR, "ref"))
import ref_carmichael as ref
import mitm_gpu2 as G
import cupy as cp

BIN = os.path.join(FR, "src", "mitm64_r7")
R = int(os.environ.get("R", "6")); NINS = int(os.environ.get("NINS", "22"))

# --- a modulus + base from the verified 148-digit incumbent (same as Stage-B 'S64' base) ---
inc = json.load(open(os.path.join(FR, "results_k64_best_global.json")))
M = int(inc["modulus"]); B0 = sorted(int(p) for p in inc["factors"]); assert len(B0) == 64
def pool(M, cap):
    out = []; q = 3
    while q < cap:
        if M % (q - 1) == 0 and M % q != 0 and ref.is_prime(q): out.append(q)
        q += 2
    return out
INS = [p for p in pool(M, 4000) if p not in set(B0)][:NINS]
P0 = 1
for p in B0: P0 *= p
topI = 1
for p in sorted(INS)[-R:]: topI *= p
ptb = 1
for p in B0[-R:]: ptb *= p
U = int(inc["n"]); Icap = max(1, min(topI, (U * ptb) // P0 + 1))
pdm = 1
for p in B0[:R]: pdm *= p
print(f"M={M} r={R} |INS|={len(INS)} Icap~{Icap.bit_length()}bits  C(64,{R})={comb(64,R):,} C({len(INS)},{R})={comb(len(INS),R):,}")

# --- CPU engine ---
inp = "\n".join([f"{M} {R} {P0 % M} {Icap} {pdm} 8 3000000", f"64 " + " ".join(map(str, B0)),
                 f"{len(INS)} " + " ".join(map(str, INS))]) + "\n"
pr = subprocess.run([BIN], input=inp, capture_output=True, text=True, timeout=900)
cpu_pairs = set()
for line in pr.stdout.splitlines():
    if line.startswith("MATCH"):
        parts = line.split("|"); D = tuple(sorted(int(x) for x in parts[1].split()[1:])); I = tuple(sorted(int(x) for x in parts[2].split()[1:]))
        cpu_pairs.add((D, I))
print("CPU:", pr.stderr.strip()); print(f"CPU matches: {len(cpu_pairs)}")

# --- GPU: D side in chunks (exercise chunking) + I side, join, decode to prime tuples ---
totD = comb(64, R); nchunks = 3; bounds = [totD * i // nchunks for i in range(nchunks + 1)]
Dk_all, Did_all = [], []
for a, b in zip(bounds, bounds[1:]):
    k, i = G.generate(B0, R, M, 1, bits=6, lo=a, hi=b); Dk_all.append(k); Did_all.append(i)
Dk = cp.concatenate(Dk_all); Did = cp.concatenate(Did_all)
k1, i1 = G.generate(B0, R, M, 1, bits=6)               # single-shot for the chunking check
assert bool(cp.all(cp.sort(Dk) == cp.sort(k1)).item()) and Dk.size == totD, "chunk union != single shot"
Ik, Iid = G.generate(INS, R, M, P0 % M, bits=8)
Dk_s, Did_s = G.sort_side(Dk, Did); Ik_s, Iid_s = G.sort_side(Ik, Iid)
pd, pi = G.join_sorted(Dk_s, Did_s, Ik_s, Iid_s)
Dt = G.decode_ids(pd, R, 6); It = G.decode_ids(pi, R, 8)
gpu_pairs = set(); bad_res = 0; best_gpu = None
for d, i in zip(Dt, It):
    D = tuple(sorted(B0[x] for x in d)); I = tuple(sorted(INS[x] for x in i))
    prodD = 1
    for p in D: prodD *= p
    prodI = 1
    for p in I: prodI *= p
    if (P0 * prodI - prodD) % M != 0: bad_res += 1
    if prodI < Icap:                                     # same product cap the CPU applies to I
        gpu_pairs.add((D, I))
        n = P0 * prodI // prodD if (P0 * prodI) % prodD == 0 else None
        if n is not None and (best_gpu is None or n < best_gpu): best_gpu = n
best_cpu = None
for D, I in cpu_pairs:
    prodD = 1
    for p in D: prodD *= p
    prodI = 1
    for p in I: prodI *= p
    n = P0 * prodI // prodD if (P0 * prodI) % prodD == 0 else None
    if n is not None and (best_cpu is None or n < best_cpu): best_cpu = n
missing = cpu_pairs - gpu_pairs
print(f"GPU pairs (Icap-filtered): {len(gpu_pairs)}  raw joined: {len(Dt)}  residue-mismatch: {bad_res}")
print(f"CPU pairs missing from GPU: {len(missing)}")
print(f"min completion  CPU={None if best_cpu is None else len(str(best_cpu))}d  GPU={None if best_gpu is None else len(str(best_gpu))}d  equal={best_cpu == best_gpu}")
ok = (bad_res == 0) and (not missing) and (len(gpu_pairs) >= len(cpu_pairs)) and (best_cpu == best_gpu)
print("VALIDATION:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
