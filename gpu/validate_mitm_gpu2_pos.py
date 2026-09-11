#!/usr/bin/env python3
"""POSITIVE-CONTROL validation of gpu/mitm_gpu2.py vs the CPU engine src/mitm64_r7.

Stage-B found the verified 148-digit incumbent N from the 'ksmall' base (the 64
smallest pool primes) at r=7 on N's own modulus. So with B0 = ksmall and an
insertion pool INS that contains N's 7 inserted primes, BOTH engines must
recover N. This exercises the real join path on genuine matches.
Checks: N recovered by CPU and by GPU; every CPU pair present in GPU pairs;
identical minimum completion; every GPU pair has an exactly equal residue.
Run: ~/erdos/.venv/bin/python gpu/validate_mitm_gpu2_pos.py
"""
import os, sys, json, subprocess
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__)); FR = os.path.join(HERE, "..")
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(FR, "ref"))
import ref_carmichael as ref
import mitm_gpu2 as G
import cupy as cp

BIN = os.path.join(FR, "src", "mitm64_r7"); NINS = int(os.environ.get("NINS", "30"))
inc = json.load(open(os.path.join(FR, "results_k64_best_global.json")))
M = int(inc["modulus"]); N = int(inc["n"]); Nf = sorted(int(p) for p in inc["factors"])
def pool(M, cap):
    out = []; q = 3
    while q < cap:
        if M % (q - 1) == 0 and M % q != 0 and ref.is_prime(q): out.append(q)
        q += 2
    return out
P = pool(M, 4000); B0 = sorted(P)[:64]; B0set = set(B0)
D_need = sorted(B0set - set(Nf)); I_need = sorted(set(Nf) - B0set)
R = len(D_need); assert R == len(I_need), (D_need, I_need)
print(f"ksmall base; N needs D={D_need} I={I_need} -> r={R}")
INS = sorted(I_need + [p for p in P if p not in B0set and p not in set(I_need)][:NINS - R])
P0 = 1
for p in B0: P0 *= p
topI = 1
for p in sorted(INS)[-R:]: topI *= p
ptb = 1
for p in B0[-R:]: ptb *= p
Icap = max(1, min(topI, (N * ptb) // P0 + 1)); pdm = 1
for p in B0[:R]: pdm *= p
print(f"M={M} r={R} |INS|={len(INS)} C(64,{R})={comb(64,R):,} C({len(INS)},{R})={comb(len(INS),R):,}")

def completion(D, I):
    prodD = 1
    for p in D: prodD *= p
    prodI = 1
    for p in I: prodI *= p
    return (P0 * prodI // prodD) if (P0 * prodI) % prodD == 0 else None, prodI, prodD

# --- CPU ---
inp = "\n".join([f"{M} {R} {P0 % M} {Icap} {pdm} 8 3000000", "64 " + " ".join(map(str, B0)),
                 f"{len(INS)} " + " ".join(map(str, INS))]) + "\n"
pr = subprocess.run([BIN], input=inp, capture_output=True, text=True, timeout=900)
cpu = {}
for line in pr.stdout.splitlines():
    if line.startswith("MATCH"):
        parts = line.split("|"); D = tuple(sorted(int(x) for x in parts[1].split()[1:])); I = tuple(sorted(int(x) for x in parts[2].split()[1:]))
        cpu[(D, I)] = completion(D, I)[0]
print("CPU:", pr.stderr.strip()); print(f"CPU matches: {len(cpu)}  N recovered: {N in cpu.values()}")

# --- GPU (D in 4 rank-chunks, I unpruned since small) ---
totD = comb(64, R); nch = 4; b = [totD * i // nch for i in range(nch + 1)]
Dk = cp.concatenate([G.generate(B0, R, M, 1, bits=6, lo=a, hi=c)[0] for a, c in zip(b, b[1:])])
Did = cp.concatenate([G.generate(B0, R, M, 1, bits=6, lo=a, hi=c)[1] for a, c in zip(b, b[1:])])
Ik, Iid = G.generate(INS, R, M, P0 % M, bits=8)
Dk_s, Did_s = G.sort_side(Dk, Did); Ik_s, Iid_s = G.sort_side(Ik, Iid)
pd, pi = G.join_sorted(Dk_s, Did_s, Ik_s, Iid_s)
gpu = {}; bad = 0
for d, i in zip(G.decode_ids(pd, R, 6), G.decode_ids(pi, R, 8)):
    D = tuple(sorted(B0[x] for x in d)); I = tuple(sorted(INS[x] for x in i))
    n, prodI, prodD = completion(D, I)
    if (P0 * prodI - prodD) % M != 0: bad += 1
    if prodI < Icap: gpu[(D, I)] = n
best_cpu = min((v for v in cpu.values() if v), default=None); best_gpu = min((v for v in gpu.values() if v), default=None)
missing = set(cpu) - set(gpu)
ok_oracle = ref.verify_certificate(Nf)[0]
print(f"GPU pairs: {len(gpu)} (raw {pd.size})  residue-mismatch: {bad}  N recovered: {N in gpu.values()}")
print(f"CPU pairs missing from GPU: {len(missing)}   min completion equal: {best_cpu == best_gpu} ({None if best_gpu is None else len(str(best_gpu))} digits)")
ok = (bad == 0) and (not missing) and (N in cpu.values()) and (N in gpu.values()) and (best_cpu == best_gpu) and ok_oracle
print("POSITIVE-CONTROL VALIDATION:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
