#!/usr/bin/env python3
"""Interval table  C_k(97) <= A_k <= S_k <= U_k  for k = 64..144.

Lower bound (certified, relaxed): one shared pass over every internally
admissible subset S of the odd primes <= y=97, extending each with its cheapest
individually-compatible tail primes; C_m(y) = min over S of prod(S)*prod_{i<=m-|S|} b_i(S).
Ignoring conflicts among the b_i only LOWERS the product, so C_m(y) <= A_m <= S_m.
Exact Python ints throughout; insufficient tail table raises (never silently excludes).
Validated: C_63 must equal the known 142-digit constant; class count 407,760.

Upper bound: our on-file incumbents U_k (results_k*.json; k=65 from the verified
Webster candidate), each re-verified here by the FROZEN oracle (distinct primes,
Korselt, product == n, exactly k factors).

Stronger certified Carmichael-exclusion lower bounds are layered on where a
complete exclusion FINISHED (k=64: none <= 10^146; k=65: none < 10^148).
Admissible and Carmichael bounds are kept separate throughout.
"""
import json, os, sys, time
from math import isqrt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ref"))
import ref_carmichael as ref

Y = 97; MMIN = 63; MMAX = 144; TAIL_LIMIT = 10_000
EXPECTED_C63 = int("1812140595092971804464959059612576188399713525008997620073689836850955613086867158799858934921557347504370550240044086773665823716701403213549")
EXPECTED_CLASSES = 407760
# completed Carmichael exclusions: k -> (bound, strict) ; strict=True means "none <= bound"
EXCLUSIONS = {64: (10**146, True), 65: (10**148, False)}   # 65: none < 10^148

def odd_primes(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b"\0\0"
    for i in range(2, isqrt(n) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [p for p in range(3, n + 1) if s[p]]

def relaxed_bounds():
    ps = odd_primes(TAIL_LIMIT)
    head = [p for p in ps if p <= Y]; tail = [p for p in ps if p > Y]
    allow = [sum(1 << i for i, q in enumerate(tail) if (q - 1) % p) for p in head]
    C = [None] * (MMAX + 1); stats = {"classes": 0, "max_tail_prime_used": 0}
    def upd(m, v):
        if C[m] is None or v < C[m]: C[m] = v
    def visit(start, chosen, P, mask):
        stats["classes"] += 1
        n0 = len(chosen); need_max = MMAX - n0
        if n0 >= MMIN: upd(n0, P)
        value, bits, cnt = P, mask, 0
        while cnt < need_max:
            if not bits: raise RuntimeError(f"insufficient tail table for class {chosen}")
            bit = bits & -bits; bits -= bit
            q = tail[bit.bit_length() - 1]; value *= q; cnt += 1
            if n0 + cnt >= MMIN:
                upd(n0 + cnt, value)
                if q > stats["max_tail_prime_used"]: stats["max_tail_prime_used"] = q
        for i in range(start, len(head)):
            q = head[i]
            if all((q - 1) % p for p in chosen):          # pairwise admissible: p !| q-1
                visit(i + 1, chosen + (q,), P * q, mask & allow[i])
    visit(0, (), 1, (1 << len(tail)) - 1)
    return C, stats

def load_U(k, root):
    # k=66..144 incumbents are published under incumbents/ (exact n + factor lists) so every
    # upper bound is independently checkable from the public repository; k=64/65 have named files.
    f = os.path.join(root, "incumbents", f"results_k{k}.json")
    if not os.path.exists(f): f = os.path.join(root, f"results_k{k}.json")
    if k == 64: f = os.path.join(root, "results_k64_best_global.json")   # the current incumbent N'
    if k == 65: f = os.path.join(root, "results_k65_webster_verified.json")
    if not os.path.exists(f): return None
    d = json.load(open(f)); n = int(d["n"]); facs = [int(x) for x in d["factors"]]
    prod = 1
    for p in facs: prod *= p
    ok = ref.verify_certificate(facs)
    return {"n": n, "factors": facs, "digits": len(str(n)), "nfactors": len(facs),
            "product_eq_n": prod == n, "oracle_ok": bool(ok[0]) and len(facs) == k,
            "largest_prime": max(facs), "source": os.path.relpath(f, root)}

def main():
    root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    t0 = time.time(); C, st = relaxed_bounds(); tC = time.time() - t0
    assert st["classes"] == EXPECTED_CLASSES, st
    assert C[63] == EXPECTED_C63, "C_63 regression mismatch"
    prim = 1; primorial_digits = {}
    for i, p in enumerate(odd_primes(2000), start=1):
        prim *= p
        if MMIN <= i <= MMAX: primorial_digits[i] = len(str(prim))
    rows = []
    for k in range(64, MMAX + 1):
        U = load_U(k, root); Ck = C[k]
        assert Ck is not None, f"C_{k} not populated"
        low_int, low_src = Ck, f"C_{k}(97)"
        if k in EXCLUSIONS:
            b, strict = EXCLUSIONS[k]; e = len(str(b)) - 1          # b == 10^e
            ex_low = b + 1 if strict else b                          # none<=10^e -> S>=10^e+1 ; none<10^e -> S>=10^e
            if ex_low > low_int:
                low_int, low_src = ex_low, (f"none ≤ 10^{e}" if strict else f"none < 10^{e}")
        rows.append({"k": k, "C_k97": str(Ck), "C_k97_digits": len(str(Ck)),
                     "primorial_digits": primorial_digits[k],
                     "lower_bound_used": str(low_int), "lower_digits": len(str(low_int)), "lower_source": low_src,
                     "U_k_digits": U["digits"] if U else None, "U_k_oracle_ok": U["oracle_ok"] if U else None,
                     "U_k_product_eq_n": U["product_eq_n"] if U else None, "U_k_largest_prime": U["largest_prime"] if U else None,
                     "U_k_source": U["source"] if U else None,
                     "U_k": str(U["n"]) if U else None,                 # exact upper endpoint (re-derive from factors)
                     "U_k_factors": U["factors"] if U else None,       # the certificate for the upper bound
                     "interval_digits": [len(str(low_int)), U["digits"]] if U else None})
    man = {"generated": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "y": Y, "tail_limit": TAIL_LIMIT,
           "head_classes": st["classes"], "max_tail_prime_used": st["max_tail_prime_used"], "relaxed_pass_seconds": round(tC, 2),
           "C63_check": "OK", "exclusions": {k: str(v[0]) for k, v in EXCLUSIONS.items()}, "rows": rows}
    json.dump(man, open(os.path.join(root, "ladder", "ladder_manifest.json"), "w"), indent=1)
    # markdown
    L = ["# Certified bounds for the least Carmichael numbers with 64–144 prime factors", "",
         "Certified lower and upper bounds C_k(97) ≤ A_k ≤ S_k ≤ U_k with exact endpoints (see `ladder_manifest.json`) and the quantified gain over the identified previous bound (the elementary primorial baseline). No priority is claimed; the bounds stand on the computations.", "",
         f"Certified relaxed lower bound C_k(97) from one shared pass over {st['classes']:,} head classes (y=97, tail primes ≤ {TAIL_LIMIT:,}, largest tail prime used {st['max_tail_prime_used']}), {tC:.1f} s; C_63 reproduces the known 142-digit constant. "
         "Every U_k re-verified by the frozen oracle. Stronger rows use a **completed** Carmichael exclusion. Admissible (C/A) and Carmichael (S) bounds are never conflated.", "",
         "| k | primorial digits | C_k(97) digits | gain | lower bound used | S_k digits ∈ | U_k digits | U_k largest p | U_k oracle |",
         "|--:|--:|--:|--:|:--|:--:|--:|--:|:--:|"]
    for r in rows:
        iv = f"[{r['interval_digits'][0]}, {r['interval_digits'][1]}]" if r["interval_digits"] else "—"
        gain = r["C_k97_digits"] - r["primorial_digits"]
        L.append(f"| {r['k']} | {r['primorial_digits']} | {r['C_k97_digits']} | +{gain} | {r['lower_source']} | **{iv}** | {r['U_k_digits'] or '—'} | {r['U_k_largest_prime'] or '—'} | {'✓' if r['U_k_oracle_ok'] else ('✗' if r['U_k_oracle_ok'] is False else '—')} |")
    open(os.path.join(root, "ladder", "LADDER_TABLE.md"), "w").write("\n".join(L) + "\n")
    bad = [r["k"] for r in rows if r["U_k_oracle_ok"] is False]
    print(f"classes={st['classes']} C63=OK pass={tC:.1f}s max_tail={st['max_tail_prime_used']} rows={len(rows)} oracle_failures={bad}")
    for k in (64, 65, 100, 137, 144):
        r = rows[k - 64]; print(f"  k={k}: C_k(97) {r['C_k97_digits']}d ({r['C_k97'][:6]}…), lower={r['lower_source']} -> S_k digits {r['interval_digits']}")

if __name__ == "__main__": main()
