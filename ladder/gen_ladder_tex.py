#!/usr/bin/env python3
"""Render paper/ladder_table.tex (longtable) from ladder/ladder_manifest.json.
Never hand-edit the .tex table: regenerate it from the manifest."""
import json, os, re
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
m = json.load(open(os.path.join(root, "ladder", "ladder_manifest.json")))

def src_tex(s):
    if s.startswith("C_"):
        mk = re.search(r"C_(\d+)", s); assert mk, s
        return f"$C_{{{mk.group(1)}}}(97)$"
    mm = re.match(r"none (≤|<) 10\^(\d+)", s); assert mm, f"unrecognised lower_source: {s}"
    op = r"\le" if mm.group(1) == "≤" else "<"
    return f"none ${op}10^{{{mm.group(2)}}}$"

hdr = (r"\toprule $k$ & primorial & $C_k(97)$ & gain & lower endpoint & digits of $S_k$ & $U_k$ & $\Pp(U_k)$\\ \midrule")
L = [r"\begin{longtable}{r r r r l c r r}",
     r"\caption{Certified digit bounds $C_k(97)\le A_k\le S_k\le U_k$, $k=64,\dots,144$. ``primorial'' = digits of the product of the first $k$ odd primes (the elementary baseline); ``gain'' = digits by which $C_k(97)$ exceeds it; ``lower endpoint'' = the certified lower bound actually used (a completed Carmichael exclusion where one exists, else $C_k(97)$); every $U_k$ re-verified by the frozen oracle. Exact endpoints are in \code{ladder/ladder\_manifest.json}.}\label{tab:ladder}\\",
     hdr + r" \endfirsthead", hdr + r" \endhead", r"\bottomrule \endfoot"]
for r in m["rows"]:
    lo, hi = r["interval_digits"]; gain = r["C_k97_digits"] - r["primorial_digits"]
    L.append(f"{r['k']} & {r['primorial_digits']} & {r['C_k97_digits']} & +{gain} & {src_tex(r['lower_source'])} & $[{lo},{hi}]$ & {r['U_k_digits']} & {r['U_k_largest_prime']}\\\\")
L.append(r"\end{longtable}")
out = os.path.join(root, "paper", "ladder_table.tex")
open(out, "w").write("\n".join(L) + "\n")
print(f"wrote {out}: {len(m['rows'])} rows, gain range +{min(r['C_k97_digits']-r['primorial_digits'] for r in m['rows'])}..+{max(r['C_k97_digits']-r['primorial_digits'] for r in m['rows'])}")
