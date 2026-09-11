#!/usr/bin/env python3
"""Coverage-aware final verdict for the partitioned k=64 exclusion below 10^147.

The exclusion was run as a union of prefix jobs that was refined four times
(SUBJOBS -> SUBJOBS2 -> SUBJOBS3 -> SUBJOBS4); each refinement REPLACED a set of
slow prefixes by all their admissible one-or-two-prime extensions (sub_gen.py).
A prefix's children cover its subtree exactly: a 64-subset with product < 10^147
has its (m)-th smallest prime p bounded by p^(65-m) < 10^147, so the next prime
after any prefix we split is far below sub_gen's 6000-prime cap, and every
extension that could carry a Carmichael number is pairwise admissible.

The script-level tallies only looked at existing .out files, which cannot detect
jobs that never ran (that produced a false ALL_NONE at 08:44Z). This verifier
instead requires: every LEAF job (all listed jobs minus the replaced prefixes)
has a verdict file containing exactly 'n=NONE'; no verdict file anywhere reports
a numeric n; and it writes DONE.verified only if both hold.
"""
import os, sys, glob
JD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "jobs_1e147")
REPLACED = {  # prefixes that were sub-partitioned (their subtrees are covered by children)
    "11", "17", "31", "19,23,29", "19,23,31", "19,29,31", "19,31,37", "19,31,41",                        # finish_1e147 (laggards)
    "19,23,31,37", "19,29,31,37", "19,29,31,41", "19,31,37,41", "19,31,41,43",                           # finish2 (quads)
    "19,23,31,37,41", "19,29,31,37,41", "19,29,31,41,43", "19,31,37,41,43", "19,31,41,43,47",             # finish3 (depth 5)
    "19,23,31,37,41,43", "19,29,31,37,41,43", "19,29,31,41,43,47", "19,31,37,41,43,47", "19,31,41,43,47,53",  # finish4 (depth 6)
}
leaves = set()
for lst in ("JOBS", "SUBJOBS", "SUBJOBS2", "SUBJOBS3", "SUBJOBS4"):
    p = os.path.join(JD, lst)
    if os.path.exists(p):
        for line in open(p):
            s = line.strip()
            if s: leaves.add(s)
leaves -= REPLACED
missing, bad, hits, none = [], [], [], 0
for pfx in sorted(leaves):
    f = os.path.join(JD, pfx.replace(",", "_") + ".out")
    if not os.path.exists(f): missing.append(pfx); continue
    txt = open(f, errors="replace").read()
    if "\nn=NONE\n" in txt or txt.startswith("n=NONE") or "\nn=NONE" in txt: none += 1
    elif any(l.startswith("n=") and l[2:3].isdigit() for l in txt.splitlines()): hits.append(pfx)
    else: bad.append(pfx)
# any hit in ANY verdict file (belt and braces, including replaced prefixes' partial runs)
extra_hits = []
for f in glob.iglob(os.path.join(JD, "*.out")):
    for l in open(f, errors="replace"):
        if l.startswith("n=") and l[2:3].isdigit(): extra_hits.append(os.path.basename(f)); break
ok = not missing and not bad and not hits and not extra_hits
rep = [f"leaf_jobs={len(leaves)} verdict_NONE={none} missing={len(missing)} no_verdict={len(bad)} hits={len(hits)} hits_anywhere={len(extra_hits)}"]
rep += [f"MISSING {m}" for m in missing[:20]] + [f"NOVERDICT {b}" for b in bad[:20]] + [f"HIT {h}" for h in hits + extra_hits]
rep.append("RESULT=ALL_NONE_VERIFIED" if ok else "RESULT=NOT_VERIFIED")
out = "\n".join(rep) + "\n"; print(out, end="")
open(os.path.join(JD, "DONE.verified" if ok else "DONE.not_verified"), "w").write(out)
sys.exit(0 if ok else 1)
