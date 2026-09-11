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
    "19,23,31,37,41,43", "19,29,31,37,41,43", "19,29,31,41,43,47", "19,31,37,41,43,47", "19,31,41,43,47,53",  # finish4 (depth 6 -> 8)
    "19,31,41,43,47,53,59,61", "19,23,31,37,41,43,53,59", "19,31,37,41,43,47,53,59",                              # finish5 (depth 8 -> 9)
    "19,31,37,41,43,47,59,61", "19,29,31,37,41,43,47,53",
    "19,23,31,37,41,43,53,59,61", "19,29,31,37,41,43,47,53,61", "19,31,37,41,43,47,53,59,61",              # finish6 (depth 9 -> 10)
    "19,31,37,41,43,47,59,61,67", "19,31,41,43,47,53,59,61,67",
    "19,23,31,37,41,43,53,59,61,67", "19,29,31,37,41,43,47,53,61,67", "19,31,37,41,43,47,53,59,61,67",  # finish7 (depth 10 -> 11)
    "19,31,37,41,43,47,59,61,67,71", "19,31,41,43,47,53,59,61,67,71",
}
leaves = set()
# SUBJOBS7 (depth 11, unpruned) was superseded: the 5 depth-10 chain prefixes are covered by the
# feasibility-pruned persistent frontier SUBJOBS8 (depth 14). SUBJOBS7's completed leaves remain on
# disk as extras (still scanned for hits) but are not part of the coverage set.
for lst in ("JOBS", "SUBJOBS", "SUBJOBS2", "SUBJOBS3", "SUBJOBS4", "SUBJOBS5", "SUBJOBS6", "SUBJOBS8"):
    p = os.path.join(JD, lst)
    if os.path.exists(p):
        for line in open(p):
            s = line.strip()
            if s: leaves.add(s)
leaves -= REPLACED
# Donation frontier (SUBJOBS9, depth 14 -> 18): the 5 hard depth-14 slices kept their single-core
# workers AND had idle cores donated to their feasibility-pruned +4 children. Either verdict covers
# the slice: its own n=NONE, or n=NONE on every one of its SUBJOBS9 children.
ALT = {}
for lst, depth in (("SUBJOBS9", 14), ("SUBJOBS10", 18)):      # donated frontiers: parent depth -> children list
    pl = os.path.join(JD, lst)
    if os.path.exists(pl):
        for line in open(pl):
            s = line.strip()
            if not s: continue
            ALT.setdefault(",".join(s.split(",")[:depth]), []).append(s)
# Secondary verdict source: duplicate runs of the same prefix/bound with a different DYNT_KMAX (a
# performance knob under the completeness invariant), written OUTSIDE jobs_1e147 so they never clash
# with a still-running primary worker's open file. A HIT anywhere always wins over NONE.
ALT_DIRS = [os.path.expanduser("~/kmax0_dup")]
def _read_verdict(f):
    if not os.path.exists(f): return "MISSING"
    txt = open(f, errors="replace").read()
    if "\nn=NONE\n" in txt or txt.startswith("n=NONE") or "\nn=NONE" in txt: return "NONE"
    if any(l.startswith("n=") and l[2:3].isdigit() for l in txt.splitlines()): return "HIT"
    return "BAD"
def verdict(pfx):
    name = pfx.replace(",", "_") + ".out"
    vs = [_read_verdict(os.path.join(JD, name))] + [_read_verdict(os.path.join(d, name)) for d in ALT_DIRS]
    if "HIT" in vs: return "HIT"
    if "NONE" in vs: return "NONE"
    return "BAD" if "BAD" in vs else "MISSING"
def covered(pfx):
    """'NONE' if pfx's own verdict is NONE or (recursively) every donated child is covered; 'HIT' if any hit; else a status."""
    v = verdict(pfx)
    if v in ("NONE", "HIT"): return v
    if pfx in ALT:
        cv = [(c, covered(c)) for c in ALT[pfx]]
        if any(x == "HIT" for _, x in cv): return "HIT"
        if all(x == "NONE" for _, x in cv): return "NONE_VIA_CHILDREN"
        return f"{v}; children unfinished {sum(1 for _, x in cv if x not in ('NONE', 'NONE_VIA_CHILDREN'))}/{len(cv)}"
    return v
missing, bad, hits, none, via_children = [], [], [], 0, 0
for pfx in sorted(leaves):
    v = covered(pfx)
    if v == "NONE": none += 1; continue
    if v == "NONE_VIA_CHILDREN": none += 1; via_children += 1; continue
    if v == "HIT": hits.append(pfx); continue
    if pfx in ALT: missing.append(f"{pfx} (own {v})"); continue
    (missing if v == "MISSING" else bad).append(pfx)
# any hit in ANY verdict file (belt and braces, including replaced prefixes' partial runs)
extra_hits = []
for f in glob.iglob(os.path.join(JD, "*.out")):
    for l in open(f, errors="replace"):
        if l.startswith("n=") and l[2:3].isdigit(): extra_hits.append(os.path.basename(f)); break
ok = not missing and not bad and not hits and not extra_hits
rep = [f"leaf_jobs={len(leaves)} verdict_NONE={none} (via_children={via_children}) missing={len(missing)} no_verdict={len(bad)} hits={len(hits)} hits_anywhere={len(extra_hits)}"]
rep += [f"MISSING {m}" for m in missing[:20]] + [f"NOVERDICT {b}" for b in bad[:20]] + [f"HIT {h}" for h in hits + extra_hits]
rep.append("RESULT=ALL_NONE_VERIFIED" if ok else "RESULT=NOT_VERIFIED")
out = "\n".join(rep) + "\n"; print(out, end="")
open(os.path.join(JD, "DONE.verified" if ok else "DONE.not_verified"), "w").write(out)
sys.exit(0 if ok else 1)
