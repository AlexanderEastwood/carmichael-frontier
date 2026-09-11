#!/usr/bin/env python3
"""Coverage check for a divisor-family run (stdlib only). Reads the census file and the records
directory written by run_mstar_family.py and accepts the family only if, for EVERY survivor modulus D:
  * a record D<D>.log exists and its header names the same D, the expected bound U, and the census rmax;
  * the '[runner] r=0' base check is present, with below_U=False;
  * the driver header confirms the same bound U (BOUND U=...) and DFILTER=1;
  * every radius 1..rmax(D) has an instance line ('dfilter' / 'host-exact' / 'nothing to search');
  * no line says INCOMPLETE and the driver ended with 'DONE complete';
  * no 'NEW BEST' (a hit) -- reported separately if present.
usage: python3 verify_family_coverage.py <census.json> <records dir> [<expected bound U>|-] [<radius cap>]
Exit 0 and 'RESULT=FAMILY_COVERED' only if all checks pass."""
import json, os, re, sys

def main() -> None:
    fam_path, rec_dir = sys.argv[1], sys.argv[2]
    fam = json.load(open(fam_path)); fam = fam["survivors"] if isinstance(fam, dict) else fam
    U_exp = int(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] != "-" else None
    rcap = int(sys.argv[4]) if len(sys.argv) > 4 else None      # discovery passes: certify coverage only up to this radius
    problems = []; hits = []; radii_total = 0
    for row in fam:
        D, rmax = int(row["D"]), int(row["rmax"]); f = os.path.join(rec_dir, f"D{D}.log")
        if rcap is not None: rmax = min(rmax, rcap)
        if not os.path.exists(f): problems.append(f"D={D}: record missing"); continue
        txt = open(f, errors="replace").read(); lines = txt.splitlines()
        m = re.search(r"\[runner\] family=(\S+) D=(\d+) U=(\d+) pool=(\d+) cap=(\d+) rmax=(\d+)", txt)
        if not m: problems.append(f"D={D}: no runner header"); continue
        if int(m.group(2)) != D: problems.append(f"D={D}: header D mismatch")
        U_rec = int(m.group(3))
        if U_exp is not None and U_rec != U_exp: problems.append(f"D={D}: bound {str(U_rec)[:12]}... != expected {str(U_exp)[:12]}...")
        if int(m.group(6)) != int(row["rmax"]): problems.append(f"D={D}: header rmax {m.group(6)} != census {row['rmax']}")
        r0 = re.search(r"\[runner\] r=0 .*carmichael=(\w+) below_U=(\w+)", txt)
        if not r0: problems.append(f"D={D}: no radius-0 base check")
        elif r0.group(2) == "True": hits.append(f"D={D}: base itself below U")
        if rmax >= 1:
            b = re.search(r"BOUND U=(\d+) \((\d+) digits\)(.*); DFILTER=(\d)", txt)
            if not b: problems.append(f"D={D}: driver header missing")
            else:
                if int(b.group(1)) != U_rec: problems.append(f"D={D}: driver bound != runner bound")
                if b.group(4) != "1": problems.append(f"D={D}: DFILTER off")
            seen = set()
            for l in lines:
                mm = re.search(r" r=(\d+) ", l)
                if mm and ("dfilter" in l or "host-exact" in l or "nothing to search" in l): seen.add(int(mm.group(1)))
            missing = [r for r in range(1, rmax + 1) if r not in seen]
            if missing: problems.append(f"D={D}: radii without a completed instance: {missing}")
            radii_total += rmax
            if "INCOMPLETE" in txt: problems.append(f"D={D}: INCOMPLETE instance(s)")
            if "DONE complete" not in txt: problems.append(f"D={D}: driver did not end with DONE complete")
        if "NEW BEST" in txt: hits.append(f"D={D}: NEW BEST")
    print(f"family={os.path.basename(fam_path)} moduli={len(fam)} exchange_instances_expected={radii_total}{f' (radii capped at {rcap}: DISCOVERY coverage only)' if rcap is not None else ''} problems={len(problems)} hits={len(hits)}")
    for p in problems[:30]: print("PROBLEM", p)
    for h in hits: print("HIT", h)
    ok = not problems and not hits
    print("RESULT=" + ("FAMILY_COVERED" if ok else ("HIT_FOUND" if hits and not problems else "NOT_COVERED")))
    sys.exit(0 if ok else 1)

if __name__ == "__main__": main()
