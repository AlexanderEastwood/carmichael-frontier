#!/usr/bin/env python3
"""OUR independent exclusion, validated against the FROZEN oracle.
- admissibility-pruned descent (pairwise Korselt: if p|(q-1) drop q after choosing p)
- R_min exact-integer product prune
- at EVERY complete leaf, the Korselt verdict is taken from ref_carmichael (frozen oracle),
  NOT from our own arithmetic -> the oracle is the arbiter.
Validates on small k (reproduce least_with_k) before the k=64 run."""
import os,sys,time
# import the frozen oracle from this repo's ref/ (portable from any checkout)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "ref"))
import ref_carmichael as ref
from math import isqrt

def odd_primes(limit):
    s=bytearray([1])*(limit+1); s[0:2]=b"\0\0"
    for p in range(2,isqrt(limit)+1):
        if s[p]: s[p*p::p]=b"\0"*len(s[p*p::p])
    return [p for p in range(3,limit+1) if s[p]]

def enumerate_k(k, Y, limit, use_oracle=True):
    ps=odd_primes(limit); m=len(ps)
    # future[i] = indices j>i admissible with ps[i] (ps[i] does NOT divide ps[j]-1)
    future=[[j for j in range(i+1,m) if (ps[j]-1)%ps[i]] for i in range(m)]
    stats={"nodes":0,"leaves":0,"carm":0,"min":None,"hits":[]}
    sys.setrecursionlimit(m+500)
    def walk(avail, chosen, P, t):
        stats["nodes"]+=1
        if t==0:
            if P>=Y: raise RuntimeError("leaf >= Y")
            stats["leaves"]+=1
            stats["min"]=P if stats["min"] is None else min(stats["min"],P)
            if use_oracle:
                # verify_certificate returns (ok, n, why); the frozen oracle is the arbiter
                res=ref.verify_certificate(list(chosen))
                is_c = res[0] is True
            else:
                n1=P-1; is_c=all(n1%(p-1)==0 for p in chosen)
            if is_c:
                stats["carm"]+=1; stats["hits"].append(P)
            return
        if len(avail)<t: return
        # R_min prune: P * product of t smallest available >= Y  -> none of this suffix works
        # build running min-product over sorted avail (avail kept sorted ascending)
        for idx in range(len(avail)-t+1):
            i=avail[idx]
            # minimal completion using t smallest available from idx on
            prodmin=P*ps[i]
            for jj in avail[idx+1:idx+t]:
                prodmin*=ps[jj]
            if prodmin>=Y: break
            # restrict remaining avail to those admissible with ps[i] AND index> i
            fut=set(future[i])
            newavail=[j for j in avail[idx+1:] if j in fut]
            walk(newavail, chosen+(ps[i],), P*ps[i], t-1)
    walk(list(range(m)), (), 1, k)
    return stats

def validate_small():
    """Reproduce the frozen oracle's brute-force least_with_k on small k.
    Prime limits are kept small so future[] (O(m^2)) and the oracle brute both stay cheap."""
    print("== small-k validation vs frozen oracle least_with_k ==", flush=True)
    ok=True
    for k,bnd,lim in [(3,2000,50),(4,60000,100),(5,10**6,100)]:
        t0=time.time()
        st=enumerate_k(k, bnd, lim, use_oracle=True)
        truth_n,_=ref.least_with_k(k,bnd)
        mymin=min(st["hits"]) if st["hits"] else None
        match=(mymin==truth_n)
        ok=ok and match
        print(f"k={k} bound={bnd} lim={lim}: our min hit={mymin} ; oracle={truth_n} ; match={match} ({time.time()-t0:.1f}s)", flush=True)
    return ok

def run64(Y_exp=145, limit=5518):
    """Oracle-gated exclusion at Y=10**Y_exp over the certified universe (odd primes <= limit).
    limit=5518 is the certified P+ bound at Y=10**145 (see carmichael_prime_bound.py)."""
    t0=time.time(); Y=10**Y_exp
    st=enumerate_k(64, Y, limit, use_oracle=True)
    print(f"[k=64 @10^{Y_exp}, universe odd primes<= {limit}] nodes={st['nodes']} "
          f"admissible_products={st['leaves']} oracle_Carmichael_hits={st['carm']} "
          f"min_product_digits={len(str(st['min'])) if st['min'] else None} ({time.time()-t0:.1f}s)", flush=True)
    print(f"[RESULT] frozen-oracle-checked exclusion -> S_64 >= 10^{Y_exp} : {st['carm']==0}", flush=True)
    return st["carm"]==0

if __name__=="__main__":
    import argparse
    ap=argparse.ArgumentParser(description="Independent, oracle-gated reproduction of S_64 >= 10^145.")
    ap.add_argument("--exponent", type=int, default=145, help="Y = 10**exponent (default 145)")
    ap.add_argument("--limit", type=int, default=5518, help="prime universe upper bound (default 5518, certified at 10^145)")
    ap.add_argument("--skip-validate", action="store_true", help="skip the small-k oracle validation")
    ap.add_argument("--skip-k64", action="store_true", help="skip the k=64 exclusion run")
    a=ap.parse_args()
    if not a.skip_validate:
        validate_small()
    if not a.skip_k64:
        run64(a.exponent, a.limit)
