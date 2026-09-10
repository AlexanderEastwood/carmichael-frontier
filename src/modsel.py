#!/usr/bin/env python3
"""
modsel.py -- Modulus-SELECTION de-risking for the k=64 MITM (Astra step-1).

Question: can a modulus/pool-selection rule seeded ONLY from KNOWN NEIGHBORS
(j in {k-2,k-1,k+1,k+2}) recover the correct search modulus lambda(N_k) WITHOUT
reading N_k's own factorization?

Why lambda(N_k) is the target: for modulus M, pool Q(M)={q prime: q-1|M, gcd(q,M)=1}
guarantees every k-subset S<=Q(M) with prod(S)==1 (mod M) is Carmichael (lambda(S)|M|(n-1)
=> Korselt). For M=lambda(N_k) the minimum such subset is exactly N_k (gate-confirmed by
mitm.cpp / gen.py). So modulus SELECTION reduces to: is lambda(N_k) in the neighbor-seeded
portfolio? -- a purely number-theoretic membership test, verified here.

Neighbor factorizations are read from ~/carmichael/README.md (oracle-verifiable);
this file NEVER uses N_k's own factorization in the RULE (only to score recovery).
Does NOT touch carmichael.cpp / carmichael_par.cpp. Verify anything Carmichael via ref oracle.
"""
import os, re, sys, itertools
REF=os.path.expanduser("~/carmichael-frontier/ref"); sys.path.insert(0,REF)
import ref_carmichael as ref
README=os.path.expanduser("~/carmichael/README.md")

def load_facts():
    F={}; curk=None
    for line in open(README):
        m=re.search(r'#+\s*\$k\s*=\s*(\d+)\$', line)
        if m: curk=int(m.group(1)); continue
        if curk is not None and '\\cdot' in line:
            F[curk]=[int(x) for x in re.findall(r'\d+',line)]; curk=None
    return F
def lam(primes):
    L={}
    for p in primes:
        x=p-1; d=2; f={}
        while d*d<=x:
            while x%d==0: f[d]=f.get(d,0)+1; x//=d
            d+=1
        if x>1: f[x]=f.get(x,0)+1
        for e,c in f.items(): L[e]=max(L.get(e,0),c)
    return L
def to_int(d):
    v=1
    for e,c in d.items(): v*=e**c
    return v

CORE_SIMPLE={2,3,5,7}
CORE_RICH={2,3,5,7,11,13}
SGPOOL=set(e for e in range(11,1000) if ref.is_prime(e) and ref.is_prime(2*e+1))  # safe-prime ell

def feasible(Lk,S,addpool,core,cb=1,db=2,ab=3):
    """True iff seed S reaches Lk by: |dexp|<=cb on each CORE prime; drop<=db existing
       non-core primes; add<=ab non-core primes drawn from addpool."""
    drops=0; adds=[]
    for e in set(Lk)|set(S):
        a=S.get(e,0); b=Lk.get(e,0)
        if e in core:
            if abs(a-b)>cb: return (False,None)
        else:
            if a==b: continue
            if b==0: drops+=1
            elif a==0:
                if e not in addpool: return (False,f"need mid-prime {e}")
                adds.append(e)
            else: return (False,f"noncore exp change {e}")
    if drops>db or len(adds)>ab: return (False,"over budget")
    return (True,f"core-shifts + drop={drops} add={sorted(adds)}")

def analyze(F, targets=(33,39,45)):
    print("="*70); print("BACKTEST: neighbor-seeded modulus recovery at cheap k"); print("="*70)
    for k in targets:
        Lk=lam(F[k]); Lki=to_int(Lk)
        ok,n,_=ref.verify_certificate(F[k]); assert ok and n%Lki==(1%Lki)  # lambda|n-1 sanity
        print(f"\n--- k={k}: lambda(N_k)={Lki}  ={dict(sorted(Lk.items()))} ---")
        # (1) simple rule: CORE_SIMPLE +-1, drop<=1, add<=1 from SG+seed mids
        # (2) rich rule:   CORE_RICH   +-1, drop<=2, add<=3 from SG+union(neighbor mids)
        nmids=set(e for j in (k-2,k-1,k+1,k+2) if j in F for e in lam(F[j]) if e not in CORE_RICH)
        for label,core,pool,cb,db,ab in [
            ("SIMPLE(+-1,drop1,add1)", CORE_SIMPLE, SGPOOL|set(), 1,1,1),
            ("RICH(+-1,drop2,add3,SG+nbr)", CORE_RICH, SGPOOL|nmids, 1,2,3)]:
            hit=None
            for j in (k-2,k-1,k+1,k+2):
                if j in F:
                    f,why=feasible(Lk,lam(F[j]),pool,core,cb,db,ab)
                    if f: hit=(j,why); break
            print(f"    {label:32s}: "+(f"RECOVERED via lambda(N_{hit[0]})  [{hit[1]}]" if hit else "NOT recovered"))

def table_scan(F,core,pool_sg,cb,db,ab,lo=20,hi=63):
    rec=0; tot=0; hard=[]
    for k in range(lo,hi+1):
        if k not in F: continue
        tot+=1; Lk=lam(F[k])
        nmids=set(e for j in (k-2,k-1,k+1,k+2) if j in F for e in lam(F[j]) if e not in core)
        pool=pool_sg|nmids
        got=any(feasible(Lk,lam(F[j]),pool,core,cb,db,ab)[0]
                for j in (k-2,k-1,k+1,k+2) if j in F)
        if got: rec+=1
        else: hard.append(k)
    return rec,tot,hard

if __name__=="__main__":
    F=load_facts()
    analyze(F)
    print("\n"+"="*70); print("FULL-TABLE CALIBRATION (k=20..63)"); print("="*70)
    r,t,h=table_scan(F,CORE_SIMPLE,set(),1,1,1)
    print(f"  SIMPLE rule (+-1 on 2/3/5/7, drop<=1, add<=1 SG): recovers {r}/{t}; misses {h}")
    r,t,h=table_scan(F,CORE_RICH,SGPOOL,1,2,3)
    print(f"  RICH rule (+-1 on 2/3/5/7/11/13, drop<=2, add<=3 SG+nbr): recovers {r}/{t}; residual {h}")
    print("\n  Residual failure mode = a non-safe-prime mid factor p=m*ell+1 (m>2), e.g.")
    print("    k=52 needs ell=107 from 643=6*107+1 (2*107+1=215 composite) -> not in any bounded neighbor pool.")
    # k=64 forward
    print("\n"+"="*70); print("FORWARD: k=64 seed geometry (neighbors 62,63,65 known; 66 unknown)"); print("="*70)
    for j in (62,63,65):
        print(f"  lambda(N_{j}) = {dict(sorted(lam(F[j]).items()))}")
    L63,L65=lam(F[63]),lam(F[65])
    d=sum(abs(L63.get(e,0)-L65.get(e,0)) for e in set(L63)|set(L65))
    print(f"  L1(lambda N63, lambda N65) = {d}  (very close; shared backbone incl. mid-prime 83)")
    print("  => richer neighbor-seed portfolio from 62/63/65 is bounded and likely contains lambda(N_64),")
    print("     UNLESS N_64 introduces a new non-safe-prime mid factor (the k=52 residual mode).")
