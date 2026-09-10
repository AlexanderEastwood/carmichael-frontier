#!/usr/bin/env python3
"""
mitm64_driver.py -- k=64 push driver.

Builds a neighbor-seeded RICH portfolio of candidate moduli M, and for each:
  * pool Q(M) = {q prime: (q-1)|M, gcd(q,M)=1}
  * guaranteed-G reduction (report) + B0 exchange base(s)
  * runs the parallel exchange-MITM (mitm64) mod M over increasing r, with
    product bounds seeded from the running incumbent U
  * every match is Carmichael-by-construction; re-verified by the FROZEN oracle
    ref_carmichael.verify_certificate; the global incumbent U (smallest verified
    64-factor Carmichael) is updated and checkpointed.

Nothing here trusts N_64's own (unknown) factorization. Neighbor facts (62,63,65)
come from ~/carmichael/README.md (all oracle-verified). Checkpoints to
results_k64.json every few minutes; hard wall cap enforced.
"""
import os, re, sys, json, time, subprocess, itertools
from math import isqrt

HOME=os.path.expanduser("~")
FR=os.path.join(HOME,"carmichael-frontier")
REF=os.path.join(FR,"ref"); sys.path.insert(0,REF)
import ref_carmichael as ref
README=os.path.join(HOME,"carmichael","README.md")
MITM=os.path.join(FR,"src","mitm64")
RESULTS=os.path.join(FR,"results_k64.json")
CKPT=os.path.join(FR,"checkpoint_k64.json")

K=64
NPROC=int(os.environ.get("NPROC","32"))
WALL_CAP=float(os.environ.get("WALL_CAP","2100"))   # seconds (~35 min); outer script also caps
R_MAX=int(os.environ.get("R_MAX","6"))
PORTFOLIO_MAX=int(os.environ.get("PORTFOLIO_MAX","250"))
INS_PRIME_CAP=int(os.environ.get("INS_PRIME_CAP","4000"))   # insertion primes bounded < this
INS_MAX=int(os.environ.get("INS_MAX","140"))                # cap insertion-pool size (bound recI)
DRES_CAP=float(os.environ.get("DRES_CAP","2e7"))            # cap C(|B0|,r) => Dres RAM (~1.6GB)
T0=time.time()

def log(*a):
    print(f"[{time.time()-T0:7.1f}s]",*a,flush=True)

# ---------- neighbor facts ----------
def load_facts():
    F={}; curk=None
    for line in open(README):
        m=re.search(r"k\s*=\s*(\d+)", line)
        if m and line.strip().startswith("####"): curk=int(m.group(1)); continue
        if curk is not None and "\\cdot" in line:
            F[curk]=[int(x) for x in re.findall(r"\d+",line)]; curk=None
    return F

def factor(x):
    f={}; d=2
    while d*d<=x:
        while x%d==0: f[d]=f.get(d,0)+1; x//=d
        d+=1 if d==2 else 2
    if x>1: f[x]=f.get(x,0)+1
    return f

def lam_f(primes):
    L={}
    for p in primes:
        for e,c in factor(p-1).items(): L[e]=max(L.get(e,0),c)
    return L

def to_int(d):
    v=1
    for e,c in d.items(): v*=e**c
    return v

CORE_RICH={2,3,5,7,11,13}
SGPOOL=set(e for e in range(11,4000) if ref.is_prime(e) and ref.is_prime(2*e+1))

# ---------- portfolio of candidate moduli ----------
def build_portfolio(F):
    """RICH neighbor-seeded moduli: from each neighbor lambda, perturb CORE exps by
    +-1, drop<=2 non-core mids, add<=3 mids from SG+neighbor-mid pool. Ordered
    best-first by L1 distance to the neighbor lambdas."""
    seeds={j:lam_f(F[j]) for j in (62,63,65)}
    seed_ints=[to_int(s) for s in seeds.values()]
    nmids=set(e for s in seeds.values() for e in s if e not in CORE_RICH)
    addpool=sorted((SGPOOL|nmids))
    cand={}   # Mint -> (dict, provenance)
    def consider(d, prov):
        v=to_int(d)
        if v.bit_length()>62: return           # keep modulus in u64 match range
        if v not in cand: cand[v]=(dict(d),prov)
    for j,base in seeds.items():
        cores=sorted(e for e in base if e in CORE_RICH)
        mids =sorted(e for e in base if e not in CORE_RICH)
        # core exponent perturbations: each core exp -1,0,+1 (bounded product of combos)
        ranges=[]
        for e in cores:
            b=base[e]; opts=[b]
            if b>1: opts=[b-1]+opts
            opts=opts+[b+1]
            ranges.append((e,opts))
        # cap core-combo explosion: only perturb up to 3 cores at once
        base_core={e:base[e] for e in cores}
        # enumerate: choose a subset (<=3) of cores to shift, each by +-1
        from itertools import combinations, product as iproduct
        core_variants=[dict(base_core)]
        for rr in (1,2,3):
            for combo in combinations(cores,rr):
                for signs in iproduct((-1,1),repeat=rr):
                    d=dict(base_core); ok=True
                    for e,s in zip(combo,signs):
                        nv=d[e]+s
                        if nv<1: ok=False; break
                        d[e]=nv
                    if ok: core_variants.append(d)
        # mid operations: drop<=2 mids, add<=3 from addpool (bounded: keep adds small)
        add_choices=[()]
        small_adds=[a for a in addpool if a not in base and a<600]   # keep adds modest
        for aa in (1,2,3):
            for combo in combinations(small_adds[:40],aa):
                add_choices.append(combo)
        drop_choices=[()]
        for dd in (1,2):
            for combo in combinations(mids,dd):
                drop_choices.append(combo)
        for cv in core_variants:
            midmap={e:base[e] for e in mids}
            for drops in drop_choices:
                for adds in add_choices:
                    d=dict(cv)
                    for e in midmap:
                        if e not in drops: d[e]=midmap[e]
                    for a in adds: d[a]=1
                    consider(d, (j,tuple(sorted(drops)),tuple(sorted(adds))))
    # rank by min L1 distance (exponent vector) to any seed lambda
    def l1(d):
        best=10**9
        for s in seeds.values():
            keys=set(d)|set(s)
            best=min(best, sum(abs(d.get(e,0)-s.get(e,0)) for e in keys))
        return best
    ranked=sorted(cand.items(), key=lambda kv:(l1(kv[1][0]), kv[0]))
    return [(v,d,prov) for v,(d,prov) in ranked][:PORTFOLIO_MAX]

# ---------- pool Q(M) ----------
def build_pool(M, qcap):
    """primes q<qcap with (q-1)|M and gcd(q,M)=1. qcap chosen from M's structure."""
    pool=[]
    q=3
    while q<qcap:
        if M % (q-1)==0 and M % q!=0 and ref.is_prime(q):
            pool.append(q)
        q+=2
    return pool

def guaranteed_G(pool, k):
    """G = prod_ell ell^{(r+1)-th largest ell-valuation across pool}, r=|pool|-k.
    Divides lambda(S) for every k-subset S of pool (necessary-condition modulus)."""
    r=len(pool)-k
    if r<0: return None
    # valuations
    vals={}
    for q in pool:
        for e,c in factor(q-1).items(): vals.setdefault(e,[]).append(c)
    G={}
    for e,lst in vals.items():
        lst.sort(reverse=True)
        # (r+1)-th largest; if fewer than r+1 entries, that prime can be fully avoided -> exp 0
        if len(lst)>=r+1:
            g=lst[r]
            if g>0: G[e]=g
    return G

# ---------- run one MITM call ----------
def run_mitm(M, r, B0, INS, P0modM, Icap, prodDmin, maxmatch=200000):
    inp=[]
    inp.append(f"{M} {r} {P0modM} {Icap} {prodDmin} {NPROC} {maxmatch}")
    inp.append(f"{len(B0)} "+" ".join(map(str,B0)))
    inp.append(f"{len(INS)} "+" ".join(map(str,INS)))
    data="\n".join(inp)+"\n"
    try:
        pr=subprocess.run([MITM],input=data,capture_output=True,text=True,timeout=600)
    except subprocess.TimeoutExpired:
        return [], "timeout"
    matches=[]
    for line in pr.stdout.splitlines():
        if line.startswith("MATCH"):
            # MATCH prodI prodD | D d.. | I i..
            parts=line.split("|")
            head=parts[0].split()
            D=[int(x) for x in parts[1].split()[1:]]
            I=[int(x) for x in parts[2].split()[1:]]
            matches.append((D,I))
    return matches, pr.stderr.strip()

# ---------- incumbent state ----------
state={"incumbent":None,"incumbent_factors":None,"incumbent_modulus":None,
       "moduli_tried":0,"matches_found":0,"portfolio_size":0,"pool_stats":[],
       "r_reached":0,"status":"running","started":time.strftime("%Y-%m-%d %H:%M:%S")}

def save():
    state["elapsed_s"]=round(time.time()-T0,1)
    tmp=CKPT+".tmp"
    json.dump(state,open(tmp,"w"),indent=1); os.replace(tmp,CKPT)

def consider_candidate(prov, factors, M):
    ok,n,why=ref.verify_certificate(sorted(factors))
    if not (ok and len(factors)==K):
        log("  REJECTED by oracle:",why,"nfac",len(factors)); return
    state["matches_found"]+=1
    if state["incumbent"] is None or n<state["incumbent"]:
        state["incumbent"]=n
        state["incumbent_factors"]=sorted(factors)
        state["incumbent_modulus"]=str(M)
        state["incumbent_prov"]=str(prov)
        state["incumbent_digits"]=len(str(n))
        log(f"  *** VERIFIED k=64 CARMICHAEL, {len(str(n))} digits, modulus prov {prov} ***")
        json.dump({"n":str(n),"factors":sorted(factors),"modulus":str(M),
                   "provenance":str(prov),"digits":len(str(n)),
                   "when":time.strftime("%Y-%m-%d %H:%M:%S")},
                  open(RESULTS,"w"),indent=1)
    save()

def main():
    assert os.path.exists(MITM), "build mitm64 first"
    F=load_facts()
    for j in (62,63,65):
        ok,_,_=ref.verify_certificate(F[j]); assert ok, f"neighbor {j} failed oracle"
    portfolio=build_portfolio(F)
    state["portfolio_size"]=len(portfolio)
    log(f"portfolio: {len(portfolio)} moduli (best-first)")
    save()

    # neighbor prime universe (for B0 / INS construction)
    Nfac={j:set(F[j]) for j in (62,63,65)}

    for idx,(M,dmod,prov) in enumerate(portfolio):
        if time.time()-T0>WALL_CAP:
            state["status"]="wall_cap_reached"; save(); log("WALL CAP hit"); break
        # choose qcap: cover neighbor primes generously
        qcap=max(3000, INS_PRIME_CAP)
        pool=build_pool(M,qcap)
        if len(pool)<K:
            continue
        G=guaranteed_G(pool,K)
        Gi=to_int(G) if G else 1
        state["moduli_tried"]+=1
        # Build candidate B0 sets (64-subsets close to neighbors), all subset of pool:
        poolset=set(pool)
        b0_variants=[]
        # (a) k-smallest of pool
        b0_variants.append(("ksmall", sorted(pool)[:K]))
        # (b) neighbor-derived: N63 factors that are in pool, padded/trimmed to 64 by size
        for j in (63,65,62):
            base=[p for p in F[j] if p in poolset]
            if len(base)>=K:
                b0=sorted(base)[:K]
            else:
                extra=[p for p in sorted(pool) if p not in set(base)]
                b0=sorted(base+extra[:K-len(base)])
            if len(b0)==K:
                b0_variants.append((f"N{j}",b0))
        # insertion pool: pool primes not in B0, bounded < INS_PRIME_CAP
        pstat={"prov":str(prov),"M":str(M),"Gdigits":len(str(Gi)),
               "pool":len(pool),"poolmax":pool[-1]}
        if idx<12 or idx%25==0:
            state["pool_stats"].append(pstat)

        for tag,B0 in b0_variants:
            if time.time()-T0>WALL_CAP:
                state["status"]="wall_cap_reached"; save(); break
            B0set=set(B0)
            INS=[p for p in pool if p not in B0set and p<INS_PRIME_CAP][:INS_MAX]
            if not INS: continue
            # P0 = product of B0; P0modM
            P0=1
            for p in B0: P0*=p
            P0modM=P0 % M
            # incumbent-seeded product bounds:
            #   ∏S = P0 * prodI / prodD  ; want ∏S < U (incumbent) if we have one
            # Icap: cap prodI so that even with prodD=1, generation stays bounded.
            # Use a generous static cap tied to insertion-prime range (r primes < INS_PRIME_CAP).
            from math import comb
            for r in range(1, R_MAX+1):
                if time.time()-T0>WALL_CAP: break
                if comb(len(B0),r) > DRES_CAP: break   # RAM guard on Dres side
                # product bounds
                # prodI <= (largest r INS primes) ; prodDmin: reach prune floor.
                # Set Icap so ∏S can still beat incumbent OR (no incumbent) a loose ceiling.
                topI=sorted(INS)[-r:] if len(INS)>=r else INS
                Icap = 1
                for p in topI: Icap*=p
                # prodDmin: require deletions worth at least the smallest r B0 primes' product
                smallD=sorted(B0)[:r]
                prodDmin=1
                for p in smallD: prodDmin*=p
                # if we have an incumbent U, tighten: need P0*prodI/prodD < U => prodD > P0*prodI/U
                # (applied per-candidate in python after match; here keep generation bounded)
                matches,stderr=run_mitm(M,r,B0,INS,P0modM,Icap,prodDmin)
                state["r_reached"]=max(state["r_reached"],r)
                if matches:
                    log(f"  M#{idx} {tag} r={r}: {len(matches)} raw matches | {stderr}")
                    for D,I in matches:
                        factors=sorted((B0set-set(D))|set(I))
                        if len(factors)!=K: continue
                        consider_candidate((prov,tag,r),factors,M)
        if idx%10==0:
            save(); log(f"progress: {idx+1}/{len(portfolio)} moduli, incumbent={state['incumbent'] is not None}")

    if state["status"]=="running":
        state["status"]="portfolio_exhausted"
    save()
    log("DONE status=",state["status"],"incumbent found:",state["incumbent"] is not None)

if __name__=="__main__":
    main()
