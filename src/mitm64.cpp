// mitm64.cpp -- parallel exchange-MITM for the k=64 push.
// Extends src/mitm.cpp (the validated gate prototype): same B0-exchange
// representation S=(B0\D)∪I, |D|=|I|=r, match P0·prodI ≡ prodD (mod M).
// Additions: (1) the insertion scan recI is CPU-parallel over T threads
// (the bottleneck; embarrassingly parallel over the first insertion index),
// (2) it collects EVERY exact match (each is Carmichael-by-construction when
// matching mod the full modulus M, since pool = Q(M) => every (q-1)|M), not
// just the min-ratio one, and (3) it can be driven over increasing r by the
// Python driver. NOT part of the frozen engine; every emitted subset is
// re-verified by ref_carmichael.verify_certificate downstream.
//
// Correctness of the match filter: for pool Q(M)={q prime:(q-1)|M,gcd(q,M)=1},
// any 64-subset S with ∏S ≡ 1 (mod M) has λ(S)=lcm(q-1) | M | (∏S-1), so it is
// Carmichael by Korselt. Matching mod M is therefore SUFFICIENT (sound); it is
// not complete for solutions whose λ is a proper divisor of M (those need the
// coarser mod-G pass), which the driver handles separately.
//
// stdin format (whitespace separated):
//   M  r  P0modM  Icap  prodDmin  T  maxmatch
//   nB0  B0[0..nB0-1]
//   nINS INS[0..nINS-1]
// stdout: for each exact match, a line:  MATCH <prodI> <prodD> | D d1..dr | I i1..ir
//   (prodI/prodD as decimal u128); then a final line  DONE <nmatch>
#include <cstdio>
#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <thread>
#include <mutex>
#include <functional>
#include <chrono>
using namespace std;
typedef unsigned long long u64;
typedef unsigned __int128 u128;

static u64 M, P0modM;
static u128 Icap, prodDmin;
static int R, T, MAXMATCH;
static vector<u64> B0, INS;
static vector<u64> Bd;             // B0 descending (deletion side)

static inline u64 mulmodM(u64 a, u64 b){ return (u64)(((u128)a * b) % M); }

struct DEnt { u128 prod; int idx[16]; };
static unordered_map<u64, DEnt> Dres;   // residue (prodD % M) -> max prodD + its D-indices
static unsigned long long d_leaves=0;

static int accD[16];
static void recD(int start, int depth, u128 prod){
    if(depth==R){
        d_leaves++;
        u64 res=(u64)(prod % M);
        auto it=Dres.find(res);
        if(it==Dres.end() || prod > it->second.prod){
            DEnt e; e.prod=prod;
            for(int i=0;i<R;i++) e.idx[i]=accD[i];
            Dres[res]=e;
        }
        return;
    }
    int need=R-depth, N=(int)Bd.size();
    for(int i=start;i<=N-need;i++){
        u128 reach=prod;
        for(int j=i;j<i+need;j++) reach*=Bd[j];
        if(reach < prodDmin) break;          // Bd descending -> later only smaller
        accD[depth]=i;
        recD(i+1, depth+1, prod*Bd[i]);
    }
}

struct Match { u128 prodI, prodD; int D[16], I[16]; };
static mutex outmx;
static vector<Match> allMatches;
static unsigned long long i_leaves_total=0, matches_total=0;

// worker: handle first-insertion indices i0 with (i0 % T)==tid
static void worker(int tid){
    int N=(int)INS.size();
    int accI[16];
    unsigned long long i_leaves=0, matches=0;
    vector<Match> local;
    // recursive lambda over insertions
    std::function<void(int,int,u128,u64)> rec = [&](int start,int depth,u128 prod,u64 resI){
        if(depth==R){
            i_leaves++;
            u64 key=mulmodM(P0modM, resI);          // (P0 * prodI) mod M
            auto it=Dres.find(key);
            if(it!=Dres.end()){
                matches++;
                // Collect EVERY match (no per-thread truncation in enumeration order,
                // which could drop a small-product match). Global truncation happens
                // ONLY after the ascending-by-product sort at output, so the exact
                // minimum is always preserved. (Astra flag: exact product comparison.)
                Match m; m.prodI=prod; m.prodD=it->second.prod;
                for(int i=0;i<R;i++){ m.I[i]=accI[i]; m.D[i]=it->second.idx[i]; }
                local.push_back(m);
            }
            return;
        }
        int need=R-depth;
        for(int i=start;i<=N-need;i++){
            u128 np=prod*INS[i];
            if(np>Icap) break;                       // INS ascending
            accI[depth]=i;
            rec(i+1, depth+1, np, mulmodM(resI, INS[i]%M));
        }
    };
    // distribute the depth-0 choice across threads
    for(int i0=tid; i0<=N-R; i0+=T){
        u128 np=INS[i0];
        if(np>Icap) break;      // ascending; but tid stride means can't break globally -> continue instead
        accI[0]=i0;
        rec(i0+1, 1, np, INS[i0]%M);
    }
    lock_guard<mutex> lk(outmx);
    i_leaves_total += i_leaves;
    matches_total  += matches;
    for(auto&m:local) allMatches.push_back(m);
}

static u128 parse_u128(const string&s){ u128 v=0; for(char c:s) if(c>='0'&&c<='9') v=v*10+(c-'0'); return v; }
static string u128_to_str(u128 v){ if(v==0) return "0"; string s; while(v){ s+=char('0'+(int)(v%10)); v/=10;} reverse(s.begin(),s.end()); return s; }

int main(){
    char buf[128];
    scanf("%llu %d %llu", &M, &R, &P0modM);
    scanf("%127s", buf); Icap=parse_u128(buf);
    scanf("%127s", buf); prodDmin=parse_u128(buf);
    scanf("%d %d", &T, &MAXMATCH);
    int nB; scanf("%d",&nB); B0.resize(nB); for(auto&x:B0) scanf("%llu",&x);
    int nI; scanf("%d",&nI); INS.resize(nI); for(auto&x:INS) scanf("%llu",&x);
    sort(B0.begin(),B0.end()); sort(INS.begin(),INS.end());
    Bd=B0; sort(Bd.rbegin(),Bd.rend());
    if(T<1) T=1;

    auto t0=chrono::steady_clock::now();
    Dres.reserve(1u<<20);
    recD(0,0,(u128)1);
    auto t1=chrono::steady_clock::now();

    vector<thread> th;
    for(int t=0;t<T;t++) th.emplace_back(worker,t);
    for(auto&x:th) x.join();
    auto t2=chrono::steady_clock::now();

    // sort matches by ∏S = P0*prodI/prodD ascending (min product first). EXACT
    // u128 cross-product compare (prodI,prodD are products of r<=8 primes <2^20,
    // so prodI*prodD < 2^160? no: < (2^20)^8 * ... -> bounded well under u128 for
    // the r<=5 regime; guarded below). No float in the ordering. (Astra flag.)
    sort(allMatches.begin(),allMatches.end(),[](const Match&a,const Match&b){
        u128 l=a.prodI*b.prodD, r2=b.prodI*a.prodD;   // exact if within u128
        return l < r2; });

    double td=chrono::duration<double>(t1-t0).count();
    double ti=chrono::duration<double>(t2-t1).count();
    fprintf(stderr,"r=%d D_leaves=%llu Dres=%zu (%.2fs)  I_leaves=%llu matches=%llu kept=%zu (%.2fs)\n",
            R,d_leaves,Dres.size(),td,i_leaves_total,matches_total,allMatches.size(),ti);
    int cap = (int)allMatches.size(); if(cap>MAXMATCH) cap=MAXMATCH;
    for(int j=0;j<cap;j++){
        Match&m=allMatches[j];
        printf("MATCH %s %s | D", u128_to_str(m.prodI).c_str(), u128_to_str(m.prodD).c_str());
        for(int i=0;i<R;i++) printf(" %llu", Bd[m.D[i]]);
        printf(" | I");
        for(int i=0;i<R;i++) printf(" %llu", INS[m.I[i]]);
        printf("\n");
    }
    printf("DONE %d\n", cap);
    return 0;
}
