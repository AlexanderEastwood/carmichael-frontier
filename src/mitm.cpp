// mitm.cpp -- Task B prototype: exchange-search MITM with EXACT modulus supplied.
// Finds the min-product k-subset of Q(lambda) with product == 1 (mod L), via
// B0-exchange at fixed r: max(prodD) per residue over deletions D<=B0, matched
// against insertions I (product <= Icap). Hot loops use unsigned __int128
// (products) and uint64 residues (mod L < 2^51). NOT part of the frozen engine;
// every Carmichael claim is re-verified by the Python oracle downstream.
#include <cstdio>
#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <chrono>
using namespace std;
typedef unsigned long long u64;
typedef unsigned __int128 u128;

static u64 L, P0modL;
static u128 Icap, prodDmin;
static int R;
static vector<u64> B0, INS;   // B0 ascending; INS ascending (bounded insertion pool)

static inline u64 mulmodL(u64 a, u64 b){ return (u64)(( (u128)a * b) % L); }

struct DEnt { u128 prod; int idx[8]; };
static unordered_map<u64, DEnt> Dres;   // residue (prodD % L) -> max prodD + D primes
static unsigned long long d_leaves=0;

// deletions: choose R from B0, descending, prod >= prodDmin (reach prune)
static vector<u64> Bd;   // B0 descending
static int accD[8];
static void recD(int start, int depth, u128 prod){
    if(depth==R){
        d_leaves++;
        u64 res = (u64)(prod % L);
        auto it = Dres.find(res);
        if(it==Dres.end() || prod > it->second.prod){
            DEnt e; e.prod=prod;
            for(int i=0;i<R;i++) e.idx[i]=accD[i];
            Dres[res]=e;
        }
        return;
    }
    int need = R-depth;
    int N=(int)Bd.size();
    for(int i=start;i<=N-need;i++){
        // reach: prod * product of the `need` largest available (Bd[i..i+need-1])
        u128 reach=prod;
        for(int j=i;j<i+need;j++) reach*=Bd[j];
        if(reach < prodDmin) break;   // descending -> only smaller later
        accD[depth]=i;
        recD(i+1, depth+1, prod*Bd[i]);
    }
}

static long double best_ratio=1e300L;
static int bestI[8], bestD[8];
static bool found=false;
static unsigned long long i_leaves=0, matches=0;
static int accI[8];
static void recI(int start, int depth, u128 prod, u64 resI){
    if(depth==R){
        i_leaves++;
        u64 key = mulmodL(P0modL, resI);   // (P0 * prodI) mod L
        auto it = Dres.find(key);
        if(it!=Dres.end()){
            matches++;
            long double ratio = (long double)prod / (long double)it->second.prod;
            if(ratio < best_ratio){
                best_ratio=ratio; found=true;
                for(int i=0;i<R;i++){ bestI[i]=accI[i]; bestD[i]=it->second.idx[i]; }
            }
        }
        return;
    }
    int need=R-depth;
    int N=(int)INS.size();
    for(int i=start;i<=N-need;i++){
        u128 np = prod * INS[i];
        if(np > Icap) break;   // ascending
        accI[depth]=i;
        recI(i+1, depth+1, np, mulmodL(resI, INS[i]%L));
    }
}

static u128 parse_u128(const string&s){ u128 v=0; for(char c:s) if(c>='0'&&c<='9') v=v*10+(c-'0'); return v; }

int main(){
    // params from stdin: L R P0modL Icap prodDmin
    //   nB0 then B0[]  ; nINS then INS[]
    string tok;
    scanf("%llu %d %llu", &L, &R, &P0modL);
    { char buf[64]; scanf("%s", buf); Icap=parse_u128(buf); scanf("%s", buf); prodDmin=parse_u128(buf); }
    int nB; scanf("%d",&nB); B0.resize(nB); for(auto&x:B0) scanf("%llu",&x);
    int nI; scanf("%d",&nI); INS.resize(nI); for(auto&x:INS) scanf("%llu",&x);
    sort(B0.begin(),B0.end()); sort(INS.begin(),INS.end());
    Bd=B0; sort(Bd.rbegin(),Bd.rend());

    auto t0=chrono::steady_clock::now();
    Dres.reserve(6000000);
    recD(0,0,(u128)1);
    auto t1=chrono::steady_clock::now();
    recI(0,0,(u128)1,1);
    auto t2=chrono::steady_clock::now();

    double td=chrono::duration<double>(t1-t0).count();
    double ti=chrono::duration<double>(t2-t1).count();
    fprintf(stderr,"D_leaves=%llu Dres_size=%zu (%.2fs)  I_leaves=%llu matches=%llu (%.2fs)\n",
            d_leaves, Dres.size(), td, i_leaves, matches, ti);
    if(found){
        // print D primes and I primes (space separated), two lines
        printf("D");
        for(int i=0;i<R;i++) printf(" %llu", Bd[bestD[i]]);
        printf("\nI");
        for(int i=0;i<R;i++) printf(" %llu", INS[bestI[i]]);
        printf("\n");
    } else {
        printf("NONE\n");
    }
    return 0;
}
