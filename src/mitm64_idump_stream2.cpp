// mitm64_idump_stream2.cpp -- streaming I-side dumper, v2 (review-driven).
//
// Same records and same stdin format as mitm64_idump_stream.cpp, with three
// additions that never change WHICH records are emitted (only how much work
// is spent on prefixes that cannot complete):
//   1. cheapest-completion bound: with t more insertions still required after
//      choosing INS[i], stop the loop when
//          partial * INS[i+1] * ... * INS[i+t]  >=  Icap        (INS ascending)
//      -- the same R_min idea as the exclusion engine; safe because that product
//      is <= any completion's product.
//   2. loop bound i <= |INS| - (remaining insertions).
//   3. optional stratification for the expanded-pool search: env IDUMP_SPLIT=<prime>
//      and IDUMP_HMAX=<h> keep only subsets with at most h insertion primes
//      > SPLIT (e.g. SPLIT=4000, HMAX=2, justified by the deletion/insertion
//      lower bound). Unset => no restriction (identical output to v1).
// Build: g++ -O3 -std=c++17 -pthread mitm64_idump_stream2.cpp -o mitm64_idump_stream2
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <thread>
#include <mutex>
typedef unsigned long long u64; typedef unsigned __int128 u128;
static u64 M, P0modM; static int R, T; static u128 Icap; static std::vector<u64> INS;
static std::vector<u128> PRE;                       // PRE[j] = prod INS[0..j-1]  (u128: <= 8 primes < 10^5)
static int SPLIT_IDX = 1 << 30; static int HMAX = 1 << 30;
static std::mutex out_mx;
static inline u64 mulmodM(u64 a, u64 b){ return (u64)(((u128)a * b) % M); }
static u128 parse_u128(const char* s){ u128 v=0; for(;*s;++s) v = v*10 + (u128)(*s-'0'); return v; }
struct Writer {
    std::vector<u64> buf; u64 n=0;
    Writer(){ buf.reserve(1<<20); }
    inline void put(u64 key,u64 id){ buf.push_back(key); buf.push_back(id); ++n; if(buf.size()>=(1u<<20)) flush(); }
    void flush(){ if(buf.empty()) return; std::lock_guard<std::mutex> g(out_mx); fwrite(buf.data(),8,buf.size(),stdout); buf.clear(); }
    ~Writer(){ flush(); }
};
// window product INS[a..a+t-1] exactly (PRE ratios need division; keep it exact by multiplying t terms)
static inline u128 window(int a, int t){ u128 v=1; for(int j=0;j<t;++j) v*=INS[a+j]; return v; }
static void rec(int start,int depth,u128 prod,u64 resI,u64 id,int nlarge,Writer& w){
    if(depth==R){ w.put(mulmodM(P0modM,resI), id); return; }
    const int n=(int)INS.size(); const int t=R-depth-1;        // insertions still needed after this one
    for(int i=start; i<=n-1-t; ++i){
        if(i>=SPLIT_IDX && nlarge+1>HMAX) break;                // all later INS are > SPLIT too
        u128 np=prod*INS[i];
        if(np>=Icap) break;
        if(t>0 && np*window(i+1,t)>=Icap) break;                // cheapest completion already too big
        rec(i+1,depth+1,np,mulmodM(resI,INS[i]%M), id | ((u64)i<<(8*depth)), nlarge + (i>=SPLIT_IDX), w);
    }
}
int main(){
    char buf[160]; u64 junk;
    if(scanf("%llu %d %llu",&M,&R,&P0modM)!=3) return 2;
    if(scanf("%159s",buf)!=1) return 2; Icap=parse_u128(buf); if(scanf("%159s",buf)!=1) return 2;
    int maxm; if(scanf("%d %d",&T,&maxm)!=2) return 2;
    int nB; if(scanf("%d",&nB)!=1) return 2; for(int i=0;i<nB;i++) if(scanf("%llu",&junk)!=1) return 2;
    int nI; if(scanf("%d",&nI)!=1) return 2; INS.resize(nI); for(auto&x:INS) if(scanf("%llu",&x)!=1) return 2;
    if(nI>256 || R>8){ fprintf(stderr,"nI<=256 and r<=8 required (8-bit packing)\n"); return 3; }
    for(int i=1;i<nI;i++) if(INS[i]<=INS[i-1]){ fprintf(stderr,"INS must be strictly ascending\n"); return 4; }
    if(const char* s=getenv("IDUMP_SPLIT")){ u64 sp=strtoull(s,nullptr,10); SPLIT_IDX=nI; for(int i=0;i<nI;i++) if(INS[i]>sp){ SPLIT_IDX=i; break; } }
    if(const char* h=getenv("IDUMP_HMAX")) HMAX=atoi(h);
    if(T<1) T=1;
    static char obuf[1<<22]; setvbuf(stdout, obuf, _IOFBF, sizeof obuf);
    std::vector<std::thread> th; std::vector<u64> counts(T,0);
    for(int t=0;t<T;t++) th.emplace_back([&,t](){
        Writer w; const int tt=R-1;
        for(int i0=t;i0<=(int)INS.size()-1-tt;i0+=T){
            if(i0>=SPLIT_IDX && 1>HMAX) break;
            u128 np=(u128)INS[i0]; if(np>=Icap) continue;                 // striped: cannot break globally
            if(tt>0 && np*window(i0+1,tt)>=Icap) continue;
            rec(i0+1,1,np,INS[i0]%M,(u64)i0,(i0>=SPLIT_IDX),w);
        }
        w.flush(); counts[t]=w.n;
    });
    for(auto&x:th) x.join();
    fflush(stdout);
    u64 tot=0; for(int t=0;t<T;t++){ tot+=counts[t]; fprintf(stderr,"t%d %llu\n",t,counts[t]); }
    fprintf(stderr,"TOTAL %llu records (r=%d, |INS|=%d, split_idx=%d, hmax=%d)\n",tot,R,nI,SPLIT_IDX,HMAX);
    return 0;
}
