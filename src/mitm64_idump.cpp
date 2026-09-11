// mitm64_idump.cpp -- I-side record DUMPER for the GPU MITM join.
//
// Same insertion enumeration as src/mitm64.cpp's recI (ascending INS, exact
// product prune prodI < Icap, depth-0 index distributed across T threads), but
// instead of probing an in-memory D-residue table it WRITES every surviving
// I-subset as a fixed 16-byte record  { u64 key = (P0 * prodI) mod M ,
// u64 id = indices into INS packed 8 bits each }  to per-thread binary files
// <prefix>.t<tid>.bin.  The GPU then sorts these and joins them against the
// (chunk-sorted, resident) D side that it generates itself -- which is what
// lifts the r<=7 CPU cap: the CPU never has to hold C(64,r) D-records.
//
// stdin (whitespace separated), matching the Stage-B engine format:
//   M  r  P0modM  Icap  (ignored)  T  (ignored)
//   nB  b1..bnB          (base, ignored here; present for format compatibility)
//   nI  i1..inI          (insertion pool, must be ascending, nI <= 256)
// argv[1] = output prefix.   stderr: per-thread and total record counts.
// Build: g++ -O3 -std=c++17 -pthread mitm64_idump.cpp -o mitm64_idump
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <thread>
typedef unsigned long long u64; typedef unsigned __int128 u128;
static u64 M, P0modM; static int R, T; static u128 Icap; static std::vector<u64> INS;
static inline u64 mulmodM(u64 a, u64 b){ return (u64)(((u128)a * b) % M); }
static u128 parse_u128(const char* s){ u128 v=0; for(;*s;++s) v = v*10 + (u128)(*s-'0'); return v; }

struct Writer {
    FILE* f; std::vector<u64> buf; u64 n=0;
    explicit Writer(const std::string& path){ f=fopen(path.c_str(),"wb"); buf.reserve(1<<20); }
    inline void put(u64 key,u64 id){ buf.push_back(key); buf.push_back(id); if(buf.size()>=(1u<<20)) flush(); ++n; }
    void flush(){ if(!buf.empty()){ fwrite(buf.data(),8,buf.size(),f); buf.clear(); } }
    ~Writer(){ flush(); fclose(f); }
};

static void rec(int start,int depth,u128 prod,u64 resI,u64 id,Writer& w){
    if(depth==R){ w.put(mulmodM(P0modM,resI), id); return; }
    const int n=(int)INS.size();
    for(int i=start;i<n;++i){
        u128 np=prod*INS[i];
        if(np>=Icap) break;                                   // INS ascending -> later only larger
        rec(i+1,depth+1,np,mulmodM(resI,INS[i]%M), id | ((u64)i<<(8*depth)), w);
    }
}
int main(int argc,char**argv){
    if(argc<2){ fprintf(stderr,"usage: %s out_prefix < stdin\n",argv[0]); return 1; }
    char buf[160]; u64 junk;
    if(scanf("%llu %d %llu",&M,&R,&P0modM)!=3) return 2;
    scanf("%159s",buf); Icap=parse_u128(buf); scanf("%159s",buf);           // Icap, prodDmin(ignored)
    int maxm; scanf("%d %d",&T,&maxm);
    int nB; scanf("%d",&nB); for(int i=0;i<nB;i++) scanf("%llu",&junk);
    int nI; scanf("%d",&nI); INS.resize(nI); for(auto&x:INS) scanf("%llu",&x);
    if(nI>256 || R>8){ fprintf(stderr,"nI<=256 and r<=8 required (8-bit packing)\n"); return 3; }
    if(T<1) T=1;
    std::vector<std::thread> th; std::vector<u64> counts(T,0);
    for(int t=0;t<T;t++) th.emplace_back([&,t](){
        Writer w(std::string(argv[1])+".t"+std::to_string(t)+".bin");
        for(int i0=t;i0<nI;i0+=T){
            u128 np=(u128)INS[i0]; if(np>=Icap) continue;               // stride: cannot break globally
            rec(i0+1,1,np,INS[i0]%M,(u64)i0,w);
        }
        counts[t]=w.n;
    });
    for(auto&x:th) x.join();
    u64 tot=0; for(int t=0;t<T;t++){ tot+=counts[t]; fprintf(stderr,"t%d %llu\n",t,counts[t]); }
    fprintf(stderr,"TOTAL %llu records (r=%d, |INS|=%d)\n",tot,R,nI);
    return 0;
}
