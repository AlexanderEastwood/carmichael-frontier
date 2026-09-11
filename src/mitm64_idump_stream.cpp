// mitm64_idump_stream.cpp -- STREAMING I-side dumper for the GPU MITM join.
//
// Identical enumeration to src/mitm64_idump.cpp (ascending INS, exact product
// prune prodI < Icap, depth-0 index striped across T threads; the file-based
// dumper is validated record-exact against the engine's I_leaves), but the
// 16-byte records {u64 key=(P0*prodI) mod M, u64 id=8-bit-packed INS indices}
// go to STDOUT as one binary stream instead of per-thread files. Threads fill
// private 8 MB buffers and flush them under one mutex, so the stream is a
// well-formed sequence of records in an arbitrary (irrelevant) order. The GPU
// driver consumes it chunk-by-chunk, so the I side never has to fit on disk;
// only a time budget bounds an instance. stderr: per-thread + TOTAL counts.
// stdin format as the engine:  M r P0modM Icap (ignored) T (ignored) / nB b.. / nI i..
// Build: g++ -O3 -std=c++17 -pthread mitm64_idump_stream.cpp -o mitm64_idump_stream
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <thread>
#include <mutex>
typedef unsigned long long u64; typedef unsigned __int128 u128;
static u64 M, P0modM; static int R, T; static u128 Icap; static std::vector<u64> INS;
static std::mutex out_mx;
static inline u64 mulmodM(u64 a, u64 b){ return (u64)(((u128)a * b) % M); }
static u128 parse_u128(const char* s){ u128 v=0; for(;*s;++s) v = v*10 + (u128)(*s-'0'); return v; }

struct Writer {
    std::vector<u64> buf; u64 n=0;
    Writer(){ buf.reserve(1u<<20); }
    inline void put(u64 key,u64 id){ buf.push_back(key); buf.push_back(id); ++n; if(buf.size()>=(1u<<20)) flush(); }
    void flush(){ if(buf.empty()) return; std::lock_guard<std::mutex> g(out_mx); fwrite(buf.data(),8,buf.size(),stdout); buf.clear(); }
    ~Writer(){ flush(); }
};
static void rec(int start,int depth,u128 prod,u64 resI,u64 id,Writer& w){
    if(depth==R){ w.put(mulmodM(P0modM,resI), id); return; }
    const int n=(int)INS.size();
    for(int i=start;i<n;++i){
        u128 np=prod*INS[i];
        if(np>=Icap) break;
        rec(i+1,depth+1,np,mulmodM(resI,INS[i]%M), id | ((u64)i<<(8*depth)), w);
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
    if(T<1) T=1;
    static char obuf[1<<22]; setvbuf(stdout, obuf, _IOFBF, sizeof obuf);
    std::vector<std::thread> th; std::vector<u64> counts(T,0);
    for(int t=0;t<T;t++) th.emplace_back([&,t](){
        Writer w;
        for(int i0=t;i0<nI;i0+=T){ u128 np=(u128)INS[i0]; if(np>=Icap) continue; rec(i0+1,1,np,INS[i0]%M,(u64)i0,w); }
        w.flush(); counts[t]=w.n;
    });
    for(auto&x:th) x.join();
    fflush(stdout);
    u64 tot=0; for(int t=0;t<T;t++){ tot+=counts[t]; fprintf(stderr,"t%d %llu\n",t,counts[t]); }
    fprintf(stderr,"TOTAL %llu records (r=%d, |INS|=%d)\n",tot,R,nI);
    return 0;
}
