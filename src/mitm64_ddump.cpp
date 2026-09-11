// mitm64_ddump.cpp -- product-filtered DELETION-side dumper for the exchange MITM.
//
// Astra (2026-09-11): apply the insertion product bound symmetrically to deletions. With base
// product P_B, insertion product P_I >= I_min,r (product of the r smallest insertion primes) and
// bound U, an improving completion n = P_B * P_I / P_D < U forces
//     P_D  >  P_B * P_I / U  >=  P_B * I_min,r / U,   i.e.   P_D >= Dmin := floor(P_B*I_min,r/U) + 1.
// At r=9 for N' this keeps ~4.1e7 of the C(64,9) = 2.75e10 deletion subsets (~670x fewer), so the
// D side fits resident on the GPU and the I side is streamed through it. The filter only drops
// deletion subsets that cannot be part of an improving completion, so joins with an Icap-bounded
// I side lose no improvement (they do lose some non-improving residue matches).
//
// Records: (key = prodD mod M, id = r packed 6-bit base indices, ascending) as u64 pairs on stdout,
// identical layout/id encoding to gpu.mitm_gpu2.generate(B0, r, M, bits=6). stderr: TOTAL line.
// stdin:  M r Dmin\n 64 <base primes ascending>\n        (Dmin as decimal, < 2^128)
// Build: g++ -O3 -std=c++17 mitm64_ddump.cpp -o mitm64_ddump
#include <cstdio>
#include <cstdlib>
#include <vector>
typedef unsigned long long u64; typedef unsigned __int128 u128;
static u64 M; static int R; static u128 Dmin; static std::vector<u64> B;   // B ascending, |B| = 64
static std::vector<u128> SUF;          // SUF[t] = product of the t largest base primes (t <= R)
static std::vector<u64> buf; static u64 total = 0;
static inline u64 mulmodM(u64 a, u64 b){ return (u64)(((u128)a * b) % M); }
static u128 parse_u128(const char* s){ u128 v=0; for(;*s;++s) v = v*10 + (u128)(*s-'0'); return v; }
static void flush(){ if(!buf.empty()){ fwrite(buf.data(),8,buf.size(),stdout); buf.clear(); } }
// choose indices in DESCENDING order (largest primes first) so the prune "partial * (t largest still
// available) < Dmin" cuts early; ids are still emitted with ascending index positions.
static void rec(int hi, int depth, u128 prod, u64 res, int idx[]){
    if(depth==R){
        u64 id=0; for(int j=R-1,pos=0;j>=0;--j,++pos) id |= ((u64)idx[j]) << (6*pos);   // ascending order
        buf.push_back(mulmodM(1,res)); buf.push_back(id); ++total; if(buf.size()>=(1u<<20)) flush(); return;
    }
    const int t = R-depth-1;                       // deletions still needed after this one
    for(int i=hi; i>=t; --i){                      // need t more indices below i
        u128 np = prod * B[i];
        // largest achievable completion from here: np * product of the t largest primes below B[i]
        u128 best = np; for(int j=0;j<t;j++) best *= B[i-1-j];
        if(best < Dmin) break;                     // smaller i only lowers it further
        idx[depth]=i; rec(i-1, depth+1, np, mulmodM(res, B[i]%M), idx);
    }
}
int main(){
    char s[160];
    if(scanf("%llu %d %159s",&M,&R,s)!=3) return 2; Dmin=parse_u128(s);
    int nB; if(scanf("%d",&nB)!=1) return 2; B.resize(nB); for(auto&x:B) if(scanf("%llu",&x)!=1) return 2;
    if(nB!=64 || R<1 || R>10){ fprintf(stderr,"need 64 base primes and 1<=r<=10\n"); return 3; }
    for(int i=1;i<nB;i++) if(B[i]<=B[i-1]){ fprintf(stderr,"base must be strictly ascending\n"); return 4; }
    static char obuf[1<<22]; setvbuf(stdout, obuf, _IOFBF, sizeof obuf); buf.reserve(1<<20);
    int idx[16]; rec(nB-1, 0, (u128)1, 1, idx); flush(); fflush(stdout);
    fprintf(stderr,"TOTAL %llu records (r=%d, filtered by prodD >= Dmin)\n", total, R);
    return 0;
}
