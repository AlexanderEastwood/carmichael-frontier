/*
Direct verification of the exclusion S65 >= 10^148.
Prime cap 2266 is supplied by the accompanying Python certificate:
  every 64-prime admissible cofactor is at least its certified global minimum A64;
  floor((10^148-1)/A64)=2266.
This file does not independently prove the cofactor bound.

All admissible complete products are visited and Korselt is tested with exact
arithmetic, without modular inverses or progression pruning. The binary
traversal enumerates the final two prime choices inside one call. Therefore
'nodes' is the count for the equivalent uncollapsed binary traversal, while
'calls' counts the actual loop entries. Overflow aborts, never certifies.

512-bit safety: running product bounds are <10^148 before multiplication by
one prime <=2266. The t=2 quotient comparison multiplies by at most 2266^2;
the 63-prime prefix was generated with room for two primes >=313 and 317.
All these quantities fit 512 bits; multiplication is also checked.
Prefix remainders times the final two-prime product fit uint64_t.

Compile: g++ -O3 -std=c++17 verify65.cpp -o verify65
Run: ./verify65 148 45
Expected: COMPLETE 1, leaves 22287138, carmichael_hits 0.
The optional last argument restricts to one smallest prime and is NOT a
standalone global exclusion. The default has no such restriction.
*/
#include <array>
#include <vector>
#include <iostream>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <algorithm>
#include <stdexcept>
#include <limits>
using U=uint64_t;
struct Big {
    std::array<U,8> x{};
    int len=1;
    Big(U a=0){x[0]=a;}
    Big times(unsigned p) const {
        Big z;z.len=len;U carry=0;
        for(int i=0;i<len;i++){
            __uint128_t v=(__uint128_t)x[i]*p+carry;
            z.x[i]=(U)v;carry=(U)(v>>64);
        }
        if(carry){if(z.len==8)throw std::overflow_error("512-bit product overflow");z.x[z.len++]=carry;}
        return z;
    }
    bool operator<(const Big& b)const{
        if(len!=b.len)return len<b.len;
        for(int i=len-1;i>=0;i--)if(x[i]!=b.x[i])return x[i]<b.x[i];
        return false;
    }
};
constexpr int MAXW=6;
using Mask=std::array<U,MAXW>;
std::vector<unsigned> primes;
std::vector<Mask> future;
std::vector<std::vector<int>> conflicts;
Big bound;
std::vector<unsigned> chosen;
U carmichael_hits=0;
U modsmall(const Big& x, unsigned m){U r=0;for(int i=x.len-1;i>=0;--i)r=(U)((((__uint128_t)r<<64)+x.x[i])%m);return r;}
int words;
U nodes=0,leaves=0,calls=0;
void add(U& x,U n){if(n>std::numeric_limits<U>::max()-x)throw std::overflow_error("counter overflow");x+=n;}
auto begun=std::chrono::steady_clock::now();
double deadline=300;
double elapsed(){return std::chrono::duration<double>(std::chrono::steady_clock::now()-begun).count();}

void walk(Mask mask,int t,const Big& P){
    while(true){
        add(nodes,1);add(calls,1);
        if((calls&1048575)==0 && elapsed()>deadline)throw std::runtime_error("TIMEOUT");
        if(t==2){
            U lo=0,hi=U(primes.back())*primes.back()+1;
            while(lo+1<hi){U mid=lo+(hi-lo)/2;if(P.times((unsigned)mid)<bound)lo=mid;else hi=mid;}
            U H=lo;
            std::vector<int> a;
            unsigned minimum=0;
            for(int w=0;w<words;w++){
                U bits=mask[w];
                while(bits){int i=64*w+__builtin_ctzll(bits);bits&=bits-1;
                    if(!minimum)minimum=primes[i];
                    if(U(minimum)*primes[i]>H)break;
                    a.push_back(i);
                }
                if(minimum && w*64<(int)primes.size() && U(minimum)*primes[std::min((int)primes.size()-1,w*64)]>H)break;
            }
            U pairs=0,r=0;int j=(int)a.size()-1;
            std::vector<std::pair<unsigned,U>> residues;
            for(unsigned p:chosen)residues.push_back({p-1,modsmall(P,p-1)});
            for(int i=0;i+1<(int)a.size() && U(primes[a[i]])*primes[a[i+1]]<=H;i++){
                ++r;
                while(j>i && U(primes[a[i]])*primes[a[j]]>H)--j;
                for(int bidx=i+1;bidx<=j;++bidx){
                    unsigned aa=primes[a[i]],bb=primes[a[bidx]];
                    if((bb-1)%aa==0)continue;
                    ++pairs; U R=U(aa)*bb;bool pass=true;
                    for(auto [m,rem]:residues)if((rem*R)%m!=1){pass=false;break;}
                    if(pass){Big n=P.times(aa).times(bb);if(modsmall(n,aa-1)==1 && modsmall(n,bb-1)==1)++carmichael_hits;}
                }
            }
            add(leaves,pairs);add(nodes,2*r+2*pairs);
            return;
        }
        if(t==1){
            // The greatest index j with P*prime[j] < bound, by exact comparison.
            int lo=0,hi=(int)primes.size();
            while(lo<hi){int mid=lo+(hi-lo)/2;if(P.times(primes[mid])<bound)lo=mid+1;else hi=mid;}
            U m=0;
            for(int w=0;w<lo/64;w++)m+=__builtin_popcountll(mask[w]);
            if(lo%64)m+=__builtin_popcountll(mask[lo/64]&((U(1)<<(lo%64))-1));
            add(leaves,m);add(nodes,2*m);
            return;
        }
        Big lower=P;int needed=t,first=-1;
        for(int w=0;w<words && needed;w++){
            U bits=mask[w];
            while(bits && needed){
                int bit=__builtin_ctzll(bits),i=64*w+bit;
                if(first<0)first=i;
                lower=lower.times(primes[i]);
                if(!(lower<bound))return;
                --needed;bits&=bits-1;
            }
        }
        if(needed)return;
        Mask child{};
        for(int w=0;w<words;w++)child[w]=mask[w]&future[first][w];
        U oldnodes=nodes,oldleaves=leaves;double before=elapsed();
        chosen.push_back(primes[first]);
        walk(child,t-1,P.times(primes[first]));
        chosen.pop_back();
        if(t==65)std::cout<<"root "<<primes[first]<<" nodes "<<nodes-oldnodes<<" leaves "<<leaves-oldleaves<<" seconds "<<elapsed()-before<<std::endl;
        mask[first/64]&=~(U(1)<<(first%64));
    }
}
int main(int argc,char**argv){
    int e=argc>1?std::atoi(argv[1]):148;
    if(e!=148)throw std::runtime_error("exponent must be 148");
    if(argc>2)deadline=std::atof(argv[2]);
    int qcap=2266;
    std::vector<bool> prime(qcap+1,true);prime[0]=prime[1]=false;
    for(int p=2;p*p<=qcap;p++)if(prime[p])for(int q=p*p;q<=qcap;q+=p)prime[q]=false;
    for(int p=3;p<=qcap;p++)if(prime[p])primes.push_back(p);
    words=(primes.size()+63)/64;if(words>MAXW)throw std::runtime_error("mask width");
    bound=Big(1);for(int i=0;i<e;i++)bound=bound.times(10);
    future.resize(primes.size());conflicts.resize(primes.size());Mask all{};
    for(int i=0;i<(int)primes.size();i++){
        all[i/64]|=U(1)<<(i%64);
        for(int j=i+1;j<(int)primes.size();j++)if((primes[j]-1)%primes[i])future[i][j/64]|=U(1)<<(j%64);else conflicts[i].push_back(j);
    }
    begun=std::chrono::steady_clock::now();bool complete=false;
    try{
        if(argc>3){
            unsigned fixed=std::atoi(argv[3]);
            auto it=std::lower_bound(primes.begin(),primes.end(),fixed);
            if(it==primes.end() || *it!=fixed)throw std::runtime_error("fixed smallest factor must be an odd prime in the universe");
            int i=it-primes.begin();chosen.push_back(fixed);walk(future[i],64,Big(fixed));chosen.pop_back();
        }else walk(all,65,Big(1));
        complete=true;
    }catch(const std::runtime_error&e){if(std::string(e.what())!="TIMEOUT")throw;}
    std::cout<<"COMPLETE "<<complete<<" exponent "<<e<<" fixed_smallest_prime "<<(argc>3?argv[3]:"ALL")<<" universe "<<primes.size()<<" nodes "<<nodes<<" leaves "<<leaves<<" carmichael_hits "<<carmichael_hits<<" calls "<<calls<<" seconds "<<elapsed()<<std::endl;
}
