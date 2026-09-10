// carmichael.cpp -- M1 correct, instrumented baseline for A006931:
//   the least Carmichael number with EXACTLY k prime factors.
//
// Preproduct-tabulation search (structure adapted from J. Webster's
// small-carmichael-numbers) with ONE deliberate difference: a FULL wheel that
// admits every odd prime (3, 5, 7, ...) as a possible factor, so it reproduces
// the TRUE minima (k=3 -> 561, not 1729). See src/README.md.
//
// Correctness over speed. All big arithmetic is GMP (mpz_class). No 64-bit
// fast paths, no mod-30 wheel, no pair-step -- those are later milestones.
//
// Build: g++ -O2 -std=c++17 carmichael.cpp -o carmichael -lgmp -lgmpxx
// Run:   ./carmichael k [bound]         (bound omitted or 0 => bound-doubling)

#include <gmpxx.h>
#include <x86intrin.h>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

// ------------------------------------------------------------ prune flags ----
// The DEFAULT build (no flags) is now M2.5F: the R_min prune + base-case R*_min
// residue-class tightening + a Fermat-on-R rejection-first leaf filter. Every
// prior behavior stays reachable behind an opt-out flag, so the A/B harnesses
// keep working on identical source:
//   (default)          M2.5F : R_min + base-case R*_min + Fermat-on-R leaf filter
//   -DNO_FERMAT        M2.5  : as default but WITHOUT the Fermat leaf filter
//                              (this is the "pre-filter M2.5" A/B baseline)
//   -DNO_PRUNE_STAR    M2    : R_min lower-bound prune only (no R*_min, no filter)
//   -DNO_PRUNE         M1    : correct instrumented baseline (no prune, no filter)
//   -DPRUNE_STAR_REC   M2.5rF: default + R*_min at every recursive child (pays one
//                              extra GMP modular inverse per surviving child)
//
// Internally the code still keys off PRUNE / PRUNE_STAR / PRUNE_STAR_REC /
// FERMAT_FILTER; the block below derives those from the opt-out flags so that
// "no flags" == the fully-tightened default engine.
#ifdef PRUNE_STAR_REC
#  ifndef NO_PRUNE_STAR
#    define PRUNE_STAR_REC_ON
#  endif
#endif
#if defined(NO_PRUNE)
   // M1: define none of PRUNE / PRUNE_STAR / FERMAT_FILTER.
#elif defined(NO_PRUNE_STAR)
#  define PRUNE                       // M2
#else
#  define PRUNE                       // M2.5 / M2.5r
#  define PRUNE_STAR
#  ifdef PRUNE_STAR_REC_ON
#    ifndef PRUNE_STAR_REC
#      define PRUNE_STAR_REC
#    endif
#  endif
#  ifndef NO_FERMAT
#    define FERMAT_FILTER             // rejection-first Fermat-on-R leaf filter
#  endif
#endif

// ------------------------------------------------------------------ sieve ----
// Odd primes only. A Carmichael number is odd: if 2 | n then n is even, but n
// must be squarefree with >=3 distinct primes, so some ODD prime p | n, and
// Korselt needs (p-1) | (n-1) with p-1 even and n-1 odd -- impossible. So 2 is
// excluded by PROOF, not by a wheel shortcut. 3 and 5 ARE admitted.
static const uint32_t SIEVE_LIMIT = 5'000'000u;
static std::vector<uint32_t> PR;  // primes 3,5,7,... < SIEVE_LIMIT

static void build_sieve() {
    std::vector<bool> comp(SIEVE_LIMIT, false);
    for (uint32_t i = 2; (uint64_t)i * i < SIEVE_LIMIT; ++i)
        if (!comp[i])
            for (uint32_t j = i * i; j < SIEVE_LIMIT; j += i) comp[j] = true;
    for (uint32_t i = 3; i < SIEVE_LIMIT; i += 2)
        if (!comp[i]) PR.push_back(i);
}

// -------------------------------------------------------------- statistics ---
struct Stats {
    unsigned long long N_E = 0;  // extensions examined (recursive descents)
    unsigned long long N_I = 0;  // terminal inversions (base-case entries)
    unsigned long long N_F = 0;  // first rejection tests (AP candidates screened)
    unsigned long long N_H = 0;  // candidates surviving the cheap screen (reach the leaf)
    unsigned long long N_G = 0;  // candidates PASSING the Fermat filter (get full factor+Korselt)
    unsigned long long cyc_first = 0;  // cycles in the cheap first-rejection screen
    unsigned long long cyc_fermat = 0; // cycles in the Fermat-on-R rejection filter
    unsigned long long cyc_deep = 0;   // cycles in full factorization + Korselt
};
static Stats S;

// ---------------------------------------------------------------- globals ----
static const mpz_class g_two(2);         // Fermat-filter base
static mpz_class g_best;                 // current bound / incumbent (search n <= g_best)
static bool g_found = false;             // any exact-k Carmichael found <= bound
static mpz_class g_best_n;               // the incumbent number
static std::vector<mpz_class> g_best_fac;// its factorization
static int g_k;                          // target factor count

// integer t-th root (floor) of x
static mpz_class iroot(const mpz_class& x, unsigned long t) {
    mpz_class r;
    mpz_root(r.get_mpz_t(), x.get_mpz_t(), t);
    return r;
}

// admissible(q): no prime p already in `chosen` divides (q-1). This keeps
// gcd(P, lambda(P)) = 1 so P is invertible mod lambda(P) in the base case.
static inline bool admissible(uint32_t q, const std::vector<uint32_t>& chosen) {
    uint32_t qm1 = q - 1;
    for (uint32_t p : chosen)
        if (qm1 % p == 0) return false;
    return true;
}

#ifdef PRUNE
// r_min: R_min = product of the t SMALLEST primes each INDIVIDUALLY admissible
// to the current preproduct P (prime, > pmax, and no prime dividing P divides
// q-1). This is a valid LOWER BOUND on the cofactor R of ANY completion of P
// into a k-factor Carmichael number:
//   * every remaining factor is a prime > pmax (generation invariant), and
//   * every remaining factor q is admissible to P -- because in ANY Carmichael
//     number no prime factor p divides (q-1) of another prime factor q: if it
//     did, then p | (q-1) | lambda(n) | (n-1) while p | n, forcing p | gcd(n,n-1)
//     = 1, impossible.
// So the t true remaining factors are t DISTINCT admissible primes > pmax, hence
// their product >= product of the t smallest admissible primes > pmax = R_min.
// (The t smallest admissible primes may be mutually incompatible -- some may
// divide another's minus-one -- so R_min itself need not be realizable; that
// only makes the true product LARGER, so R_min remains a valid lower bound.)
// Returns false if the sieve is exhausted before t admissible primes are found
// (then we conservatively do NOT prune).
static bool r_min(uint32_t pmax, int t, const std::vector<uint32_t>& chosen,
                  mpz_class& out) {
    out = 1;
    int cnt = 0;
    for (uint32_t q : PR) {
        if (q <= pmax) continue;
        if (!admissible(q, chosen)) continue;
        out *= q;
        if (++cnt == t) return true;
    }
    return false;  // sieve exhausted -> cannot form the bound -> do not prune
}
#endif

// factor_R: is R a product of EXACTLY t DISTINCT primes, all > pmax?
// (Factors must be > pmax so every Carmichael number is generated exactly once,
//  by the unique prefix at which the base case fires.) Returns primes in `out`.
static bool factor_R(mpz_class R, uint32_t pmax, int t, std::vector<mpz_class>& out) {
    out.clear();
    // trial-divide by sieve primes > pmax
    for (uint32_t p : PR) {
        if (p <= pmax) continue;
        mpz_class mp(p);
        if (mp * mp > R) break;
        if (mpz_divisible_ui_p(R.get_mpz_t(), p)) {
            R /= p;
            if (mpz_divisible_ui_p(R.get_mpz_t(), p)) return false;  // square -> not squarefree
            out.push_back(mp);
            if ((int)out.size() > t) return false;
        }
    }
    if (R > 1) {
        // remaining cofactor must be a single prime > pmax
        if (R <= pmax) return false;
        if (mpz_probab_prime_p(R.get_mpz_t(), 30) == 0) return false;  // composite w/ large factors
        out.push_back(R);
    }
    return (int)out.size() == t;
}

// cheap first-rejection screen: does R have a prime factor <= pmax?
// (Such an R belongs to a different prefix; reject fast.) True => reject.
static inline bool has_small_factor(const mpz_class& R, uint32_t pmax) {
    for (uint32_t p : PR) {
        if (p > pmax) break;
        if (mpz_divisible_ui_p(R.get_mpz_t(), p)) return true;
    }
    return false;
}

// verify Korselt for the full prime list of n = P*R.
static bool korselt(const mpz_class& n, const std::vector<uint32_t>& chosen,
                    const std::vector<mpz_class>& rf) {
    mpz_class nm1 = n - 1;
    for (uint32_t p : chosen)
        if (!mpz_divisible_ui_p(nm1.get_mpz_t(), p - 1)) return false;
    for (const mpz_class& f : rf) {
        mpz_class fm1 = f - 1;
        if (!mpz_divisible_p(nm1.get_mpz_t(), fm1.get_mpz_t())) return false;
    }
    return true;
}

// base case: P*lambda(P) > bound (or t<=1). All valid completions have
// R = P^{-1} mod lambda(P) + j*lambda(P), R <= bound/P.
static void base_case(const mpz_class& P, const mpz_class& lam, uint32_t pmax,
                      int t, std::vector<uint32_t>& chosen) {
#ifdef PRUNE
    // R_min prune: if even the t smallest admissible primes > pmax already push
    // P*R_min past the bound, no completion here can be <= bound -> skip this
    // base-case entry entirely, BEFORE spending a modular inverse (N_I) and the
    // whole R = inv + j*lambda AP scan (N_F) / factorizations (N_H).
    mpz_class rmin;
    bool have_rmin = r_min(pmax, t, chosen, rmin);
    if (have_rmin && P * rmin > g_best) return;
#endif
    ++S.N_I;
    mpz_class inv;
    if (mpz_invert(inv.get_mpz_t(), P.get_mpz_t(), lam.get_mpz_t()) == 0) return;  // gcd!=1 (shouldn't happen)
    mpz_class r = inv;
    if (r == 0) r = lam;                 // want r in [1, lam]
#ifdef PRUNE_STAR
    // R*_min tightening (M2.5). Every completion has cofactor R with
    //   R >= R_min  AND  R ≡ a (mod lam),  a = P^{-1} mod lam = `inv`.
    // So the least feasible R is R*_min = R_min + ((a - R_min) mod lam): the least
    // value >= R_min in the forced residue class. `inv` is already in hand here
    // (base case computes it anyway) so this adds NO modular inverse over M2. Two
    // gains: (i) start the AP scan at R*_min, skipping every term below R_min that
    // M2 still screens (N_F) / factors (N_H); (ii) a strictly stronger whole-node
    // prune. Cannot drop the true minimum: any real completion has R >= R*_min.
    if (have_rmin) {
        mpz_class m;
        mpz_class d = inv - rmin;
        mpz_mod(m.get_mpz_t(), d.get_mpz_t(), lam.get_mpz_t());  // m in [0,lam)
        mpz_class rstar = rmin + m;                              // least >=rmin, ≡inv (mod lam)
        if (P * rstar > g_best) return;                         // stronger node prune
        r = rstar;                                              // advance AP start
    }
#endif
    mpz_class limit = g_best / P;        // R <= bound/P
    for (; r <= limit; r += lam) {
        ++S.N_F;
        unsigned long long c0 = __rdtsc();
        bool small = (r <= 1) || has_small_factor(r, pmax);
        S.cyc_first += __rdtsc() - c0;
        if (small) continue;

        ++S.N_H;
#ifdef FERMAT_FILTER
        // Rejection-first Fermat filter on R (the leaf's main cost-saver). If
        // n = P*R is Carmichael then 2^(n-1) ≡ 1 (mod n), hence 2^(n-1) ≡ 1
        // (mod R). So any candidate with 2^(P*R-1) mod R != 1 CANNOT be a
        // Carmichael number -> reject it before the expensive factor_R trial
        // division. This is a NECESSARY condition only: it may REJECT but must
        // NEVER accept, so every survivor still runs the full exact factor_R +
        // Korselt check below (a Fermat pass never short-circuits acceptance).
        // Correct for all t = k - w(P): when t > 1, R is EXPECTED composite, so
        // we deliberately do NOT use a compositeness/BPSW test here -- the
        // Fermat-on-R test is the right filter for every t. Exponent is the
        // full n-1 = P*R-1 (NOT R-1, which would be a primality test).
        {
            unsigned long long f0 = __rdtsc();
            mpz_class nm1 = P * r - 1;                 // n - 1 = P*R - 1
            mpz_class fres;
            mpz_powm(fres.get_mpz_t(), g_two.get_mpz_t(),
                     nm1.get_mpz_t(), r.get_mpz_t());  // 2^(n-1) mod R
            bool pass = (fres == 1);
            S.cyc_fermat += __rdtsc() - f0;
            if (!pass) continue;                       // definitive rejection
        }
        ++S.N_G;                                       // survivor -> full check
#endif
        unsigned long long d0 = __rdtsc();
        std::vector<mpz_class> rf;
        bool ok = factor_R(r, pmax, t, rf);
        mpz_class n;
        if (ok) {
            n = P * r;
            ok = korselt(n, chosen, rf);
        }
        S.cyc_deep += __rdtsc() - d0;

        if (ok && n < g_best) {
            g_best = n;
            g_best_n = n;
            g_found = true;
            g_best_fac.clear();
            for (uint32_t p : chosen) g_best_fac.push_back(mpz_class(p));
            for (const mpz_class& f : rf) g_best_fac.push_back(f);
        }
    }
}

// recursive case: extend P by one admissible prime q in (pmax, (bound/P)^(1/t)].
static void rec(const mpz_class& P, const mpz_class& lam, uint32_t pmax, int t,
                std::vector<uint32_t>& chosen) {
    if (t <= 1 || P * lam > g_best) {
        base_case(P, lam, pmax, t, chosen);
        return;
    }
    mpz_class hi = iroot(g_best / P, (unsigned long)t);
    if (hi > PR.back()) {
        std::cerr << "FATAL: needed prime interval exceeds sieve limit ("
                  << hi << " > " << PR.back() << "); raise SIEVE_LIMIT.\n";
        std::exit(2);
    }
    for (uint32_t q : PR) {
        if (q <= pmax) continue;
        if (q > hi) break;
        if (!admissible(q, chosen)) continue;
        mpz_class Pq = P * q;
#ifdef PRUNE
        // R_min prune on the CHILD (preproduct Pq, largest prime q, t-1 factors
        // left): skip the whole descent when Pq * R_min(child) > bound. Done
        // BEFORE counting the extension, so N_E measures descents actually taken.
        {
            chosen.push_back(q);            // admissibility of child is w.r.t. Pq
            mpz_class rmin;
            bool have_rmin = r_min(q, t - 1, chosen, rmin);
            bool skip = have_rmin && Pq * rmin > g_best;
#ifdef PRUNE_STAR_REC
            // Expensive variant (M2.5r): tighten the child prune to R*_min. The
            // residue constraint R ≡ (Pq)^{-1} (mod lam(Pq)) holds at EVERY node,
            // not just base cases, so R*_min is a valid stronger lower bound here
            // too. But unlike the base case, a recursive node has no inverse in
            // hand -- this pays one extra GMP modular inverse (and the lcm early)
            // per surviving child that plain R_min never needs. Measured, not
            // assumed to be a win.
            if (have_rmin && !skip) {
                mpz_class lamq;
                mpz_lcm_ui(lamq.get_mpz_t(), lam.get_mpz_t(), q - 1);
                mpz_class a;
                if (mpz_invert(a.get_mpz_t(), Pq.get_mpz_t(), lamq.get_mpz_t())) {
                    mpz_class m, d = a - rmin;
                    mpz_mod(m.get_mpz_t(), d.get_mpz_t(), lamq.get_mpz_t());
                    mpz_class rstar = rmin + m;
                    if (Pq * rstar > g_best) skip = true;
                }
            }
#endif
            chosen.pop_back();
            if (skip) continue;
        }
#endif
        ++S.N_E;
        mpz_class lamq;
        mpz_lcm_ui(lamq.get_mpz_t(), lam.get_mpz_t(), q - 1);
        chosen.push_back(q);
        rec(Pq, lamq, q, t - 1, chosen);
        chosen.pop_back();
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " k [bound]\n";
        return 1;
    }
    g_k = std::stoi(argv[1]);
    if (g_k < 3) { std::cerr << "k must be >= 3\n"; return 1; }
    mpz_class user_bound(0);
    if (argc >= 3) user_bound = mpz_class(argv[2]);

    build_sieve();

    // calibrate rdtsc ticks/sec
    auto cal0 = std::chrono::steady_clock::now();
    unsigned long long tsc0 = __rdtsc();
    while (std::chrono::duration<double>(std::chrono::steady_clock::now() - cal0).count() < 0.05) {}
    unsigned long long tsc1 = __rdtsc();
    auto cal1 = std::chrono::steady_clock::now();
    double tps = (tsc1 - tsc0) /
                 std::chrono::duration<double>(cal1 - cal0).count();

    auto wall0 = std::chrono::steady_clock::now();

    // bound-doubling if no bound given: independent of any known table.
    std::vector<mpz_class> bounds;
    if (user_bound > 0) {
        bounds.push_back(user_bound);
    } else {
        mpz_class b = 1000;
        for (int i = 0; i < 130 && !g_found; ++i) {
            bounds.push_back(b);
            b *= 4;
        }
    }

    for (const mpz_class& B : bounds) {
        if (g_found) break;
        g_best = B;
        std::vector<uint32_t> chosen;
        mpz_class P(1), lam(1);
        rec(P, lam, 2, g_k, chosen);  // pmax=2 so first prime is >=3
    }

    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();
    double t_first = S.cyc_first / tps;
    double t_fermat = S.cyc_fermat / tps;
    double t_deep = S.cyc_deep / tps;
    double t_enum = wall - t_first - t_fermat - t_deep;
    if (t_enum < 0) t_enum = 0;

    std::cout << "k=" << g_k << "\n";
    if (g_found) {
        std::cout << "n=" << g_best_n << "\n";
        std::cout << "FACTORS:";
        for (const mpz_class& f : g_best_fac) std::cout << " " << f;
        std::cout << "\n";
        std::cout << "nfactors=" << g_best_fac.size() << "\n";
    } else {
        std::cout << "n=NONE\nFACTORS:\nnfactors=0\n";
    }
    std::cout << "N_E=" << S.N_E << " N_I=" << S.N_I
              << " N_F=" << S.N_F << " N_H=" << S.N_H << " N_G=" << S.N_G << "\n";
    // survivor fraction s = fraction of leaf candidates (N_H, past the cheap
    // screen) that pass the Fermat filter and reach the full factor+Korselt.
    std::cout << "survivor_frac=" << (S.N_H > 0 ? (double)S.N_G / S.N_H : 0) << "\n";
    std::cout << "wall_s=" << wall
              << " t_enum_s=" << t_enum
              << " t_first_s=" << t_first
              << " t_fermat_s=" << t_fermat
              << " t_deep_s=" << t_deep << "\n";
    std::cout << "frac_enum=" << (wall > 0 ? t_enum / wall : 0)
              << " frac_first=" << (wall > 0 ? t_first / wall : 0)
              << " frac_fermat=" << (wall > 0 ? t_fermat / wall : 0)
              << " frac_deep=" << (wall > 0 ? t_deep / wall : 0) << "\n";
    return 0;
}
