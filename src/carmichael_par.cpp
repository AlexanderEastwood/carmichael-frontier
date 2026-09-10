// carmichael_par.cpp -- M3: single-box PARALLEL search for A006931 (least
// Carmichael number with EXACTLY k prime factors), with a SHARED GLOBAL
// INCUMBENT bound. Built beside the M1/M2 single-proc engine (carmichael.cpp);
// the search math (sieve, admissibility, R_min prune, factor_R, Korselt) is
// identical -- this file only adds (1) coarse prefix-subtree job partitioning
// and (2) a shared, cheaply-read incumbent so every worker prunes against the
// current global best without the deep-copy contention of Jon's dws.
//
// Design (contrast with dws, which scales NEGATIVELY):
//   * PARTITION: enumerate all preproduct prefixes down to a small split depth
//     -> a frontier cut of the search tree = many independent subtree jobs.
//     Running rec() on every frontier node covers the whole single-proc tree
//     exactly once (completeness), so the parallel minimum == the M2 minimum.
//   * DYNAMIC QUEUE: one atomic job index; workers pull coarse jobs. No
//     fine-grained work-stealing, no per-steal GMP deep copies.
//   * SHARED INCUMBENT: essentially ONE integer (+ its winning factorization)
//     behind a mutex. The hot path NEVER locks: each worker keeps a thread-local
//     copy of the bound and a version stamp; per node it does ONE relaxed atomic
//     load of a global version counter, and only re-copies the mpz (under the
//     mutex) on the rare event that some worker lowered the bound. Writes happen
//     only when a strictly-better n is found. A stale local bound is always
//     >= the true shared bound, so it only ever prunes LESS -> never removes the
//     true minimum. Correctness is preserved for any worker count / timing.
//
// Build: g++ -O2 -std=c++17 carmichael_par.cpp -o carmichael_par -lgmp -lgmpxx -pthread
// Run:   ./carmichael_par k [bound] [workers]   (bound 0/omitted => bound-doubling)

#include <gmpxx.h>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

// ------------------------------------------------------------------ sieve ----
static const uint32_t SIEVE_LIMIT = 5'000'000u;
static std::vector<uint32_t> PR;  // odd primes 3,5,7,... < SIEVE_LIMIT

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
    unsigned long long N_E = 0, N_I = 0, N_F = 0, N_H = 0, N_G = 0;
    void add(const Stats& o) { N_E += o.N_E; N_I += o.N_I; N_F += o.N_F; N_H += o.N_H; N_G += o.N_G; }
};

static const mpz_class g_two(2);         // Fermat-filter base

// ------------------------------------------------ shared incumbent state -----
static std::mutex g_mtx;                     // guards the four below
static mpz_class  g_best;                    // shared bound / incumbent
static mpz_class  g_best_n;                  // incumbent number
static std::vector<mpz_class> g_best_fac;    // its factorization
static bool       g_found = false;
static std::atomic<uint64_t> g_ver{0};       // bumped on every improvement
static int        g_k;

// integer t-th root (floor)
static mpz_class iroot(const mpz_class& x, unsigned long t) {
    mpz_class r; mpz_root(r.get_mpz_t(), x.get_mpz_t(), t); return r;
}

static inline bool admissible(uint32_t q, const std::vector<uint32_t>& chosen) {
    uint32_t qm1 = q - 1;
    for (uint32_t p : chosen) if (qm1 % p == 0) return false;
    return true;
}

// R_min lower bound (identical to carmichael.cpp M2). Product of the t smallest
// primes individually admissible to P (prime, > pmax, none dividing P divides q-1).
static bool r_min(uint32_t pmax, int t, const std::vector<uint32_t>& chosen, mpz_class& out) {
    out = 1; int cnt = 0;
    for (uint32_t q : PR) {
        if (q <= pmax) continue;
        if (!admissible(q, chosen)) continue;
        out *= q;
        if (++cnt == t) return true;
    }
    return false;  // sieve exhausted -> do not prune
}

static bool factor_R(mpz_class R, uint32_t pmax, int t, std::vector<mpz_class>& out) {
    out.clear();
    for (uint32_t p : PR) {
        if (p <= pmax) continue;
        mpz_class mp(p);
        if (mp * mp > R) break;
        if (mpz_divisible_ui_p(R.get_mpz_t(), p)) {
            R /= p;
            if (mpz_divisible_ui_p(R.get_mpz_t(), p)) return false;  // not squarefree
            out.push_back(mp);
            if ((int)out.size() > t) return false;
        }
    }
    if (R > 1) {
        if (R <= pmax) return false;
        if (mpz_probab_prime_p(R.get_mpz_t(), 30) == 0) return false;
        out.push_back(R);
    }
    return (int)out.size() == t;
}

static inline bool has_small_factor(const mpz_class& R, uint32_t pmax) {
    for (uint32_t p : PR) {
        if (p > pmax) break;
        if (mpz_divisible_ui_p(R.get_mpz_t(), p)) return true;
    }
    return false;
}

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

// --------------------------------------------------------------- worker -------
struct Worker {
    mpz_class local_best;   // thread-local copy of the shared bound
    uint64_t  local_ver;    // version stamp of local_best
    Stats S;

    void init() {  // sync from shared under lock (call once per job start is enough)
        std::lock_guard<std::mutex> lk(g_mtx);
        local_best = g_best;
        local_ver  = g_ver.load(std::memory_order_relaxed);
    }
    // Hot path: one relaxed/acquire atomic load; lock+copy only on a real change.
    inline void refresh() {
        uint64_t v = g_ver.load(std::memory_order_acquire);
        if (v != local_ver) {
            std::lock_guard<std::mutex> lk(g_mtx);
            local_best = g_best;
            local_ver  = g_ver.load(std::memory_order_relaxed);
        }
    }
    // Rare: a candidate beat our local bound. Commit under the mutex iff it also
    // beats the true shared incumbent; then refresh local from shared.
    void submit(const mpz_class& n, const std::vector<uint32_t>& chosen,
                const std::vector<mpz_class>& rf) {
        std::lock_guard<std::mutex> lk(g_mtx);
        if (!g_found || n < g_best) {
            g_best = n; g_best_n = n; g_found = true;
            g_best_fac.clear();
            for (uint32_t p : chosen) g_best_fac.push_back(mpz_class(p));
            for (const mpz_class& f : rf) g_best_fac.push_back(f);
            g_ver.fetch_add(1, std::memory_order_release);
        }
        local_best = g_best;
        local_ver  = g_ver.load(std::memory_order_relaxed);
    }
};

// base case: all completions have R = P^{-1} mod lam + j*lam, R <= bound/P.
static void base_case(const mpz_class& P, const mpz_class& lam, uint32_t pmax,
                      int t, std::vector<uint32_t>& chosen, Worker& w) {
    w.refresh();
    // R_min prune before spending the inverse / AP scan.
    mpz_class rmin;
    bool have_rmin = r_min(pmax, t, chosen, rmin);
    if (have_rmin && P * rmin > w.local_best) return;
    ++w.S.N_I;
    mpz_class inv;
    if (mpz_invert(inv.get_mpz_t(), P.get_mpz_t(), lam.get_mpz_t()) == 0) return;
    mpz_class r = inv;
    if (r == 0) r = lam;
    // M2.5 base-case R*_min tightening: every completion has cofactor R with
    // R >= R_min AND R ≡ inv (mod lam). The least feasible R is therefore
    // R*_min = R_min + ((inv - R_min) mod lam) -- the least value >= R_min in the
    // forced residue class. `inv` is already in hand, so this adds NO extra
    // inverse over M2. It advances the AP start past every term below R_min and
    // gives a strictly stronger whole-node prune. Cannot drop the true minimum:
    // any real completion already has R >= R*_min.
    if (have_rmin) {
        mpz_class m, d = inv - rmin;
        mpz_mod(m.get_mpz_t(), d.get_mpz_t(), lam.get_mpz_t());  // m in [0,lam)
        mpz_class rstar = rmin + m;
        if (P * rstar > w.local_best) return;                    // stronger node prune
        r = rstar;                                               // advance AP start
    }
    mpz_class limit = w.local_best / P;
    for (; r <= limit; r += lam) {
        ++w.S.N_F;
        if (r <= 1 || has_small_factor(r, pmax)) continue;
        ++w.S.N_H;
        // Rejection-first Fermat filter on R (necessary condition; may only
        // REJECT). If n = P*R is Carmichael then 2^(n-1) ≡ 1 (mod R); reject any
        // candidate failing that before the expensive factor_R. Survivors STILL
        // run the full exact factor_R + Korselt below -- never short-circuited.
        // Correct for all t: when t > 1, R is expected composite, so this is
        // deliberately a Fermat test, not a compositeness/BPSW test. Exponent is
        // the full n-1 = P*R-1 (NOT R-1).
        {
            mpz_class nm1 = P * r - 1;
            mpz_class fres;
            mpz_powm(fres.get_mpz_t(), g_two.get_mpz_t(),
                     nm1.get_mpz_t(), r.get_mpz_t());
            if (fres != 1) continue;   // definitive rejection
        }
        ++w.S.N_G;
        std::vector<mpz_class> rf;
        bool ok = factor_R(r, pmax, t, rf);
        mpz_class n;
        if (ok) { n = P * r; ok = korselt(n, chosen, rf); }
        if (ok && n < w.local_best) {
            w.submit(n, chosen, rf);
            limit = w.local_best / P;  // tighten our own scan immediately
        }
    }
}

// recursive case: extend P by one admissible prime q in (pmax, (bound/P)^(1/t)].
static void rec(const mpz_class& P, const mpz_class& lam, uint32_t pmax, int t,
                std::vector<uint32_t>& chosen, Worker& w) {
    w.refresh();
    if (t <= 1 || P * lam > w.local_best) { base_case(P, lam, pmax, t, chosen, w); return; }
    mpz_class hi = iroot(w.local_best / P, (unsigned long)t);
    if (hi > PR.back()) {
        std::cerr << "FATAL: needed prime interval exceeds sieve limit (" << hi
                  << " > " << PR.back() << "); raise SIEVE_LIMIT.\n";
        std::exit(2);
    }
    for (uint32_t q : PR) {
        if (q <= pmax) continue;
        if (q > hi) break;
        if (!admissible(q, chosen)) continue;
        mpz_class Pq = P * q;
        {   // R_min prune on the child
            chosen.push_back(q);
            mpz_class rmin;
            bool skip = r_min(q, t - 1, chosen, rmin) && Pq * rmin > w.local_best;
            chosen.pop_back();
            if (skip) continue;
        }
        ++w.S.N_E;
        mpz_class lamq;
        mpz_lcm_ui(lamq.get_mpz_t(), lam.get_mpz_t(), q - 1);
        chosen.push_back(q);
        rec(Pq, lamq, q, t - 1, chosen, w);
        chosen.pop_back();
    }
}

// ----------------------------------------------------- job partitioning ------
struct Job {
    std::vector<uint32_t> chosen;
    mpz_class P, lam;
    uint32_t pmax;
    int t;
};

// Enumerate a FRONTIER CUT of the tree using the loosest bound B (no tightening,
// so it prunes the least -> only ever prunes subtrees whose every completion
// already exceeds B, i.e. cannot be the answer). A node becomes a job when it
// reaches split depth OR is already a base case. Uses the SAME child-generation
// and prune as rec(), so the jobs' rec() runs reconstruct the whole search.
static void gen_jobs(const mpz_class& B, const mpz_class& P, const mpz_class& lam,
                     uint32_t pmax, int t, int depth, int dsplit,
                     std::vector<uint32_t>& chosen, std::vector<Job>& jobs) {
    if (t <= 1 || P * lam > B || depth >= dsplit) {
        jobs.push_back(Job{chosen, P, lam, pmax, t});
        return;
    }
    mpz_class hi = iroot(B / P, (unsigned long)t);
    if (hi > PR.back()) {
        std::cerr << "FATAL: sieve too small in gen_jobs\n"; std::exit(2);
    }
    for (uint32_t q : PR) {
        if (q <= pmax) continue;
        if (q > hi) break;
        if (!admissible(q, chosen)) continue;
        mpz_class Pq = P * q;
        {
            chosen.push_back(q);
            mpz_class rmin;
            bool skip = r_min(q, t - 1, chosen, rmin) && Pq * rmin > B;
            chosen.pop_back();
            if (skip) continue;
        }
        mpz_class lamq;
        mpz_lcm_ui(lamq.get_mpz_t(), lam.get_mpz_t(), q - 1);
        chosen.push_back(q);
        gen_jobs(B, Pq, lamq, q, t - 1, depth + 1, dsplit, chosen, jobs);
        chosen.pop_back();
    }
}

// pick the smallest split depth that yields >= target jobs (cap depth).
static std::vector<Job> build_jobs(const mpz_class& B, int k, size_t target) {
    std::vector<Job> jobs;
    for (int dsplit = 1; dsplit <= 40; ++dsplit) {
        jobs.clear();
        std::vector<uint32_t> chosen;
        gen_jobs(B, mpz_class(1), mpz_class(1), 2, k, 0, dsplit, chosen, jobs);
        if (jobs.size() >= target) break;
    }
    return jobs;
}

static size_t g_tpw = 64;      // target jobs per worker (tunable via argv[4])
static long long g_maxjob_us = 0;  // wall of the single heaviest job (load-balance probe)
static double g_enum_s = 0;        // serial time spent building the job list (Amdahl probe)

// run one parallel pass over bound B with `nw` workers; accumulate stats.
static void run_pass(const mpz_class& B, int k, int nw, Stats& out, size_t& njobs) {
    {   // set shared bound = B for this pass
        std::lock_guard<std::mutex> lk(g_mtx);
        g_best = B;
        g_ver.fetch_add(1, std::memory_order_release);
    }
    auto e0 = std::chrono::steady_clock::now();
    std::vector<Job> jobs = build_jobs(B, k, std::max((size_t)256, g_tpw * (size_t)nw));
    g_enum_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - e0).count();
    njobs = jobs.size();
    std::atomic<size_t> idx{0};
    std::atomic<long long> maxjob_us{0};
    std::vector<Worker> workers(nw);
    std::vector<std::thread> pool;
    for (int wi = 0; wi < nw; ++wi) {
        pool.emplace_back([&, wi]() {
            Worker& w = workers[wi];
            w.init();
            for (;;) {
                size_t i = idx.fetch_add(1, std::memory_order_relaxed);
                if (i >= jobs.size()) break;
                w.init();  // resync bound at each job boundary
                Job& j = jobs[i];
                std::vector<uint32_t> chosen = j.chosen;  // per-worker mutable copy
                auto j0 = std::chrono::steady_clock::now();
                rec(j.P, j.lam, j.pmax, j.t, chosen, w);
                long long us = std::chrono::duration_cast<std::chrono::microseconds>(
                    std::chrono::steady_clock::now() - j0).count();
                long long prev = maxjob_us.load(std::memory_order_relaxed);
                while (us > prev && !maxjob_us.compare_exchange_weak(prev, us)) {}
            }
        });
    }
    for (auto& th : pool) th.join();
    for (auto& w : workers) out.add(w.S);
    g_maxjob_us = maxjob_us.load();
}

int main(int argc, char** argv) {
    if (argc < 2) { std::cerr << "usage: " << argv[0] << " k [bound] [workers]\n"; return 1; }
    g_k = std::stoi(argv[1]);
    if (g_k < 3) { std::cerr << "k must be >= 3\n"; return 1; }
    mpz_class user_bound(0);
    if (argc >= 3) user_bound = mpz_class(argv[2]);
    int nw = (argc >= 4) ? std::stoi(argv[3]) : (int)std::thread::hardware_concurrency();
    if (nw < 1) nw = 1;
    if (argc >= 5) g_tpw = (size_t)std::stoul(argv[4]);

    build_sieve();
    auto wall0 = std::chrono::steady_clock::now();

    Stats total;
    size_t njobs_last = 0;
    if (user_bound > 0) {
        run_pass(user_bound, g_k, nw, total, njobs_last);
    } else {
        mpz_class b = 1000;
        for (int i = 0; i < 130; ++i) {
            g_found = false;                 // fresh pass
            run_pass(b, g_k, nw, total, njobs_last);
            if (g_found) break;
            b *= 4;
        }
    }
    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();

    std::cout << "k=" << g_k << "\n";
    if (g_found) {
        std::cout << "n=" << g_best_n << "\nFACTORS:";
        for (const mpz_class& f : g_best_fac) std::cout << " " << f;
        std::cout << "\nnfactors=" << g_best_fac.size() << "\n";
    } else {
        std::cout << "n=NONE\nFACTORS:\nnfactors=0\n";
    }
    std::cout << "workers=" << nw << " jobs=" << njobs_last << " tpw=" << g_tpw << "\n";
    std::cout << "N_E=" << total.N_E << " N_I=" << total.N_I
              << " N_F=" << total.N_F << " N_H=" << total.N_H << " N_G=" << total.N_G << "\n";
    std::cout << "survivor_frac=" << (total.N_H > 0 ? (double)total.N_G / total.N_H : 0) << "\n";
    std::cout << "wall_s=" << wall << " maxjob_s=" << (g_maxjob_us / 1e6)
              << " enum_s=" << g_enum_s << "\n";
    return 0;
}
