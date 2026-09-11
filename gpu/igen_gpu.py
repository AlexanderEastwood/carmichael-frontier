#!/usr/bin/env python3
"""igen_gpu.py -- GPU generator for the FILTERED insertion side of the exchange MITM.

Emits the same records as src/mitm64_idump_stream2 with IDUMP_RANK=1: one (key, id) per r-subset I of
the ascending insertion pool INS with  prod(I) < Icap, at most HMAX primes >= SPLIT, where
    key = (P0 * prod(I)) mod M      (u64, M < 2^62, exact via 128-bit mulmod)
    id  = colex rank = sum_j C(i_j, j+1)   (u64; the caller guarantees C(n, r) < 2^64)
Pruning uses float64 log2 sums with a conservative margin: a branch is cut only when the bound is
exceeded by more than EPS, so the GPU set is a SUPERSET of the exact-integer set (extras have exact
product >= Icap and are discarded by the exact re-check of every join match). Thread per (i,j,k)
prefix on a 2-D grid (blockIdx.x=i, blockIdx.y=j, threadIdx.x=k; non-increasing triples exit at once),
then an in-thread iterative DFS over the remaining r-3 positions with the same break rules as the
C++ dumper (ascending pool => a failed bound ends the level). Two passes: count, then write.
Supports 3 <= r <= 8 and n <= 1024.
"""
import math
import numpy as np
import cupy as cp

EPS = 1e-7

_SRC = r"""
extern "C" __global__ void igen(
    const double* __restrict__ lp,          // log2 INS[x]
    const double* __restrict__ ps,          // prefix sums of lp, length n+1
    const unsigned long long* __restrict__ pm,   // INS[x] mod M
    const unsigned long long* __restrict__ binom,// (n+1)*(R+1) row-major C(a,b)
    const int n, const int R, const unsigned long long M, const unsigned long long P0m,
    const double LB, const int split_idx, const int hmax, const int mode,
    const int i_lo, const int i_hi,           // only prefixes with i in [i_lo, i_hi)
    unsigned long long* __restrict__ counter, // mode 0: per-i counters (length n); mode 1: single write cursor
    const unsigned long long cap,
    unsigned long long* __restrict__ out_keys, unsigned long long* __restrict__ out_ids)
{
    const int i = blockIdx.x, j = blockIdx.y, k = blockIdx.z * blockDim.x + threadIdx.x;
    if (i < i_lo || i >= i_hi) return;
    if (!(i < j && j < k) || k >= n) return;
    const int stride = R + 1;
    int idx[9]; double part[10]; unsigned long long acc[10], rid[10]; int nl[10];
    part[0] = 0.0; acc[0] = P0m % M; rid[0] = 0ULL; nl[0] = 0;
    // fixed prefix i, j, k with the same rules as the DFS (return instead of break)
    int pre[3] = {i, j, k};
    for (int d = 0; d < 3; ++d) {
        const int x = pre[d]; const int t = R - d - 1;
        if (x > n - 1 - t) return;
        const double np = part[d] + lp[x];
        if (np > LB) return;
        if (t > 0 && np + (ps[x + 1 + t] - ps[x + 1]) > LB) return;
        if (x >= split_idx && nl[d] + 1 > hmax) return;
        idx[d] = x; part[d + 1] = np; nl[d + 1] = nl[d] + (x >= split_idx);
        unsigned __int128 tt = (unsigned __int128)acc[d] * pm[x]; acc[d + 1] = (unsigned long long)(tt % M);
        rid[d + 1] = rid[d] + binom[(unsigned long long)x * stride + (d + 1)];
    }
    #define EMIT(KEY, ID) { if (mode == 0) atomicAdd(counter + i, 1ULL); else { unsigned long long pos = atomicAdd(counter, 1ULL); if (pos < cap) { out_keys[pos] = (KEY); out_ids[pos] = (ID); } } }
    if (R == 3) { EMIT(acc[3], rid[3]); return; }
    int d = 3; idx[3] = k + 1;
    while (true) {
        const int t = R - d - 1; const int x = idx[d];
        bool cut = false;
        if (x > n - 1 - t) cut = true;
        else {
            const double np = part[d] + lp[x];
            if (np > LB) cut = true;
            else if (t > 0 && np + (ps[x + 1 + t] - ps[x + 1]) > LB) cut = true;
            else if (x >= split_idx && nl[d] + 1 > hmax) cut = true;
            else {
                part[d + 1] = np; nl[d + 1] = nl[d] + (x >= split_idx);
                unsigned __int128 tt = (unsigned __int128)acc[d] * pm[x]; acc[d + 1] = (unsigned long long)(tt % M);
                rid[d + 1] = rid[d] + binom[(unsigned long long)x * stride + (d + 1)];
                if (d + 1 == R) { EMIT(acc[d + 1], rid[d + 1]); idx[d] = x + 1; continue; }
                d += 1; idx[d] = x + 1; continue;
            }
        }
        // this level is exhausted (ascending pool): back up
        d -= 1; if (d < 3) return; idx[d] += 1;
    }
}
"""
_kernel = cp.RawKernel(_SRC, "igen", options=("--device-int128",))

def _binom_table(n, r):
    B = np.zeros((n + 1, r + 1), dtype=object); B[:, 0] = 1
    for a in range(1, n + 1):
        for b in range(1, min(a, r) + 1): B[a, b] = B[a - 1, b - 1] + B[a - 1, b]
    assert B[n, r] < 2 ** 64, "rank ids would overflow u64"
    return B.astype(np.uint64)

def generate_chunks(INS, r, M, P0, Icap, split=None, hmax=None, max_records=200_000_000):
    """Yield (keys, ids) cupy uint64 chunks (each <= max_records unless a single first index exceeds it)
    covering exactly the filtered r-subsets of INS. Total count is available as the generator's return
    value via StopIteration.value; simpler: use generate() for small instances."""
    n = len(INS); assert 3 <= r <= 8 and n <= 1024 and M < (1 << 62)
    assert all(INS[i] < INS[i + 1] for i in range(n - 1)), "INS must be strictly ascending"
    if Icap <= 1: return 0
    lp = np.array([math.log2(p) for p in INS], dtype=np.float64)
    ps = np.concatenate([[0.0], np.cumsum(lp)]).astype(np.float64)
    LB = math.log2(Icap) + EPS
    split_idx = n if split is None else next((i for i, p in enumerate(INS) if p > split), n)
    hmax = r if hmax is None else hmax
    d_lp = cp.asarray(lp); d_ps = cp.asarray(ps); d_pm = cp.asarray(np.array([p % M for p in INS], dtype=np.uint64))
    d_b = cp.asarray(_binom_table(n, r).reshape(-1)); dummy = cp.empty(1, dtype=cp.uint64)
    block = (256,); grid = (n, n, (n + 255) // 256)
    base = [d_lp, d_ps, d_pm, d_b, np.int32(n), np.int32(r), np.uint64(M), np.uint64(P0 % M), np.float64(LB), np.int32(split_idx), np.int32(hmax)]
    per_i = cp.zeros(n, dtype=cp.uint64)
    _kernel(grid, block, tuple(base + [np.int32(0), np.int32(0), np.int32(n), per_i, np.uint64(0), dummy, dummy]))
    counts = per_i.get().astype(np.int64); total = int(counts.sum())
    a = 0
    while a < n:
        b = a + 1; acc = int(counts[a])
        while b < n and acc + int(counts[b]) <= max_records: acc += int(counts[b]); b += 1
        if acc > 0:
            keys = cp.empty(acc, dtype=cp.uint64); ids = cp.empty(acc, dtype=cp.uint64); cur = cp.zeros(1, dtype=cp.uint64)
            _kernel(grid, block, tuple(base + [np.int32(1), np.int32(a), np.int32(b), cur, np.uint64(acc), keys, ids]))
            assert int(cur.get()[0]) == acc, (a, b, int(cur.get()[0]), acc)
            yield keys, ids
        a = b
    return total

def generate(INS, r, M, P0, Icap, split=None, hmax=None):
    """Single-shot convenience: concatenated (keys, ids, count). Use generate_chunks for big instances."""
    ks, is_ = [], []
    for k, i in generate_chunks(INS, r, M, P0, Icap, split, hmax, max_records=1 << 62): ks.append(k); is_.append(i)
    if not ks: return cp.empty(0, cp.uint64), cp.empty(0, cp.uint64), 0
    k = cp.concatenate(ks); i = cp.concatenate(is_); return k, i, int(k.size)
