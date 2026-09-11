#!/usr/bin/env python3
"""igen_gpu.py -- GPU generator for the FILTERED insertion side of the exchange MITM (exact arithmetic).

Emits exactly the records of src/mitm64_idump_stream2 with IDUMP_RANK=1: one (key, id) per r-subset I of
the ascending insertion pool INS = a_0 < ... < a_{n-1} with prod(I) < Icap and at most HMAX primes > SPLIT:
    key = (P0 * prod(I)) mod M      (u64, M < 2^64, exact via 128-bit mulmod)
    id  = colex rank = sum_j C(i_j, j+1)   (u64; the caller guarantees C(n, r) < 2^64)

Exact pruning (Astra, 2026-09-11) -- no floating point. With t positions still to fill after choosing a_x
and A the exact partial product before that choice, the host table
    H[t][x] = floor((Icap - 1) / (a_x * a_{x+1} * ... * a_{x+t}))        (x + t < n)
combines the product cap and the cheapest-completion bound: the consecutive window is the smallest
completion starting at x and grows with x, so  if A > H[t][x]: end this level;  otherwise
A * a_x <= A * a_x ... a_{x+t} <= Icap - 1, hence the next product cannot overflow 128 bits.
Requires Icap - 1 < 2^128 (checked; wider instances must use the CPU dumper).

Work decomposition: one thread per (i, j, k) prefix (grid (b-a) x n x ceil(n/256) for first indices
i in [a, b)); the thread runs an iterative DFS over the remaining r-3 positions. Pass 1 counts records
per prefix (no atomics, no residue/rank work); the host scans counts into offsets; pass 2 writes each
prefix's records into its own interval, and count/write agreement is checked per prefix.
Supports 3 <= r <= 8, n <= 1024, M < 2^64.
"""
import numpy as np
import cupy as cp

_SRC = r"""
typedef unsigned long long u64; typedef unsigned __int128 u128;
__device__ __forceinline__ u128 ld128(const u64* p, u64 idx) { return ((u128)p[2*idx+1] << 64) | (u128)p[2*idx]; }

extern "C" __global__ void igen(
    const u64* __restrict__ primes,  // INS[x]
    const u64* __restrict__ pm,      // INS[x] mod M
    const u64* __restrict__ H,       // R*n entries of u128 as (lo,hi) pairs: H[t*n + x], t = positions left after x
    const u64* __restrict__ binom,   // (n+1)*(R+1) row-major C(a,b)
    const int i0, const int n, const int R, const u64 M, const u64 P0m, const int split_idx, const int hmax,
    const int mode,                  // 0: count into cnt[prefix]; 1: write from off[prefix], record count in wrote[prefix]
    u64* __restrict__ cnt, const u64* __restrict__ off, u64* __restrict__ wrote,
    u64* __restrict__ out_keys, u64* __restrict__ out_ids)
{
    const int i = blockIdx.x + i0, j = blockIdx.y, k = blockIdx.z * blockDim.x + threadIdx.x;
    if (!(i < j && j < k) || k >= n) return;
    const u64 pid = ((u64)(i - i0) * n + j) * n + k;    // prefix id relative to the block's first index
    const int stride = R + 1;
    int idx[9]; u128 A[10]; u64 acc[10], rid[10]; int nl[10];
    A[0] = 1; acc[0] = P0m % M; rid[0] = 0ULL; nl[0] = 0;
    int pre[3] = {i, j, k};
    for (int d = 0; d < 3; ++d) {
        const int x = pre[d]; const int t = R - d - 1;
        if (x > n - 1 - t) return;
        if (A[d] > ld128(H, (u64)t * n + x)) return;
        if (x >= split_idx && nl[d] + 1 > hmax) return;
        idx[d] = x; A[d + 1] = A[d] * (u128)primes[x]; nl[d + 1] = nl[d] + (x >= split_idx);
        u128 tt = (u128)acc[d] * pm[x]; acc[d + 1] = (u64)(tt % M);
        rid[d + 1] = rid[d] + binom[(u64)x * stride + (d + 1)];
    }
    u64 c = 0; u64 pos = (mode == 1) ? off[pid] : 0;
    #define EMIT(KEY, ID) { if (mode == 1) { out_keys[pos] = (KEY); out_ids[pos] = (ID); ++pos; } ++c; }
    if (R == 3) { EMIT(acc[3], rid[3]); }
    else {
        int d = 3; idx[3] = k + 1;
        while (true) {
            const int t = R - d - 1; const int x = idx[d];
            const bool cut = (x > n - 1 - t) || (A[d] > ld128(H, (u64)t * n + x)) || (x >= split_idx && nl[d] + 1 > hmax);
            if (!cut) {
                A[d + 1] = A[d] * (u128)primes[x]; nl[d + 1] = nl[d] + (x >= split_idx);
                u128 tt = (u128)acc[d] * pm[x]; acc[d + 1] = (u64)(tt % M);
                rid[d + 1] = rid[d] + binom[(u64)x * stride + (d + 1)];
                if (d + 1 == R) { EMIT(acc[d + 1], rid[d + 1]); idx[d] = x + 1; continue; }
                d += 1; idx[d] = x + 1; continue;
            }
            d -= 1; if (d < 3) break; idx[d] += 1;      // level exhausted (ascending pool): back up
        }
    }
    if (mode == 0) cnt[pid] = c; else wrote[pid] = c;
}
"""
_kernel = cp.RawKernel(_SRC, "igen", options=("--device-int128",))

def _binom_table(n, r):
    B = np.zeros((n + 1, r + 1), dtype=object); B[:, 0] = 1
    for a in range(1, n + 1):
        for b in range(1, min(a, r) + 1): B[a, b] = B[a - 1, b - 1] + B[a - 1, b]
    assert B[n, r] < 2 ** 64, "rank ids would overflow u64"
    return B.astype(np.uint64)

def _u128_pairs(vals):
    a = np.zeros(2 * len(vals), dtype=np.uint64); mask = (1 << 64) - 1
    for i, v in enumerate(vals): a[2 * i] = v & mask; a[2 * i + 1] = v >> 64
    return a

def quotient_table(INS, r, Icap):
    """H[t][x] = floor((Icap-1) / (a_x ... a_{x+t})) for 0 <= t <= r-1 and x + t < n (0 otherwise); exact."""
    n = len(INS); top = Icap - 1; H = [0] * (r * n)
    for t in range(r):
        w = 1
        for x in range(t + 1): w *= INS[x]
        for x in range(n - t):
            if x > 0: w = (w // INS[x - 1]) * INS[x + t]
            H[t * n + x] = top // w
    return H

def _prepare(INS, r, M, P0, Icap, split, hmax):
    n = len(INS); assert 3 <= r <= 8 and n <= 1024 and 1 <= M < (1 << 64)
    assert all(INS[i] < INS[i + 1] for i in range(n - 1)), "INS must be strictly ascending"
    assert Icap - 1 < (1 << 128), "Icap needs more than 128 bits: use the CPU dumper"
    split_idx = n if split is None else next((i for i, p in enumerate(INS) if p > split), n)   # primes > SPLIT, as the C++ dumper
    hmax = r if hmax is None else hmax
    d_primes = cp.asarray(np.array(INS, dtype=np.uint64)); d_pm = cp.asarray(np.array([p % M for p in INS], dtype=np.uint64))
    d_H = cp.asarray(_u128_pairs(quotient_table(INS, r, Icap))); d_b = cp.asarray(_binom_table(n, r).reshape(-1))
    base = [d_primes, d_pm, d_H, d_b]; tail = [np.int32(n), np.int32(r), np.uint64(M), np.uint64(P0 % M), np.int32(split_idx), np.int32(hmax)]
    return n, base, tail

def generate_chunks(INS, r, M, P0, Icap, split=None, hmax=None, max_records=200_000_000, block_bytes=512 << 20):
    """Yield (keys, ids) cupy uint64 chunks covering exactly the filtered r-subsets of INS; the generator's
    return value (StopIteration.value) is the total record count. First indices are processed in blocks
    [a, b) sized so the per-prefix count/offset arrays ((b-a)*n*n u64) stay under block_bytes; within a block
    the records are written in one pass (a block may exceed max_records if a single first index does)."""
    if Icap <= 1: return 0
    n, base, tail = _prepare(INS, r, M, P0, Icap, split, hmax)
    blk = (256,); zb = (n + 255) // 256; dummy = cp.empty(1, dtype=cp.uint64)
    ib = max(1, min(n, block_bytes // (8 * n * n))); total = 0; a = 0
    while a < n:
        b = min(n, a + ib); m = (b - a) * n * n
        cnt = cp.zeros(m, dtype=cp.uint64)
        _kernel((b - a, n, zb), blk, tuple(base + [np.int32(a)] + tail + [np.int32(0), cnt, dummy, dummy, dummy, dummy]))
        # prefix ids in the kernel are global ((i*n+j)*n+k); shift so the block's arrays start at 0
        per_i = cnt.reshape(b - a, n * n).sum(axis=1).get().astype(np.int64)
        # sub-split the block by cumulative record count so that no single write exceeds max_records
        # (a single first index may still exceed it; then it is written alone)
        a2 = a
        while a2 < b:
            b2 = a2 + 1; acc = int(per_i[a2 - a])
            while b2 < b and acc + int(per_i[b2 - a]) <= max_records: acc += int(per_i[b2 - a]); b2 += 1
            if acc:
                sl = slice((a2 - a) * n * n, (b2 - a) * n * n); csl = cnt[sl]
                off = cp.zeros(csl.size, dtype=cp.uint64); off[:] = cp.cumsum(csl) - csl; wrote = cp.zeros(csl.size, dtype=cp.uint64)
                keys = cp.empty(acc, dtype=cp.uint64); ids = cp.empty(acc, dtype=cp.uint64)
                _kernel((b2 - a2, n, zb), blk, tuple(base + [np.int32(a2)] + tail + [np.int32(1), dummy, off, wrote, keys, ids]))
                assert bool((wrote == csl).all().get()), ("per-prefix count/write mismatch", a2, b2)
                del off, wrote; total += acc
                yield keys, ids
            a2 = b2
        del cnt
        a = b
    return total

def generate(INS, r, M, P0, Icap, split=None, hmax=None):
    """Single-shot convenience: concatenated (keys, ids, count)."""
    ks, is_ = [], []
    for k, i in generate_chunks(INS, r, M, P0, Icap, split, hmax, max_records=1 << 62): ks.append(k); is_.append(i)
    if not ks: return cp.empty(0, cp.uint64), cp.empty(0, cp.uint64), 0
    k = cp.concatenate(ks); i = cp.concatenate(is_); return k, i, int(k.size)
