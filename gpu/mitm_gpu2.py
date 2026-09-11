#!/usr/bin/env python3
"""
mitm_gpu2.py -- r-unbounded, chunkable generalisation of mitm_gpu.py.

Same exact-residue MITM join (key = mult * prod(subset) mod M as a u64, M < 2^62,
so equal key <=> equal residue), with two changes that lift the r<=6 cap:

  * per-side index width `bits`: ids pack r indices of `bits` bits each, so the
    constraint is r*bits <= 64 (D side over 64 base primes: bits=6 -> r<=10;
    I side over <=256 pool primes: bits=8 -> r<=8). Nothing else changes.
  * rank-range chunking: generate(..., lo, hi) emits only combinations with
    combinatorial-number-system rank in [lo, hi). C(64,8) ~ 4.4e9 records will
    not sort in one allocation, so the caller sorts rank-chunks independently and
    joins against each. Chunks partition the rank space exactly (no overlap, no
    gap), so the union of chunk joins is the full join.

Exact product / disjointness / Korselt are still re-checked in host big-int by
the caller against the frozen oracle. Nothing here is trusted without it.
"""
import numpy as np
import cupy as cp
from math import comb

_GEN_SRC = r"""
extern "C" __global__ void gen_records(
    const unsigned long long* __restrict__ primes, // pool primes, length n
    const unsigned long long* __restrict__ binom,  // (n+1)*(R+1) row-major C(a,b)
    const int n, const int R, const int bits,
    const unsigned long long M,
    const unsigned long long mult,   // 1 for D-side, (P0 mod M) for I-side
    const unsigned long long lo,     // first rank in this chunk
    const unsigned long long count,  // ranks in this chunk
    unsigned long long* __restrict__ out_keys,
    unsigned long long* __restrict__ out_ids)
{
    unsigned long long tid = blockIdx.x * (unsigned long long)blockDim.x + threadIdx.x;
    if (tid >= count) return;
    unsigned long long rank = lo + tid;

    // --- combinatorial-number-system unranking (lexicographic, increasing) ---
    int idx[10];
    int start = 0;
    const int stride = R + 1;
    for (int pos = 0; pos < R; ++pos) {
        int remaining = R - pos - 1;
        for (int c = start; c <= n - 1; ++c) {
            int a = n - 1 - c;
            unsigned long long cnt = binom[(unsigned long long)a * stride + remaining];
            if (rank < cnt) { idx[pos] = c; start = c + 1; break; }
            rank -= cnt;
        }
    }

    // --- exact residue: mult * prod(primes[idx]) mod M, via __int128 mulmod ---
    unsigned long long acc = mult % M;
    unsigned long long id = 0ULL;
    for (int pos = 0; pos < R; ++pos) {
        unsigned long long p = primes[idx[pos]];
        unsigned __int128 t = (unsigned __int128)acc * (p % M);
        acc = (unsigned long long)(t % M);
        id |= ((unsigned long long)idx[pos]) << (pos * bits);
    }
    out_keys[tid] = acc;
    out_ids[tid]  = id;
}
"""
_gen_kernel = cp.RawKernel(_GEN_SRC, "gen_records", options=("--device-int128",))


def _binom_table(n: int, r: int) -> np.ndarray:
    B = np.zeros((n + 1, r + 1), dtype=np.uint64)
    B[:, 0] = 1
    for a in range(1, n + 1):
        for b in range(1, min(a, r) + 1):
            B[a, b] = B[a - 1, b - 1] + B[a - 1, b]
    return B


def generate(primes, r: int, M: int, mult: int = 1, bits: int = 10, lo: int = 0, hi=None):
    """Residue records for the r-combinations of `primes` with CNS rank in [lo, hi).
    Returns cupy uint64 (keys, ids); ids pack r indices of `bits` bits each."""
    n = len(primes)
    assert 1 <= r <= 10 and r * bits <= 64, f"r={r}, bits={bits}: need r*bits<=64, r<=10"
    assert n <= (1 << bits), f"pool {n} does not fit {bits}-bit indices"
    assert M < (1 << 62), "M must stay < 2^62 so the residue key fits u64 uniquely"
    total = comb(n, r)
    if hi is None: hi = total
    assert 0 <= lo <= hi <= total, (lo, hi, total)
    count = hi - lo
    keys = cp.empty(count, dtype=cp.uint64); ids = cp.empty(count, dtype=cp.uint64)
    if count == 0: return keys, ids
    d_primes = cp.asarray(np.asarray(primes, dtype=np.uint64))
    d_binom = cp.asarray(_binom_table(n, r).reshape(-1))
    threads = 256; blocks = (count + threads - 1) // threads
    _gen_kernel((blocks,), (threads,),
                (d_primes, d_binom, np.int32(n), np.int32(r), np.int32(bits),
                 np.uint64(M), np.uint64(mult), np.uint64(lo), np.uint64(count), keys, ids))
    return keys, ids


def sort_side(keys, ids):
    """Sort one side by key (CUB radix via argsort). Returns (sorted_keys, sorted_ids)."""
    o = cp.argsort(keys, kind="stable")
    return keys[o], ids[o]


def join_sorted(Dk, Did, Ik, Iid):
    """Sort-merge join of two ALREADY-SORTED sides. Emits EVERY equal-key (D,I) pair."""
    nI = Ik.size
    lo = cp.searchsorted(Dk, Ik, side="left"); hi = cp.searchsorted(Dk, Ik, side="right")
    counts = (hi - lo).astype(cp.int64)
    total = int(counts.sum().item())
    if total == 0:
        return cp.empty(0, dtype=cp.uint64), cp.empty(0, dtype=cp.uint64)
    cum = cp.zeros(nI + 1, dtype=cp.int64); cp.cumsum(counts, out=cum[1:])
    Irow = cp.repeat(cp.arange(nI, dtype=cp.int64), counts)
    within = cp.arange(total, dtype=cp.int64) - cum[Irow]
    Dpos = lo[Irow] + within
    return Did[Dpos], Iid[Irow]


def decode_ids(ids, r: int, bits: int):
    """Unpack ids -> list of index tuples (host)."""
    if isinstance(ids, cp.ndarray): ids = ids.get()
    mask = (1 << bits) - 1
    return [tuple((int(v) >> (p * bits)) & mask for p in range(r)) for v in ids]
