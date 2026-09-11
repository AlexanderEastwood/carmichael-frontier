# A reproducible exclusion: S65 >= 10^148

This is a computational lower bound, not a determination of S65 and not
a claim of novelty. It was obtained while calibrating the next search.

## Python reproduction

Python 3.10+, standard library only:

```sh
python s65_exclusion.py
```

It first recomputes the certified cofactor bound and all 127,092 admissible
64-prime products below 10^145. Their global minimum A64 is a lower bound
for every 64-prime cofactor of a 65-factor Carmichael. Thus any such number
below 10^148 has largest prime <=2266, a universe of 334 odd primes.

An exhaustive ordered-prefix search with exact product pruning stops
descending when P*lambda(P)>10^148 and checks the at-most-one cofactor in
the required residue class. It enforces exactly 65 distinct factors and
tests full Korselt with exact primality from the sieve. Expected output:

- 64,355 search nodes;
- 5,012 modular inversions;
- 66 progression candidates;
- maximum remaining factor count at the progression step: 10;
- no Carmichael numbers;
- `CERTIFIED: S65 >= 10**148`.

Initial timing was about 0.35 seconds of search, excluding the roughly
9-second cofactor setup. Timings vary by machine and contention.

The bundled s65_transfer_certificate.py also independently runs its earlier
3.17*10^147 exclusion if invoked directly. Its functions supply the cofactor
setup here; the setup never assumes that a cofactor is itself Carmichael.

## Direct native cross-check

```sh
g++ -O3 -std=c++17 verify65.cpp -o verify65
./verify65 148 45
```

This uses the certified prime cap 2266, an independently generated sieve,
checked 512-bit products, and binary admissible-product enumeration. It
does not use modular inverses or the progression shortcut. At two remaining
primes it visits each compatible pair individually. For each complete
product it tests Korselt, stopping the per-prime checks at the first failure.

Expected final output includes:

- `COMPLETE 1`;
- `leaves 22287138`;
- `carmichael_hits 0`;
- `calls 1804229`;
- `nodes 51969795` (the equivalent uncollapsed binary-tree count, not the
  actual number of C++ function calls).

Measured search time: about 1.73 seconds. This cross-check relies on the
Python certificate for the prime bound; it is not a separate end-to-end
derivation of that bound. All verification is computational, not a formal
proof-assistant development.

## Budgeted next-decade calibration

```sh
python s65_exclusion.py --exponent 149 --seconds 30
```

A 30-second sample at 10^149 was INCOMPLETE after 4,341,760 nodes,
956,560 inversions and 13,950 progression candidates. No Carmichael was
found in that partial run. This does not prove a bound at 10^149 and does
not give the total runtime. The full exact minimum remains undetermined.
