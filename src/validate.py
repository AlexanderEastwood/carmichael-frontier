#!/usr/bin/env python3
"""Thin driver over the FROZEN oracle ref/ref_carmichael.py. Does NOT reimplement
anything -- it imports and calls the oracle so every engine result is checked by
independent code.

  validate.py verify P1 P2 ...     -> prints "OK <count>" or "FAIL <reason>"
  validate.py brute K BOUND        -> prints oracle least_with_k(K, BOUND)
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ref"))
import ref_carmichael as ref  # frozen oracle


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: validate.py verify P1 P2 ... | brute K BOUND", file=sys.stderr)
        return 1
    cmd = sys.argv[1]
    if cmd == "verify":
        primes = [int(x) for x in sys.argv[2:]]
        ok, n, why = ref.verify_certificate(primes)
        if ok:
            print(f"OK {len(primes)} n={n} :: {why}")
            return 0
        print(f"FAIL {why}")
        return 3
    if cmd == "brute":
        k = int(sys.argv[2])
        bound = int(sys.argv[3])
        n, fs = ref.least_with_k(k, bound)
        if n is None:
            print("NONE")
        else:
            print(f"{n} {' '.join(map(str, fs))}")
        return 0
    print(f"unknown command {cmd}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
