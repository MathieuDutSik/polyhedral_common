#!/usr/bin/env python3
"""
The data of the tests of GENUS_GlueIsometry factors (RunCharPolyTests.sh).

For every genus LAT with an oracle table charpoly/oracle_LAT.txt (lines
"sig mass" from CharPolyOracle.g, sig = d1^b1,d2^b2,...) and a Gram matrix
charpoly/LAT.txt, write
- charpoly/LAT_sig.spec: the factors, "orth 1 b" for Phi_1^b, "orth -1 b"
  for Phi_2^b, "file charpoly/herm_d{d}_b{b}_N{N}_E{E}.txt" otherwise;
- charpoly/jobs.txt: the piece files to make with
  julia HermitianPieces.jl --jobs charpoly/jobs.txt
- charpoly_tests.txt: the lines "LAT sig mass".
For an odd lattice the pieces are those of its even sublattice (of the same
characteristic polynomial, its determinant 4 det).

The bounds of a piece P = Phi_d^b, rank b phi(d), of a lattice L of rank n
and determinant det: the glue with the rest is killed by R = prod_e
|Res(Phi_d, Phi_e)| over the other factors e, so det P divides N = det
R^min(rank P, n - rank P) and the exponent of D(P) divides E = e R, e the
exponent of D(L) (P^# is the projection of L^#, and the projection of L is
P up to the glue). The resultants are from the closed form (Apostol): for
a > b, |Res(Phi_a, Phi_b)| = p^phi(b) when a / b is a power of the prime
p, Phi_a(1) for b = 1 (p for a = p^k, 1 otherwise), 1 otherwise.

Usage (from CI_tests/15B_GenusEnumeration): python3 MakeCharPolyTests.py LAT...

The whole data of a genus LAT, from charpoly/LAT.txt (a Gram matrix) and
charpoly/genus_LAT.g (its classes and generators of their automorphism
groups, e.g. from GENUS_Enumerate and LATT_Automorphism):
  1. the oracle, in GAP:
       echo 'InFile:="charpoly/genus_LAT.g";; OutFile:="charpoly/oracle_LAT.txt";;
             Read("CharPolyOracle.g");' | gap -q
  2. python3 MakeCharPolyTests.py LAT ... (all the genera at once, as it
     writes charpoly_tests.txt and charpoly/jobs.txt)
  3. julia HermitianPieces.jl --jobs charpoly/jobs.txt (Hecke.jl; a file
     already made can be left out of the jobs)
  4. RunCharPolyTests.sh charpoly_tests.txt, and charpoly_tests_ci.txt the
     fast cases.
"""
import sys
from fractions import Fraction
from math import gcd


def phi(n):
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


def prime_of_power(n):
    if n < 2:
        return 0
    p = 2
    while n % p != 0:
        p += 1
    while n % p == 0:
        n //= p
    return p if n == 1 else 0


def res_cyclotomic(a, b):
    if a == b:
        raise ValueError("equal factors")
    if a < b:
        a, b = b, a
    if b == 1:
        p = prime_of_power(a)
        return p if p else 1
    if a % b == 0 and prime_of_power(a // b):
        return prime_of_power(a // b) ** phi(b)
    return 1


def det_int(M):
    M = [[Fraction(x) for x in row] for row in M]
    n = len(M)
    det = Fraction(1)
    for i in range(n):
        piv = next((k for k in range(i, n) if M[k][i] != 0), None)
        if piv is None:
            return 0
        if piv != i:
            M[i], M[piv] = M[piv], M[i]
            det = -det
        det *= M[i][i]
        for k in range(i + 1, n):
            f = M[k][i] / M[i][i]
            for j in range(i, n):
                M[k][j] -= f * M[i][j]
    return int(det)


def exponent_discriminant(M):
    """The exponent of D(L): the lcm of the denominators of M^-1."""
    n = len(M)
    aug = [[Fraction(x) for x in row] + [Fraction(int(i == j)) for j in range(n)]
           for i, row in enumerate(M)]
    for i in range(n):
        piv = next(k for k in range(i, n) if aug[k][i] != 0)
        aug[i], aug[piv] = aug[piv], aug[i]
        p = aug[i][i]
        aug[i] = [x / p for x in aug[i]]
        for k in range(n):
            if k != i and aug[k][i] != 0:
                f = aug[k][i]
                aug[k] = [a - f * b for a, b in zip(aug[k], aug[i])]
    e = 1
    for row in aug:
        for x in row[n:]:
            e = e * x.denominator // gcd(e, x.denominator)
    return e


def even_sublattice(A):
    """The Gram matrix of the even sublattice of an odd lattice, in the basis
    of EvenSublatticeBasis (IsometryGlueFactors.h); A itself when even."""
    n = len(A)
    odd = [i for i in range(n) if A[i][i] % 2 != 0]
    if not odd:
        return A
    i0 = odd[0]
    B = [[int(i == j) for j in range(n)] for i in range(n)]
    for j in odd:
        if j == i0:
            B[j][j] = 2
        else:
            B[j][i0] = 1
    BA = [[sum(B[i][k] * A[k][j] for k in range(n)) for j in range(n)]
          for i in range(n)]
    return [[sum(BA[i][k] * B[j][k] for k in range(n)) for j in range(n)]
            for i in range(n)]


def read_gram(path):
    words = open(path).read().split()
    n = int(words[0])
    vals = list(map(int, words[2:2 + n * n]))
    return [vals[i * n:(i + 1) * n] for i in range(n)]


def main():
    jobs = {}
    tests = []
    for lat in sys.argv[1:]:
        # an odd lattice: GENUS_GlueIsometry glues the pieces of its even
        # sublattice, so the bounds are those of the even sublattice
        A = even_sublattice(read_gram(f"charpoly/{lat}.txt"))
        n = len(A)
        det = det_int(A)
        expo = exponent_discriminant(A)
        for line in open(f"charpoly/oracle_{lat}.txt"):
            sig, mass = line.split()
            factors = [tuple(map(int, f.split("^"))) for f in sig.split(",")]
            spec = []
            for d, b in factors:
                if d == 1:
                    spec.append(f"orth 1 {b}")
                elif d == 2:
                    spec.append(f"orth -1 {b}")
                else:
                    R = 1
                    for e, _ in factors:
                        if e != d:
                            R *= res_cyclotomic(d, e)
                    rank = b * phi(d)
                    N = det * R ** min(rank, n - rank)
                    E = expo * R
                    name = f"charpoly/herm_d{d}_b{b}_N{N}_E{E}.txt"
                    jobs[name] = (d, b, N, E)
                    spec.append(f"file {name}")
            sig_name = sig.replace("^", "p").replace(",", "_")
            with open(f"charpoly/{lat}_{sig_name}.spec", "w") as f:
                f.write("\n".join(spec) + "\n")
            tests.append(f"{lat} {sig_name} {mass}")
    with open("charpoly/jobs.txt", "w") as f:
        for name, (d, b, N, E) in sorted(jobs.items()):
            f.write(f"{d} {b} {N} {E} {name}\n")
    with open("charpoly_tests.txt", "w") as f:
        f.write("# LAT sig mass: the oracle of CharPolyOracle.g\n")
        f.write("\n".join(tests) + "\n")
    print(f"{len(tests)} tests, {len(jobs)} piece files")


if __name__ == "__main__":
    main()
