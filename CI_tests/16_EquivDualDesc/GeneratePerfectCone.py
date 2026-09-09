#!/usr/bin/env python3
"""Generate the EXT files of the perfect cones of the root lattices E7 and E8.

The perfect domain (Voronoi cone) of a perfect form Q is the cone spanned by
the rank one forms v v^T for v running over the minimal vectors of Q. As Q is
perfect those rank one forms span the whole space of symmetric matrices, so the
cone is full dimensional in Sym^n(R), of dimension n(n+1)/2.

All the rank one forms satisfy <Q, v v^T> = min(Q). Therefore the rays lie in
an affine hyperplane and the cone is the homogenization of a polytope. We use
that to emit the usual "first column equal to 1" EXT format: the coordinate
X(0,0) is traded for the constant functional <Q, X> / min(Q).

Coordinates on Sym^n(R) follow the convention of SymmetricMatrixToVector in
basic_common_cpp/src_matrix/MAT_MatrixFund.h, that is the pairs (i,j) with
j <= i ordered lexicographically by (i,j).

Usage: GeneratePerfectCone.py E7 PerfectE7/PerfectE7.ext
       GeneratePerfectCone.py E8 PerfectE8/PerfectE8.ext
"""
import sys
from fractions import Fraction
from itertools import combinations, product


def e8_roots():
    """The 240 roots of E8, of norm 2, in the standard coordinates."""
    roots = []
    for i, j in combinations(range(8), 2):
        for si, sj in product([1, -1], repeat=2):
            v = [Fraction(0)] * 8
            v[i] = Fraction(si)
            v[j] = Fraction(sj)
            roots.append(tuple(v))
    for signs in product([1, -1], repeat=8):
        if signs.count(-1) % 2 == 0:
            roots.append(tuple(Fraction(s, 2) for s in signs))
    assert len(roots) == 240
    return roots


def e8_simple_roots():
    """The Bourbaki simple roots of E8, a Z-basis of the root lattice."""
    def unit(i):
        v = [Fraction(0)] * 8
        v[i] = Fraction(1)
        return v

    def comb(coefs):
        return tuple(sum((Fraction(c) * u[k] for c, u in coefs), Fraction(0))
                     for k in range(8))

    a1 = tuple(Fraction(s, 2) for s in [1, -1, -1, -1, -1, -1, -1, 1])
    a2 = comb([(1, unit(0)), (1, unit(1))])
    rest = [comb([(1, unit(k)), (-1, unit(k - 1))]) for k in range(1, 7)]
    return [a1, a2] + rest


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def solve(basis, v):
    """Coordinates of v on the given basis, both expressed as Gram data."""
    n = len(basis)
    gram = [[dot(basis[i], basis[j]) for j in range(n)] for i in range(n)]
    rhs = [dot(basis[i], v) for i in range(n)]
    return gauss(gram, rhs)


def gauss(mat, rhs):
    n = len(mat)
    a = [row[:] + [rhs[i]] for i, row in enumerate(mat)]
    for c in range(n):
        piv = next(r for r in range(c, n) if a[r][c] != 0)
        a[c], a[piv] = a[piv], a[c]
        inv = Fraction(1) / a[c][c]
        a[c] = [x * inv for x in a[c]]
        for r in range(n):
            if r != c and a[r][c] != 0:
                f = a[r][c]
                a[r] = [x - f * y for x, y in zip(a[r], a[c])]
    return [a[r][n] for r in range(n)]


def determinant(mat):
    n = len(mat)
    a = [row[:] for row in mat]
    det = Fraction(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if a[r][c] != 0), None)
        if piv is None:
            return Fraction(0)
        if piv != c:
            a[c], a[piv] = a[piv], a[c]
            det = -det
        det *= a[c][c]
        inv = Fraction(1) / a[c][c]
        for r in range(c + 1, n):
            f = a[r][c] * inv
            if f != 0:
                a[r] = [x - f * y for x, y in zip(a[r], a[c])]
    return det


def lattice_data(name):
    """Return (Gram matrix of the chosen Z-basis, minimal vectors in it)."""
    roots = e8_roots()
    simple = e8_simple_roots()
    coords = [solve(simple, r) for r in roots]
    for c in coords:
        assert all(x.denominator == 1 for x in c), "non integral expansion"
    if name == "E8":
        sel = coords
        dim = 8
    elif name == "E7":
        sel = [c[:7] for c in coords if c[7] == 0]
        dim = 7
    else:
        raise ValueError("unknown lattice " + name)
    gram = [[dot(simple[i], simple[j]) for j in range(dim)] for i in range(dim)]
    assert len(sel) == {7: 126, 8: 240}[dim]
    assert determinant(gram) == {7: 2, 8: 1}[dim]
    return gram, [tuple(int(x) for x in c) for c in sel]


def pair_representatives(vects):
    """One vector out of each antipodal pair, the first nonzero entry > 0."""
    out = []
    for v in vects:
        nz = next(x for x in v if x != 0)
        if nz > 0:
            out.append(v)
    assert 2 * len(out) == len(vects)
    return out


def build_ext(gram, vects):
    """The EXT matrix of the perfect cone, in homogeneous polytope format."""
    n = len(gram)
    minimum = min(sum(gram[i][j] * v[i] * v[j] for i in range(n)
                      for j in range(n)) for v in vects)
    idx = [(i, j) for i in range(n) for j in range(i + 1) if (i, j) != (0, 0)]
    rows = []
    for v in vects:
        norm = sum(gram[i][j] * v[i] * v[j] for i in range(n) for j in range(n))
        assert norm == minimum
        rows.append([1] + [v[i] * v[j] for (i, j) in idx])
    assert len(rows[0]) == (n * (n + 1)) // 2
    return rows


def main():
    name, outfile = sys.argv[1], sys.argv[2]
    gram, vects = lattice_data(name)
    rows = build_ext(gram, pair_representatives(vects))
    with open(outfile, "w") as f:
        f.write("%d %d\n" % (len(rows), len(rows[0])))
        for row in rows:
            f.write(" ".join(str(x) for x in row) + "\n")
    print("%s: %d vertices, %d columns" % (name, len(rows), len(rows[0])))


if __name__ == "__main__":
    main()
