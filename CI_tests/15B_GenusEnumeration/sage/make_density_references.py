"""Reference values for the local density code of src_genera, from Sage.

    sage -python make_density_references.py <task> <output.g>

writes a GAP file  return [ [inputs..., value], ... ];  that the GAP test
TestLocalDensity.g compares the C++ program GENUS_LocalDensity against.
The values are computed here with Sage's own special values (zeta__exact,
quadratic_L_function__exact) and symbolic simplification, sharing no code
with the C++ side, which carries the transcendental parts exactly instead.

Tasks:
  generic_part   alpha_oo(S,T) times the product of alpha_p(S,T) over the
                 primes not dividing 2 det S det T, for S of rank m and
                 determinant detS, T of rank n and determinant detT
                 (conventions of the C++ header LocalDensity.h).
  siegel_series  Katsurada's polynomial F_p(S/2, X) for even Gram matrices S
                 (root lattices and random even forms) at p = 2, 3, 5, 7,
                 from the vendored package external/siegel_series.
  alpha_unimodular
                 alpha_p(S, T) for T unimodular at p of rank 2k and
                 character xi_T, assembled from that polynomial.
  rank_one_density
                 alpha_p(a, T) from Sage's QuadraticForm.local_density, for
                 random even forms T with structure at p, p = 2, 3, 5, 7.
  unit_chain     alpha_p(S, T) at an odd p dividing det T by splitting off
                 the unit diagonal entries of S (Sage's local_normal_form and
                 local_density), or fail when S has two or more non-unit
                 entries at p.
  certified_mass_ratio <reference_dir>
                 M_S / mass from completely enumerated genera: reads the
                 genus representatives and the certified representation
                 masses of the strongly perfect lattice project
                 (results/density_calib/reps and results/reference_MR).
  genus_mass     the mass of the genus of positive definite even forms
                 (root lattices and random forms with structure at one or
                 two primes, ranks 2 to 14), from Sage's Genus.mass, the
                 Conway-Sloane formula on Sage's own genus symbols, with
                 the canonical 2-adic symbol of each form.

The Euler factors of the good-prime density are those of LocalDensity.h;
they are checked independently, against the orders of the orthogonal
groups over F_p, by the GAP test itself.
"""

import os
import random
import sys

from sage.all import (QQ, SR, ZZ, CartanMatrix, block_diagonal_matrix,
                      fundamental_discriminant, gamma, identity_matrix,
                      kronecker_symbol, matrix, pi, prod,
                      PolynomialRing)
from sage.quadratic_forms.special_values import (quadratic_L_function__exact,
                                                 zeta__exact)


def euler_terms(n, m, detS, detT):
    if n % 2 == 0:
        r = n // 2
        terms = [(r, fundamental_discriminant((-1)**r * detT), 1)]
        if m % 2 == 0:
            rp = (n - m) // 2
            DW = fundamental_discriminant((-1)**rp * detS * detT)
            terms += [(2 * i, 1, 1) for i in range(rp, r)]
            terms += [(rp, DW, -1)]
        else:
            rp = (n - m - 1) // 2
            terms += [(2 * i, 1, 1) for i in range(rp + 1, r)]
        return terms
    r = (n - 1) // 2
    if m % 2 == 1:
        rp = (n - m) // 2
        DW = fundamental_discriminant((-1)**rp * detS * detT)
        terms = [(2 * i, 1, 1) for i in range(rp, r + 1)]
        terms += [(rp, DW, -1)]
    else:
        rp = (n - m - 1) // 2
        terms = [(2 * i, 1, 1) for i in range(rp + 1, r + 1)]
    return terms


def euler_value(terms, p):
    v = QQ(1)
    for (s, D, e) in terms:
        chi = 1 if D == 1 else kronecker_symbol(D, p)
        v *= (1 - QQ(chi) / QQ(p)**s)**e
    return v


def generic_part(n, m, detS, detT):
    terms = euler_terms(n, m, detS, detT)
    num = SR(pi)**(QQ(m * (2 * n - m + 1)) / 4)
    den = prod(gamma(QQ(j) / 2) for j in range(n - m + 1, n + 1))
    v = (num / den * SR(detS)**(QQ(n - m - 1) / 2)
         / SR(detT)**(QQ(m) / 2))
    for (s, D, e) in terms:
        L = zeta__exact(s) if D == 1 else quadratic_L_function__exact(s, D)
        v *= SR(L)**(-e)
    for p in sorted(set(ZZ(2 * detS * detT).prime_divisors())):
        v /= euler_value(terms, p)
    return QQ(v.canonicalize_radical())


def grid_generic_part():
    out = []
    for n in (8, 13, 14, 16):
        for m in (1, 2, 3, 4, 5, 7, 10):
            if n - m < 3:
                continue
            for detS in (1, 2, 3, 4, 5, 8, 12):
                for detT in (1, 3, 180, 216, 243, 351, 507):
                    out.append((n, m, detS, detT,
                                generic_part(n, m, detS, detT)))
    return out


def cartan(letter, rank):
    C = CartanMatrix([letter, rank])
    return matrix(ZZ, C.nrows(), C.ncols(), list(C))


def ade_types(max_rank):
    """All root lattice types of rank <= max_rank, as lists of (letter, rank)."""
    comps = [("A", r) for r in range(1, max_rank + 1)]
    comps += [("D", r) for r in range(4, max_rank + 1)]
    comps += [("E", r) for r in (6, 7, 8) if r <= max_rank]
    out = []

    def rec(start, rank_left, current):
        if current:
            out.append(list(current))
        for i in range(start, len(comps)):
            if comps[i][1] <= rank_left:
                current.append(comps[i])
                rec(i, rank_left - comps[i][1], current)
                current.pop()
    rec(0, max_rank, [])
    return out


def ade_gram(typ):
    return block_diagonal_matrix([cartan(l, r) for (l, r) in typ],
                                 subdivide=False)


def test_matrices():
    """Root lattices, and random even positive definite forms R^t E R."""
    mats = [ade_gram(t) for t in ade_types(6)]
    extra = [[("A", 7)], [("A", 8)], [("D", 7)], [("D", 8)], [("E", 7)],
             [("E", 8)], [("A", 1), ("E", 7)], [("A", 2), ("A", 1), ("D", 4)],
             [("A", 1)] * 7, [("A", 2)] * 4]
    mats += [ade_gram(t) for t in extra]
    rng = random.Random(20261026)
    for i in range(150):
        m = rng.randint(1, 6)
        kind = rng.choice(["2I", "A", "2A", "6I", "D", "10I"])
        if kind == "2I":
            E = 2 * identity_matrix(ZZ, m)
        elif kind == "6I":
            E = 6 * identity_matrix(ZZ, m)
        elif kind == "10I":
            E = 10 * identity_matrix(ZZ, m)
        elif kind == "A":
            E = cartan("A", m)
        elif kind == "2A":
            E = 2 * cartan("A", m)
        else:
            E = cartan("D", m) if m >= 4 else cartan("A", m)
        while True:
            R = matrix(ZZ, m, m, [rng.randint(-2, 2) for _ in range(m * m)])
            if R.det() != 0:
                break
        mats.append(R.transpose() * E * R)
    return mats


def load_siegel_series():
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.join(here, "external"))
    from siegel_series.impl import siegel_series_polynomial
    return siegel_series_polynomial


def poly_coeffs(F):
    R = PolynomialRing(QQ, "x")
    G = R.fraction_field()(F)
    num, den = R(G.numerator()), R(G.denominator())
    assert den.degree() == 0, "the Siegel series is not a polynomial"
    P = num / den.constant_coefficient()
    return [QQ(c) for c in R(P).list()]


def grid_siegel_series():
    ssp = load_siegel_series()
    out = []
    for S in test_matrices():
        for p in (2, 3, 5, 7):
            F = ssp(S / 2, ZZ(p))
            out.append((p, S.nrows(), list(S.list()), poly_coeffs(F)))
    return out


def alpha_unimodular(ssp, S, p, k, xi_T):
    m = S.nrows()
    p = ZZ(p)
    X = QQ(xi_T) / p**k
    coeffs = poly_coeffs(ssp(S / 2, p))
    Fv = sum(c * X**i for i, c in enumerate(coeffs))
    pre = (1 - X) * prod(1 - p**(2 * i) * X**2 for i in range(1, m // 2 + 1))
    if m % 2 == 0:
        D = fundamental_discriminant((-1)**(m // 2) * S.det())
        pre /= (1 - p**(m // 2) * kronecker_symbol(D, p) * X)
    if p == 2:
        pre *= ZZ(2)**(m - 1)
    return pre * Fv


def grid_alpha_unimodular():
    ssp = load_siegel_series()
    out = []
    cases = [(2, 4, 1), (2, 7, 1), (2, 7, -1), (2, 8, 1), (3, 7, 1),
             (3, 7, -1), (5, 7, 1), (5, 7, -1)]
    for t in ade_types(6):
        S = ade_gram(t)
        for (p, k, xi) in cases:
            if 2 * k - S.nrows() >= 3:
                out.append((p, k, xi, S.nrows(), list(S.list()),
                            alpha_unimodular(ssp, S, p, k, xi)))
    return out


def random_even_form(rng, n, p):
    """R^t E R with E an even form carrying structure at p."""
    blocks = []
    rest = n
    while rest > 0:
        kind = rng.choice(["A1", "A2", "pA1", "pA2", "p2A1", "A3"])
        if kind in ("A2", "pA2") and rest < 2:
            kind = "A1"
        if kind == "A3" and rest < 3:
            kind = "A1"
        if kind == "A1":
            blocks.append(cartan("A", 1)); rest -= 1
        elif kind == "A2":
            blocks.append(cartan("A", 2)); rest -= 2
        elif kind == "A3":
            blocks.append(cartan("A", 3)); rest -= 3
        elif kind == "pA1":
            blocks.append(p * cartan("A", 1)); rest -= 1
        elif kind == "pA2":
            blocks.append(p * cartan("A", 2)); rest -= 2
        else:
            blocks.append(p * p * cartan("A", 1)); rest -= 1
    E = block_diagonal_matrix(blocks, subdivide=False)
    while True:
        R = matrix(ZZ, n, n, [rng.randint(-1, 1) for _ in range(n * n)])
        if abs(R.det()) == 1:
            return R.transpose() * E * R


def grid_rank_one_density():
    from sage.all import QuadraticForm, gcd, valuation
    rng = random.Random(20261027)
    out = []
    for p in (2, 3, 5, 7):
        for i in range(25):
            n = rng.randint(3, 8)
            T = random_even_form(rng, n, p)
            Q = QuadraticForm(ZZ, T)
            # Sage's local_density divides the form and the target by the
            # content of the form, which for Q = p^v Q' divides the density
            # by p^v; the definition (and the C++ code) does not, so the
            # factor is restored. (Checked by exhaustive counting modulo p^k
            # on forms of content 2 and 3.)
            content = gcd(Q.coefficients())
            fix = QQ(p)**valuation(content, p)
            for a in (1, 2, 3, 4, 5, 6, 8, 9, 12, 25, 27, 49):
                out.append((p, n, list(T.list()), a,
                            fix * Q.local_density(p, a)))
    # larger primes, where counting modulo p^k is out of reach and the C++
    # code uses its recursion only
    rng = random.Random(20261026)
    for p in (11, 13, 397, 1061):
        for i in range(10):
            n = rng.randint(3, 8)
            T = random_even_form(rng, n, p)
            Q = QuadraticForm(ZZ, T)
            content = gcd(Q.coefficients())
            fix = QQ(p)**valuation(content, p)
            for a in (1, 2, 3, p, 2 * p, 3 * p, p * p, 2 * p * p, p**3):
                out.append((p, n, list(T.list()), a,
                            fix * Q.local_density(p, a)))
    return out


def true_local_density(Q, p, a):
    """The local density of the definition: Sage's local_density divides the
    form and the target by the content of the form, which for an imprimitive
    Q = p^v Q' divides the density by p^v (see grid_rank_one_density)."""
    from sage.all import gcd, valuation
    return QQ(p)**valuation(gcd(Q.coefficients()), p) * Q.local_density(p, a)


def unit_chain_sage(S, T, p):
    from sage.all import QuadraticForm, valuation
    def diag(M):
        N = QuadraticForm(ZZ, M).local_normal_form(p)
        return [ZZ(N[i, i]) for i in range(N.dim())]
    def qf(coeffs):
        n = len(coeffs)
        return QuadraticForm(ZZ, n, [coeffs[i] if i == j else 0
                                     for i in range(n) for j in range(i, n)])
    D = diag(T)
    DS = diag(S)
    units = [s for s in DS if valuation(s, p) == 0]
    nonunits = [s for s in DS if valuation(s, p) > 0]
    if len(nonunits) > 1:
        return "fail"
    chain = units + nonunits
    alpha = QQ(1)
    for k, s in enumerate(chain):
        alpha *= true_local_density(qf(D), p, s)
        if k == len(chain) - 1:
            break
        uni = [c for c in D if valuation(c, p) == 0]
        rest = [c for c in D if valuation(c, p) > 0]
        if len(uni) < 2:
            # no unimodular part left in the complement
            D = rest
            continue
        prod_uni = prod(QQ(c) for c in uni)
        keep = uni[:-2]
        last = prod_uni / (QQ(s) * prod(QQ(c) for c in keep))
        D = keep + [ZZ(last.numerator() * last.denominator())] + rest
    return alpha


def grid_unit_chain():
    rng = random.Random(20261028)
    out = []
    types = [t for t in ade_types(4)]
    for p in (3, 5, 7):
        for i in range(12):
            n = rng.choice([6, 8, 10])
            T = random_even_form(rng, n, p)
            if T.det() % p != 0:
                continue
            for t in types:
                S = ade_gram(t)
                if n - S.nrows() < 3:
                    continue
                out.append((p, n, list(T.list()), S.nrows(), list(S.list()),
                            unit_chain_sage(S, T, p)))
    return out


def grid_certified_mass_ratio(ref_dir):
    genera = [("243", "0"), ("351", "0"), ("375", "0")]
    out = []
    for (d, g) in genera:
        tok = open(os.path.join(ref_dir, "density_calib", "reps",
                                "d%sg%s" % (d, g), "l.txt")).read().split()
        n = int(tok[1])
        T = [int(x) for x in tok[3:3 + n * n]]
        path = os.path.join(ref_dir, "reference_MR",
                            "det%s_genus%s_MR.txt" % (d, g))
        for line in open(path):
            if line.startswith("#") or not line.strip():
                continue
            name, _, ratio = line.split()[:3]
            if name == "0":
                continue
            typ = [(c[0], int(c[1:])) for c in name.split("+")]
            S = ade_gram(typ)
            if S.nrows() > 6:
                continue
            out.append((n, T, S.nrows(), list(S.list()), QQ(ratio)))
    return out


def grid_genus_mass():
    from sage.quadratic_forms.genera.genus import Genus
    rng = random.Random(20261026)
    forms = [ade_gram(t) for t in ade_types(14) if len(t) <= 2]
    for n in range(2, 15):
        for p in (2, 3, 5):
            forms.append(random_even_form(rng, n, p))
        for (p, q) in ((2, 3), (2, 5), (3, 5)):
            a = rng.randint(1, n - 1)
            forms.append(block_diagonal_matrix(
                [random_even_form(rng, a, p), random_even_form(rng, n - a, q)],
                subdivide=False))
    out = []
    for T in forms:
        G = Genus(T)
        sym2 = [s for s in G.local_symbols() if s.prime() == 2][0]
        canonical = [[int(x) for x in c] for c in sym2.canonical_symbol()]
        out.append([T.nrows(), T.list(), G.mass(), canonical])
    return out


def gap_string(x):
    if isinstance(x, (list, tuple)):
        return "[" + ", ".join(gap_string(y) for y in x) + "]"
    return str(x)


def main():
    task, out_file = sys.argv[1], sys.argv[2]
    if task == "certified_mass_ratio" and len(sys.argv) < 4:
        raise ValueError("certified_mass_ratio needs the reference directory")
    if task == "generic_part":
        rows = grid_generic_part()
    elif task == "siegel_series":
        rows = grid_siegel_series()
    elif task == "alpha_unimodular":
        rows = grid_alpha_unimodular()
    elif task == "rank_one_density":
        rows = grid_rank_one_density()
    elif task == "unit_chain":
        rows = grid_unit_chain()
    elif task == "certified_mass_ratio":
        rows = grid_certified_mass_ratio(sys.argv[3])
    elif task == "genus_mass":
        rows = grid_genus_mass()
    else:
        raise ValueError("unknown task " + task)
    with open(out_file, "w") as fh:
        fh.write("# Generated by sage/make_density_references.py %s\n" % task)
        fh.write("return [\n")
        fh.write(",\n".join(gap_string(r) for r in rows))
        fh.write("];\n")
    print("%d records written to %s" % (len(rows), out_file))


if __name__ == "__main__":
    main()
