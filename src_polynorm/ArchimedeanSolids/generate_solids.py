"""Exact vertex coordinates of the Platonic and Archimedean solids, written
for polyhedral_common: rationals, Qsqrt2 / Qsqrt5 (x = sqrt d), or a real
algebraic field with generator x (RealAlgebraic=FileDesc)."""
import itertools, sys
import sympy as sp
from sympy import Rational as R, sqrt, Poly, QQ, symbols, invert, rem
import mpmath

g = symbols('x')

def cyc(v):
    a, b, c = v
    return [(a, b, c), (b, c, a), (c, a, b)]

def allperm(v):
    # The distinct permutations, in a deterministic order (a set of sympy
    # expressions would iterate in hash order, which changes between runs).
    out = []
    for p in itertools.permutations(v):
        if p not in out:
            out.append(p)
    return out

def signs(v):
    out = []
    for s in itertools.product([1, -1], repeat=3):
        out.append(tuple(si * vi for si, vi in zip(s, v)))
    return out

def dedup(L):
    seen = []
    for p in L:
        q = tuple(sp.expand(c) for c in p)
        if all(any(sp.simplify(a - b) != 0 for a, b in zip(q, r)) for r in seen):
            seen.append(q)
    return seen

phi = (1 + sqrt(5)) / 2

def solids():
    S = {}
    S['tetrahedron'] = ('Q', [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)], 3)
    S['cube'] = ('Q', list(itertools.product([1, -1], repeat=3)), 3)
    S['octahedron'] = ('Q', [p for v in allperm((1, 0, 0)) for p in signs(v)], 4)
    S['dodecahedron'] = ('Qsqrt5', list(itertools.product([1, -1], repeat=3)) +
                         [p for v in cyc((0, 1 / phi, phi)) for p in signs(v)], 3)
    S['icosahedron'] = ('Qsqrt5', [p for v in cyc((0, 1, phi)) for p in signs(v)], 5)
    S['truncated_tetrahedron'] = ('Q', [p for v in allperm((3, 1, 1)) for p in signs(v)
                                        if sum(1 for c in p if c < 0) % 2 == 0], 3)
    S['cuboctahedron'] = ('Q', [p for v in allperm((1, 1, 0)) for p in signs(v)], 4)
    S['truncated_cube'] = ('Qsqrt2', [p for v in allperm((sqrt(2) - 1, 1, 1)) for p in signs(v)], 3)
    S['truncated_octahedron'] = ('Q', [p for v in allperm((0, 1, 2)) for p in signs(v)], 3)
    S['rhombicuboctahedron'] = ('Qsqrt2', [p for v in allperm((1, 1, 1 + sqrt(2))) for p in signs(v)], 4)
    S['truncated_cuboctahedron'] = ('Qsqrt2', [p for v in allperm((1, 1 + sqrt(2), 1 + 2 * sqrt(2))) for p in signs(v)], 3)
    S['icosidodecahedron'] = ('Qsqrt5', [p for v in cyc((0, 0, phi)) for p in signs(v)] +
                              [p for v in cyc((R(1, 2), phi / 2, phi**2 / 2)) for p in signs(v)], 4)
    S['truncated_dodecahedron'] = ('Qsqrt5', [p for base in [(0, 1 / phi, 2 + phi), (1 / phi, phi, 2 * phi), (phi, 2, phi + 1)]
                                              for v in cyc(base) for p in signs(v)], 3)
    S['truncated_icosahedron'] = ('Qsqrt5', [p for base in [(0, 1, 3 * phi), (1, 2 + phi, 2 * phi), (phi, 2, phi**3)]
                                             for v in cyc(base) for p in signs(v)], 3)
    S['rhombicosidodecahedron'] = ('Qsqrt5', [p for base in [(1, 1, phi**3), (phi**2, phi, 2 * phi), (2 + phi, 0, phi**2)]
                                              for v in cyc(base) for p in signs(v)], 4)
    S['truncated_icosidodecahedron'] = ('Qsqrt5', [p for base in [(1 / phi, 1 / phi, 3 + phi), (2 / phi, phi, 1 + 2 * phi),
                                                                  (1 / phi, phi**2, -1 + 3 * phi), (2 * phi - 1, 2, 2 + phi), (phi, 3, 2 * phi)]
                                                   for v in cyc(base) for p in signs(v)], 3)
    # Snub cube: t the tribonacci constant, generator x = t.
    t = g
    sc = []
    for perm in itertools.permutations(range(3)):
        even = sp.combinatorics.Permutation(list(perm)).is_even
        base = (1, 1 / t, t)
        v = tuple(base[perm[i]] for i in range(3))
        for s in itertools.product([1, -1], repeat=3):
            nplus = sum(1 for si in s if si > 0)
            if (even and nplus % 2 == 0) or ((not even) and nplus % 2 == 1):
                sc.append(tuple(si * vi for si, vi in zip(s, v)))
    S['snub_cube'] = (('Alg', [-1, -1, -1, 1]), sc, 5)
    # Snub dodecahedron: xi^3 - 2 xi = phi, generator x = xi, phi = x^3 - 2x.
    xi = g
    ph = xi**3 - 2 * xi
    al = xi - 1 / xi
    be = xi * ph + ph**2 + ph / xi
    bases = [(2 * al, 2, 2 * be),
             (al + be / ph + ph, -al * ph + be + 1 / ph, al / ph + be * ph - 1),
             (-al / ph + be * ph + 1, -al + be / ph - ph, al * ph + be - 1 / ph),
             (-al / ph + be * ph - 1, al - be / ph - ph, al * ph + be + 1 / ph),
             (al + be / ph - ph, al * ph - be + 1 / ph, al / ph + be * ph + 1)]
    sd = []
    for base in bases:
        for v in cyc(base):
            for s in itertools.product([1, -1], repeat=3):
                if sum(1 for si in s if si > 0) % 2 == 0:
                    sd.append(tuple(si * vi for si, vi in zip(s, v)))
    mp_sd = Poly(sp.expand((xi**3 - 2 * xi)**2 - (xi**3 - 2 * xi) - 1), xi)
    S['snub_dodecahedron'] = (('Alg', [int(c) for c in reversed(mp_sd.all_coeffs())]), sd, 5)
    return S

def root_value(coeffs):
    # The real root used: tribonacci > 1 for the snub cube; for the snub
    # dodecahedron the root of xi^3 - 2 xi = phi near 1.7155615.
    mpmath.mp.dps = 3000
    p = coeffs[::-1]
    roots = mpmath.polyroots(p, maxsteps=500, extraprec=3000)
    real = [r.real for r in roots if abs(r.imag) < mpmath.mpf(10)**-1000]
    return real

def to_poly(expr, coeffs):
    """expr a rational function in x: reduce modulo the minimal polynomial."""
    mp = Poly(list(reversed(coeffs)), g, domain=QQ)
    expr = sp.together(sp.sympify(expr))
    num, den = sp.fraction(expr)
    num = Poly(num, g, domain=QQ)
    den = Poly(den, g, domain=QQ)
    inv = invert(den, mp)
    res = (num * inv).rem(mp)
    return res

def fmt_poly(P):
    terms = []
    for (k,), c in sorted(P.terms(), key=lambda t: t[0][0]):
        if c == 0:
            continue
        cs = str(c)
        if k == 0:
            terms.append(cs)
        else:
            mon = 'x' if k == 1 else f'x^{k}'
            if c == 1:
                terms.append(mon)
            elif c == -1:
                terms.append('-' + mon)
            else:
                terms.append(cs + '*' + mon)
    if not terms:
        return '0'
    s = terms[0]
    for t in terms[1:]:
        s += t if t.startswith('-') else '+' + t
    return s

def fmt_quad(c, d):
    c = sp.nsimplify(sp.expand(c))
    a = sp.Rational(sp.expand(c).subs(sqrt(d), 0)) if c.has(sqrt(d)) else sp.Rational(c)
    b = sp.Rational(sp.expand((c - a) / sqrt(d)))
    if b == 0:
        return str(a)
    bs = 'x' if b == 1 else ('-x' if b == -1 else f'{b}*x')
    if a == 0:
        return bs
    return f'{a}' + (bs if bs.startswith('-') else '+' + bs)

def numeric(c, coeffs=None, root=None):
    if coeffs is None:
        return float(sp.N(c, 30))
    return float(sp.N(sp.sympify(c).subs(g, root), 30))

def write_desc(coeffs, val, fname):
    mpmath.mp.dps = 3000
    with open(fname, 'w') as f:
        f.write(f'{len(coeffs) - 1}\n')
        f.write(' '.join(str(c) for c in coeffs) + '\n')
        f.write(repr(float(val)) + '\n')
        n_expo = 100
        f.write(f'{n_expo}\n')
        # continued fraction convergents of val
        terms = []
        y = val
        for _ in range(5 * n_expo + 2):
            a = mpmath.floor(y)
            terms.append(int(a))
            y = 1 / (y - a)
        def conv(n):
            h1, h0, k1, k0 = 1, 0, 0, 1
            for a in terms[:n]:
                h1, h0 = a * h1 + h0, h1
                k1, k0 = a * k1 + k0, k1
            return sp.Rational(h1, k1)
        for i in range(1, n_expo + 1):
            a1, a2 = conv(5 * i), conv(5 * i + 1)
            lo, hi = (a1, a2) if a1 < a2 else (a2, a1)
            assert lo <= val <= hi or abs(lo - val) < mpmath.mpf(10)**-2000
            f.write(f'{lo} {hi}\n')

def main(outdir):
    S = solids()
    for name, (field, verts, degree) in S.items():
        if isinstance(field, tuple):
            coeffs = field[1]
            reals = root_value(coeffs)
            if name == 'snub_cube':
                root = [r for r in reals if r > 1][0]
            else:
                root = [r for r in reals if abs(r - mpmath.mpf('1.7155615')) < 1e-3][0]
            num = [tuple(numeric(c, coeffs, root) for c in v) for v in verts]
        else:
            num = [tuple(numeric(c) for c in v) for v in verts]
        # deduplicate numerically
        uniq, idx = [], []
        for i, p in enumerate(num):
            if all(max(abs(a - b) for a, b in zip(p, q)) > 1e-9 for q in uniq):
                uniq.append(p)
                idx.append(i)
        verts = [verts[i] for i in idx]
        # checks: equal edges, vertex degree
        import math
        d2 = [[sum((a - b)**2 for a, b in zip(p, q)) for q in uniq] for p in uniq]
        dmin = min(d2[i][j] for i in range(len(uniq)) for j in range(len(uniq)) if i != j)
        degs = set(sum(1 for j in range(len(uniq)) if j != i and abs(d2[i][j] - dmin) < 1e-7 * dmin)
                   for i in range(len(uniq)))
        print(f'{name}: {len(uniq)} vertices, field {field if not isinstance(field, tuple) else "alg"}, degrees {degs} (expected {degree})')
        if degs != {degree}:
            print('  CHECK FAILED')
        fname = f'{outdir}/{name}.ext'
        with open(fname, 'w') as f:
            f.write(f'{len(verts)} 4\n')
            for v in verts:
                if field == 'Q':
                    s = ' '.join(str(sp.Rational(c)) for c in v)
                elif field == 'Qsqrt2':
                    s = ' '.join(fmt_quad(c, 2) for c in v)
                elif field == 'Qsqrt5':
                    s = ' '.join(fmt_quad(c, 5) for c in v)
                else:
                    s = ' '.join(fmt_poly(to_poly(c, field[1])) for c in v)
                f.write('1 ' + s + '\n')
        if isinstance(field, tuple):
            write_desc(field[1], root, f'{outdir}/{name}.field')

main(sys.argv[1])
