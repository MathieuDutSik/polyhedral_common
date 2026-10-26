from sage.all import (cached_function, ZZ, QQ, Integer, matrix, identity_matrix,
                      valuation, least_quadratic_nonresidue, legendre_symbol,
                      mul, hilbert_symbol, cached_method)
from degree2.utils import list_group_by
from itertools import groupby

_hy_det = {'h': -ZZ(1) / ZZ(4), 'y': ZZ(3) / ZZ(4)}
_hy_rational_diags = {'h': [ZZ(1), ZZ(-1)], 'y': [ZZ(1), ZZ(3)]}

two = ZZ(2)

# h = matrix([[0, 1/2], [1/2, 0]])
# y = matrix([[1, 1/2], [1/2, 1]])
# u1, u2 :units
# q1 <-> 2^m * diag(u1, u2), h or y.

mat_dict = {"h": matrix([[QQ(0), QQ(1) / QQ(2)],
                         [QQ(1) / QQ(2), QQ(0)]]),
            "y": matrix([[QQ(1), QQ(1) / QQ(2)],
                         [QQ(1) / QQ(2), QQ(1)]])}


class JordanBlock2(object):

    def __init__(self, units_or_hy, m):
        '''
        units_or_hy is a list of two units or "h" or "y".
        m is an integer.
        '''
        self._m = m
        if units_or_hy in ("h", "y"):
            self._type = units_or_hy
            self._mat_prim = mat_dict[units_or_hy]
        else:
            self._type = "u"
            self._mat_prim = matrix([[units_or_hy[0], 0],
                                     [0, units_or_hy[1]]])

    @property
    def gram_mat(self):
        return self._mat_prim * two ** self.m

    @property
    def type(self):
        return self._type

    @property
    def m(self):
        return self._m


class JordanBlocks(object):

    def __init__(self, blocks, p):
        self._blocks = blocks
        self._p = ZZ(p)

    @property
    def blocks(self):
        return self._blocks

    @property
    def p(self):
        return self._p

    @cached_method
    def dim(self):
        p = self.p
        if p != 2:
            return len(self.blocks)
        else:
            return sum(2 if isinstance(b, str) else 1 for _, b in self.blocks)

    @cached_method
    def hasse_invariant__OMeara(self):
        p = self.p
        if p != 2:
            rational_diags = [p ** a * b for a, b in self.blocks]
        else:
            rational_diags = []
            for a, b in self.blocks:
                if isinstance(b, str):
                    rational_diags.extend(
                        [p ** a * x for x in _hy_rational_diags[b]])
                else:
                    rational_diags.append(p ** a * b)
        n = len(rational_diags)
        res = ZZ(1)
        for i in range(n):
            for j in range(i + 1, n):
                res *= hilbert_symbol(rational_diags[i], rational_diags[j], p)
        return res * hilbert_symbol(self.Gram_det(), -1, p)

    def content_order(self):
        return min([a for a, _ in self.blocks])

    @cached_method
    def Gram_det(self):
        p = self.p
        if p != 2:
            return mul(p ** a * b for a, b in self.blocks)
        else:
            res = ZZ(1)
            for a, b in self.blocks:
                # b is 'h' or 'y'
                if isinstance(b, str):
                    res *= p ** (2 * a) * _hy_det[b]
                else:
                    res *= p ** a * b
            return res

    def __radd__(self, other):
        if isinstance(other, JordanBlock2):
            if other.type in ('h', 'y'):
                blcs = [(other.m, other.type)] + self.blocks
            else:
                blcs = [(other.m, other._mat_prim[i, i])
                        for i in range(2)] + self.blocks
            return JordanBlocks(blcs, two)
        elif isinstance(other, JordanBlocks):
            return JordanBlocks(other.blocks + self.blocks, self.p)
        else:
            raise NotImplementedError

    def __add__(self, other):
        if isinstance(other, JordanBlocks):
            return JordanBlocks(other.blocks + self.blocks, self.p)
        else:
            raise NotImplementedError

    def __hash__(self):
        return hash(('JordanBlocks', tuple(self.blocks), self.p))

    def __eq__(self, other):
        if isinstance(other, JordanBlocks):
            return self.blocks == other.blocks and self.p == other.p
        else:
            return False


def bracket_action(a, b):
    return b.transpose() * a * b


def perm_mat(i, j, n):
    l = identity_matrix(QQ, n).columns()
    l[i], l[j] = l[j], l[i]
    return matrix(l)


def perm_mat2(i, j, n):
    # (i, j) => (0, 1)
    l = identity_matrix(QQ, n).columns()
    rem = [a for a in range(n) if a not in (i, j)]
    l1 = [l[i], l[j]] + [l[a] for a in rem]
    return matrix(QQ, l1).transpose()


def find_min_ord_elet(S, p):
    '''
    S: ((1+delta_ij)/2 s_ij) half integral matrix
    Retruns (i, j) such that s_ij != 0
    and ord((1+delta_ij)/2 s_ij) is min.
    If p = 2, then index (i, j) (i != j) is preferred.
    If p is odd, index (i, j) (i == j) is preferred.
    '''
    p = Integer(p)
    n = len(S.columns())
    elts_with_idx = [(S[(i, j)], (i, j)) if i == j else
                     (S[(i, j)] * 2, (i, j))
                     for i in range(n) for j in range(n)
                     if i <= j and S[(i, j)] != 0]
    val_with_idx = [(valuation(e, p), i) for e, i in elts_with_idx]
    min_val = min([v for v, _ in val_with_idx])
    min_val_idcs = [i for v, i in val_with_idx if v == min_val]
    if p != 2:
        for i, j in min_val_idcs:
            if i == j:
                return (i, j)
        return min_val_idcs[0]
    else:
        for i, j in min_val_idcs:
            if i != j:
                return (i, j)
        return min_val_idcs[0]


def _jordan_decomposition_odd_p(S, p):
    '''
    Input:
      S: half integral matrix
      p: prime
    Output:
      a list l = [p^n1 * u1, p^n2 * u2, ... p^nk * uk]
      such that n1 <= n2 <= ... nk and
      diag(l) is Z_p equivalent to S.
    '''
    p = Integer(p)
    acc = []
    n = len(S.columns())
    while True:
        if n == 1:
            return acc + [S[0, 0]]
        i0, j0 = find_min_ord_elet(S, p)
        if i0 != j0:
            u = identity_matrix(QQ, n)
            u[(j0, i0)] = 1
            S = bracket_action(S, u)
            j0 = i0

        S = _jordan_dcomp_diag(i0, n, S)
        acc.append(S[(0, 0)])
        S = S.submatrix(row=1, col=1)
        n -= 1


def _jordan_dcomp_diag(i0, n, S):
    u = perm_mat(0, i0, n)
    S = bracket_action(S, u)

    # (0, 0) elment is the min order.
    u = identity_matrix(QQ, n)
    a = S[(0, 0)]
    for j in range(1, n):
        u[(0, j)] = -S[(0, j)] / a
    return bracket_action(S, u)


def _jordan_decomposition_2(S):
    '''
    Input:
      S: half integral matrix
    Output:
      list of tuples (a, b)
      Here a is an integer which is an exponent.
      b is equal to an element of [1, 3, 5, 7] or 'h' or 'y'.
    '''
    n = len(S.columns())
    acc = []
    while True:
        i0, j0 = find_min_ord_elet(S, 2)
        if (n == 1) or (n == 2 and i0 != j0):
            mat_lst = acc + [S]
            break

        if i0 == j0:
            S = _jordan_dcomp_diag(i0, n, S)
            acc.append(S.submatrix(nrows=1, ncols=1))
            S = S.submatrix(row=1, col=1)
            n -= 1
        else:
            # (0, 1) element has the minimal order.
            S = bracket_action(S, perm_mat2(i0, j0, n))
            u = identity_matrix(QQ, n)
            for j in range(2, n):
                s00, s01, s11 = S[(0, 0)], S[(0, 1)], S[(1, 1)]
                d = s00 * s11 - s01 ** 2
                s0j = S[(0, j)]
                s1j = S[(1, j)]
                a0 = (-s11 * s0j + s01 * s1j) / d
                a1 = (s01 * s0j - s00 * s1j) / d
                u[(0, j)] = a0
                u[(1, j)] = a1
            S = bracket_action(S, u)
            acc.append(S.submatrix(nrows=2, ncols=2))
            S = S.submatrix(row=2, col=2)
            n -= 2
    res = []
    for m in mat_lst:
        if m.ncols() == 1:
            a = m[0, 0]
            e = valuation(a, two)
            u = (a / two ** e) % 8
            res.append((e, u))
        else:
            e = valuation(m[0, 1], two) + 1
            m = m / two ** (e - 1)
            if m.det() % 8 == 3:
                h_or_y = 'y'
            else:
                h_or_y = 'h'
            res.append((e, h_or_y))
    return res


@cached_function
def does_2adically_rep_zero(a, b, c):
    '''
    Should a, b, c in [1, 3, 5, 7].
    Returns True if ax^2 + by^2 +cz^2 is 2-adically represents zero
    otherwise False.
    '''
    for x in [0, 1, 4]:
        for y in [0, 1, 4]:
            for z in [0, 1, 4]:
                if (any((w % 2 for w in [x, y, z])) and
                        (a * x + b * y + c * z) % 8 == 0):
                    return True
    return False


def _trans_jordan_dec_2(us):
    '''
    us = (u1, ... , un)
    Returns a jordan decomp of diag(us) so that the number of diag elements
    is less that or equal to 2.
    For example,
    x1^2 + x2^2 - x3^2 ~ x1^2 + 2*x1*x2.
    x1^2 + x2^2 + x3^2 ~ 3*x1^2 + 2(x2^2 + x2*x3 + x3^2).
    '''
    hyp2_name = "h"
    y2_name = "y"
    if len(us) <= 2:
        return us
    u1, u2, u3 = us[:3]
    u = u1 * u2 * u3
    l = _trans_jordan_dec_2(us[3:])
    if does_2adically_rep_zero(u1, u2, u3):
        l.extend([(-u) % 8, hyp2_name])
    else:
        l.extend([(3 * u) % 8, y2_name])
    l_diag = [a for a in l if a not in [hyp2_name, y2_name]]
    if len(l_diag) <= 2:
        return l
    else:
        l_non_diag = [a for a in l if a in [hyp2_name, y2_name]]
        l1 = _trans_jordan_dec_2(l_diag)
        return l1 + l_non_diag


def jordan_blocks_odd(S, p):
    '''
    p: odd prime,
    S: half integral matrix.
    Let S ~ sum p**a_i * b_i * x_i^2 be a jordan decomposition.
    Returns the instance of JordanBlocks attached to
    a list of (a_i, b_i) sorted so that a_0 >= a_1 >= ...
    '''
    p = Integer(p)
    res = []
    u = least_quadratic_nonresidue(p)
    for e, ls in groupby(_jordan_decomposition_odd_p(S, p),
                         key=lambda x: valuation(x, p)):
        ls = [x // p ** e for x in ls]
        part_res = [(e, ZZ(1)) for _ in ls]
        if legendre_symbol(mul(ls), p) == -1:
            part_res[-1] = (e, u)
        res.extend(part_res)
    return JordanBlocks(list(reversed(res)), p)


def jordan_blocks_2(S):
    '''
    S: half integral matrix
    Returns the instance of JordanBlocks attached to
    a list of tuples (a, b)
    `a' is a integer which represents an exponent.
    `b' corresponds a quadratic form with size <= 2,
    that is, `b' is in [1, 3, 5, 7] or is equal to
    strings 'h' or 'y'.
    If `b' is equal to 'h', then the gram matrix of
    the corresponding quadratic form is equal to
    matrix([[0, 1/2], [1/2, 0]]).
    If `b' is equal to 'y', then the gram matrix of
    the corresponding quadratic form is equal to
    matrix([[1, 1/2], [1/2, 1]]).
    If `b' is in [1, 3, 5, 7], then the corresponding
    gram matrix is equal to matrix([[b]]).
    These tuples are sorted so that assumptions of Theorem 4.1 or Theorem 4.2
    hold.
    '''
    ls = _jordan_decomposition_2(S)
    # ls is a list of jordan blocks but it may contain diagonal entries of
    # units of length >= 3.
    # We take another jordan blocks by _trans_jordan_dec_2.
    res = [(a, b) for a, b in ls if isinstance(b, str)]
    diag_elts = [(a, b) for a, b in ls if not isinstance(b, str)]
    for a, us_w_idx in list_group_by(diag_elts, lambda x: x[0]):
        l = _trans_jordan_dec_2([_b for _, _b in us_w_idx])
        for b in l:
            if isinstance(b, str):
                res.append((a + 1, b))
            else:
                res.append((a, b))
    # Sort the result so that assumptions of Theorem 4.1 and 4.2 hold.
    non_diags = ("h", "y")
    res1 = []
    for _, blcs in sorted(list_group_by(res, lambda x: x[0]),
                          key=lambda x: -x[0]):
        unit_diags = [(a, qf) for a, qf in blcs if qf not in non_diags]
        non_unit_diags = [(a, qf) for a, qf in blcs if qf in non_diags]
        blcs = unit_diags + non_unit_diags
        res1.extend(blcs)
    return JordanBlocks(res1, two)
