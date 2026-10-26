# -*- coding: utf-8; mode: sage -*-
'''
Compute Siegel series by using recursive equations given by
Katsurada, `An explicit formula for Siegel series,
American jornal of Mathematics 121 (199),415 - 452.'
'''

from sage.all import (ZZ, fundamental_discriminant, valuation,
                      kronecker_symbol, cached_function)

from .rec_coeff import cbb2_0, cbb2_1, _pol_ring, cbb2_dict
from .jordan_block import (jordan_blocks_odd, jordan_blocks_2,
                           JordanBlock2, JordanBlocks)

X = _pol_ring().gens()[0]


def siegel_series_polynomial(B, p):
    '''
    B is a half integral matrix and p is a rational prime.
    It returns the polynomial part of the local Siegel series.
    '''
    p = ZZ(p)
    if p == 2:
        blcs = jordan_blocks_2(B)
        return _siegel_series_polynomial_2(blcs)
    else:
        blcs = jordan_blocks_odd(B, p)
        return _siegel_series_polynomial_odd(blcs, p)


def swap_ls(l, i, j):
    l[i], l[j] = l[j], l[i]


@cached_function
def _siegel_series_polynomial_2(q):
    non_diags = ("h", "y")
    max_expt = q.blocks[0][0]
    unit_diags_first_blc = [qf for expt, qf in q.blocks
                            if expt == max_expt and qf not in non_diags]
    # Use known formulas if dim(q) in (1, 2).
    if q.dim() == 1:
        return siegel_series_dim1(q, ZZ(2))
    elif q.dim() == 2:
        return siegel_series_dim2(q, ZZ(2))

    if len(unit_diags_first_blc) == 1:
        # Use the recursive equation given in Theorem 4.1.

        blcs_q2 = q.blocks[1:]
        q2 = JordanBlocks(blcs_q2, 2)
        pol = _siegel_series_polynomial_2(q2)
        return (cbb2_1(q, q2) * pol.subs({X: ZZ(2) * X}) +
                cbb2_0(q, q2) * pol)

    # Else use the recursive equation given in Theorem 4.2.

    blcs = q.blocks
    if blcs[0][1] in non_diags:
        m, units_or_hy = blcs[0]
        blcs_q2 = blcs[1:]
    else:
        (m, u1), (_, u2) = blcs[0], blcs[1]
        units_or_hy = [u1, u2]
        blcs_q2 = blcs[2:]

    q2 = JordanBlocks(blcs_q2, 2)
    b1 = JordanBlock2(units_or_hy, m)
    coeff_dict = cbb2_dict(b1, q2)
    c11 = coeff_dict['11']
    c10 = coeff_dict['10']
    c21 = coeff_dict['21']
    c20 = coeff_dict['20']
    times2 = {X: ZZ(2) * X}
    times4 = {X: ZZ(4) * X}
    coeffs = [c11 * c21.subs(times2),
              c11 * c20.subs(times2) + c10 * c21,
              c10 * c20]
    pol = _siegel_series_polynomial_2(q2)
    pols = [pol.subs(times4), pol.subs(times2), pol]
    return sum((a * b for a, b in zip(coeffs, pols)))


@cached_function
def _siegel_series_polynomial_odd(q, p):
    # Use known formulas when dim(q) = 1 or = 2.
    if q.dim() == 1:
        return siegel_series_dim1(q, p)
    elif q.dim() == 2:
        return siegel_series_dim2(q, p)
    # Else use the recursive equation given in Theorem 4.1
    q2 = JordanBlocks(q.blocks[1:], p)
    pol = _siegel_series_polynomial_odd(q2, p)
    return (cbb2_1(q, q2) * pol.subs({X: p * X}) +
            cbb2_0(q, q2) * pol)


def siegel_series_dim1(q, p):
    return sum(((p * X) ** a for a in range(q.content_order() + 1)))


def siegel_series_dim2(q, p):
    det_4 = q.Gram_det() * ZZ(4)
    c = q.content_order()
    fd = fundamental_discriminant(-det_4)
    f = (valuation(det_4, p) - valuation(fd, p)) / ZZ(2)
    return (_siegel_series_dim2(p, c, f + 1) -
            kronecker_symbol(fd, p) * p * X * _siegel_series_dim2(p, c, f))


def _siegel_series_dim2(p, a, b):
    if b == 0:
        return 0
    a = min(a, b - 1)
    r1 = (1 - (p ** 2 * X) ** (a + 1)) / (1 - p ** 2 * X)
    rn2 = (p ** 3 * X ** 2) ** b * p * X - p ** b * (p * X) ** (2 * b - a)
    rd2 = p * X - 1
    return (r1 - rn2 / rd2) / (1 - p ** 3 * X ** 2)
