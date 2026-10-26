# -*- coding: utf-8; mode: sage -*-
'''
Local invariants from Katsurada's paper "an explicit formula for
Siegel series" section 3.
'''

from sage.all import (valuation, ZZ, kronecker_symbol, hilbert_symbol)


def _chi_p_odd(a, p):
    p = ZZ(p)
    r = valuation(a, p)
    if r % 2 == 0:
        c = a // (p ** r)
        return kronecker_symbol(c, p)
    else:
        return 0


def _chi_2(a):
    two = ZZ(2)
    r = valuation(a, two)
    if r % 2 == 1:
        return 0
    else:
        c = a // (two ** r)
        rem = c % 8
        if rem == 1:
            return 1
        elif rem == 5:
            return -1
        else:
            return 0


def chi_p(a, p):
    if p == 2:
        return _chi_2(a)
    else:
        return _chi_p_odd(a, p)


def eta_p(q):
    '''
    q: an instance of jordan_block.JordanBlocks.
    '''
    p = q.p
    hasse_inv = q.hasse_invariant__OMeara()
    detb = q.Gram_det()
    n = q.dim()
    h1 = hilbert_symbol(detb, detb * (-1) ** ((n - 1) // 2), p)
    h2 = hilbert_symbol(-1, -1, p)
    return hasse_inv * h1 * h2 ** ((n ** 2 - 1) // 8)


def xi_p(q):
    '''
    q: an instance of jordan_block.JordanBlocks
    '''
    n = q.dim()
    p = q.p
    return ZZ(chi_p((-1) ** (n // 2) * q.Gram_det(), p))


def xi_prime_p(q):
    xi = xi_p(q)
    return 1 + xi - xi ** 2


def xi_to_xi_dash(xi):
    return 1 + xi - xi ** 2


def small_d(q):
    '''
    d(B) given in Thoerem 3.2.
    '''
    two = ZZ(2)
    n = q.dim()
    p = q.p
    larged = two ** (2 * (n // 2)) * q.Gram_det()
    return valuation(larged, p)


def delta_p(q):
    r'''
    \delta(B) given in Theorem 3.2.
    '''
    n = q.dim()
    p = q.p
    db = small_d(q)
    if n % 2 == 0:
        dl = 1 if p == 2 else 0
        return 2 * ((db + 1 - dl) // 2)
    else:
        return db


def e_p(q):
    '''
    e(B) given in Theorem 3.2.
    '''
    n = q.dim()
    dlt = delta_p(q)
    if n % 2:
        return dlt
    else:
        return dlt - ZZ(2) + ZZ(2) * xi_p(q) ** 2


def zeta_p(q, p):
    r'''
    \zeta(B) given in Theorem 3.2.
    '''
    if q.dim() % 2:
        return eta_p(q)
    else:
        return ZZ(1)
