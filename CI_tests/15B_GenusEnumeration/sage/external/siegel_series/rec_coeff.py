# -*- coding: utf-8; mode: sage -*-
from .local_invariants import (xi_p, eta_p,
                               delta_p, small_d, xi_to_xi_dash)
from sage.all import PolynomialRing, QQ, cached_function, ZZ, hilbert_symbol
from .jordan_block import JordanBlocks


@cached_function
def _pol_ring():
    R = PolynomialRing(QQ, names="x")
    return R


def cbb2_1(q, q2):
    '''
    A polynomial defined before Theorem 4.1 in the Katsurada's paper.
    q, q2: instances of jordan_block.JordanBlocks.
    q: degree n
    q2: degree n - 1
    '''
    p = q.p
    X = _pol_ring().gens()[0]
    n = q.dim()
    if n % 2 == 0:
        xi = xi_p(q)
        return (1 - p ** (n // 2) * xi * X) / (1 - p ** (n + 1) * X ** 2)
    else:
        if n > 1:
            xi_tilde = xi_p(q2)
        else:
            xi_tilde = 1
        return ZZ(1) / (1 - p ** ((n + 1) // 2) * xi_tilde * X)


def _invariants_1_even(q, q2):
    xi = xi_p(q)
    xi_dash = xi_to_xi_dash(xi)
    eta_tilde = eta_p(q2)
    delta = delta_p(q)
    delta_tilde = delta_p(q2)
    return {"xi_dash": xi_dash,
            "eta_tilde": eta_tilde,
            "xi": xi,
            "delta": delta,
            "delta_tilde": delta_tilde}


def _invariants_1_odd(q, q2):
    n = q.dim()
    eta = eta_p(q)
    delta = delta_p(q)
    if n > 1:
        xi_tilde = xi_p(q2)
        xi_tilde_dash = xi_to_xi_dash(xi_tilde)
        delta_tilde = delta_p(q2)
    else:
        xi_tilde = xi_tilde_dash = 1
        delta_tilde = 0
    return {"eta": eta,
            "delta": delta,
            "xi_tilde": xi_tilde,
            "xi_tilde_dash": xi_tilde_dash,
            "delta_tilde": delta_tilde}


def _rat_func_1_even(p, n,
                     xi_dash=None, eta_tilde=None, xi=None,
                     delta=None, delta_tilde=None):
    X = _pol_ring().gens()[0]
    num = (ZZ(-1) ** (xi + 1) * xi_dash * eta_tilde *
           (1 - p ** (n // 2 + 1) * xi * X) *
           (p ** (n // 2) * X) ** (delta - delta_tilde + xi ** 2) *
           p ** ZZ(QQ(delta) / 2))
    denom = 1 - p ** (n + 1) * X ** 2
    return num / denom


def _rat_func_1_odd(p, n,
                    delta_tilde=None,
                    eta=None,
                    xi_tilde_dash=None,
                    xi_tilde=None,
                    delta=None):
    X = _pol_ring().gens()[0]
    num = (ZZ(-1) ** xi_tilde * xi_tilde_dash * eta *
           (p ** ((n - 1) // 2) * X) ** (delta - delta_tilde + 2 - xi_tilde ** 2) *
           p ** ZZ(QQ(2 * delta - delta_tilde + 2) / 2))
    denom = 1 - p ** ZZ(QQ(n + 1) / 2) * X * xi_tilde
    return num / denom


def cbb2_0(q, q2):
    '''
    A polynomial defined before Theorem 4.1 in the Katsurada's paper.
    q: degree n
    q2: degree n - 1
    '''
    p = q.p
    n = q.dim()
    if n % 2 == 0:
        d = _invariants_1_even(q, q2)
        return _rat_func_1_even(p, n, **d)
    else:
        d = _invariants_1_odd(q, q2)
        return _rat_func_1_odd(p, n, **d)


two = ZZ(2)


def _invariants_2_common(b1, q2):
    '''
    b1 is an instance of JordanBlock2.
    q2: an instance of jordan_block.JordanBlocks of dim n - 2
    '''
    n = q2.dim() + 2
    m = b1.m
    q = b1 + q2
    q3 = JordanBlocks([(m, 1)], two) + q2
    delta = delta_p(q)
    delta_tilde = delta_p(q3)
    delta_hat = delta_p(q2)
    # Definition of sigma
    if ((n % 2 == 0 and b1.type == 'u' and small_d(q) % 2 == 1)
            or
            (n % 2 == 0 and b1.type != 'u' and xi_p(q2) == 0)):
        sigma = (2 * delta_tilde - delta - delta_hat + 2) / two
    elif n % 2 == 1 and b1.type != 'u' and small_d(q3) % 2 == 0:
        sigma = ZZ(2)
    else:
        sigma = ZZ(0)
    return {'sigma': sigma,
            'delta': delta,
            'delta_tilde': delta_tilde,
            'delta_hat': delta_hat}


def _invariants_2_even(b1, q2):
    n = q2.dim() + 2
    m = b1.m
    q = b1 + q2
    xi = xi_p(q)
    xi_dash = xi_to_xi_dash(xi)
    xi_hat = xi_p(q2)
    xi_hat_dash = xi_to_xi_dash(xi_hat)
    # Definition of eta_tilde
    if b1.type == 'u' and small_d(q2) % 2 == 0:
        _q = JordanBlocks([(m, b1._mat_prim[(1, 1)])], two)
        eta_tilde = eta_p(_q + q2)
    elif b1.type != 'u' and xi_hat != 0:
        eta_tilde = (ZZ(-1) ** (((n - 1) ** 2 - 1) // 8) * q2.hasse_invariant__OMeara()
                     * hilbert_symbol(two ** m,
                                      ZZ(-1) ** (n // 2 - 1) * q2.Gram_det(),
                                      2))
    else:
        eta_tilde = ZZ(1)

    res = {"xi": xi,
           "xi_dash": xi_dash,
           "xi_hat": xi_hat,
           "xi_hat_dash": xi_hat_dash,
           "eta_tilde": eta_tilde}
    res.update(_invariants_2_common(b1, q2))
    return res


def _invariants_2_odd(b1, q2):
    m = b1.m
    q = b1 + q2
    eta = eta_p(q)
    eta_hat = eta_p(q2)
    q3 = JordanBlocks([(m, 1)], two) + q2
    if b1.type != 'u' and small_d(q3) % 2 == 0:
        xi_tilde = 1
    else:
        xi_tilde = 0
    res = {'eta': eta,
           'eta_hat': eta_hat,
           'xi_tilde': xi_tilde}
    res.update(_invariants_2_common(b1, q2))
    return res


def _rat_funcs_even(n, xi=None, xi_dash=None, xi_hat=None,
                    xi_hat_dash=None, eta_tilde=None, sigma=None,
                    delta=None, delta_tilde=None, delta_hat=None):
    res = {}
    X = _pol_ring().gens()[0]
    # 11
    num = 1 - two ** (n // 2) * xi * X
    denom = 1 - two ** (n + 1) * X ** 2
    res['11'] = num / denom
    # 10
    expt = delta - delta_tilde + xi ** 2 + sigma
    num = (ZZ(-1) ** (xi + 1) * xi_dash * eta_tilde *
           (1 - two ** (n // 2 + 1) * X * xi) *
           X ** expt * 2 ** (delta / two + (n // 2) * expt))
    denom = 1 - two ** (n + 1) * X ** 2
    res['10'] = num / denom
    # 21
    num = ZZ(1)
    denom = 1 - 2 ** (n // 2) * xi_hat * X
    res['21'] = num / denom
    # 20
    expt = delta_tilde - delta_hat + 2 - xi_hat ** 2 - sigma
    num = (ZZ(-1) ** xi_hat * xi_hat_dash * eta_tilde *
           X ** expt *
           2 ** (((n - 2) // 2) * expt +
                 (2 * delta_tilde - delta_hat + 2 - 2 * sigma) / two))
    denom = 1 - 2 ** (n // 2) * X * xi_hat
    res['20'] = num / denom
    return res


def _rat_funcs_odd(n, sigma=None, delta=None, delta_tilde=None,
                   delta_hat=None, eta=None, eta_hat=None, xi_tilde=None):
    res = {}
    X = _pol_ring().gens()[0]
    # 11
    num = ZZ(1)
    denom = 1 - 2 ** ((n + 1) // 2) * xi_tilde * X
    res['11'] = num / denom
    # 10
    expt = delta - delta_tilde + 2 - xi_tilde ** 2 + sigma
    num = (ZZ(-1) ** xi_tilde * eta *
           X ** expt *
           2 ** ((n - 1) // 2 * expt + (2 * delta - delta_tilde + 2 + sigma) / two))
    denom = 1 - 2 ** ((n + 1) // 2) * X * xi_tilde
    res['10'] = num / denom
    # 21
    num = 1 - 2 ** ((n - 1) // 2) * xi_tilde * X
    denom = 1 - 2 ** n * X ** 2
    res['21'] = num / denom
    # 20
    expt = delta_tilde - delta_hat + xi_tilde ** 2 - sigma
    num = (ZZ(-1) ** (xi_tilde + 1) * eta_hat *
           (1 - 2 ** ((n + 1) // 2) * X * xi_tilde) *
           X ** expt * 2 ** ((n - 1) // 2 * expt + (delta_tilde - sigma) / two))
    denom = 1 - 2 ** n * X ** 2
    res['20'] = num / denom
    return res


def cbb2_dict(b1, q2):
    n = q2.dim() + 2
    if n % 2 == 0:
        d = _invariants_2_even(b1, q2)
        return _rat_funcs_even(n, **d)

    if n % 2 == 1:
        d = _invariants_2_odd(b1, q2)
        return _rat_funcs_odd(n, **d)
