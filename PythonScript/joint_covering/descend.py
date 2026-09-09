# Joint local descent on (Q, C) for Z^5 + {0, C_1, ..., C_{m-1}}: outer loop
# of re-tessellation + SLSQP on (cholesky(Q), C) under the circumradius
# constraints of all cell classes, each round verified against a fresh
# tessellation so a stale cell list cannot fake an improvement.
import numpy as np, math, sys
from scipy.optimize import minimize
from evaluator import covering_density, as_cosets

N, MU2T = 5, 0.0625
il = np.tril_indices(N)

def batch_R2(LAT, COS, Q, C):
    P = LAT + C[COS]
    V = P[:, 1:, :] - P[:, :1, :]
    G = V @ Q @ V.transpose(0, 2, 1)
    q = np.diagonal(G, axis1=1, axis2=2)
    a = np.linalg.solve(G, q[:, :, None])[:, :, 0]
    return np.einsum('ij,ij->i', q, a) / 4.0

def descend(Q, c, rounds=12, verbose=True, freeze_c_rounds=0, max_rng=4):
    C = as_cosets(c)
    m = len(C)
    best = None; stall = 0
    for it in range(rounds):
        th, mu2, cells = covering_density(Q, C, max_rng=max_rng)
        Q = Q * (MU2T / mu2)
        LAT = np.array([cl[0] for cl in cells], dtype=float)
        COS = np.array([cl[1] for cl in cells], dtype=int)
        L0 = np.linalg.cholesky(Q)
        x0 = np.concatenate([L0[il], C[1:].ravel()])
        frozen = it < freeze_c_rounds
        def unpack(x):
            L = np.zeros((N, N)); L[il] = x[:15]
            Cx = C if frozen else np.vstack([np.zeros(N),
                                             x[15:].reshape(m - 1, N)])
            return L @ L.T, Cx
        def obj(x):
            L = np.zeros((N, N)); L[il] = x[:15]
            return -2 * np.sum(np.log(np.abs(np.diag(L)) + 1e-300))
        def cons(x):
            Qx, Cx = unpack(x)
            return MU2T - batch_R2(LAT, COS, Qx, Cx)
        res = minimize(obj, x0, method="SLSQP",
                       constraints=[{"type": "ineq", "fun": cons}],
                       options={"maxiter": 300, "ftol": 1e-14})
        Qn, Cn = unpack(res.x)
        thn, mu2n, _ = covering_density(Qn, Cn, max_rng=max_rng)
        if verbose:
            print(f"  round {it}: SLSQP {res.status}, verified Theta = {thn:.12f}")
        sys.stdout.flush()
        if best is None or thn < best[0] - 1e-11:
            best = (thn, Qn.copy(), Cn.copy()); stall = 0
        else:
            stall += 1
            if it >= freeze_c_rounds and stall >= 2:
                return best
        Q, C = Qn, Cn
    return best
