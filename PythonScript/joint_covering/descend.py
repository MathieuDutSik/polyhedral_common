# Joint local descent on (Q, c) for the covering density of Z^5 + {0, c}.
#
# Outer loop: tessellate (scipy periodic Delaunay), then with the cell list
# fixed run SLSQP on (cholesky(Q), c):
#     minimize -log det Q   s.t.  R^2_S(Q, c) <= mu2_target  for every class S.
# Re-tessellate and repeat until the verified density stops improving. The
# verification at each round is a fresh tessellation, so a stale cell list
# (whose max would under-estimate mu) cannot produce a fake improvement.
import numpy as np, math, sys
from scipy.optimize import minimize
from evaluator import cell_classes, circumradius2, positions, covering_density

N, KAPPA5, MU2T = 5, 8*math.pi**2/15, 0.0625
il = np.tril_indices(N)

def batch_R2(LAT, COS, Q, c):
    P = LAT + COS[:, :, None] * c
    V = P[:, 1:, :] - P[:, :1, :]
    G = V @ Q @ V.transpose(0, 2, 1)
    q = np.diagonal(G, axis1=1, axis2=2)
    a = np.linalg.solve(G, q[:, :, None])[:, :, 0]
    return np.einsum('ij,ij->i', q, a) / 4.0

def descend(Q, c, rounds=12, verbose=True, freeze_c_rounds=0, max_rng=4):
    """freeze_c_rounds: rounds during which c is held fixed, so Q first
    reaches the per-domain (fixed-c) optimum before the joint phase."""
    best = None; stall = 0
    for it in range(rounds):
        th, mu2, cells = covering_density(Q, c, max_rng=max_rng)
        Q = Q * (MU2T / mu2)                     # normalize the scale
        LAT = np.array([cl[0] for cl in cells], dtype=float)
        COS = np.array([cl[1] for cl in cells], dtype=float)
        L0 = np.linalg.cholesky(Q)
        x0 = np.concatenate([L0[il], c])
        def unpack(x):
            L = np.zeros((N, N)); L[il] = x[:15]
            return L @ L.T, x[15:]
        def obj(x):
            L = np.zeros((N, N)); L[il] = x[:15]
            return -2 * np.sum(np.log(np.abs(np.diag(L)) + 1e-300))
        frozen = it < freeze_c_rounds
        def cons(x):
            Qx, cx = unpack(x)
            if frozen:
                cx = c
            return MU2T - batch_R2(LAT, COS, Qx, cx)
        res = minimize(obj, x0, method="SLSQP",
                       constraints=[{"type": "ineq", "fun": cons}],
                       options={"maxiter": 300, "ftol": 1e-14})
        Qn, cn = unpack(res.x)
        if frozen:
            cn = c
        thn, mu2n, _ = covering_density(Qn, cn, max_rng=max_rng)  # verified
        if verbose:
            print(f"  round {it}: SLSQP {res.status}, verified Theta = {thn:.12f}")
        sys.stdout.flush()
        if best is None or thn < best[0] - 1e-11:
            best = (thn, Qn.copy(), cn.copy()); stall = 0
        else:
            stall += 1
            # never stop inside the frozen phase, and give the joint phase
            # two non-improving rounds before giving up
            if it >= freeze_c_rounds and stall >= 2:
                return best
        Q, c = Qn, cn
    return best

if __name__ == "__main__":
    r2 = math.sqrt(2)
    Q0 = np.array([
        [(9+2*r2)/80, (-9+6*r2)/80, 0, r2/10, 9/80],
        [(-9+6*r2)/80, (27+18*r2)/80, -9/80, 3*r2/10, -9/80],
        [0, -9/80, 9/40, 9/80, 9/80],
        [r2/10, 3*r2/10, 9/80, (9+16*r2)/40, 0],
        [9/80, -9/80, 9/80, 0, 9/40]])
    c0 = np.array([0.75, 0.25, 0.5, 0.25, 0.0])
    print("joint descent from the certified optimum (2.160060765053):")
    th, Q, c = descend(Q0, c0)
    print(f"final verified Theta = {th:.12f}")
    print("final c =", np.round(c, 8))
