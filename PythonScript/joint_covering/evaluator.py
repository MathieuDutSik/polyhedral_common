# Numerical evaluator for the covering density of Z^5 + {0, c} at a form Q,
# via a genuine periodic Delaunay triangulation (scipy/qhull), plus the
# extraction of the translation classes of cells used by the descent.
import numpy as np
from scipy.spatial import Delaunay
import itertools, math

N = 5
KAPPA5 = 8 * math.pi**2 / 15

def lll_unimodular(Q, delta=0.75):
    """U integral unimodular with U Q U^T reduced (float LLL on the basis)."""
    Lc = np.linalg.cholesky(Q)
    B = Lc.copy()                      # rows = lattice basis vectors
    U = np.eye(N, dtype=int)
    def gso(B):
        Bs = B.copy(); mu = np.zeros((N, N))
        for i in range(N):
            for j in range(i):
                mu[i, j] = B[i] @ Bs[j] / (Bs[j] @ Bs[j])
                Bs[i] -= mu[i, j] * Bs[j]
        return Bs, mu
    k = 1
    while k < N:
        Bs, mu = gso(B)
        for j in range(k - 1, -1, -1):
            r = round(mu[k, j])
            if r:
                B[k] -= r * B[j]; U[k] -= r * U[j]
                Bs, mu = gso(B)
        if (Bs[k] @ Bs[k]) >= (delta - mu[k, k-1]**2) * (Bs[k-1] @ Bs[k-1]):
            k += 1
        else:
            B[[k-1, k]] = B[[k, k-1]]; U[[k-1, k]] = U[[k, k-1]]
            k = max(k - 1, 1)
    return U

def periodic_points(c, rng=2, U=None):
    """Lattice + coset points; returns (positions_frac, labels) where a
    label is (lattice_vector, coset_index)."""
    if U is None:
        U = np.eye(N, dtype=int)
    pts, labs = [], []
    for b in itertools.product(range(-rng, rng + 1), repeat=N):
        a = np.array(b, dtype=int) @ U      # reduced-basis box, original coords
        af = a.astype(float)
        binf = max(abs(x) for x in b)       # for the boundary-touch test
        pts.append(af); labs.append((tuple(int(x) for x in a), 0, binf))
        pts.append(af + c); labs.append((tuple(int(x) for x in a), 1, binf))
    return np.array(pts), labs

def circumradius2(P, Q):
    """Squared circumradius (metric Q) of the simplex with rows P (6 x 5)."""
    V = P[1:] - P[0]
    G = V @ Q @ V.T
    q = np.diag(G).copy()
    a = np.linalg.solve(G, q)
    return float(q @ a) / 4.0

def simplex_volume(P):
    V = P[1:] - P[0]
    return abs(np.linalg.det(V))

def cell_classes(Q, c, rng=2, max_rng=4):
    """Translation classes of Delone cells of Z^5+{0,c} at Q.
    Each class: (6 x 5 int lattice parts, 6 coset indices).

    If a kept cell touches the shell of the enumeration box, the box was too
    small for the star and qhull filled the region with spurious cells whose
    circumradii over-estimate mu; the box then grows and everything is
    redone -- iteratively, releasing the previous triangulation first, since
    a 5-dimensional qhull run on a large box costs gigabytes. Beyond max_rng
    a RuntimeError is raised and the caller treats the form as unusable.
    """
    U = lll_unimodular(Q)
    L = np.linalg.cholesky(Q)
    while True:
        pts, labs = periodic_points(c, rng, U=U)
        tri = Delaunay(pts @ L, qhull_options="QJ")
        classes = {}
        touches_shell = False
        for simp in tri.simplices:
            lat = np.array([labs[i][0] for i in simp], dtype=int)
            cos = np.array([labs[i][1] for i in simp], dtype=int)
            if simplex_volume(lat.astype(float) + np.outer(cos, c)) < 1e-9:
                continue
            if not np.any(np.all(lat == 0, axis=1) & (cos == 0)):
                continue
            if max(labs[i][2] for i in simp) >= rng:
                touches_shell = True
                break
            i0 = int(np.where(np.all(lat == 0, axis=1) & (cos == 0))[0][0])
            lat0 = lat - lat[i0]
            order = np.lexsort(np.c_[lat0, cos].T[::-1])
            key = tuple(map(tuple, np.c_[lat0, cos][order]))
            classes[key] = (lat0[order], cos[order])
        del tri, pts, labs
        if not touches_shell:
            return list(classes.values())
        if rng >= max_rng:
            raise RuntimeError("enumeration box insufficient at max_rng=%d" % max_rng)
        rng += 1

def positions(cls, c):
    lat, cos = cls
    return lat.astype(float) + np.outer(cos, c)

def covering_density(Q, c, cells=None, rng=2, max_rng=4):
    if cells is None:
        cells = cell_classes(Q, c, rng, max_rng)
    mu2 = max(circumradius2(positions(cl, c), Q) for cl in cells)
    det = np.linalg.det(Q)
    return 2 * KAPPA5 * mu2**2.5 / math.sqrt(det), mu2, cells

if __name__ == "__main__":
    # the certified optimum: Q in Q(sqrt2), c = (3/4,1/4,1/2,1/4,0)
    r2 = math.sqrt(2)
    Q0 = np.array([
        [(9+2*r2)/80, (-9+6*r2)/80, 0, r2/10, 9/80],
        [(-9+6*r2)/80, (27+18*r2)/80, -9/80, 3*r2/10, -9/80],
        [0, -9/80, 9/40, 9/80, 9/80],
        [r2/10, 3*r2/10, 9/80, (9+16*r2)/40, 0],
        [9/80, -9/80, 9/80, 0, 9/40]])
    c0 = np.array([0.75, 0.25, 0.5, 0.25, 0.0])
    th, mu2, cells = covering_density(Q0, c0)
    print(f"cell classes: {len(cells)}")
    print(f"mu^2 = {mu2:.12f}   (certified: 1.0)")
    print(f"Theta = {th:.12f}   (certified: 2.160060765053)")
    r2s = sorted(circumradius2(positions(cl, c0), Q0) for cl in cells)
    print("largest five R^2:", [f"{v:.6f}" for v in r2s[-5:]])
