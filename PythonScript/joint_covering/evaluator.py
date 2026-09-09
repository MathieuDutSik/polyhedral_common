# Covering density of the periodic set Z^5 + {C_0=0, C_1, ..., C_{m-1}} at a
# form Q, via a genuine periodic Delaunay triangulation (scipy/qhull).
#
# Soundness notes.
# * The enumeration box is LLL-reduced and grows when a kept cell touches
#   its shell: a too-small box makes qhull fill the star with spurious cells
#   that over-estimate mu. The growth is iterative and releases the previous
#   triangulation (a recursive version held every level alive, 14 GB).
# * A cell is kept when ANY vertex has lattice part 0 -- anchoring only on
#   coset-0 vertices would drop the classes of cells whose vertices all lie
#   in other cosets, which exist since the set is not translation-symmetric
#   by a coset. The class key is the minimum over all vertex anchors, so a
#   class is counted once however it meets the central copy.
import numpy as np
from scipy.spatial import Delaunay
import itertools, math

N = 5
KAPPA5 = 8 * math.pi**2 / 15

def as_cosets(c):
    """Accept a single coset vector (m=2) or an (m-1) x 5 / m x 5 matrix."""
    C = np.atleast_2d(np.asarray(c, dtype=float))
    if not np.allclose(C[0], 0):
        C = np.vstack([np.zeros(N), C])
    return C

def lll_unimodular(Q, delta=0.75):
    Lc = np.linalg.cholesky(Q)
    B = Lc.copy(); U = np.eye(N, dtype=int)
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

def points_in_ball(Q, C, R, U):
    """All points a + C_t of the periodic set with Q-norm <= R, enumerated
    through the LLL-reduced basis. Returns (positions, labels), a label
    being (lattice_vector, coset_index)."""
    Qr = U @ Q @ U.T
    sig_min = math.sqrt(np.linalg.eigvalsh(Qr)[0])
    B0 = int(math.ceil(R / sig_min)) + 1
    axes = [np.arange(-B0, B0 + 1)] * N
    Bm = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, N)
    A = Bm @ U
    pts, labs = [], []
    for t in range(len(C)):
        V = A + C[t]
        keep = np.einsum('ij,jk,ik->i', V, Q, V) <= R * R
        for a, v in zip(A[keep], V[keep]):
            pts.append(v)
            labs.append((tuple(int(x) for x in a), t))
    return np.array(pts), labs

def simplex_volume(P):
    return abs(np.linalg.det(P[1:] - P[0]))

def circumradius2(P, Q):
    V = P[1:] - P[0]
    G = V @ Q @ V.T
    q = np.diag(G).copy()
    a = np.linalg.solve(G, q)
    return float(q @ a) / 4.0

def circumcenter_R2(P, Q):
    """Circumcenter and squared circumradius of the simplex with rows P."""
    V = P[1:] - P[0]
    G = V @ Q @ V.T
    q = np.diag(G).copy()
    a = np.linalg.solve(G, q)
    z = P[0] + 0.5 * a @ V
    return z, float(q @ a) / 4.0

def cell_classes(Q, c, rng=None, max_rng=None, growth=1.35, max_iter=6):
    """Translation classes of Delone cells of Z^5 + cosets at Q.

    The point cloud is the ball |x|_Q <= R around the origin, so a cell is
    provably a Delone cell of the full periodic set once its circumsphere
    lies inside the ball: |center|_Q + radius <= R. Cells with a lattice-0
    vertex failing that test mean R was too small, and R grows. rng/max_rng
    are accepted and ignored for compatibility with older callers.
    """
    C = as_cosets(c)
    U = lll_unimodular(Q)
    L = np.linalg.cholesky(Q)
    # covering radius of the sublattice Z^5 bounds the one of the union
    Qr = U @ Q @ U.T
    mu_bound = 0.5 * math.sqrt(np.trace(Qr))
    R = 2.5 * mu_bound
    for _ in range(max_iter):
        pts, labs = points_in_ball(Q, C, R, U)
        tri = Delaunay(pts @ L, qhull_options="QJ")
        classes = {}
        ok = True
        for simp in tri.simplices:
            lat = np.array([labs[i][0] for i in simp], dtype=int)
            cos = np.array([labs[i][1] for i in simp], dtype=int)
            if not np.any(np.all(lat == 0, axis=1)):
                continue
            P = lat.astype(float) + C[cos]
            if simplex_volume(P) < 1e-9:
                continue
            z, R2 = circumcenter_R2(P, Q)
            if math.sqrt(z @ Q @ z) + math.sqrt(max(R2, 0)) > R:
                ok = False
                break
            best_key = None
            for i0 in range(len(simp)):
                lat0 = lat - lat[i0]
                order = np.lexsort(np.c_[lat0, cos].T[::-1])
                key = tuple(map(tuple, np.c_[lat0, cos][order]))
                if best_key is None or key < best_key:
                    best_key = key; best_val = (lat0[order], cos[order])
            classes[best_key] = best_val
        del tri, pts, labs
        if ok:
            return list(classes.values())
        R *= growth
    raise RuntimeError("ball radius did not stabilize")

def positions(cl, c):
    C = as_cosets(c)
    lat, cos = cl
    return lat.astype(float) + C[cos]

def covering_density(Q, c, cells=None, rng=2, max_rng=4):
    C = as_cosets(c)
    if cells is None:
        cells = cell_classes(Q, C, rng, max_rng)
    mu2 = max(circumradius2(positions(cl, C), Q) for cl in cells)
    det = np.linalg.det(Q)
    return len(C) * KAPPA5 * mu2**2.5 / math.sqrt(det), mu2, cells

if __name__ == "__main__":
    r2 = math.sqrt(2)
    Q0 = np.array([
        [(9+2*r2)/80, (-9+6*r2)/80, 0, r2/10, 9/80],
        [(-9+6*r2)/80, (27+18*r2)/80, -9/80, 3*r2/10, -9/80],
        [0, -9/80, 9/40, 9/80, 9/80],
        [r2/10, 3*r2/10, 9/80, (9+16*r2)/40, 0],
        [9/80, -9/80, 9/80, 0, 9/40]])
    c0 = np.array([0.75, 0.25, 0.5, 0.25, 0.0])
    th, mu2, cells = covering_density(Q0, c0)
    print(f"m=2 regression: cells={len(cells)}  Theta={th:.12f}  (certified 2.160060765053)")
