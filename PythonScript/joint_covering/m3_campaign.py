# The m=3 campaign in the continuum. Three seed families:
#   mode walk : from the m3N6 flip-walk optimum (2.3766)
#   mode add  : the m=2 optimum (2.16006) plus a random third coset -- can a
#               third point pay for its 3/2 density factor?
#   mode rand : random (Q, c1, c2)
import numpy as np, re, sys, math
from fractions import Fraction
from descend import descend

S = "/private/tmp/claude-501/-Users-mathieudutoursikiric-GITall-GITnon-lattice-covering-polyhedral-common/5558689f-f195-4111-b432-bbf333af0ed0/scratchpad"
mode, lane = sys.argv[1], int(sys.argv[2])
rng = np.random.default_rng(770000 + lane * 3391)
r2 = math.sqrt(2)
Q2 = np.array([
    [(9+2*r2)/80, (-9+6*r2)/80, 0, r2/10, 9/80],
    [(-9+6*r2)/80, (27+18*r2)/80, -9/80, 3*r2/10, -9/80],
    [0, -9/80, 9/40, 9/80, 9/80],
    [r2/10, 3*r2/10, 9/80, (9+16*r2)/40, 0],
    [9/80, -9/80, 9/80, 0, 9/40]])
c2opt = np.array([0.75, 0.25, 0.5, 0.25, 0.0])
out = open(f"{S}/joint/m3_{mode}{lane}.txt", "a", buffering=1)
best = 1e9
for start in range(10000):
    if mode == "walk":
        t = open(f"{S}/rewalk3/out.g").read()
        mm = re.search(r'best_gram:=\[(.*?)\]\);', t, re.S)
        Q = np.array([[float(x) for x in r.split(',')]
                      for r in re.findall(r'\[([^\[\]]*)\]', mm.group(1))])
        lines = open(f"{S}/rewalk3/c.txt").read().strip().split('\n')
        C = np.array([[float(Fraction(x)) for x in l.split()] for l in lines[1:]])
        if start > 0:   # after the exact seed, jitter it
            C = C + np.vstack([np.zeros(5), rng.normal(0, 0.04, size=(2, 5))])
    elif mode == "add":
        Q = Q2.copy()
        C = np.vstack([np.zeros(5), c2opt, rng.uniform(0, 1, size=5)])
    else:
        A = rng.normal(size=(5, 5))
        Q = A @ A.T + 0.05 * np.eye(5)
        C = np.vstack([np.zeros(5), rng.uniform(0, 1, size=(2, 5))])
    try:
        th, Qf, Cf = descend(Q, C, rounds=15, verbose=False)
    except Exception as e:
        out.write(f"{mode}{lane} {start} FAIL {type(e).__name__}\n"); continue
    out.write(f"{mode}{lane} {start} {th:.12f}\n")
    if th < best:
        best = th
        np.savez(f"{S}/joint/m3_best_{mode}{lane}.npz", Q=Qf, C=Cf, theta=th)
