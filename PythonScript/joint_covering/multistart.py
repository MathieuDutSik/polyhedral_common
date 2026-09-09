# Multistart joint descent: random (Q, c), descend, log. The c coordinate is
# continuous -- no denominator, no -I condition -- so this searches the full
# 2-point family, which the exact flipping machinery could not.
import numpy as np, math, sys
from descend import descend

lane = int(sys.argv[1])
rng = np.random.default_rng(900000 + lane * 1231)
out = open(f"ms_lane{lane}.txt", "a", buffering=1)
best = 1e9
for start in range(10000):
    A = rng.normal(size=(5, 5))
    Q = A @ A.T + 0.05 * np.eye(5)
    c = rng.uniform(0, 1, size=5)
    try:
        th, Qf, cf = descend(Q, c, rounds=15, verbose=False, max_rng=3)
    except Exception as e:
        out.write(f"{lane} {start} FAIL {type(e).__name__}\n"); continue
    out.write(f"{lane} {start} {th:.12f} c={np.round(cf,6).tolist()}\n")
    if th < best:
        best = th
        np.savez(f"ms_best_lane{lane}.npz", Q=Qf, c=cf, theta=th)
