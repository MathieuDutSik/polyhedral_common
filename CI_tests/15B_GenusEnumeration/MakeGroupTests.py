#!/usr/bin/env python3
"""
The data of the tests of GENUS_GlueIsometry factors for a finite group G,
abelian or not (RunCharPolyTests.sh with the table group_tests.txt).

For every genus LAT and group G with an oracle table
group/oracle_LAT_G.txt (GroupOracle.g: a comment line with the degrees of
Irr(G), then lines "m1,...,mk mass", the multiplicities of the irreducible
characters in the order of GAP) and a description group/group_G.txt
(GroupIrreps.g: a line per irreducible character, "orth s_1 ... s_k",
"herm d e_1 ... e_k", "conj i", "irrep FILE m" or "unsupported ..."),
write
- charpoly/LAT_gG_m1_..._mk.spec: "group k", then a factor per irreducible
  character of nonzero multiplicity m: "orth s_1 ... s_k m" for a linear
  character of values s_j = +-1 on the generators, "irrep FILE m" for an
  absolutely irreducible character of higher degree, FILE holding a
  G-lattice in its rational representation, of Schur index m: c = the
  multiplicity / m (IsotypicPieces.h); "herm d b e_1 ... e_k" for a
  linear character of values zeta_d^e_j and its Galois conjugates, all of
  multiplicity b; no spec when a character of nonzero multiplicity is
  unsupported;
- group_tests.txt: the lines "LAT gG_m1_..._mk mass".

Usage (from CI_tests/15B_GenusEnumeration), the generators of G as
permutations, e.g. S3 = [(1,2,3),(1,2)], D8 = [(1,2,3,4),(1,3)],
S4 = [(1,2,3,4),(1,2)]:
  1. the description, in GAP:
       GroupGens:=[(1,2,3),(1,2)];; OutFile:="group/group_S3.txt";;
       Prefix:="group/irrep_S3";; Read("GroupIrreps.g");
  2. the oracle, in GAP, for each genus LAT:
       InFile:="charpoly/genus_LAT.g";; OutFile:="group/oracle_LAT_S3.txt";;
       GroupGens:=[(1,2,3),(1,2)];; Read("GroupOracle.g");
  3. python3 MakeGroupTests.py LAT_G ...
  4. RunCharPolyTests.sh group_tests.txt
"""
import sys


def read_group(grp):
    irr = []
    with open(f"group/group_{grp}.txt") as f:
        for line in f:
            w = line.split()
            if not w:
                continue
            if w[0] == "orth":
                irr.append(("orth", [int(x) for x in w[1:]]))
            elif w[0] == "herm":
                irr.append(("herm", [int(x) for x in w[1:]]))
            elif w[0] == "conj":
                irr.append(("conj", int(w[1]) - 1))
            elif w[0] == "irrep":
                irr.append(("irrep", (w[1], int(w[2]) if len(w) > 2 else 1)))
            else:
                irr.append(("unsupported", " ".join(w[1:])))
    n_gen = max((len(d) for k, d in irr if k == "orth"), default=0)
    return {"n_gen": n_gen, "irr": irr}


def main(args):
    lines = []
    for arg in args:
        lat, grp = arg.rsplit("_", 1)
        g = read_group(grp)
        with open(f"group/oracle_{lat}_{grp}.txt") as f:
            for line in f:
                line = line.strip()
                if line.startswith("# degrees"):
                    degs = [int(x) for x in line.split()[2].split(",")]
                    if len(degs) != len(g["irr"]):
                        sys.exit(f"{arg}: {len(degs)} characters in the "
                                 f"oracle, {len(g['irr'])} in the group")
                    continue
                if not line or line.startswith("#"):
                    continue
                typ, mass = line.split()
                mult = [int(x) for x in typ.split(",")]
                sig = "g" + grp + "_" + "_".join(str(m) for m in mult)
                spec = [f"group {g['n_gen']}"]
                if any(m != 0 and kind == "unsupported"
                       for (kind, data), m in zip(g["irr"], mult)):
                    continue
                ok = True
                for (kind, data), m in zip(g["irr"], mult):
                    if kind == "conj" and m != mult[data]:
                        ok = False
                    if kind == "irrep" and m % data[1] != 0:
                        ok = False
                if not ok:
                    sys.exit(f"{arg} {typ}: multiplicities inconsistent "
                             "with the Galois orbits or Schur indices")
                for (kind, data), m in zip(g["irr"], mult):
                    if m == 0 or kind == "conj":
                        continue
                    if kind == "orth":
                        spec.append("orth " + " ".join(str(s) for s in data)
                                    + f" {m}")
                    elif kind == "herm":
                        spec.append(f"herm {data[0]} {m} "
                                    + " ".join(str(e) for e in data[1:]))
                    else:
                        spec.append(f"irrep {data[0]} {m // data[1]}")
                with open(f"charpoly/{lat}_{sig}.spec", "w") as f2:
                    f2.write("\n".join(spec) + "\n")
                lines.append(f"{lat} {sig} {mass}")
    with open("group_tests.txt", "w") as f:
        f.write("# LAT sig mass: the oracle of GroupOracle.g\n")
        for line in lines:
            f.write(line + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
