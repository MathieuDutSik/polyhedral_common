// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "PeriodicDelaunay.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

/*
  The subgroup of the transformations preserving a periodic point set is
  checked against the direct computation: a group of transformations is
  built by hand, all its elements are enumerated, and the elements that
  preserve the point set (as decided by PeriodicCosetPermutation, itself
  checked against a membership test in TEST_PeriodicCVP) must be exactly
  the elements of the subgroup returned.
 */

using T = mpq_class;
using Tring = mpz_class;
using Tidx = uint32_t;
using Telt = permutalib::SingleSidedPerm<Tidx>;
using Tint_grp = mpz_class;
using Tgroup = permutalib::Group<Telt, Tint_grp>;
using Ttrans = MyMatrix<T>;

void check(bool test, std::string const &context) {
  if (!test) {
    std::cerr << "TEST_PeriodicCosetSubgroup failed: " << context << "\n";
    throw TerminalException{1};
  }
}

int main() {
  try {
    std::ostream &os = std::cerr;
    int n = 2;
    // The point set Z^2 + {(0,0), (1/3,0)}, of denominator 3. It is
    // genuinely periodic: (1/3,0) + (1/3,0) is (2/3,0), not a coset, so
    // the cosets are not a group and the set is not a lattice.
    MyMatrix<T> Cosets(2, n);
    Cosets(0, 0) = 0;
    Cosets(0, 1) = 0;
    Cosets(1, 0) = T(1) / T(3);
    Cosets(1, 1) = 0;
    PeriodicPointSet<Tring> pps = PeriodicPointSetFromRational<Tring>(Cosets);
    check(pps.N == 3, "denominator of the point set");
    check(!IsPeriodicPointSetLattice(pps), "the point set is not a lattice");
    //
    // Two transformations of denominator 2, acting on a set of four
    // formal points so that the permutation determines the transformation:
    //   g1: preserves the set,
    //   g2: does not.
    // g1 negates the second coordinate, which the cosets, all of second
    // coordinate zero, do not see: it preserves the set. g2 negates the
    // first, sending the coset 1/3 to -1/3, that is 2/3, which is not a
    // coset: it does not. Both are involutions and they commute.
    // As homogeneous integral affine matrices, no translation.
    Ttrans g1 = IdentityMat<T>(n + 1);
    g1(2, 2) = -1;
    Ttrans g2 = IdentityMat<T>(n + 1);
    g2(1, 1) = -1;
    check(IsPeriodicPointSetPreserved<Tring>(pps, g1), "g1 preserves the set");
    check(!IsPeriodicPointSetPreserved<Tring>(pps, g2), "g2 does not preserve it");
    //
    // The group is realized on four points, the permutation of the points
    // determining which product of g1 and g2 is meant. Both generators are
    // involutions and they commute, so the group is (Z/2)^2 acting on the
    // four pairs of exponents.
    //   point i encodes (a, b) with i = a + 2 b, meaning g1^a g2^b.
    auto get_perm = [&](int i_gen) -> Telt {
      std::vector<Tidx> eList(4);
      for (int i = 0; i < 4; i++) {
        int a = i % 2, b = i / 2;
        if (i_gen == 0) {
          a = 1 - a;
        } else {
          b = 1 - b;
        }
        eList[i] = static_cast<Tidx>(a + 2 * b);
      }
      return Telt(eList);
    };
    std::vector<Telt> LGen{get_perm(0), get_perm(1)};
    Tgroup GRP(LGen, 4);
    check(GRP.size() == 4, "the group has order 4");
    // The transformation realized by a permutation, read off the image of
    // the identity point.
    auto f_trans = [&](Telt const &eElt) -> Ttrans {
      int img = OnPoints(0, eElt);
      int a = img % 2, b = img / 2;
      Ttrans ret = IdentityMat<T>(n + 1);
      if (a == 1) {
        ret = ret * g1;
      }
      if (b == 1) {
        ret = ret * g2;
      }
      return ret;
    };
    Tgroup GRPsub =
        PeriodicCosetPreservingSubgroup<Tring, Tgroup, decltype(f_trans)>(
            pps, GRP, f_trans, os);
    // Direct computation: the elements whose transformation preserves the
    // point set. The subgroup acts on the same four points, so the two can
    // be compared element by element.
    size_t n_pres = 0;
    for (auto &eElt : GRP.get_all_element()) {
      if (IsPeriodicPointSetPreserved<Tring>(pps, f_trans(eElt))) {
        n_pres++;
      }
    }
    std::cerr << "|GRP|=" << GRP.size() << " |GRPsub|=" << GRPsub.size()
              << " n_preserving=" << n_pres << "\n";
    check(UniversalScalarConversion<size_t, Tint_grp>(GRPsub.size()) == n_pres,
          "the subgroup has the number of preserving elements");
    // Every element of the subgroup preserves the point set.
    for (auto &eElt : GRPsub.get_all_element()) {
      check(IsPeriodicPointSetPreserved<Tring>(pps, f_trans(eElt)),
            "an element of the subgroup preserves the point set");
    }
    check(n_pres == 2, "exactly the two elements without g2");
    //
    // The enumeration of the Delaunay cells, and with it every piece of
    // the periodic path: the initial cell and the adjacent ones from the
    // geometry templated on the solver, the stabilizers, the orbits of
    // facets from the dual description and the recognition of an already
    // found cell from the periodic equivalence.
    //
    // The properties checked are those defining the answer rather than
    // values read off a previous run: each cell is a Delaunay cell of the
    // point set, its stabilizer is the one the stabilizer computation
    // gives on its own, each adjacency transformation preserves the point
    // set, and cells of different orbits are not equivalent.
    {
      MyMatrix<T> GramMat = IdentityMat<T>(n);
      MyMatrix<T> SHV(4, n);
      SHV(0, 0) = 1;  SHV(0, 1) = 0;
      SHV(1, 0) = -1; SHV(1, 1) = 0;
      SHV(2, 0) = 0;  SHV(2, 1) = 1;
      SHV(3, 0) = 0;  SHV(3, 1) = -1;
      MyMatrix<Tring> Graver(4, n);
      Graver(0, 0) = 1;  Graver(0, 1) = 0;
      Graver(1, 0) = -1; Graver(1, 1) = 0;
      Graver(2, 0) = 0;  Graver(2, 1) = 1;
      Graver(3, 0) = 0;  Graver(3, 1) = -1;
      PolyHeuristicSerial<Tint_grp> AllArr =
          AllStandardHeuristicSerial<T, Tint_grp>(n + 1, os);
      PeriodicDataDelaunay<T, Tring, Tgroup> data =
          GetPeriodicDataDelaunay<T, Tring, Tgroup>(GramMat, pps, SHV, Graver,
                                                    AllArr, os);
      PeriodicDataDelaunayFunc<T, Tring, Tgroup> data_func{data};
      auto f_incorrect = [&]([[maybe_unused]] PeriodicDelaunay_Obj<
                             Tring, Tgroup> const &x) -> bool {
        return false;
      };
      std::optional<std::vector<DatabaseEntry_Serial<
          PeriodicDelaunay_Obj<Tring, Tgroup>,
          PeriodicDelaunay_AdjO<Tring>>>>
          opt = EnumerateAndStore_Serial(data_func, f_incorrect, 0);
      check(opt.has_value(), "the enumeration terminated");
      std::vector<DatabaseEntry_Serial<PeriodicDelaunay_Obj<Tring, Tgroup>,
                                       PeriodicDelaunay_AdjO<Tring>>> const
          &l_ent = *opt;
      std::cerr << "|orbits of periodic Delaunay cells|=" << l_ent.size()
                << "\n";
      check(l_ent.size() > 0, "the enumeration found something");
      PeriodicCVPSolver<T, Tring> psolver(GramMat, pps, os);
      for (auto &eEnt : l_ent) {
        MyMatrix<Tring> const &EXT = eEnt.x.EXT;
        int nbVert = EXT.rows();
        std::cerr << "   |vertices|=" << nbVert
                  << " |stabilizer|=" << eEnt.x.GRP.size()
                  << " |adjacencies|=" << eEnt.ListAdj.size() << "\n";
        check(nbVert >= n + 1, "the cell is full dimensional");
        // Its vertices are points of the point set.
        for (int i = 0; i < nbVert; i++) {
          MyVector<Tring> u(n);
          for (int j2 = 0; j2 < n; j2++) {
            u(j2) = EXT(i, j2 + 1);
          }
          check(GetCosetIndex(pps, u).has_value(),
                "a vertex of a cell is a point of the point set");
        }
        // It is a Delaunay cell: nothing of the set is strictly closer to
        // its circumcenter than its vertices, and they are all at the
        // circumradius.
        MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tring>(EXT);
        CP<T> cp = CenterRadiusDelaunayPolytopeGeneral<T>(GramMat, EXT_T);
        MyVector<T> Cent(n);
        for (int j2 = 0; j2 < n; j2++) {
          Cent(j2) = cp.eCent(j2 + 1);
        }
        resultCVP<T, Tring> res = psolver.nearest_vectors(Cent);
        check(res.TheNorm == cp.SquareRadius,
              "nothing of the set is closer than the circumradius");
        check(res.ListVect.rows() == nbVert,
              "the closest points are exactly the vertices");
        // The stabilizer of the enumeration is the stabilizer.
        Tgroup GRPdirect = PeriodicDelaunay_Stabilizer<T, Tring, Tgroup>(
            GramMat, pps, SHV, EXT, os);
        check(GRPdirect.size() == eEnt.x.GRP.size(),
              "the stabilizer of the enumeration is the direct one");
        // Every adjacency transformation preserves the point set.
        for (auto &eAdj : eEnt.ListAdj) {
          check(IsPeriodicPointSetPreserved<Tring>(pps, eAdj.x.eBigMat),
                "an adjacency transformation preserves the point set");
        }
      }
      // A cell is equivalent to its translate by a period of the set, and
      // cells of different orbits are not equivalent to each other.
      {
        MyMatrix<Tring> const &EXT0 = l_ent[0].x.EXT;
        MyMatrix<Tring> EXT_tr = EXT0;
        for (int i = 0; i < EXT_tr.rows(); i++) {
          EXT_tr(i, 1) += pps.N;
        }
        check(PeriodicDelaunay_TestEquivalence<T, Tring, Tgroup>(
                  GramMat, pps, SHV, EXT0, EXT_tr, os)
                  .has_value(),
              "a cell is equivalent to its translate by a period");
        for (size_t i1 = 0; i1 < l_ent.size(); i1++) {
          for (size_t i2 = i1 + 1; i2 < l_ent.size(); i2++) {
            check(!PeriodicDelaunay_TestEquivalence<T, Tring, Tgroup>(
                       GramMat, pps, SHV, l_ent[i1].x.EXT, l_ent[i2].x.EXT, os)
                       .has_value(),
                  "cells of different orbits are not equivalent");
          }
        }
      }
    }
    //
    // The tessellation in the form the iso-Delaunay machinery consumes.
    // It needs a GENERIC form, whose Delaunay cells are simplices: the
    // identity form above gives rectangles, which are co-spherical only
    // on a wall and where the inequalities are not defined. The form is
    // perturbed off the diagonal to reach the interior of a domain.
    {
      MyMatrix<T> GramMat(n, n);
      GramMat(0, 0) = 1;
      GramMat(0, 1) = T(1) / T(7);
      GramMat(1, 0) = T(1) / T(7);
      GramMat(1, 1) = 1;
      MyMatrix<T> SHV(4, n);
      SHV(0, 0) = 1;  SHV(0, 1) = 0;
      SHV(1, 0) = -1; SHV(1, 1) = 0;
      SHV(2, 0) = 0;  SHV(2, 1) = 1;
      SHV(3, 0) = 0;  SHV(3, 1) = -1;
      MyMatrix<Tring> Graver(4, n);
      Graver(0, 0) = 1;  Graver(0, 1) = 0;
      Graver(1, 0) = -1; Graver(1, 1) = 0;
      Graver(2, 0) = 0;  Graver(2, 1) = 1;
      Graver(3, 0) = 0;  Graver(3, 1) = -1;
      PolyHeuristicSerial<Tint_grp> AllArr =
          AllStandardHeuristicSerial<T, Tint_grp>(n + 1, os);
      PeriodicDataDelaunay<T, Tring, Tgroup> data =
          GetPeriodicDataDelaunay<T, Tring, Tgroup>(GramMat, pps, SHV, Graver,
                                                    AllArr, os);
      PeriodicDataDelaunayFunc<T, Tring, Tgroup> data_func{data};
      auto f_inc = [&]([[maybe_unused]] PeriodicDelaunay_Obj<Tring, Tgroup>
                       const &x) -> bool { return false; };
      auto opt = EnumerateAndStore_Serial(data_func, f_inc, 0);
      check(opt.has_value(), "the enumeration at the generic form finished");
      std::cerr << "|orbits at the generic form|=" << opt->size() << "\n";
      // Generic means simplicial: n+1 vertices per cell.
      for (auto &eEnt : *opt) {
        std::cerr << "   |vertices|=" << eEnt.x.EXT.rows() << "\n";
        check(eEnt.x.EXT.rows() == n + 1,
              "the form is generic, its cells are simplices");
      }
      // The T-space of all the symmetric matrices of order n.
      std::vector<std::vector<T>> ListGram;
      for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
          MyMatrix<T> eMat = ZeroMatrix<T>(n, n);
          eMat(i, j) = 1;
          eMat(j, i) = 1;
          ListGram.push_back(GetLineVector(eMat));
        }
      }
      DelaunayTesselationIneq<T, Tgroup> DTI =
          BuildPeriodicDelaunayTesselationIneq<T, Tring, Tgroup>(
              *opt, pps, ListGram, os);
      check(DTI.l_dels.size() == opt->size(),
            "the tessellation has the orbits of the enumeration");
      // Consistency: across each adjacency the image of the neighbour
      // meets the cell exactly in the vertices of eInc.
      for (auto &eDel : DTI.l_dels) {
        MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tring>(eDel.EXT);
        ContainerMatrix<T> cont(EXT_T);
        for (auto &eAdj : eDel.ListAdj) {
          MyMatrix<T> EXTadj =
              UniversalMatrixConversion<T, Tring>(DTI.l_dels[eAdj.iOrb].EXT) *
              UniversalMatrixConversion<T, Tring>(eAdj.eBigMat);
          Face f_att(EXT_T.rows());
          for (int u = 0; u < EXTadj.rows(); u++) {
            MyVector<T> V = GetMatrixRow(EXTadj, u);
            std::optional<size_t> o = cont.GetIdx_v(V);
            if (o) {
              f_att[*o] = 1;
            }
          }
          check(f_att == eAdj.eInc,
                "the adjacency meets the cell in the vertices of eInc");
          check(eAdj.eIneq.size() > 0, "the adjacency carries an inequality");
        }
      }
      std::cerr << "periodic tessellation built and consistent\n";
    }
    std::cerr << "Normal termination of TEST_PeriodicCosetSubgroup\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of TEST_PeriodicCosetSubgroup\n";
    exit(e.eVal);
  }
}
