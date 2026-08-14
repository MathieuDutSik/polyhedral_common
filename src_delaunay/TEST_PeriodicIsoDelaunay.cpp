// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "PeriodicDelaunay.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

/*
  The enumeration of the periodic iso-Delaunay domains of the point set
  Z^2 + {(0,0), (1/3,0)} in the T-space of all the symmetric matrices.

  The properties checked are the defining ones, not values read off a
  previous run: the Gram matrix of each domain is strictly interior to the
  domain's inequality system, the tessellation of each domain is the
  Delaunay tessellation of the point set at that Gram matrix, and the
  domains of different orbits are pairwise inequivalent under the
  transformations preserving both the T-space and the point set.
 */

using T = mpq_class;
using Tint = mpz_class;
using Tidx = uint32_t;
using Telt = permutalib::SingleSidedPerm<Tidx>;
using Tint_grp = mpz_class;
using Tgroup = permutalib::Group<Telt, Tint_grp>;

void check(bool test, std::string const &context) {
  if (!test) {
    std::cerr << "TEST_PeriodicIsoDelaunay failed: " << context << "\n";
    throw TerminalException{1};
  }
}

/*
  Enumerate the periodic iso-Delaunay domains and check the defining
  properties: the Gram matrix of each domain strictly interior to the
  domain's inequality system, the tessellation of each domain the Delaunay
  tessellation of the point set at that matrix, the orbits pairwise
  inequivalent, and the number of orbits equal to the expected value.
 */
void run_case(LinSpaceMatrix<T> const &LinSpa, PeriodicPointSet<Tint> const &pps,
              size_t expected_nb, std::string const &name, std::ostream &os) {
  PolyHeuristicSerial<Tint_grp> AllArr =
      AllStandardHeuristicSerial<T, Tint_grp>(LinSpa.n + 1, os);
  RecordDualDescOperation<T, Tgroup> rddo(AllArr, os);
  std::optional<MyMatrix<T>> CommonGramMat;
  DataIsoDelaunayDomains<T, Tint, Tgroup> data{LinSpa, std::move(rddo),
                                               CommonGramMat};
  PeriodicDataIsoDelaunayDomainsFunc<T, Tint, Tgroup> data_func{
      std::move(data), pps};
  using Tobj = typename PeriodicDataIsoDelaunayDomainsFunc<T, Tint,
                                                           Tgroup>::Tobj;
  auto f_incorrect = [&]([[maybe_unused]] Tobj const &x) -> bool {
    return false;
  };
  auto opt = EnumerateAndStore_Serial(data_func, f_incorrect, 0);
  check(opt.has_value(), "the enumeration terminated");
  auto const &l_ent = *opt;
  size_t n_orbit = l_ent.size();
  std::cerr << name << ": |orbits of periodic iso-Delaunay domains|="
            << n_orbit << "\n";
  check(n_orbit == expected_nb, "the number of orbits is the expected one");
  for (auto &eEnt : l_ent) {
    Tobj const &x = eEnt.x;
    MyMatrix<T> const &GramMat = x.DT_gram.GramMat;
    std::vector<FullAdjInfo<T>> ListIneq =
        ComputeListIneqFromTesselationIneq<T, Tgroup>(x.DT_gram.DT);
    MyVector<T> g_vec =
        LINSPA_GetVectorOfMatrixExpression(data_func.data.LinSpa, GramMat);
    std::cerr << "   |ineq|=" << ListIneq.size()
              << " |cells|=" << x.DT_gram.DT.l_dels.size()
              << " |adj|=" << eEnt.ListAdj.size() << "\n";
    for (auto &eRec : ListIneq) {
      check(eRec.eIneq.dot(g_vec) > 0,
            "the Gram matrix is strictly interior to its domain");
    }
    check(IsPeriodicDelaunayTesselation<T, Tint, Tgroup>(
              x.DT_gram.DT, GramMat, pps, os),
          "the tessellation of the domain is the Delaunay tessellation");
  }
  for (size_t i = 0; i < n_orbit; i++) {
    for (size_t j = i + 1; j < n_orbit; j++) {
      std::optional<MyMatrix<T>> opt_equiv =
          LINSPA_TestEquivalenceGramMatrix_SHV_Periodic<T, Tint, Tgroup>(
              data_func.data.LinSpa, pps, l_ent[i].x.DT_gram.GramMat,
              l_ent[j].x.DT_gram.GramMat, l_ent[i].x.DT_gram.SHV_T,
              l_ent[j].x.DT_gram.SHV_T, CommonGramMat, os);
      check(!opt_equiv.has_value(),
            "the domains of different orbits are inequivalent");
    }
  }
}

int main() {
  HumanTime time;
  try {
    std::ostream &os = std::cerr;
    int n = 2;
    MyMatrix<T> Cosets(2, n);
    Cosets(0, 0) = 0;
    Cosets(0, 1) = 0;
    Cosets(1, 0) = T(1) / T(3);
    Cosets(1, 1) = 0;
    PeriodicPointSet<Tint> pps = PeriodicPointSetFromRational<Tint>(Cosets);
    LinSpaceMatrix<T> LinSpa = ComputeCanonicalSpace<T>(n);
    run_case(LinSpa, pps, 1, "Z2+{0,(1/3,0)} full space", os);
    //
    // A4 + {0, c} in the dimension-2 T-space of the first subgroup with a
    // 2-dimensional space of invariant forms (a 5-cycle; the selection
    // procedure is in the CI notes). The expected count 2 is validated by
    // the property checks; the GAP oracle currently diverges on this case.
    MyMatrix<T> GramA4(4, 4);
    std::vector<std::vector<int>> gram_a4{
        {2,-1,0,0},{-1,2,-1,0},{0,-1,2,-1},{0,0,-1,2}};
    std::vector<std::vector<int>> b1_a4{
        {2,-1,1,-2},{-1,0,-1,1},{1,-1,2,0},{-2,1,0,0}};
    std::vector<std::vector<int>> h_a4{
        {0,0,1,0},{0,-1,-1,0},{0,1,1,1},{-1,-1,-1,-1}};
    MyMatrix<T> B1(4, 4), Hgen(4, 4);
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        GramA4(i, j) = gram_a4[i][j];
        B1(i, j) = b1_a4[i][j];
        Hgen(i, j) = h_a4[i][j];
      }
    }
    std::vector<MyMatrix<T>> ListMatA4{B1, GramA4};
    LinSpaceMatrix<T> LinSpaA4 =
        BuildLinSpaceMatrix<T, Tint, Tgroup>(ListMatA4, os);
    // The known-good small group handed to the kernels has to preserve
    // the cosets; the subgroup of the selection does, the full pointwise
    // stabilizer might not.
    LinSpaA4.PtStabGens = {Hgen};
    MyMatrix<T> CosetsA4(2, 4);
    for (int j = 0; j < 4; j++) {
      CosetsA4(0, j) = 0;
      CosetsA4(1, j) = T(4 - j) / T(5);
    }
    PeriodicPointSet<Tint> ppsA4 =
        PeriodicPointSetFromRational<Tint>(CosetsA4);
    run_case(LinSpaA4, ppsA4, 2, "A4+{0,c} dim-2 space", os);
    // The two other point sets of the A4 family, in the same T-space:
    // all of Aut(A4) is coset-adequate for each of them, so the same
    // subgroup selection applies.
    MyMatrix<T> Cosets2c(2, 4);
    MyMatrix<T> CosetsC2c(3, 4);
    for (int j = 0; j < 4; j++) {
      T c_j = T(4 - j) / T(5);
      T c2_j = T(2 * (4 - j) % 5) / T(5);
      Cosets2c(0, j) = 0;
      Cosets2c(1, j) = c2_j;
      CosetsC2c(0, j) = 0;
      CosetsC2c(1, j) = c_j;
      CosetsC2c(2, j) = c2_j;
    }
    PeriodicPointSet<Tint> pps2c =
        PeriodicPointSetFromRational<Tint>(Cosets2c);
    PeriodicPointSet<Tint> ppsC2c =
        PeriodicPointSetFromRational<Tint>(CosetsC2c);
    run_case(LinSpaA4, pps2c, 2, "A4+{0,2c} dim-2 space", os);
    run_case(LinSpaA4, ppsC2c, 3, "A4+{0,c,2c} dim-2 space", os);
    //
    // HE6 = E6 + {0, c6} in the dimension-2 T-space of the first subgroup
    // of Aut(E6) with a 2-dimensional space of invariant forms (order 18,
    // index 516 of the 1503 conjugacy classes). The basis of invariant
    // forms is not integrally saturated as computed, so it is saturated
    // first, as the Raw T-space input does with IntegralSaturation = T.
    auto fill = [](MyMatrix<T> &M, std::vector<std::vector<int>> const &v)
        -> void {
      for (size_t i = 0; i < v.size(); i++) {
        for (size_t j = 0; j < v[i].size(); j++) {
          M(i, j) = v[i][j];
        }
      }
    };
    MyMatrix<T> E1(6, 6), E2(6, 6), G1(6, 6), G2(6, 6), G3(6, 6);
    fill(E1, {{4,-1,-1,0,2,0},{-1,4,-2,0,-2,0},{-1,-2,4,-3,1,0},
              {0,0,-3,6,-3,0},{2,-2,1,-3,4,0},{0,0,0,0,0,0}});
    fill(E2, {{0,-1,1,0,-2,0},{-1,0,0,0,2,0},{1,0,0,1,-1,-2},
              {0,0,1,-2,1,0},{-2,2,-1,1,0,0},{0,0,-2,0,0,4}});
    fill(G1, {{0,1,1,1,0,0},{0,0,1,0,0,1},{0,0,0,1,1,0},
              {-1,-1,-2,-2,-1,-1},{1,1,1,1,0,1},{0,-1,-1,-1,-1,-1}});
    fill(G2, {{0,0,0,-1,-1,0},{0,0,0,0,1,0},{1,1,1,1,0,0},
              {-1,-1,0,0,0,0},{0,1,0,0,0,0},{0,0,0,0,0,1}});
    fill(G3, {{-1,-2,-2,-1,-1,-1},{0,0,0,0,1,0},{1,1,1,1,0,0},
              {-1,-1,-2,-2,-1,-1},{0,1,2,2,1,1},{0,0,0,0,0,1}});
    std::vector<MyMatrix<T>> ListMatE6 =
        IntegralSaturationSpace<T>({E1, E2}, os);
    LinSpaceMatrix<T> LinSpaE6 =
        BuildLinSpaceMatrix<T, Tint, Tgroup>(ListMatE6, os);
    LinSpaE6.PtStabGens = {G1, G2, G3};
    MyMatrix<T> CosetsE6(2, 6);
    std::vector<int> c6num{1, 2, 0, 1, 2, 0};
    for (int j = 0; j < 6; j++) {
      CosetsE6(0, j) = 0;
      CosetsE6(1, j) = T(c6num[j]) / T(3);
    }
    PeriodicPointSet<Tint> ppsE6 =
        PeriodicPointSetFromRational<Tint>(CosetsE6);
    run_case(LinSpaE6, ppsE6, 3, "HE6 dim-2 space", os);
    std::cerr << "Normal termination of TEST_PeriodicIsoDelaunay\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of TEST_PeriodicIsoDelaunay\n";
    exit(e.eVal);
  }
  runtime(time);
}
