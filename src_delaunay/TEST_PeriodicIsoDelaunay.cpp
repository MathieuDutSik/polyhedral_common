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
    PolyHeuristicSerial<Tint_grp> AllArr =
        AllStandardHeuristicSerial<T, Tint_grp>(n + 1, os);
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
    std::cerr << "|orbits of periodic iso-Delaunay domains|=" << n_orbit
              << "\n";
    check(n_orbit > 0, "the enumeration found something");
    for (auto &eEnt : l_ent) {
      Tobj const &x = eEnt.x;
      MyMatrix<T> const &GramMat = x.DT_gram.GramMat;
      // The Gram matrix is strictly interior to its own domain.
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
      // The tessellation is the Delaunay tessellation of the point set.
      check(IsPeriodicDelaunayTesselation<T, Tint, Tgroup>(
                x.DT_gram.DT, GramMat, pps, os),
            "the tessellation of the domain is the Delaunay tessellation");
    }
    // The domains of different orbits are pairwise inequivalent.
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
    std::cerr << "Normal termination of TEST_PeriodicIsoDelaunay\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of TEST_PeriodicIsoDelaunay\n";
    exit(e.eVal);
  }
  runtime(time);
}
