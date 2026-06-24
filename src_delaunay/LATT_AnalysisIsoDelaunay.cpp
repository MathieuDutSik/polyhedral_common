// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
#include "IsoDelaunayDomains.h"
#include "Permutation.h"
#include "Group.h"
#include <boost/archive/text_iarchive.hpp>
// clang-format on

/*
  Reads an iso-Delaunay domain (boost text_oarchive of IsoDelaunayDomain<T,
  Tint, Tgroup> as written by LATT_SerialComputeDelaunay's
  FileIsoDelaunayDomain QUERIES option) and reports:
    * |ListIneq|   — full set of defining inequalities
                     (ComputeDefiningIneqIsoDelaunayDomain)
    * |ListIrred|  — irredundant inequalities (get_non_redundant_indices)
    * For each extreme ray r of the cone, the rank of
        sum_u r_u * LinSpa.ListMat[u]
      and a tally of how many rays fall in each rank bucket.
 */

template <typename T, typename Tint>
void process(std::string const &FileIso, std::ostream &os_out) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  std::ostream &os = std::cerr;
  //
  IsoDelaunayDomain<T, Tint, Tgroup> x;
  {
    std::ifstream ifs(FileIso);
    boost::archive::text_iarchive ia(ifs);
    ia >> x;
  }
  int n = x.GramMat.rows();
  int sym_dim = (n * (n + 1)) / 2;
  LinSpaceMatrix<T> LinSpa = ComputeCanonicalSpace<T>(n);
  os << "ANA: n=" << n << " sym_dim=" << sym_dim
     << " |DT|=" << x.DT.l_dels.size() << "\n";
  //
  HumanTime t1;
  std::vector<FullAdjInfo<T>> ListIneq =
      ComputeDefiningIneqIsoDelaunayDomain<T, Tgroup>(
          x.DT, LinSpa.ListLineMat, os);
  int n_ineq = ListIneq.size();
  os << "ANA: ComputeDefiningIneqIsoDelaunayDomain done n_ineq=" << n_ineq
     << " |elapsed|=" << t1 << "\n";
  //
  // Direct dual description on the raw FAC. lrs (the default backend
  // chosen by get_dual_desc_method) handles redundant rows in the input
  // itself; doing it this way avoids the Clarkson shooting loop in
  // get_non_redundant_indices, which is O(n_ineq^2) LPs in rational
  // arithmetic and is the bottleneck on dim-6 primitive L-types.
  HumanTime t2;
  MyMatrix<T> FAC = GetFACineq(ListIneq);
  MyMatrix<T> EXT = DirectDualDescription_mat(FAC, os);
  int n_ray = EXT.rows();
  os << "ANA: DirectDualDescription_mat (FAC->EXT) done n_ray=" << n_ray
     << " |elapsed|=" << t2 << "\n";
  //
  // Irredundant facet count: come back EXT -> FAC. Each row of the result
  // is a defining inequality (a facet); redundant rows of the original FAC
  // are dropped.
  HumanTime t3;
  MyMatrix<T> FACred = DirectDualDescription_mat(EXT, os);
  int n_irred = FACred.rows();
  os << "ANA: DirectDualDescription_mat (EXT->FAC) done n_irred=" << n_irred
     << " |elapsed|=" << t3 << "\n";
  //
  std::map<int, int> rank_tally;
  int dimSpace = LinSpa.ListMat.size();
  for (int i_row = 0; i_row < n_ray; i_row++) {
    MyMatrix<T> RayMat = ZeroMatrix<T>(n, n);
    for (int u = 0; u < dimSpace; u++) {
      RayMat += EXT(i_row, u) * LinSpa.ListMat[u];
    }
    int rk = RankMat(RayMat);
    rank_tally[rk]++;
  }
  //
  os_out << "n=" << n << " sym_dim=" << sym_dim << "\n";
  os_out << "n_ineq=" << n_ineq << "\n";
  os_out << "n_irred=" << n_irred << "\n";
  os_out << "n_ray=" << n_ray << "\n";
  os_out << "rank_tally=[";
  bool first = true;
  for (auto &kv : rank_tally) {
    if (!first) {
      os_out << ", ";
    }
    first = false;
    os_out << "(rank=" << kv.first << ", count=" << kv.second << ")";
  }
  os_out << "]\n";
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3 && argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_AnalysisIsoDelaunay [arith] [file.iso] [OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_AnalysisIsoDelaunay [arith] [file.iso]\n";
      std::cerr << "arith   : gmp | gmp_boost | multi_boost | safe\n";
      std::cerr << "file.iso: a boost text-archive file produced by "
                   "LATT_SerialComputeDelaunay's FileIsoDelaunayDomain "
                   "QUERIES option\n";
      std::cerr << "OutFile : path to write the analysis (defaults to stderr)\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileIso = argv[2];
    std::string OutFile = "stderr";
    if (argc == 4) {
      OutFile = argv[3];
    }
    auto f = [&](std::ostream &os_out) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return process<T, Tint>(FileIso, os_out);
      }
      if (arith == "gmp_boost") {
        using T = boost::multiprecision::mpq_rational;
        using Tint = boost::multiprecision::mpz_int;
        return process<T, Tint>(FileIso, os_out);
      }
      if (arith == "multi_boost") {
        using T = boost::multiprecision::cpp_rational;
        using Tint = boost::multiprecision::cpp_int;
        return process<T, Tint>(FileIso, os_out);
      }
      if (arith == "safe") {
        using T = Rational<SafeInt64>;
        using Tint = SafeInt64;
        return process<T, Tint>(FileIso, os_out);
      }
      std::cerr << "LATT_AnalysisIsoDelaunay: Failed to find a matching type "
                   "for arith="
                << arith << "\n";
      std::cerr << "Available types: gmp, gmp_boost, multi_boost, safe\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(OutFile, f);
    std::cerr << "Normal termination of LATT_AnalysisIsoDelaunay\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_AnalysisIsoDelaunay\n";
    exit(e.eVal);
  }
  runtime(time);
}
