// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#include "Iminimum.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

#ifdef DEBUG
#define DEBUG_SERIAL_COMPUTE_IMINIMUM
#endif

/*
  For an i-dimensional affine space X and a lattice L, the distance
  d(X, L) is the minimal distance between X and the points of L.
  The i-covering radius of L is the supremum of d(X, L) over all the
  i-dimensional affine spaces X.
  ---
  For a fixed direction space V = span(L_i) with L_i an i-dimensional
  sublattice of L, moving along V is free, so the maximum over the
  affine spaces X of direction V of d(X, L) is the covering radius of
  the projection of L onto the orthogonal complement of V.
  ---
  This program enumerates the orbits of the i-dimensional sublattices
  L_i of determinant at most max_det and computes for each the covering
  radius of the projected lattice, that is the profile g_i(d) of the
  maxima per determinant d <= max_det. The maximum of those covering
  radii is a lower bound for the i-minimum which is certified by
  LATT_SerialComputeIminimum.
 */

template <typename T> struct IminimumEntry {
  T TheDet;
  T CovSqr;
  size_t orbit_size;
  std::string str_latt;
};

template <typename Tgroup, typename T, typename Tint>
void compute_i_minimum(int const &i, std::string const &strMaxDet,
                       std::string const &FileGram,
                       std::string const &OutFormat, std::ostream &os_out) {
  std::ostream &os = std::cerr;
  MyMatrix<T> GramMat = ReadMatrixFile<T>(FileGram);
  T MaxDet = ParseScalar<T>(strMaxDet);
  int n = GramMat.rows();
  if (i < 1 || i >= n) {
    std::cerr << "We should have 1 <= i <= n-1 with n = " << n
              << " but i = " << i << "\n";
    throw TerminalException{1};
  }
  std::vector<SublatticeOrbit<Tint>> l_orbit =
      Rankin_k_level_orbits<T, Tint, Tgroup>(GramMat, i, MaxDet, os);
#ifdef DEBUG_SERIAL_COMPUTE_IMINIMUM
  os << "IMINIMUM: |l_orbit|=" << l_orbit.size() << "\n";
#endif
  T MaxCovSqr(0);
  std::vector<IminimumEntry<T>> l_entry;
  for (auto &eOrbit : l_orbit) {
    MyMatrix<Tint> const &eLatt = eOrbit.representative;
    MyMatrix<T> eLatt_T = UniversalMatrixConversion<T, Tint>(eLatt);
    MyMatrix<T> eGram = eLatt_T * GramMat * eLatt_T.transpose();
    T TheDet = DeterminantMat(eGram);
    MyMatrix<T> RedGram = GetProjectedGramMatrix(GramMat, eLatt);
    T CovSqr = ComputeCoveringRadiusSquared<T, Tint, Tgroup>(RedGram, os);
    if (CovSqr > MaxCovSqr) {
      MaxCovSqr = CovSqr;
    }
#ifdef DEBUG_SERIAL_COMPUTE_IMINIMUM
    os << "IMINIMUM: TheDet=" << TheDet << " CovSqr=" << CovSqr
       << " OrbitSize=" << eOrbit.orbit_size << "\n";
#endif
    IminimumEntry<T> entry{std::move(TheDet), std::move(CovSqr),
                           eOrbit.orbit_size, StringMatrixGAP_line(eLatt)};
    l_entry.push_back(std::move(entry));
  }
  if (OutFormat == "GAP") {
    os_out << "return rec(n:=" << n << ", i:=" << i << ", MaxDet:=" << MaxDet
           << ",\nListOrbit:=[";
    bool IsFirst = true;
    for (auto &entry : l_entry) {
      if (!IsFirst) {
        os_out << ",\n";
      }
      IsFirst = false;
      os_out << "rec(TheLattice:=" << entry.str_latt
             << ", TheDet:=" << entry.TheDet << ", CovSqr:=" << entry.CovSqr
             << ", OrbitSize:=" << entry.orbit_size << ")";
    }
    os_out << "],\nMaxCovSqr:=" << MaxCovSqr << ");\n";
    return;
  }
  if (OutFormat == "TXT") {
    os_out << "table of (determinant, squared covering radius):\n";
    for (auto &entry : l_entry) {
      os_out << entry.TheDet << " " << entry.CovSqr << "\n";
    }
    os_out << "MaxCovSqr = " << MaxCovSqr << "\n";
    return;
  }
  std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat
            << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 5 && argc != 7) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_SerialComputeIminimumGiD [arith] [i] [max_det] "
                   "[FileGram]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_SerialComputeIminimumGiD [arith] [i] [max_det] "
                   "[FileGram] [OutFormat] [OutFile]\n";
      std::cerr << "\n";
      std::cerr << "arith: gmp, gmp_boost, multi_boost\n";
      std::cerr << "i: The dimension of the sublattices being considered\n";
      std::cerr << "max_det: The maximum determinant of the i-dimensional\n";
      std::cerr << "    sublattices being enumerated\n";
      std::cerr << "FileGram: The Gram matrix on input\n";
      std::cerr << "OutFormat: GAP (default) or TXT\n";
      std::cerr << "OutFile: The file for the output, stderr and stdout\n";
      std::cerr << "    being special values\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string strI = argv[2];
    std::string strMaxDet = argv[3];
    std::string FileGram = argv[4];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 7) {
      OutFormat = argv[5];
      OutFile = argv[6];
    }
    int i = ParseScalar<int>(strI);
    //
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint_grp = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint_grp>;
    auto f = [&](std::ostream &os_out) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return compute_i_minimum<Tgroup, T, Tint>(i, strMaxDet, FileGram,
                                                  OutFormat, os_out);
      }
      if (arith == "gmp_boost") {
        using T = boost::multiprecision::mpq_rational;
        using Tint = boost::multiprecision::mpz_int;
        return compute_i_minimum<Tgroup, T, Tint>(i, strMaxDet, FileGram,
                                                  OutFormat, os_out);
      }
      if (arith == "multi_boost") {
        using T = boost::multiprecision::cpp_rational;
        using Tint = boost::multiprecision::cpp_int;
        return compute_i_minimum<Tgroup, T, Tint>(i, strMaxDet, FileGram,
                                                  OutFormat, os_out);
      }
      std::cerr << "Failed to find a matching entry for arith=" << arith
                << "\n";
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    std::cerr << "Normal termination of LATT_SerialComputeIminimumGiD\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_SerialComputeIminimumGiD\n";
    exit(e.eVal);
  }
  runtime(time);
}
