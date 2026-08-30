// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#include "Iminimum.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

/*
  Computes the i-minimum mu_i(L)^2 of the lattice with certificate:
  the orbits of the i-dimensional sublattices are enumerated up to a
  determinant bound, the covering radii of the projections give the
  candidate value, and the certification algorithm BOUND of
  Iminimum.h proves that no direction exceeds it, raising the
  enumeration bound as needed. See Iminimum.h for the description of
  the algorithm; LATT_SerialComputeIminimumGiD computes the profile
  per determinant for a determinant bound given by the user.
 */

template <typename Tgroup, typename T, typename Tint>
void compute_i_minimum_certified(int const &i, std::string const &FileGram,
                                 std::string const &OutFormat,
                                 std::ostream &os_out) {
  std::ostream &os = std::cerr;
  MyMatrix<T> GramMat = ReadMatrixFile<T>(FileGram);
  int n = GramMat.rows();
  if (i < 1 || i >= n) {
    std::cerr << "We should have 1 <= i <= n-1 with n = " << n
              << " but i = " << i << "\n";
    throw TerminalException{1};
  }
  IminimumResult<T, Tint> result =
      ComputeIminimumCertified<T, Tint, Tgroup>(GramMat, i, os);
  if (OutFormat == "GAP") {
    os_out << "return rec(n:=" << n << ", i:=" << result.i;
    if (result.certified) {
      os_out << ", certified:=true";
    } else {
      os_out << ", certified:=false, failure:=\"" << result.failure_reason
             << "\"";
    }
    os_out << ", IminimumSqr:=" << result.IminimumSqr
           << ", MaxDetEnumerated:=" << result.MaxDetEnumerated
           << ",\nExceptionalDets:=[";
    bool IsFirst = true;
    for (auto &det : result.exceptional_dets) {
      if (!IsFirst) {
        os_out << ",";
      }
      IsFirst = false;
      os_out << det;
    }
    os_out << "],\nListOrbit:=[";
    IsFirst = true;
    for (auto &entry : result.l_orbit) {
      if (!IsFirst) {
        os_out << ",\n";
      }
      IsFirst = false;
      os_out << "rec(TheLattice:=" << StringMatrixGAP_line(entry.representative)
             << ", TheDet:=" << entry.det << ", CovSqr:=" << entry.covsqr
             << ", OrbitSize:=" << entry.orbit_size << ")";
    }
    os_out << "]);\n";
    return;
  }
  if (OutFormat == "TXT") {
    if (result.certified) {
      os_out << "mu_" << result.i << "^2 = " << result.IminimumSqr
             << " (certified)\n";
    } else {
      os_out << "mu_" << result.i << "^2 >= " << result.IminimumSqr
             << " (NOT certified: " << result.failure_reason << ")\n";
    }
    os_out << "enumerated up to determinant " << result.MaxDetEnumerated
           << "\n";
    os_out << "table of (determinant, squared covering radius, orbit "
              "size):\n";
    for (auto &entry : result.l_orbit) {
      os_out << entry.det << " " << entry.covsqr << " " << entry.orbit_size
             << "\n";
    }
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
    if (argc != 4 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_SerialComputeIminimum [arith] [i] [FileGram]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_SerialComputeIminimum [arith] [i] [FileGram] "
                   "[OutFormat] [OutFile]\n";
      std::cerr << "\n";
      std::cerr << "arith: gmp, gmp_boost, multi_boost\n";
      std::cerr << "i: The dimension of the direction spaces\n";
      std::cerr << "FileGram: The Gram matrix on input\n";
      std::cerr << "OutFormat: GAP (default) or TXT\n";
      std::cerr << "OutFile: The file for the output, stderr and stdout\n";
      std::cerr << "    being special values\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string strI = argv[2];
    std::string FileGram = argv[3];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      OutFile = argv[5];
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
        return compute_i_minimum_certified<Tgroup, T, Tint>(i, FileGram,
                                                            OutFormat, os_out);
      }
      if (arith == "gmp_boost") {
        using T = boost::multiprecision::mpq_rational;
        using Tint = boost::multiprecision::mpz_int;
        return compute_i_minimum_certified<Tgroup, T, Tint>(i, FileGram,
                                                            OutFormat, os_out);
      }
      if (arith == "multi_boost") {
        using T = boost::multiprecision::cpp_rational;
        using Tint = boost::multiprecision::cpp_int;
        return compute_i_minimum_certified<Tgroup, T, Tint>(i, FileGram,
                                                            OutFormat, os_out);
      }
      std::cerr << "Failed to find a matching entry for arith=" << arith
                << "\n";
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    std::cerr << "Normal termination of LATT_SerialComputeIminimum\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_SerialComputeIminimum\n";
    exit(e.eVal);
  }
  runtime(time);
}
