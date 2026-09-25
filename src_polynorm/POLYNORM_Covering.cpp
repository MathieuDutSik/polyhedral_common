// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Group.h"
#include "Permutation.h"
#include "PolyNorm_Covering.h"
// clang-format on

template <typename T>
void ComputeCovering(std::string const &FileEXT, std::string const &OutFormat,
                     std::ostream &os_out) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  MyMatrix<T> EXT = ReadMatrixFile<T>(FileEXT);
  PolyNormData<T, Tint, Tgroup> data =
      PolyNorm_BuildData<T, Tint, Tgroup>(EXT, std::cerr);
  PolyNormCovering<T, Tint> res =
      ComputePolyNormCovering<T, Tint, Tgroup>(data, std::cerr);
  // The last covered point in the coordinates of the input.
  MyVector<T> p_input = res.p + data.center;
  if (OutFormat == "text") {
    os_out << "Dimension n=" << data.n << " with " << EXT.rows()
           << " vertices and " << data.Lmat.rows() << " facets\n";
    os_out << "Symmetry group order |{A in GL_n(Z) : PA = P}|="
           << data.GRP.size() << "\n";
    os_out << "Upper bound mu0=" << res.mu0 << "\n";
    os_out << "Search: " << res.n_node << " nodes, " << res.n_lp
           << " linear programs, " << res.n_skip_canonical
           << " nodes skipped by symmetry\n";
    os_out << "Covering radius mu=" << res.mu << "\n";
    os_out << "Last covered point p=" << StringVector(p_input) << "\n";
    os_out << "Lattice points on the boundary of p - mu P, with the facet "
           << "of P they touch (" << res.ListTight.size() << "):\n";
    for (auto &ePair : res.ListTight) {
      os_out << "  z=" << StringVector(ePair.z) << " facet " << ePair.i
             << "\n";
    }
    return;
  }
  if (OutFormat == "GAP") {
    os_out << "return rec(mu:=" << res.mu << ", GRPorder:=" << data.GRP.size()
           << ", p:=";
    WriteVectorGAP(os_out, p_input);
    os_out << ", ListTight:=[";
    bool IsFirst = true;
    for (auto &ePair : res.ListTight) {
      if (!IsFirst) {
        os_out << ",";
      }
      IsFirst = false;
      os_out << "rec(z:=";
      WriteVectorGAP(os_out, ePair.z);
      os_out << ", facet:=" << (ePair.i + 1) << ")";
    }
    os_out << "]);\n";
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
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLYNORM_Covering [arith] [FileEXT] [OutFormat] [OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "POLYNORM_Covering [arith] [FileEXT]\n";
      std::cerr << "\n";
      std::cerr << "arith values:\n";
      std::cerr << "  gmp   : mpq_class / mpz_class\n";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "  flint : fmpq_class / fmpz_class\n";
#endif
      std::cerr << "FileEXT   : the vertices of the polytope, one per row\n";
      std::cerr << "            in homogeneous coordinates (1, v)\n";
      std::cerr << "OutFormat : text (default) or GAP\n";
      std::cerr << "OutFile   : the output file, or stderr / stdout\n";
      std::cerr << "\n";
      std::cerr << "Computes the smallest mu such that the translates\n";
      std::cerr << "mu P + z, z in Z^n, cover R^n, with a last covered point\n";
      std::cerr << "and the lattice points pinning its empty translate\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileEXT = argv[2];
    std::string OutFormat = "text";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    auto prt = [&](std::ostream &os_out) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        return ComputeCovering<T>(FileEXT, OutFormat, os_out);
      }
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "flint") {
        using T = fmpq_class;
        return ComputeCovering<T>(FileEXT, OutFormat, os_out);
      }
#endif
      std::cerr << "Failed to find a matching entry for arith=" << arith
                << "\n";
      std::cerr << "Available possibilities: gmp";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << ", flint";
#else
      std::cerr << " (build with ENABLE_FLINT_SUPPORT=1 for flint)";
#endif
      std::cerr << "\n";
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(OutFile, prt);
    std::cerr << "Normal termination of POLYNORM_Covering\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLYNORM_Covering\n";
    exit(e.eVal);
  }
  runtime(time);
}
