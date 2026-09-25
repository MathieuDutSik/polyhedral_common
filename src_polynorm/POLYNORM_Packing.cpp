// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Group.h"
#include "Permutation.h"
#include "PolyNorm_Packing.h"
// clang-format on

template <typename T>
void ComputePacking(std::string const &FileEXT, std::string const &OutFormat,
                    std::ostream &os_out) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  MyMatrix<T> EXT = ReadMatrixFile<T>(FileEXT);
  PolyNormData<T, Tint, Tgroup> data =
      PolyNorm_BuildData<T, Tint, Tgroup>(EXT, std::cerr);
  PolyNormPacking<T, Tint> res =
      ComputePolyNormPacking<T, Tint>(data.LmatDiff, std::cerr);
  if (OutFormat == "text") {
    os_out << "Dimension n=" << data.n << " with " << EXT.rows()
           << " vertices and " << data.Lmat.rows() << " facets\n";
    os_out << "Difference body with " << data.LmatDiff.rows()
           << " facets\n";
    os_out << "Symmetry group order |{A in GL_n(Z) : PA = P}|="
           << data.GRP.size() << "\n";
    os_out << "Packing scalar alpha=" << res.alpha << "\n";
    os_out << "Contact vectors (" << res.ListContact.size() << "):\n";
    for (auto &z : res.ListContact) {
      os_out << "  " << StringVector(z) << "\n";
    }
    return;
  }
  if (OutFormat == "GAP") {
    os_out << "return rec(alpha:=" << res.alpha
           << ", GRPorder:=" << data.GRP.size() << ", ListContact:=[";
    bool IsFirst = true;
    for (auto &z : res.ListContact) {
      if (!IsFirst) {
        os_out << ",";
      }
      IsFirst = false;
      WriteVectorGAP(os_out, z);
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
    if (argc != 2 && argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLYNORM_Packing [FileEXT] [OutFormat] [OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "POLYNORM_Packing [FileEXT]\n";
      std::cerr << "\n";
      std::cerr << "FileEXT   : the vertices of the polytope, one per row\n";
      std::cerr << "            in homogeneous coordinates (1, v)\n";
      std::cerr << "OutFormat : text (default) or GAP\n";
      std::cerr << "OutFile   : the output file, or stderr / stdout\n";
      std::cerr << "\n";
      std::cerr << "Computes the largest alpha such that the translates\n";
      std::cerr << "alpha P + z, z in Z^n, have disjoint interiors, with the\n";
      std::cerr << "lattice vectors z for which alpha P and alpha P + z touch\n";
      return -1;
    }
    std::string FileEXT = argv[1];
    std::string OutFormat = "text";
    std::string OutFile = "stderr";
    if (argc == 4) {
      OutFormat = argv[2];
      OutFile = argv[3];
    }
    auto prt = [&](std::ostream &os_out) -> void {
      using T = mpq_class;
      return ComputePacking<T>(FileEXT, OutFormat, os_out);
    };
    FILE_PrintStderrStdoutFile(OutFile, prt);
    std::cerr << "Normal termination of POLYNORM_Packing\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLYNORM_Packing\n";
    exit(e.eVal);
  }
  runtime(time);
}
