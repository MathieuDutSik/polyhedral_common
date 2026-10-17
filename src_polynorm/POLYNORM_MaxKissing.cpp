// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Group.h"
#include "Permutation.h"
#include "PolyNorm_Kissing.h"
// clang-format on

template <typename T>
void ComputeKissing(std::string const &FileEXT, std::string const &OutFormat,
                    std::ostream &os_out) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  // The budget of the admissibility branch and bound of one configuration.
  int64_t adm_budget_per_leaf = 20000;
  MyMatrix<T> EXT = ReadMatrixFile<T>(FileEXT);
  KissingResult<T> res =
      ComputeMaxKissing<T, Tgroup>(EXT, adm_budget_per_leaf, std::cerr);
  if (OutFormat == "text") {
    os_out << "Maximal kissing number kiss=" << res.n_contact << "\n";
    if (res.n_undecided > 0) {
      os_out << "WARNING: " << res.n_undecided << " configurations exhausted "
             << "the admissibility budget, the kissing number is in ["
             << res.n_contact << ", " << res.n_contact_possible << "]\n";
    }
    if (res.witness) {
      KissingWitness<T> const &wit = *res.witness;
      os_out << "Witness lattice (class " << wit.i_class << ", det="
             << wit.det << "), basis:\n";
      WriteMatrix(os_out, wit.LattBasis);
      os_out << "Contact vectors (" << wit.ListContact.size() << "):\n";
      for (auto &x : wit.ListContact) {
        os_out << "  " << StringVector(x) << "\n";
      }
    }
    return;
  }
  if (OutFormat == "GAP") {
    os_out << "return rec(kiss:=" << res.n_contact
           << ", kiss_upper:=" << res.n_contact_possible
           << ", n_undecided:=" << res.n_undecided;
    if (res.witness) {
      KissingWitness<T> const &wit = *res.witness;
      os_out << ", det:=" << wit.det << ", basis:=";
      WriteMatrixGAP(os_out, wit.LattBasis);
      os_out << ", ListContact:=[";
      bool IsFirst = true;
      for (auto &x : wit.ListContact) {
        if (!IsFirst) {
          os_out << ",";
        }
        IsFirst = false;
        WriteVectorGAP(os_out, x);
      }
      os_out << "]";
    }
    os_out << ");\n";
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
      std::cerr << "POLYNORM_MaxKissing [arith] [FileEXT] [OutFormat] "
                   "[OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "POLYNORM_MaxKissing [arith] [FileEXT]\n";
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
      std::cerr << "Computes the maximal kissing number of the polytope P:\n";
      std::cerr << "the maximum over the lattices L with P + L a packing of\n";
      std::cerr << "the number of translates P + v, v in L, touching P, with\n";
      std::cerr << "a witness lattice and its contact vectors\n";
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
        return ComputeKissing<T>(FileEXT, OutFormat, os_out);
      }
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "flint") {
        using T = fmpq_class;
        return ComputeKissing<T>(FileEXT, OutFormat, os_out);
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
    std::cerr << "Normal termination of POLYNORM_MaxKissing\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLYNORM_MaxKissing\n";
    exit(e.eVal);
  }
  runtime(time);
}
