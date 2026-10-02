// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Group.h"
#include "Permutation.h"
#include "PolyNorm_LatticeOpt.h"
// clang-format on

template <typename T>
void PrintCertificate(PolyNormCertificate const &cert, std::ostream &os_out) {
  os_out << "  certified_local_optimum=" << cert.local_optimum << " ("
         << cert.message << ")\n";
  os_out << "  contacts/holes=" << cert.n_contact
         << " constraints=" << cert.n_constraint
         << " positive_multipliers=" << cert.n_positive
         << " dim_tangent=" << cert.dim_tangent << "\n";
}

template <typename T>
void ComputeCertify(std::string const &FileEXT, std::string const &FileLattice,
                    std::ostream &os_out) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  MyMatrix<T> EXT = ReadMatrixFile<T>(FileEXT);
  MyMatrix<T> W = ReadMatrixFile<T>(FileLattice);
  PolyNormData<T, Tint, Tgroup> data =
      PolyNorm_BuildData<T, Tint, Tgroup>(EXT, std::cerr);
  int n = data.n;
  if (W.rows() != n || W.cols() != n) {
    std::cerr << "The lattice basis should be a square matrix of size " << n
              << "\n";
    throw TerminalException{1};
  }
  if (DeterminantMat(W) == 0) {
    std::cerr << "The lattice basis is singular\n";
    throw TerminalException{1};
  }
  if (DeterminantMat(W) < 0) {
    for (int b = 0; b < n; b++) {
      W(0, b) = -W(0, b);
    }
  }
  T vol = rev_search::Kernel_VolumePolytope(data.EXT);
  os_out << "Dimension n=" << n << ", vol(P)=" << vol
         << ", det(L)=" << DeterminantMat(W) << "\n";
  // Packing.
  PolyNormLatticePackingEval<T, Tint> evalP =
      PolyNorm_EvaluatePacking<T, Tint>(data.LmatDiff, W, std::cerr);
  T pw(1);
  for (int i = 0; i < n; i++) {
    pw *= evalP.alpha;
  }
  T dens_pack = vol * pw / evalP.det;
  os_out << "Packing: alpha=" << evalP.alpha << " density=" << dens_pack
         << " ~ " << PolyNorm_ToDouble(dens_pack) << "\n";
  PolyNormCertificate certP = PolyNorm_CertifyPacking<T, Tint>(
      data.LmatDiff, W / evalP.alpha, 100000, std::cerr);
  PrintCertificate<T>(certP, os_out);
  // Covering.
  PolyNormLatticeCoveringEval<T, Tint> evalC =
      PolyNorm_EvaluateCovering<T, Tint, Tgroup>(data, W, {}, std::cerr);
  T pwc(1);
  for (int i = 0; i < n; i++) {
    pwc *= evalC.mu;
  }
  T dens_cov = vol * pwc / evalC.det;
  os_out << "Covering: mu=" << evalC.mu << " density=" << dens_cov << " ~ "
         << PolyNorm_ToDouble(dens_cov) << "\n";
  PolyNormCertificate certC = PolyNorm_CertifyCovering<T, Tint, Tgroup>(
      data, W / evalC.mu, 20000, std::cerr);
  PrintCertificate<T>(certC, os_out);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 4 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLYNORM_LatticeCertify [arith] [FileEXT] [FileLattice] "
                   "[OutFile]\n";
      std::cerr << "\n";
      std::cerr << "FileEXT     : the vertices of the polytope, one per row\n";
      std::cerr << "              in homogeneous coordinates (1, v)\n";
      std::cerr << "FileLattice : a basis of the lattice, one vector per row\n";
      std::cerr << "OutFile     : the output file, or stderr / stdout\n";
      std::cerr << "\n";
      std::cerr << "Computes the packing and covering densities of the\n";
      std::cerr << "polytope for the lattice and whether the lattice is an\n";
      std::cerr << "exactly certified strict local optimum of each problem.\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileEXT = argv[2];
    std::string FileLattice = argv[3];
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFile = argv[4];
    }
    auto prt = [&](std::ostream &os_out) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        return ComputeCertify<T>(FileEXT, FileLattice, os_out);
      }
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "flint") {
        using T = fmpq_class;
        return ComputeCertify<T>(FileEXT, FileLattice, os_out);
      }
#endif
      std::cerr << "Failed to find a matching entry for arith=" << arith
                << "\n";
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(OutFile, prt);
    std::cerr << "Normal termination of POLYNORM_LatticeCertify\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLYNORM_LatticeCertify\n";
    exit(e.eVal);
  }
  runtime(time);
}
