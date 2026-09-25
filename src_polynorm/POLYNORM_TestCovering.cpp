// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Group.h"
#include "Permutation.h"
#include "PolyNorm_Packing.h"
#include "PolyNorm_Covering.h"
// clang-format on

/*
  Cross-check of the packing and covering computations on one polytope:
  --- the branch and bound covering radius against the brute force
      enumeration of the paper (when the enumeration is small enough),
  --- both against expected values when they are given.
  Exits with an error when a mismatch is found.
 */
template <typename T>
void TestPolytope(std::string const &FileEXT,
                  std::optional<T> const &alpha_expected,
                  std::optional<T> const &mu_expected,
                  size_t max_brute_force, std::ostream &os) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  MyMatrix<T> EXT = ReadMatrixFile<T>(FileEXT);
  PolyNormData<T, Tint, Tgroup> data =
      PolyNorm_BuildData<T, Tint, Tgroup>(EXT, os);
  os << "TEST: n=" << data.n << " n_vert=" << EXT.rows()
     << " m=" << data.Lmat.rows() << " |GRP|=" << data.GRP.size() << "\n";
  //
  PolyNormPacking<T, Tint> pack =
      ComputePolyNormPacking<T, Tint>(data.LmatDiff, os);
  os << "TEST: alpha=" << pack.alpha << " with " << pack.ListContact.size()
     << " contact vectors\n";
  if (alpha_expected && *alpha_expected != pack.alpha) {
    std::cerr << "TEST: alpha=" << pack.alpha << " but expected "
              << *alpha_expected << "\n";
    throw TerminalException{1};
  }
  //
  PolyNormCovering<T, Tint> cov =
      ComputePolyNormCovering<T, Tint, Tgroup>(data, os);
  os << "TEST: mu=" << cov.mu << " p=" << StringVector(cov.p)
     << " n_node=" << cov.n_node << " n_lp=" << cov.n_lp
     << " n_skip_canonical=" << cov.n_skip_canonical << " mu0=" << cov.mu0
     << " |tight|=" << cov.ListTight.size() << "\n";
  if (mu_expected && *mu_expected != cov.mu) {
    std::cerr << "TEST: mu=" << cov.mu << " but expected " << *mu_expected
              << "\n";
    throw TerminalException{1};
  }
  // The brute force enumerates m C(|L| m - 1, n) systems for the lattice
  // points L of mu0 (P - P).
  size_t n_latt =
      PolyNorm_LatticePoints<T, Tint>(data.LmatDiff, cov.mu0, os).size();
  size_t n_pair = n_latt * data.Lmat.rows();
  size_t n_system = data.Lmat.rows();
  for (int r = 0; r < data.n; r++) {
    n_system = n_system * (n_pair - 1 - r) / (r + 1);
  }
  if (n_system <= max_brute_force) {
    std::pair<T, MyVector<T>> brute = PolyNorm_CoveringBruteForce(data, os);
    os << "TEST: brute force mu=" << brute.first
       << " p=" << StringVector(brute.second) << " (" << n_system
       << " systems)\n";
    if (brute.first != cov.mu) {
      std::cerr << "TEST: branch and bound mu=" << cov.mu
                << " but brute force mu=" << brute.first << "\n";
      throw TerminalException{1};
    }
  } else {
    os << "TEST: brute force skipped, " << n_system << " systems\n";
  }
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 3 && argc != 5 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLYNORM_TestCovering [arith] [FileEXT]\n";
      std::cerr << "    or\n";
      std::cerr << "POLYNORM_TestCovering [arith] [FileEXT] [alpha] [mu]\n";
      std::cerr << "    or\n";
      std::cerr << "POLYNORM_TestCovering [arith] [FileEXT] [alpha] [mu] "
                << "[max_brute_force]\n";
      std::cerr << "\n";
      std::cerr << "arith   : gmp";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << " or flint";
#endif
      std::cerr << "\n";
      std::cerr << "FileEXT : the vertices in homogeneous coordinates\n";
      std::cerr << "alpha   : the expected packing scalar, or none\n";
      std::cerr << "mu      : the expected covering radius, or none\n";
      std::cerr << "max_brute_force : the largest number of systems for\n";
      std::cerr << "          which the brute force check is run\n";
      std::cerr << "          (default 2000000)\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileEXT = argv[2];
    std::string s_alpha = "none";
    std::string s_mu = "none";
    size_t max_brute_force = 2000000;
    if (argc >= 5) {
      s_alpha = argv[3];
      s_mu = argv[4];
    }
    if (argc == 6) {
      max_brute_force = ParseScalar<size_t>(argv[5]);
    }
    auto run = [&]<typename T>() -> void {
      std::optional<T> alpha_expected;
      std::optional<T> mu_expected;
      if (s_alpha != "none") {
        alpha_expected = ParseScalar<T>(s_alpha);
      }
      if (s_mu != "none") {
        mu_expected = ParseScalar<T>(s_mu);
      }
      TestPolytope<T>(FileEXT, alpha_expected, mu_expected, max_brute_force,
                      std::cerr);
    };
    if (arith == "gmp") {
      run.template operator()<mpq_class>();
#ifdef ENABLE_FLINT_SUPPORT
    } else if (arith == "flint") {
      run.template operator()<fmpq_class>();
#endif
    } else {
      std::cerr << "Failed to find a matching entry for arith=" << arith
                << "\n";
      throw TerminalException{1};
    }
    std::cerr << "Normal termination of POLYNORM_TestCovering\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLYNORM_TestCovering\n";
    exit(e.eVal);
  }
  runtime(time);
}
