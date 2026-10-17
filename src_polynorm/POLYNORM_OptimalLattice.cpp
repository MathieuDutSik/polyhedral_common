// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Group.h"
#include "Permutation.h"
#include "PolyNorm_LatticeOpt.h"
#include "PolyNorm_Packing2D.h"
#include <iomanip>
// clang-format on

template <typename T>
void PrintLatticeResult(PolyNormLatticeResult<T> const &res,
                        std::ostream &os_out) {
  os_out << "  density=" << res.density_d;
  if (res.polished) {
    os_out << " exact=" << res.density;
  }
  if (res.density_refined) {
    os_out << " refined_floating_point=" << std::setprecision(15)
           << *res.density_refined << std::setprecision(6)
           << " (critical point not rational)";
  }
  os_out << "\n";
  if (res.polished && res.density == 1) {
    // A density 1 lattice is a tiling: optimal for both problems.
    os_out << "  global optimum: density 1, a lattice tiling\n";
  }
  os_out << "  polished=" << res.polished
         << " certified_local_optimum=" << res.cert.local_optimum << " ("
         << res.cert.message << ")\n";
  os_out << "  contacts/holes=" << res.cert.n_contact
         << " constraints=" << res.cert.n_constraint
         << " positive_multipliers=" << res.cert.n_positive
         << " dim_tangent=" << res.cert.dim_tangent << "\n";
  os_out << "  basis (rows, radius normalized to 1):\n";
  for (int k = 0; k < res.W.rows(); k++) {
    os_out << "   " << StringVector(GetMatrixRow(res.W, k)) << "\n";
  }
}

template <typename T>
void ComputeOptimalLattice(std::string const &FileEXT,
                           std::string const &problem, int n_start,
                           uint64_t seed, std::string const &FileStart,
                           std::ostream &os_out) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  if (problem != "packing" && problem != "covering") {
    std::cerr << "The problem should be packing or covering, not " << problem
              << "\n";
    throw TerminalException{1};
  }
  MyMatrix<T> EXT = ReadMatrixFile<T>(FileEXT);
  PolyNormData<T, Tint, Tgroup> data =
      PolyNorm_BuildData<T, Tint, Tgroup>(EXT, std::cerr);
  T vol = DirectVolumePolytope(data.EXT, TriangulationProgram::lrs, std::cerr);
  os_out << "Dimension n=" << data.n << " with " << data.EXT.rows()
         << " vertices, " << data.Lmat.rows() << " facets, vol(P)=" << vol
         << "\n";
  os_out << "Difference body with " << data.LmatDiff.rows() << " facets\n";
  if (problem == "packing" && data.n == 2) {
    PolyNormPacking2D<T> res2 =
        PolyNorm_ExactPacking2D<T, Tint>(data.LmatDiff, std::cerr);
    os_out << "Exact 2-dimensional algorithm:\n";
    os_out << "  critical determinant of P - P = " << res2.Delta << "\n";
    T dens = vol / res2.Delta;
    os_out << "  lattice packing density = " << dens << " ~ "
           << PolyNorm_ToDouble(dens) << "\n";
    os_out << "  basis (rows):\n";
    for (int k = 0; k < 2; k++) {
      os_out << "   " << StringVector(GetMatrixRow(res2.W, k)) << "\n";
    }
  }
  PolyNormLatticeOptParam param;
  param.problem = problem;
  param.n_start = n_start;
  param.seed = seed;
  param.max_iter = 400;
  param.r_init = 0.2;
  param.r_min = 1.0e-9;
  param.near_ratio = 0.9;
  param.max_system = 20000;
  param.max_branch = 100000;
  param.r_pattern_min = 1.0e-6;
  param.max_round = 20;
  std::vector<MyMatrix<T>> ListStart;
  if (FileStart != "none") {
    ListStart.push_back(ReadMatrixFile<T>(FileStart));
  }
  std::vector<PolyNormLatticeResult<T>> ListRes =
      PolyNorm_LatticeOptimize(data, param, ListStart, std::cerr);
  bool is_packing = problem == "packing";
  std::optional<size_t> i_best;
  for (size_t i = 0; i < ListRes.size(); i++) {
    os_out << "Start " << i << ":\n";
    PrintLatticeResult(ListRes[i], os_out);
    if (!i_best) {
      i_best = i;
    } else {
      bool better = is_packing ? ListRes[i].density > ListRes[*i_best].density
                               : ListRes[i].density < ListRes[*i_best].density;
      if (better) {
        i_best = i;
      }
    }
  }
  if (i_best) {
    os_out << "Best lattice " << problem << ":\n";
    PrintLatticeResult(ListRes[*i_best], os_out);
  }
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 6 && argc != 7 && argc != 8) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLYNORM_OptimalLattice [arith] [FileEXT] [problem] "
                   "[n_start] [seed] [FileStart] [OutFile]\n";
      std::cerr << "\n";
      std::cerr << "arith values:\n";
      std::cerr << "  gmp   : mpq_class / mpz_class\n";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "  flint : fmpq_class / fmpz_class\n";
#endif
      std::cerr << "FileEXT   : the vertices of the polytope, one per row\n";
      std::cerr << "            in homogeneous coordinates (1, v)\n";
      std::cerr << "problem   : packing or covering\n";
      std::cerr << "n_start   : the number of random starting lattices\n";
      std::cerr << "seed      : the seed of the random starts\n";
      std::cerr << "FileStart : a starting basis (rows) or none (default)\n";
      std::cerr << "OutFile   : the output file, or stderr / stdout\n";
      std::cerr << "\n";
      std::cerr << "Searches the lattice of densest packing or thinnest\n";
      std::cerr << "covering of R^n by translates of the polytope, by local\n";
      std::cerr << "optimization from several starts, and certifies the\n";
      std::cerr << "local optimality of the results exactly. In dimension\n";
      std::cerr << "2 the packing problem is also solved exactly.\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileEXT = argv[2];
    std::string problem = argv[3];
    int n_start = ParseScalar<int>(argv[4]);
    uint64_t seed = ParseScalar<uint64_t>(argv[5]);
    std::string FileStart = "none";
    std::string OutFile = "stderr";
    if (argc >= 7) {
      FileStart = argv[6];
    }
    if (argc == 8) {
      OutFile = argv[7];
    }
    auto prt = [&](std::ostream &os_out) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        return ComputeOptimalLattice<T>(FileEXT, problem, n_start, seed,
                                        FileStart, os_out);
      }
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "flint") {
        using T = fmpq_class;
        return ComputeOptimalLattice<T>(FileEXT, problem, n_start, seed,
                                        FileStart, os_out);
      }
#endif
      std::cerr << "Failed to find a matching entry for arith=" << arith
                << "\n";
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(OutFile, prt);
    std::cerr << "Normal termination of POLYNORM_OptimalLattice\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLYNORM_OptimalLattice\n";
    exit(e.eVal);
  }
  runtime(time);
}
