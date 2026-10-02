// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryQuadField.h"
#include "NumberTheoryRealField.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Group.h"
#include "Permutation.h"
#include "PolyNorm_BetkeHenk.h"
#include <iomanip>
// clang-format on

template <typename T>
void ComputeBetkeHenk(std::string const &FileEXT, int n_thread,
                      std::ostream &os_out) {
  MyMatrix<T> EXT = ReadMatrixFile<T>(FileEXT);
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tgroup = permutalib::Group<Telt, mpz_class>;
  BetkeHenkBody<T> body = BetkeHenk_BuildBody(EXT, std::cerr);
  BetkeHenk_ComputeOrbits<T, Tgroup>(body, std::cerr);
  double vol_P = UniversalScalarConversion<double, T>(body.vol_P);
  os_out << std::setprecision(12);
  os_out << "P: " << EXT.rows() << " vertices, vol(P)=" << body.vol_P << " ~ "
         << vol_P << "\n";
  os_out << "D = P - P: " << body.m << " facets in " << body.n_orbit
         << " orbits under its linear symmetry group\n";
  double det_bound = std::numeric_limits<double>::max();
  BetkeHenkSolver<T> solver(body, det_bound, std::cerr);
  for (auto &eCase : BetkeHenk_Cases()) {
    HumanTime time;
    size_t n_before = solver.ListCand.size();
    BetkeHenkStat stat_before = solver.stat;
    BetkeHenk_EnumerateParallel(solver, eCase, n_thread);
    os_out << "Case " << eCase.case_id << ": assignments="
           << solver.stat.n_assignment - stat_before.n_assignment
           << " feasible="
           << solver.stat.n_lp_feasible - stat_before.n_lp_feasible
           << " partial_lp="
           << solver.stat.n_partial_lp - stat_before.n_partial_lp
           << " rank_ok=" << solver.stat.n_rank_ok - stat_before.n_rank_ok
           << " real_roots=" << solver.stat.n_root - stat_before.n_root
           << " candidates=" << solver.ListCand.size() - n_before
           << " time=" << time << "\n";
  }
  os_out << "Time: linear programs " << solver.stat.t_lp / 1.0e6
         << " s, homotopies " << solver.stat.t_homotopy / 1.0e6
         << " s, admissibility " << solver.stat.t_admissible / 1.0e6 << " s\n";
  os_out << "Threads: " << n_thread << "\n";
  os_out << "Triple tests: " << solver.triples.memo.size() << " ("
         << solver.triples.n_lp << " linear programs, "
         << solver.triples.n_box_reject << " rejected by the boxes)\n";
  os_out << "Degenerate case IV choices (treated through the edges and "
            "vertices of F_6): "
         << solver.stat.n_case4_degenerate << "\n";
  os_out << "Homotopy paths stopped in the middle (t <= 0.9): "
         << solver.stat.n_path_stalled_mid
         << ", in the endgame of a singular endpoint: "
         << solver.stat.n_path_stalled_end << "\n";
  os_out << "Singular critical points: " << solver.stat.n_singular
         << " (components of dimension >= 3 not resolved: "
         << solver.stat.n_singular_unresolved << ")\n";
  std::vector<BetkeHenkCandidate> &L = solver.ListCand;
  std::sort(L.begin(), L.end(),
            [](BetkeHenkCandidate const &a, BetkeHenkCandidate const &b) {
              return a.det < b.det;
            });
  double tol = 1.0e-7;
  // The candidates are scanned by increasing determinant; the first
  // admissible one is the critical lattice. A candidate before it is not
  // admissible; that is conclusive except for:
  // --- singular critical points (positive dimensional components, on
  //     which det is constant and another point could be admissible). In
  //     the cases I and II, the test points of the point found are on the
  //     boundary of D, and by the remark after Lemma 4.9 of Betke and Henk
  //     the case can then be discarded (a critical lattice of that
  //     component would also be found in case III). In the cases III and
  //     IV the test set on the boundary implies admissibility (Lemma 2.3),
  //     so a non-admissible point there indicates a numerical problem.
  // --- candidates whose admissibility could not be computed.
  std::optional<size_t> i_best;
  size_t n_unresolved = 0;
  size_t n_discarded_remark = 0;
  std::vector<size_t> ListUnresolved;
  for (size_t i = 0; i < L.size(); i++) {
    if (L[i].min_gauge >= 1 - tol) {
      i_best = i;
      break;
    }
    if (L[i].min_gauge < 0 || (L[i].singular && L[i].case_id >= 3)) {
      n_unresolved++;
      ListUnresolved.push_back(i);
    } else if (L[i].singular) {
      n_discarded_remark++;
    }
  }
  if (!i_best) {
    os_out << "No admissible candidate found\n";
    return;
  }
  BetkeHenkCandidate const &best = L[*i_best];
  os_out << "Critical determinant Delta(P - P) ~ " << best.det << "\n";
  os_out << "Densest lattice packing density ~ " << vol_P / best.det << "\n";
  os_out << "Found in case " << best.case_id << " with facets";
  for (auto &j : best.facets) {
    os_out << " " << j;
  }
  os_out << (best.singular ? " (singular critical point)" : "") << "\n";
  os_out << "Basis (rows):\n";
  for (int a = 0; a < 3; a++) {
    os_out << "  " << best.W(a, 0) << " " << best.W(a, 1) << " " << best.W(a, 2)
           << "\n";
  }
  os_out << "Singular candidates of smaller determinant discarded by the "
            "remark after Lemma 4.9: "
         << n_discarded_remark << "\n";
  for (auto &i : ListUnresolved) {
    os_out << "  unresolved: det=" << L[i].det
           << " min_gauge=" << L[i].min_gauge << " case=" << L[i].case_id
           << " facets";
    for (auto &j : L[i].facets) {
      os_out << " " << j;
    }
    os_out << "\n";
  }
  if (solver.stat.n_path_stalled_mid > 0) {
    os_out << "WARNING: " << solver.stat.n_path_stalled_mid
           << " homotopy paths stopped in the middle, a critical point "
              "could be missing\n";
  } else if (n_unresolved > 0) {
    os_out << "WARNING: " << n_unresolved
           << " candidates of smaller determinant were not resolved\n";
  } else {
    os_out << "Every candidate of smaller determinant is resolved: the "
              "result is the critical determinant (up to the numerical "
              "tolerances)\n";
  }
  os_out << "Smallest candidates (det, admissibility min ||zW||_D, case):\n";
  for (size_t i = 0; i < L.size() && i < 8; i++) {
    os_out << "  " << L[i].det << " " << L[i].min_gauge << " " << L[i].case_id
           << (L[i].singular ? " singular" : "") << "\n";
  }
}

int main(int argc, char *argv[]) {
  maybe_install_tls_gmp_pool();
  HumanTime time;
  try {
    if (argc != 3 && argc != 4 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr
          << "POLYNORM_BetkeHenk [arith] [FileEXT] [OutFile] [n_thread]\n";
      std::cerr << "\n";
      std::cerr << "arith values:\n";
      std::cerr << "  gmp                : the rationals\n";
      std::cerr << "  Qsqrt2, Qsqrt5     : the fields Q(sqrt 2), Q(sqrt 5)\n";
      std::cerr << "  RealAlgebraic=File : the real algebraic field of File\n";
      std::cerr << "FileEXT : the vertices of the 3-polytope, one per row\n";
      std::cerr << "          in homogeneous coordinates (1, v)\n";
      std::cerr << "OutFile : the output file, or stderr / stdout\n";
      std::cerr << "n_thread: the number of threads (default 1)\n";
      std::cerr << "\n";
      std::cerr << "Computes the density of the densest lattice packing of\n";
      std::cerr << "the polytope by the algorithm of Betke and Henk.\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileEXT = argv[2];
    std::string OutFile = "stderr";
    if (argc >= 4) {
      OutFile = argv[3];
    }
    int n_thread = 1;
    if (argc == 5) {
      n_thread = ParseScalar<int>(argv[4]);
    }
    auto prt = [&](std::ostream &os_out) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        return ComputeBetkeHenk<T>(FileEXT, n_thread, os_out);
      }
      if (arith == "Qsqrt5") {
        using T = QuadField<Trat_quad_field, 5>;
        return ComputeBetkeHenk<T>(FileEXT, n_thread, os_out);
      }
      if (arith == "Qsqrt2") {
        using T = QuadField<Trat_quad_field, 2>;
        return ComputeBetkeHenk<T>(FileEXT, n_thread, os_out);
      }
      std::optional<std::string> opt_realalgebraic =
          get_postfix(arith, "RealAlgebraic=");
      if (opt_realalgebraic) {
        std::string FileAlgebraicField = *opt_realalgebraic;
        if (!FILE_IsExistingFile(FileAlgebraicField)) {
          std::cerr << "FileAlgebraicField=" << FileAlgebraicField
                    << " is missing\n";
          throw TerminalException{1};
        }
        HelperClassRealField<Trat_real_field> hcrf(FileAlgebraicField);
        int const idx_real_algebraic_field = 1;
        insert_helper_real_algebraic_field(idx_real_algebraic_field, hcrf);
        using T = RealField<idx_real_algebraic_field>;
        return ComputeBetkeHenk<T>(FileEXT, n_thread, os_out);
      }
      std::cerr << "Failed to find a matching entry for arith=" << arith
                << "\n";
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(OutFile, prt);
    std::cerr << "Normal termination of POLYNORM_BetkeHenk\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLYNORM_BetkeHenk\n";
    exit(e.eVal);
  }
  runtime(time);
}
