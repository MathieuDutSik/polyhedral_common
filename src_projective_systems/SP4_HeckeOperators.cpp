// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "sp4_hecke.h"
#include <fstream>
#include <string>
#include <vector>
// clang-format on

/*
  The Hecke operators of Sp_4(Z) on the homology of the subgroups
  Gamma0_A(N) (the Siegel parabolic subgroup modulo N) and Gamma0_B(N)
  (the Klingen parabolic subgroup), as computed by the driver
  HeckeOperatorSp4/TheBigComputation.g with the GAP code of
  ProjectiveSystem.g:
  - the cell complex of Sp_4(Z) of the projective systems,
  - the chain map of the double coset of the Hecke matrix on the levels
    of the single vectors up to the level eDimMax of GAP,
  - the homology of the well rounded cells (the cells at infinity
    eliminated) of the quotient by each subgroup, with the Hecke matrix.
  The Hecke matrices are those of SYMPL_GetHeckeOperatorCandidate_A/B/C:
  A = diag(p, p, 1/p, 1/p), B = diag(p, p, 1, 1), C = diag(p, p^2, p, 1),
  multiplied by a scalar to be integral.
 */

namespace {

MyMatrix<sp4::T> get_hecke_matrix(std::string const &choice, int const &p) {
  std::vector<sp4::T> l_diag;
  if (choice == "A") {
    l_diag = {sp4::T(p * p), sp4::T(p * p), sp4::T(1), sp4::T(1)};
  } else if (choice == "B") {
    l_diag = {sp4::T(p), sp4::T(p), sp4::T(1), sp4::T(1)};
  } else if (choice == "C") {
    l_diag = {sp4::T(p), sp4::T(p * p), sp4::T(p), sp4::T(1)};
  } else {
    std::cerr << "SP4_HeckeOperators: the choice should be A, B or C\n";
    throw TerminalException{1};
  }
  MyMatrix<sp4::T> x = ZeroMatrix<sp4::T>(4, 4);
  for (int i = 0; i < 4; i++) {
    x(i, i) = l_diag[i];
  }
  return x;
}

std::vector<int> parse_list(std::string const &str) {
  std::vector<int> l;
  for (auto &s : STRING_Split(str, ",")) {
    l.push_back(ParseScalar<int>(s));
  }
  return l;
}

}  // namespace

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 8) {
      std::cerr << "SP4_HeckeOperators [choice] [p] [eDimMax] [subgroups] [levels] [prime] [prefix]\n";
      std::cerr << "  choice: A, B or C (the Hecke matrix), p: its prime\n";
      std::cerr << "  eDimMax: the levels of GAP up to which the chain map is computed (4 in the driver)\n";
      std::cerr << "  subgroups: comma separated list of Gamma0_A (Siegel), Gamma0_B (Klingen)\n";
      std::cerr << "  levels: comma separated list of the levels N of the subgroups\n";
      std::cerr << "  prime: 0 for the rational homology, otherwise the homology modulo it\n";
      std::cerr << "  prefix: the GAP results are written to prefix_[subgroup]_[N].g\n";
      return -1;
    }
    std::string choice = argv[1];
    int p = ParseScalar<int>(argv[2]);
    int eDimMax = ParseScalar<int>(argv[3]);
    std::vector<std::string> l_sub = STRING_Split(argv[4], ",");
    std::vector<int> l_level = parse_list(argv[5]);
    int64_t prime = ParseScalar<int64_t>(argv[6]);
    std::string prefix = argv[7];
    std::ostream &os = std::cerr;
    sp4::Sp4CellComplex cx = sp4::compute_sp4_cell_complex(os);
    int n_levels = cx.levels.size();
    std::cerr << "Sp4 complex: |top cells|=" << cx.l_top.size() << " |generators|=" << cx.l_gens_G.size() << "\n";
    for (int index = 0; index < n_levels; index++) {
      std::cerr << "  index=" << index << " (GAP level " << n_levels - index << ") orbits:";
      for (auto &orbit : cx.levels[index]) {
        int rnk = RankMat(orbit.EXT);
        bool orientable = true;
        for (auto &g : orbit.l_stab_gens) {
          if (sp4::orbit_orientation(orbit, g) == -1) {
            orientable = false;
          }
        }
        std::cerr << " [|EXT|=" << orbit.EXT.rows() << " rank=" << rnk << " |Stab on EXT|=" << orbit.perm_grp.size()
                  << " orientable=" << orientable << " |faces|=" << orbit.l_bound.size() << "]";
      }
      std::cerr << "\n";
    }
    MyMatrix<sp4::T> x = get_hecke_matrix(choice, p);
    int index_min = n_levels - eDimMax;
    HeckeChainMap<sp4::T, sp4::Tint> hcm = sp4::compute_sp4_hecke_chain_map(cx, x, index_min, os);
    std::cerr << "Chain map computed, |Hecke cosets|=" << hcm.cosets.size() << "\n";
    for (auto &sub : l_sub) {
      std::string type;
      if (sub == "Gamma0_A") {
        type = "Gamma0Siegel";
      } else if (sub == "Gamma0_B") {
        // In Sp_4(Z) the last row congruent to (0, 0, 0, *) implies the
        // other conditions of Gamma0_B of the driver.
        type = "Gamma0";
      } else {
        std::cerr << "SP4_HeckeOperators: the subgroups are Gamma0_A and Gamma0_B, not " << sub << "\n";
        throw TerminalException{1};
      }
      for (auto &N : l_level) {
        FiniteIndexSubgroup<sp4::Tint> gamma{};
        gamma.type = type;
        gamma.N = sp4::Tint(N);
        gamma.n = 4;
        HeckeHomologyResult<sp4::T> result =
            compute_hecke_homology(sp4::Sp4ComplexView{cx}, hcm, gamma, true, false, false, prime, os);
        std::string file = prefix + "_" + sub + "_" + std::to_string(N) + ".g";
        std::ofstream os_out(file);
        os_out << "return ";
        WriteHeckeHomologyResultGAP(os_out, result);
        os_out << ";\n";
        std::cerr << sub << " N=" << N << " index=" << result.n_coset_gamma << " |Hecke cosets of Gamma|="
                  << result.n_hecke_coset << " dimensions:";
        for (auto &level : result.l_level) {
          std::cerr << " H(index " << level.index << ")=" << level.harmonic.rows();
        }
        std::cerr << " -> " << file << "\n";
      }
    }
    std::cerr << "Normal termination of SP4_HeckeOperators\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of SP4_HeckeOperators\n";
    exit(e.eVal);
  } catch (HeckeFailure const &e) {
    std::cerr << "SP4_HeckeOperators: clean failure: " << e.message << "\n";
    exit(1);
  }
  runtime(time);
}
