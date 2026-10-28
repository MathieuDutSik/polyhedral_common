// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "projective_hecke.h"
#include <fstream>
#include <string>
#include <vector>
// clang-format on

/*
  The Hecke operators on the homology of finite index subgroups of the
  group of a projective system, as computed by the driver
  HeckeOperatorSp4/TheBigComputation.g with the GAP code of
  ProjectiveSystem.g, for any projective system given by a description
  file (see README.md):
  - the cell complex of the projective system,
  - the chain map of the double coset of the Hecke matrix on the levels
    up to the level eDimMax of GAP (the levels of the cells of eDimMax
    vectors and less for the configurations of GAP),
  - the homology of the well rounded cells (the cells at infinity
    eliminated) of the quotient by each subgroup, with the Hecke matrix.
 */

namespace {

// A matrix of rationals: the numbers of rows and columns, then the
// entries (a, -a, a/b).
MyMatrix<psys::T> read_rational_matrix(std::string const &file) {
  psys::TokenReader tr{read_description_tokens(file)};
  int n_row = tr.next_int();
  int n_col = tr.next_int();
  MyMatrix<psys::T> M(n_row, n_col);
  for (int i = 0; i < n_row; i++) {
    for (int j = 0; j < n_col; j++) {
      M(i, j) = ParseScalar<psys::T>(tr.next());
    }
  }
  return M;
}

using Level = RingElementName;

}  // namespace

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 8) {
      std::cerr << "PROJ_HeckeOperators [description] [hecke] [eDimMax] [subgroups] [levels] [prime] [prefix]\n";
      std::cerr << "  description: the description file of the projective system\n";
      std::cerr << "  hecke: the file of the Hecke matrix (realified, rational entries)\n";
      std::cerr << "  eDimMax: the number of levels, from the single vectors, on which the chain map is computed\n";
      std::cerr << "  subgroups: comma separated list of subgroup types of src_perfect (Gamma0, Gamma0Siegel, ...)\n";
      std::cerr << "  levels: comma separated list of the levels N of the subgroups (over O_K: a+bw, w = omega)\n";
      std::cerr << "  prime: 0 for the rational homology, otherwise the homology modulo it; comma separated for several\n";
      std::cerr << "  prefix: the GAP results are written to prefix_[subgroup]_[N].g (prefix_[prime]_[subgroup]_[N].g for several primes)\n";
      return -1;
    }
    std::string file_desc = argv[1];
    std::string file_hecke = argv[2];
    int eDimMax = ParseScalar<int>(argv[3]);
    std::vector<std::string> l_sub = STRING_Split(argv[4], ",");
    std::vector<Level> l_level = parse_ring_elements(argv[5]);
    std::vector<int64_t> l_prime;
    for (auto &str : STRING_Split(argv[6], ",")) {
      l_prime.push_back(ParseScalar<int64_t>(str));
    }
    std::string prefix = argv[7];
    std::ostream &os = std::cerr;
    psys::ProjectiveDescription desc = psys::read_description(file_desc);
    psys::ProjectiveCellComplex cx = psys::compute_projective_cell_complex(desc, os);
    psys::ProjectiveRing const &ring = cx.ring();
    int n_levels = cx.levels.size();
    std::cerr << "Complex: |top cells|=" << cx.l_top.size() << " |generators|=" << cx.l_gens_G.size() << "\n";
    for (int index = 0; index < n_levels; index++) {
      std::cerr << "  index=" << index << " orbits:";
      for (auto &orbit : cx.levels[index]) {
        int rnk = psys::ring_rank(ring, orbit.EXT);
        bool orientable = true;
        for (auto &g : orbit.l_stab_gens) {
          if (psys::orbit_orientation(ring, orbit, g) == -1) {
            orientable = false;
          }
        }
        std::cerr << " [|EXT|=" << orbit.EXT.rows() << " rank=" << rnk << " |Stab|=" << orbit.stab_order << " |Stab on EXT|=" << orbit.perm_grp.size()
                  << " orientable=" << orientable << " |faces|=" << orbit.l_bound.size() << "]";
      }
      std::cerr << "\n";
    }
    int index_min = n_levels - eDimMax;
    HeckeChainMap<psys::T, psys::Tint> hcm;
    if (file_hecke == "identity") {
      // The identity: the homology only, the chain map being the identity.
      hcm = psys::identity_chain_map(cx, os);
    } else {
      MyMatrix<psys::T> x = read_rational_matrix(file_hecke);
      if (x.rows() != ring.n() || x.cols() != ring.n()) {
        std::cerr << "PROJ_HeckeOperators: the Hecke matrix should be of dimension " << ring.n() << "\n";
        throw TerminalException{1};
      }
      hcm = psys::compute_projective_hecke_chain_map(cx, x, index_min, os);
    }
    std::cerr << "Chain map computed, |Hecke cosets|=" << hcm.cosets.size() << "\n";
    for (auto &sub : l_sub) {
      for (auto &N : l_level) {
        FiniteIndexSubgroup<psys::Tint> gamma{};
        gamma.type = sub;
        gamma.N = psys::Tint(static_cast<long>(N.a));
        gamma.n = ring.n();
        if (ring.d == 2) {
          gamma.set_field_ideal(ring.m, ring.s, ring.p, {{N.a, N.b}});
        } else if (N.b != 0) {
          std::cerr << "PROJ_HeckeOperators: a level a+bw needs a ring of integers\n";
          throw TerminalException{1};
        }
        for (auto &prime : l_prime) {
          HeckeHomologyResult<psys::T> result =
              compute_hecke_homology(psys::ProjectiveComplexView{cx}, hcm, gamma, true, false, false, prime, os);
          std::string file = prefix + "_" + sub + "_" + N.name + ".g";
          if (l_prime.size() > 1) {
            file = prefix + "_" + std::to_string(prime) + "_" + sub + "_" + N.name + ".g";
          }
          std::ofstream os_out(file);
          os_out << "return ";
          WriteHeckeHomologyResultGAP(os_out, result);
          os_out << ";\n";
          std::cerr << sub << " N=" << N.name << " index=" << result.n_coset_gamma
                    << " |Hecke cosets of Gamma|=" << result.n_hecke_coset << " dimensions:";
          for (auto &level : result.l_level) {
            std::cerr << " H(index " << level.index << ")=" << level.harmonic.rows();
          }
          std::cerr << " -> " << file << "\n";
        }
      }
    }
    std::cerr << "Normal termination of PROJ_HeckeOperators\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of PROJ_HeckeOperators\n";
    exit(e.eVal);
  } catch (HeckeFailure const &e) {
    std::cerr << "PROJ_HeckeOperators: clean failure: " << e.message << "\n";
    exit(1);
  }
  runtime(time);
}
