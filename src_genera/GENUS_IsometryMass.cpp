// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Group.h"
#include "Permutation.h"
#include "LatticeStabEquiCan.h"
#include "IsometryMass.h"
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
// clang-format on

/*
  The mass of the genus of a pair (L, g), g an isometry of finite order of
  the even positive definite lattice L, see IsometryMass.h.

  The check mode is a brute-force test on small lattices: it lists the
  elements of O(L), keeps those tau with the same characteristic
  polynomial as g, computes the mass of the genus of (L, tau) for each, and
  groups them by the mass and the genera of the components. When the genus
  of L has a single class, the pairs over L are all the pairs of these
  genera, and sum_{tau in a genus of pairs X} 1 / |O(L)| = mass(X), so for
  every group the count divided by |O(L)| mass is the number of genera of
  pairs in it: an integer.
 */

template <typename T>
std::map<int, T> ReadHermitianMasses(std::string const &file) {
  std::map<int, T> res;
  if (file == "none") {
    return res;
  }
  std::ifstream is(file);
  if (!is.good()) {
    std::cerr << "GENUS_IsometryMass: cannot read " << file << "\n";
    throw TerminalException{1};
  }
  int d;
  while (is >> d) {
    T num, den;
    is >> num >> den;
    res[d] = num / den;
  }
  return res;
}

template <typename T, typename Tint>
void ProcessMass(std::string const &FileGram, std::string const &FileIsom,
                 std::string const &FileHerm, std::ostream &os) {
  MyMatrix<T> A = ReadMatrixFile<T>(FileGram);
  MyMatrix<T> g = ReadMatrixFile<T>(FileIsom);
  std::map<int, T> herm = ReadHermitianMasses<T>(FileHerm);
  IsometryMassBounds bounds;
  IsometryMassResult<T> res =
      IsometryGenusMass<T, Tint>(A, g, herm, bounds, std::cerr);
  WriteIsometryMassResult(os, res);
}

template <typename T>
std::vector<int64_t> FlattenMatrix(MyMatrix<T> const &M) {
  std::vector<int64_t> v;
  for (int i = 0; i < M.rows(); i++) {
    for (int j = 0; j < M.cols(); j++) {
      v.push_back(UniversalScalarConversion<int64_t, T>(M(i, j)));
    }
  }
  return v;
}

template <typename T, typename Tint>
std::vector<int> KernelRanks(MyMatrix<T> const &P, int m) {
  std::vector<int> res;
  for (int d = 1; d <= m; d++) {
    if (m % d != 0) {
      continue;
    }
    MyMatrix<T> Phi = EvaluatePolynomialMatrix(CyclotomicPolynomial(d), P);
    res.push_back(P.rows() - RankMat(Phi));
  }
  return res;
}

template <typename T, typename Tint, typename Tgroup>
void ProcessCheck(std::string const &FileGram, std::string const &FileIsom,
                  std::string const &FileHerm, std::ostream &os) {
  MyMatrix<T> A = ReadMatrixFile<T>(FileGram);
  MyMatrix<T> g = ReadMatrixFile<T>(FileIsom);
  int n = A.rows();
  MyMatrix<T> P0 = g.transpose();
  int m = IsometryOrder(P0);
  std::vector<int> ranks0 = KernelRanks<T, Tint>(P0, m);
  IntAutomorphismGroup<Tint, Tgroup> aut =
      IteratedShortestAutomorphismGroup<T, Tint, Tgroup>(A, std::cerr);
  std::vector<MyMatrix<T>> gens;
  for (auto const &eGen : aut.ListGen) {
    MyMatrix<T> M = UniversalMatrixConversion<T, Tint>(eGen);
    if (M * A * M.transpose() == A) {
      gens.push_back(M);
    } else if (M.transpose() * A * M == A) {
      gens.push_back(M.transpose());
    } else {
      std::cerr << "GENUS_IsometryMass: a generator is not an isometry\n";
      throw TerminalException{1};
    }
  }
  size_t max_elt = 2000000;
  std::set<std::vector<int64_t>> seen;
  std::vector<MyMatrix<T>> elts{IdentityMat<T>(n)};
  seen.insert(FlattenMatrix(elts[0]));
  for (size_t pos = 0; pos < elts.size(); pos++) {
    for (auto const &h : gens) {
      MyMatrix<T> x = elts[pos] * h;
      if (seen.insert(FlattenMatrix(x)).second) {
        elts.push_back(x);
        if (elts.size() > max_elt) {
          std::cerr << "GENUS_IsometryMass: O(L) is too large for the check\n";
          throw TerminalException{1};
        }
      }
    }
  }
  T order_T = UniversalScalarConversion<T, typename Tgroup::Tint>(aut.order);
  if (ScalarFromInt64<T>(static_cast<int64_t>(elts.size())) != order_T) {
    std::cerr << "GENUS_IsometryMass: listed " << elts.size()
              << " elements but |O(L)| = " << aut.order << "\n";
    throw TerminalException{1};
  }
  std::map<int, T> herm = ReadHermitianMasses<T>(FileHerm);
  IsometryMassBounds bounds;
  std::map<std::string, std::pair<int64_t, T>> groups;
  for (auto const &P : elts) {
    if (IsometryOrder(P) != m || KernelRanks<T, Tint>(P, m) != ranks0) {
      continue;
    }
    IsometryMassResult<T> res =
        IsometryGenusMass<T, Tint>(A, P.transpose(), herm, bounds, std::cerr);
    if (!res.mass) {
      std::cerr << "GENUS_IsometryMass: the check needs the masses of all "
                   "the components\n";
      throw TerminalException{1};
    }
    std::ostringstream key;
    key << "mass=" << *res.mass;
    for (auto const &cm : res.component_masses) {
      key << " [d=" << cm.d << " rank=" << cm.rank << " det=" << cm.det
          << " mass=" << *cm.mass << "]";
    }
    key << " glue=" << res.glue_factor;
    if (res.image_assumed) {
      key << " image_assumed";
    }
    auto iter = groups.find(key.str());
    if (iter == groups.end()) {
      groups[key.str()] = {1, *res.mass};
    } else {
      iter->second.first++;
    }
  }
  os << "|O(L)| = " << elts.size() << "\n";
  bool all_integral = true;
  for (auto const &kv : groups) {
    T ratio = ScalarFromInt64<T>(kv.second.first) / (order_T * kv.second.second);
    os << "count=" << kv.second.first << " ratio=" << ratio << " " << kv.first
       << "\n";
    if (!IsInteger(ratio)) {
      all_integral = false;
    }
  }
  os << "check: " << (all_integral ? "every ratio is an integer"
                                   : "FAILED, a ratio is not an integer")
     << "\n";
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    using T = mpq_class;
    using Tint = mpz_class;
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using TintGroup = mpz_class;
    using Tgroup = permutalib::Group<Telt, TintGroup>;
    if (argc != 4 && argc != 5 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GENUS_IsometryMass mass [Gram] [Isometry] [HermMass] "
                   "[OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "GENUS_IsometryMass check [Gram] [Isometry] [HermMass] "
                   "[OutFile]\n";
      std::cerr << "\n";
      std::cerr << "Gram     (input) : the Gram matrix A of an even positive\n";
      std::cerr << "                   definite lattice L\n";
      std::cerr << "Isometry (input) : g in GL_n(Z) of finite order with\n";
      std::cerr << "                   g^T A g = A\n";
      std::cerr << "HermMass (input) : 'none', or lines 'd num den', the\n";
      std::cerr << "                   masses of the Hermitian components\n";
      std::cerr << "                   L cap ker Phi_d(g), d >= 3 (from\n";
      std::cerr << "                   Hecke); without them the program\n";
      std::cerr << "                   gives the glue factor and the\n";
      std::cerr << "                   orthogonal masses only\n";
      std::cerr << "OutFile          : stderr (default), stdout or a file\n";
      std::cerr << "\n";
      std::cerr << "mass  : the mass of the genus of the pair (L, g)\n";
      std::cerr << "check : the brute-force test over the elements of O(L),\n";
      std::cerr << "        for small lattices whose genus has one class;\n";
      std::cerr << "        HermMass is then used for every element tau,\n";
      std::cerr << "        so it is only meaningful when the Hermitian\n";
      std::cerr << "        components all have the same mass (rank one\n";
      std::cerr << "        over Z[zeta_3] or Z[i]: 1/6 or 1/4)\n";
      return -1;
    }
    std::string mode = argv[1];
    std::string FileGram = argv[2];
    std::string FileIsom = argv[3];
    if (mode == "mass") {
      std::string FileHerm = (argc >= 5) ? argv[4] : "none";
      std::string OutFile = (argc >= 6) ? argv[5] : "stderr";
      auto f = [&](std::ostream &os) -> void {
        ProcessMass<T, Tint>(FileGram, FileIsom, FileHerm, os);
      };
      FILE_PrintStderrStdoutFile(OutFile, f);
    } else if (mode == "check") {
      std::string FileHerm = (argc >= 5) ? argv[4] : "none";
      std::string OutFile = (argc >= 6) ? argv[5] : "stderr";
      auto f = [&](std::ostream &os) -> void {
        ProcessCheck<T, Tint, Tgroup>(FileGram, FileIsom, FileHerm, os);
      };
      FILE_PrintStderrStdoutFile(OutFile, f);
    } else {
      std::cerr << "GENUS_IsometryMass: the mode must be mass or check\n";
      throw TerminalException{1};
    }
    std::cerr << "Normal termination of GENUS_IsometryMass\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GENUS_IsometryMass\n";
    exit(e.eVal);
  }
  runtime(time);
}
