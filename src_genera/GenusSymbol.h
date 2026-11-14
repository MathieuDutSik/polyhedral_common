// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_GENUSSYMBOL_H_
#define SRC_GENERA_GENUSSYMBOL_H_

// clang-format off
#include "GeneraArithmetic.h"
#include "LocalDensity.h"
#include "LocalDensityOverlattice.h"
#include "ClassicLLL.h"
#include "MAT_MatrixInt.h"
#include "ResidueArithmetic.h"
#include <array>
#include <cstdint>
#include <map>
#include <optional>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_GENUS_SYMBOL
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_GENUS_SYMBOL
#endif

/*
  The genus of a positive definite integral lattice as a symbol, and a
  lattice in a genus given by its symbol.

  The genus of L is determined by its rank, its parity and its Jordan
  decompositions at the primes p dividing 2 det(L) (Conway-Sloane, Sphere
  Packings, Lattices and Groups, Chapter 15):
  * at an odd p, a component q^{eps n} per scale q = p^k: its rank n and
    the Legendre symbol eps of the product of its units;
  * at p = 2, a component q^{eps n}_t or q^{eps n}_II per scale q = 2^k:
    eps the Jacobi symbol (2 / unit determinant), and for an odd component
    (type I) its oddity t, the trace of a diagonalization mod 8. Only the
    sums of the oddities over a compartment (a run of odd components of
    consecutive scales) are invariants, and the signs can be walked along
    a train; CanonicalTwoAdicReduction picks one symbol per class.

  The text form is that of the labels of DefLatDB (Brandhorst-Hofmann),
  "rank.n_plus.det.level." followed by the symbol, the prefix being
  optional:
      II{1}^{8}_{II}{2}^{2}_{2}{5}^{2}
      I[{1}^{-8}{2}^{2}]_{4}{3}^{1}
  "II" for an even lattice, "I" for an odd one; {q}^{n} a component, the
  sign of n being eps; at p = 2 the subscript _{t} or _{II} gives the type,
  and [ ... ]_{t} a compartment with its total oddity. The unimodular
  components at the odd primes can be omitted: their rank and sign follow
  from the others and the determinant.

  A lattice of the genus is built as in Hecke.jl and SageMath
  (representative), but with search steps in place of the p-adic normal
  forms:
  1. a local Gram matrix G_p per prime p | 2d, block diagonal with the
     components of the symbol (TwoAdicLocalGram, OddPrimeLocalGram);
  2. the rational quadratic space V: a positive diagonal form of
     determinant d with the Hasse invariants of the G_p at p | 2d and 1
     elsewhere (RationalDiagonalWithInvariants); it exists iff the
     product formula holds, which is the oddity formula;
  3. an even lattice of V, made maximal at every prime of its determinant
     by isotropic steps L -> L + Z y/p (IsotropicRadicalVector). The even
     lattices that are maximal at p are all isometric at p, which fails
     for the integral ones at p = 2 (Z^8 and E_8 span one space), so an
     odd genus is built as 2 L, of Gram 4 G and even, and halved at the
     end, as in SageMath. Outside 2d the maximal lattice is unimodular
     since V has a unit determinant and Hasse invariant 1;
  4. at every p | 2d, the descent from the maximal lattice to one locally
     isometric to G_p by index p sublattices (DescendToLocalSymbol),
     following backwards the local symbols of a chain going up from G_p to
     a maximal lattice: every lattice of V_p lies in a maximal one and the
     maximal even lattices of V_p are isometric, so the target is reached;
  5. an LLL reduction, and the check that the genus is the given one.
 */

// The symbol of a genus. The 2-adic part is in canonical form: one
// (scale exponent, rank, sign, type, oddity) per Jordan component, sign
// +-1, type 1 for odd components, the oddity of a compartment on its first
// component. The odd primes are those dividing the determinant, with
// (scale exponent, rank, sign) per component, the unimodular one included.
struct GenusSymbol {
  int rank;
  bool even;
  std::vector<std::array<int, 5>> sym2;
  std::map<int64_t, std::vector<std::array<int, 3>>> sym_odd;
};

inline bool operator==(GenusSymbol const &a, GenusSymbol const &b) {
  return a.rank == b.rank && a.even == b.even && a.sym2 == b.sym2 &&
         a.sym_odd == b.sym_odd;
}

inline bool operator!=(GenusSymbol const &a, GenusSymbol const &b) {
  return !(a == b);
}

template <typename Tint> Tint GenusSymbolDeterminant(GenusSymbol const &s) {
  Tint det(1);
  for (auto const &c : s.sym2) {
    det *= MyPow<Tint>(Tint(2), c[0] * c[1]);
  }
  for (auto const &[p, comps] : s.sym_odd) {
    Tint p_T(static_cast<long>(p));
    for (auto const &c : comps) {
      det *= MyPow<Tint>(p_T, c[0] * c[1]);
    }
  }
  return det;
}

inline std::string GenusSymbolString(GenusSymbol const &s) {
  std::string str = s.even ? "II" : "I";
  auto comp = [](int64_t q, int rank, int sign) -> std::string {
    return "{" + std::to_string(q) + "}^{" +
           std::to_string(sign * rank) + "}";
  };
  auto power = [](int64_t p, int k) -> int64_t {
    int64_t q = MyPow<int64_t>(p, k);
    return q;
  };
  size_t i = 0;
  while (i < s.sym2.size()) {
    auto const &c = s.sym2[i];
    if (c[3] == 0) {
      str += comp(power(2, c[0]), c[1], c[2]) + "_{II}";
      i++;
      continue;
    }
    // the compartment starting at i
    size_t j = i + 1;
    while (j < s.sym2.size() && s.sym2[j][3] == 1 &&
           s.sym2[j][0] == s.sym2[j - 1][0] + 1) {
      j++;
    }
    int oddity = 0;
    std::string inner;
    for (size_t k = i; k < j; k++) {
      inner += comp(power(2, s.sym2[k][0]), s.sym2[k][1], s.sym2[k][2]);
      oddity += s.sym2[k][4];
    }
    oddity = static_cast<int>(ResMod(oddity, 8));
    if (j == i + 1) {
      str += inner + "_{" + std::to_string(oddity) + "}";
    } else {
      str += "[" + inner + "]_{" + std::to_string(oddity) + "}";
    }
    i = j;
  }
  for (auto const &[p, comps] : s.sym_odd) {
    for (auto const &c : comps) {
      if (c[0] > 0) {
        str += comp(power(p, c[0]), c[1], c[2]);
      }
    }
  }
  return str;
}

// The Jordan components (scale exponent, rank, sign) of G at the odd p.
template <typename T, typename Tint>
std::vector<std::array<int, 3>> OddPrimeSymbolOfGram(MyMatrix<T> const &G,
                                                     int64_t p) {
  std::map<int, std::pair<int, int64_t>> levels;
  for (auto const &x : JordanDiagonalOdd<T, Tint>(G, p)) {
    PadicClass c = PadicFromRational<T, Tint>(x, p);
    auto iter = levels.find(c.val);
    if (iter == levels.end()) {
      iter = levels.emplace(c.val, std::make_pair(0, int64_t(1))).first;
    }
    iter->second.first++;
    iter->second.second = ResMul(iter->second.second, c.unit, p);
  }
  std::vector<std::array<int, 3>> res;
  for (auto const &[k, rp] : levels) {
    res.push_back({k, rp.first, kronecker_symbol<int64_t>(rp.second, p)});
  }
  return res;
}

template <typename T, typename Tint>
GenusSymbol GenusSymbolOfGram(MyMatrix<T> const &G) {
  GenusSymbol s;
  int n = G.rows();
  s.rank = n;
  s.even = true;
  for (int i = 0; i < n; i++) {
    if (!IsInteger(G(i, i) / T(2))) {
      s.even = false;
    }
  }
  // the Jordan splitting of G, as TwoAdicSymbolRaw reads it from 2 G
  std::vector<std::array<int, 5>> raw = TwoAdicSymbolRaw<T, Tint>(T(2) * G);
  for (auto &c : raw) {
    c[0] -= 1;
  }
  s.sym2 = CanonicalTwoAdicReduction(raw);
  T det = DeterminantMat(G);
  Tint det_i = UniversalScalarConversion<Tint, T>(det);
  for (int64_t p : OddPrimeDivisors<Tint>(det_i)) {
    s.sym_odd[p] = OddPrimeSymbolOfGram<T, Tint>(G, p);
  }
  return s;
}

//
// Parsing
//

template <typename Tint> GenusSymbol ParseGenusSymbol(std::string const &input) {
  std::string str;
  for (char c : input) {
    if (c != ' ' && c != '\t' && c != '\n') {
      str += c;
    }
  }
  // the optional prefix rank.n_plus.det.level.
  std::vector<std::string> fields;
  {
    size_t start = 0;
    while (true) {
      size_t pos = str.find('.', start);
      if (pos == std::string::npos) {
        fields.push_back(str.substr(start));
        break;
      }
      fields.push_back(str.substr(start, pos - start));
      start = pos + 1;
    }
  }
  std::optional<int> prefix_rank, prefix_nplus;
  std::optional<Tint> prefix_det;
  if (fields.size() == 5) {
    prefix_rank = std::stoi(fields[0]);
    prefix_nplus = std::stoi(fields[1]);
    prefix_det = Tint(fields[2]);
    str = fields[4];
  } else if (fields.size() != 1) {
    std::cerr << "GSYM: the genus symbol " << input
              << " is neither a symbol nor rank.n_plus.det.level.symbol\n";
    throw TerminalException{1};
  }
  size_t pos = 0;
  auto fail = [&](std::string const &what) -> void {
    std::cerr << "GSYM: cannot parse the genus symbol " << input << ": "
              << what << " at position " << pos << "\n";
    throw TerminalException{1};
  };
  GenusSymbol s;
  if (str.compare(0, 2, "II") == 0) {
    s.even = true;
    pos = 2;
  } else if (str.compare(0, 1, "I") == 0) {
    s.even = false;
    pos = 1;
  } else {
    fail("no I or II prefix");
  }
  // an item, braced or not: digits with a sign, or the letters II
  auto read_item = [&]() -> std::string {
    std::string item;
    if (pos < str.size() && str[pos] == '{') {
      size_t end = str.find('}', pos);
      if (end == std::string::npos) {
        fail("an unclosed brace");
      }
      item = str.substr(pos + 1, end - pos - 1);
      pos = end + 1;
    } else {
      while (pos < str.size() &&
             (isdigit(str[pos]) || str[pos] == '-' || str[pos] == '+' ||
              str[pos] == 'I')) {
        item += str[pos];
        pos++;
      }
    }
    if (item.empty()) {
      fail("an empty item");
    }
    return item;
  };
  // (scale exponent, rank, sign, type, oddity), oddity -1 when unknown
  std::vector<std::array<int, 5>> comps2;
  std::map<int64_t, std::vector<std::array<int, 3>>> comps_odd;
  // reads {q}^{n}, returns (p, k, rank, sign)
  auto read_component = [&]() -> std::array<int64_t, 4> {
    int64_t q = std::stoll(read_item());
    if (pos >= str.size() || str[pos] != '^') {
      fail("a missing ^");
    }
    pos++;
    int e = std::stoi(read_item());
    auto [p, k] = PrimePowerDecomposition(q);
    return {p, k, std::abs(e), (e < 0) ? -1 : 1};
  };
  while (pos < str.size()) {
    if (str[pos] == '[') {
      pos++;
      std::vector<std::array<int, 5>> compart;
      while (pos < str.size() && str[pos] != ']') {
        auto [p, k, r, sign] = read_component();
        if (p != 2) {
          fail("an odd prime in a compartment");
        }
        compart.push_back({static_cast<int>(k), static_cast<int>(r),
                           static_cast<int>(sign), 1, 0});
      }
      if (pos >= str.size()) {
        fail("an unclosed compartment");
      }
      pos++;
      if (pos >= str.size() || str[pos] != '_') {
        fail("a compartment without its oddity");
      }
      pos++;
      int t = std::stoi(read_item());
      compart[0][4] = static_cast<int>(ResMod(t, 8));
      comps2.insert(comps2.end(), compart.begin(), compart.end());
      continue;
    }
    auto [p, k, r, sign] = read_component();
    if (p == 2) {
      if (pos >= str.size() || str[pos] != '_') {
        fail("a 2-adic component without its type");
      }
      pos++;
      std::string sub = read_item();
      if (sub == "II") {
        comps2.push_back({static_cast<int>(k), static_cast<int>(r),
                          static_cast<int>(sign), 0, 0});
      } else {
        comps2.push_back({static_cast<int>(k), static_cast<int>(r),
                          static_cast<int>(sign), 1,
                          static_cast<int>(ResMod(std::stoi(sub), 8))});
      }
    } else {
      if (pos < str.size() && str[pos] == '_') {
        fail("a subscript on an odd prime");
      }
      comps_odd[p].push_back({static_cast<int>(k), static_cast<int>(r),
                              static_cast<int>(sign)});
    }
  }
  // the 2-adic components, by scale; those of rank 0 are dropped
  std::sort(comps2.begin(), comps2.end(),
            [](auto const &a, auto const &b) { return a[0] < b[0]; });
  std::vector<std::array<int, 5>> raw2;
  int rank = 0;
  for (auto const &c : comps2) {
    if (!raw2.empty() && raw2.back()[0] == c[0]) {
      fail("two 2-adic components of the same scale");
    }
    if (c[1] > 0) {
      raw2.push_back(c);
      rank += c[1];
    }
  }
  if (rank == 0) {
    fail("no 2-adic component");
  }
  s.rank = rank;
  if (prefix_rank && *prefix_rank != rank) {
    std::cerr << "GSYM: the 2-adic components of " << input << " have rank "
              << rank << ", not " << *prefix_rank << "\n";
    throw TerminalException{1};
  }
  if (prefix_nplus && *prefix_nplus != rank) {
    std::cerr << "GSYM: the genus " << input
              << " is not positive definite, only those are handled\n";
    throw TerminalException{1};
  }
  for (auto const &c : raw2) {
    if (c[3] == 0 && c[1] % 2 != 0) {
      fail("an even 2-adic component of odd rank");
    }
  }
  bool has_odd_unimodular = !raw2.empty() && raw2[0][0] == 0 && raw2[0][3] == 1;
  if (s.even == has_odd_unimodular) {
    std::cerr << "GSYM: in " << input
              << ", the prefix I / II does not match the type of the "
                 "unimodular 2-adic component\n";
    throw TerminalException{1};
  }
  s.sym2 = CanonicalTwoAdicReduction(raw2);
  // the determinant, from the components
  Tint det(1);
  for (auto const &c : raw2) {
    det *= MyPow<Tint>(Tint(2), c[0] * c[1]);
  }
  for (auto const &[p, comps] : comps_odd) {
    Tint p_T(static_cast<long>(p));
    for (auto const &c : comps) {
      det *= MyPow<Tint>(p_T, c[0] * c[1]);
    }
  }
  if (prefix_det && *prefix_det != det) {
    std::cerr << "GSYM: the components of " << input << " have determinant "
              << det << ", not " << *prefix_det << "\n";
    throw TerminalException{1};
  }
  // the sign of the 2-adic symbol is that of the determinant
  {
    Tint u = det;
    while (ResInt(u, Tint(2)) == 0) {
      u /= Tint(2);
    }
    int64_t u8 = UniversalScalarConversion<int64_t, Tint>(ResInt(u, Tint(8)));
    int sign = (u8 == 1 || u8 == 7) ? 1 : -1;
    int prod = 1;
    for (auto const &c : raw2) {
      prod *= c[2];
    }
    if (prod != sign) {
      std::cerr << "GSYM: the 2-adic signs of " << input
                << " do not match the determinant " << det << "\n";
      throw TerminalException{1};
    }
  }
  // the odd primes, with their unimodular components
  for (auto &[p, comps] : comps_odd) {
    std::sort(comps.begin(), comps.end());
    Tint p_T(static_cast<long>(p));
    Tint u = det;
    while (ResInt(u, p_T) == 0) {
      u /= p_T;
    }
    int64_t u_p = UniversalScalarConversion<int64_t, Tint>(ResInt(u, p_T));
    int sign_det = kronecker_symbol<int64_t>(u_p, p);
    int prod = 1, r_other = 0;
    bool has_unimodular = false;
    std::vector<std::array<int, 3>> kept;
    for (size_t i = 0; i < comps.size(); i++) {
      if (i > 0 && comps[i][0] == comps[i - 1][0]) {
        fail("two components of the same scale");
      }
      if (comps[i][0] == 0) {
        has_unimodular = true;
      } else {
        r_other += comps[i][1];
      }
      prod *= comps[i][2];
      if (comps[i][1] > 0) {
        kept.push_back(comps[i]);
      }
    }
    if (r_other > rank) {
      std::cerr << "GSYM: the components at " << p << " of " << input
                << " have a rank above " << rank << "\n";
      throw TerminalException{1};
    }
    if (!has_unimodular) {
      int r0 = rank - r_other;
      if (r0 > 0) {
        kept.insert(kept.begin(), {0, r0, prod * sign_det});
      } else if (prod != sign_det) {
        std::cerr << "GSYM: the signs at " << p << " of " << input
                  << " do not match the determinant\n";
        throw TerminalException{1};
      }
    } else {
      int r_total = 0;
      for (auto const &c : kept) {
        r_total += c[1];
      }
      if (r_total != rank || prod != sign_det) {
        std::cerr << "GSYM: the components at " << p << " of " << input
                  << " do not match the rank and the determinant\n";
        throw TerminalException{1};
      }
    }
    if (kept.size() > 1 || (kept.size() == 1 && kept[0][0] > 0)) {
      s.sym_odd[p] = kept;
    }
  }
  return s;
}

//
// The local Gram matrices
//

// Odd p: p^k diag(1, ..., 1, delta) per component, delta a non-residue for
// the sign -1.
template <typename T>
MyMatrix<T> OddPrimeLocalGram(int64_t p,
                              std::vector<std::array<int, 3>> const &comps,
                              int n) {
  MyMatrix<T> G = ZeroMatrix<T>(n, n);
  int pos = 0;
  for (auto const &c : comps) {
    T q = MyPow<T>(UniversalScalarConversion<T, int64_t>(p), c[0]);
    for (int i = 0; i < c[1]; i++) {
      G(pos + i, pos + i) = q;
    }
    if (c[2] == -1) {
      G(pos + c[1] - 1, pos + c[1] - 1) *= T(static_cast<long>(LeastNonResidue(p)));
    }
    pos += c[1];
  }
  return G;
}

// Units u_1, ..., u_r of Z_2 (mod 8) with product of sign `sign` (+1 for
// +-1 mod 8) and sum t mod 8, by a dynamic programming over (product, sum).
inline std::optional<std::vector<int64_t>> OddUnitDiagonal(int r, int sign,
                                                           int t) {
  // states: product in {1, 3, 5, 7} (index u / 2) and sum mod 8
  std::vector<std::array<int, 32>> pred(r + 1);
  std::vector<std::array<bool, 32>> reach(r + 1);
  for (auto &a : reach) {
    a.fill(false);
  }
  auto idx = [](int64_t prod, int64_t sum) -> int {
    return static_cast<int>((prod / 2) * 8 + sum);
  };
  reach[0][idx(1, 0)] = true;
  for (int i = 0; i < r; i++) {
    for (int64_t prod : {1, 3, 5, 7}) {
      for (int64_t sum = 0; sum < 8; sum++) {
        if (!reach[i][idx(prod, sum)]) {
          continue;
        }
        for (int64_t u : {1, 3, 5, 7}) {
          int j = idx((prod * u) % 8, (sum + u) % 8);
          if (!reach[i + 1][j]) {
            reach[i + 1][j] = true;
            pred[i + 1][j] = idx(prod, sum) * 8 + static_cast<int>(u);
          }
        }
      }
    }
  }
  for (int64_t prod : {1, 3, 5, 7}) {
    int s = (prod == 1 || prod == 7) ? 1 : -1;
    if (s != sign || !reach[r][idx(prod, t)]) {
      continue;
    }
    std::vector<int64_t> units;
    int j = idx(prod, t);
    for (int i = r; i > 0; i--) {
      units.push_back(pred[i][j] % 8);
      j = pred[i][j] / 8;
    }
    return units;
  }
  return {};
}

// A Gram matrix of the 2-adic class: the canonical symbol fixes only the
// sums of the oddities over the compartments and the products of the signs
// over the trains, so the individual (sign, oddity) of the components are
// searched among those whose canonical reduction is the given one.
template <typename T>
MyMatrix<T> TwoAdicLocalGram(std::vector<std::array<int, 5>> const &sym2,
                             int n) {
  size_t r = sym2.size();
  // per component, the possible (sign, oddity)
  std::vector<std::vector<std::pair<int, int>>> options(r);
  for (size_t i = 0; i < r; i++) {
    for (int sign : {1, -1}) {
      if (sym2[i][3] == 0) {
        options[i].push_back({sign, 0});
        continue;
      }
      for (int t = 0; t < 8; t++) {
        if (OddUnitDiagonal(sym2[i][1], sign, t)) {
          options[i].push_back({sign, t});
        }
      }
    }
  }
  std::vector<size_t> choice(r, 0);
  size_t n_tried = 0;
  while (true) {
    std::vector<std::array<int, 5>> raw = sym2;
    for (size_t i = 0; i < r; i++) {
      raw[i][2] = options[i][choice[i]].first;
      raw[i][4] = options[i][choice[i]].second;
    }
    if (CanonicalTwoAdicReduction(raw) == sym2) {
      MyMatrix<T> G = ZeroMatrix<T>(n, n);
      int pos = 0;
      for (auto const &c : raw) {
        T q = MyPow<T>(UniversalScalarConversion<T, int64_t>(2), c[0]);
        if (c[3] == 1) {
          std::vector<int64_t> units = *OddUnitDiagonal(c[1], c[2], c[4]);
          for (int i = 0; i < c[1]; i++) {
            G(pos + i, pos + i) = q * T(static_cast<long>(units[i]));
          }
        } else {
          // H = [[0, 1], [1, 0]] (sign +1), V = [[2, 1], [1, 2]] (sign -1)
          for (int i = 0; i < c[1]; i += 2) {
            G(pos + i, pos + i + 1) = q;
            G(pos + i + 1, pos + i) = q;
          }
          if (c[2] == -1) {
            G(pos, pos) = T(2) * q;
            G(pos + 1, pos + 1) = T(2) * q;
          }
        }
        pos += c[1];
      }
      return G;
    }
    n_tried++;
    size_t i = 0;
    while (i < r) {
      choice[i]++;
      if (choice[i] < options[i].size()) {
        break;
      }
      choice[i] = 0;
      i++;
    }
    if (i == r || n_tried > 4000000) {
      break;
    }
  }
  std::cerr << "GSYM: no 2-adic lattice has the 2-adic symbol of the genus\n";
  throw TerminalException{1};
}

//
// Rational invariants
//

// The diagonal of a diagonalization over Q of the non-degenerate S.
template <typename T> std::vector<T> RationalDiagonalForm(MyMatrix<T> S) {
  int n = S.rows();
  std::vector<T> diag;
  for (int k = 0; k < n; k++) {
    if (S(k, k) == 0) {
      int j_diag = -1, j_off = -1;
      for (int j = k + 1; j < n; j++) {
        if (S(j, j) != 0 && j_diag == -1) {
          j_diag = j;
        }
        if (S(k, j) != 0 && j_off == -1) {
          j_off = j;
        }
      }
      if (j_diag != -1) {
        S.row(k).swap(S.row(j_diag));
        S.col(k).swap(S.col(j_diag));
      } else if (j_off != -1) {
        S.row(k) += S.row(j_off);
        S.col(k) += S.col(j_off);
      } else {
        std::cerr << "GSYM: RationalDiagonalForm of a degenerate matrix\n";
        throw TerminalException{1};
      }
    }
    T a = S(k, k);
    diag.push_back(a);
    for (int i = k + 1; i < n; i++) {
      T f = S(i, k) / a;
      for (int j = k + 1; j < n; j++) {
        S(i, j) -= f * S(k, j);
      }
    }
  }
  return diag;
}

template <typename T, typename Tint>
int HilbertSymbolRational(T const &x, T const &y, int64_t p) {
  return HilbertSymbolClass(PadicFromRational<T, Tint>(x, p),
                            PadicFromRational<T, Tint>(y, p), p);
}

// The Hasse invariant prod_{i < j} (a_i, a_j)_p of diag(a).
template <typename T, typename Tint>
int HasseInvariantDiagonal(std::vector<T> const &diag, int64_t p) {
  int h = 1;
  T prod(1);
  for (size_t i = 0; i < diag.size(); i++) {
    if (i > 0) {
      h *= HilbertSymbolRational<T, Tint>(prod, diag[i], p);
    }
    prod *= diag[i];
  }
  return h;
}

/*
  A positive diagonal form over Q of rank n and determinant d with the
  Hasse invariant hasse[p] at the primes of `hasse` and 1 at the others:
  (d) for n = 1, (x, d x) for n = 2, (1, ..., 1, x, y, d x y) for n >= 3,
  with squarefree x, y found by search. The invariants at a prime outside
  those of 2 d x y are 1 since all the entries are units there.
 */
template <typename T, typename Tint>
std::vector<T> RationalDiagonalWithInvariants(int n, Tint const &d,
                                              std::map<int64_t, int> const &hasse) {
  T d_T = UniversalScalarConversion<T, Tint>(d);
  auto required = [&](int64_t p) -> int {
    auto iter = hasse.find(p);
    return (iter == hasse.end()) ? 1 : iter->second;
  };
  auto matches = [&](std::vector<T> const &diag,
                     std::vector<int64_t> const &extra) -> bool {
    std::set<int64_t> primes;
    for (auto const &kv : hasse) {
      primes.insert(kv.first);
    }
    primes.insert(extra.begin(), extra.end());
    for (int64_t p : primes) {
      if (HasseInvariantDiagonal<T, Tint>(diag, p) != required(p)) {
        return false;
      }
    }
    return true;
  };
  if (n == 1) {
    std::vector<T> diag{d_T};
    if (matches(diag, {})) {
      return diag;
    }
  } else if (n == 2) {
    for (int64_t x = 1; x <= 1000000; x++) {
      if (!IsSquareFree(x)) {
        continue;
      }
      T x_T(static_cast<long>(x));
      std::vector<T> diag{x_T, d_T * x_T};
      if (matches(diag, SmallPrimeDivisors(x))) {
        return diag;
      }
    }
  } else {
    for (int64_t sum = 2; sum <= 4000; sum++) {
      for (int64_t x = 1; x < sum; x++) {
        int64_t y = sum - x;
        if (!IsSquareFree(x) || !IsSquareFree(y)) {
          continue;
        }
        T x_T(static_cast<long>(x)), y_T(static_cast<long>(y));
        std::vector<T> three{x_T, y_T, d_T * x_T * y_T};
        std::vector<int64_t> extra = SmallPrimeDivisors(x);
        for (int64_t q : SmallPrimeDivisors(y)) {
          extra.push_back(q);
        }
        if (matches(three, extra)) {
          std::vector<T> diag(n - 3, T(1));
          diag.insert(diag.end(), three.begin(), three.end());
          return diag;
        }
      }
    }
  }
  std::cerr << "GSYM: no rational quadratic form has the invariants of the "
               "genus (the oddity formula fails?)\n";
  throw TerminalException{1};
}

//
// Lattices of the rational space
//

// A lattice of V by its basis B in the coordinates of V and its Gram G.
template <typename T> struct SymbolLattice {
  MyMatrix<T> B;
  MyMatrix<T> G;
};

template <typename T>
void ApplyBasisChange(SymbolLattice<T> &L, MyMatrix<T> const &E) {
  L.B = E * L.B;
  L.G = E * L.G * E.transpose();
}

// A non-zero isotropic vector of the quadratic form c M c^T over F_p, p
// odd, M symmetric, by a diagonalization W M W^T = diag(delta).
inline std::optional<std::vector<int64_t>>
IsotropicVectorFp(std::vector<std::vector<int64_t>> M, int64_t p) {
  int m = M.size();
  std::vector<std::vector<int64_t>> W(m, std::vector<int64_t>(m, 0));
  for (int i = 0; i < m; i++) {
    W[i][i] = 1;
  }
  auto add_row = [&](int k, int j, int64_t f) -> void {
    // e_k <- e_k + f e_j
    for (int c = 0; c < m; c++) {
      W[k][c] = ResMod(W[k][c] + f * W[j][c], p);
    }
    for (int c = 0; c < m; c++) {
      M[k][c] = ResMod(M[k][c] + f * M[j][c], p);
    }
    for (int c = 0; c < m; c++) {
      M[c][k] = M[k][c];
    }
    M[k][k] = ResMod(M[k][k] + f * M[j][k], p);
  };
  std::vector<int64_t> delta;
  for (int k = 0; k < m; k++) {
    if (M[k][k] == 0) {
      int j_off = -1;
      for (int j = k + 1; j < m; j++) {
        if (M[k][j] != 0) {
          j_off = j;
          break;
        }
      }
      if (j_off == -1) {
        // e_k is in the radical of the remaining form, so isotropic
        return W[k];
      }
      if (M[j_off][j_off] == 0) {
        add_row(k, j_off, 1);
      } else {
        std::swap(W[k], W[j_off]);
        std::swap(M[k], M[j_off]);
        for (int c = 0; c < m; c++) {
          std::swap(M[c][k], M[c][j_off]);
        }
      }
    }
    int64_t inv = ResInv(M[k][k], p);
    for (int i = k + 1; i < m; i++) {
      if (M[i][k] != 0) {
        add_row(i, k, ResMod(-ResMul(M[i][k], inv, p), p));
      }
    }
    delta.push_back(M[k][k]);
  }
  auto combine = [&](std::vector<std::pair<int, int64_t>> const &terms)
      -> std::vector<int64_t> {
    std::vector<int64_t> v(m, 0);
    for (auto const &[i, f] : terms) {
      for (int c = 0; c < m; c++) {
        v[c] = ResMod(v[c] + ResMul(f, W[i][c], p), p);
      }
    }
    return v;
  };
  for (int i = 0; i < m; i++) {
    for (int j = i + 1; j < m; j++) {
      int64_t r = ResMul(p - delta[i], ResInv(delta[j], p), p);
      if (kronecker_symbol<int64_t>(r, p) == 1) {
        return combine({{i, 1}, {j, SqrtModPrime(r, p)}});
      }
    }
  }
  if (m >= 3) {
    // delta_0 + delta_1 y^2 + delta_2 z^2 = 0 has a solution
    int64_t inv2 = ResInv(delta[2], p);
    for (int64_t y = 0; y < p; y++) {
      int64_t val = ResMod(delta[0] + ResMul(delta[1], ResMul(y, y, p), p), p);
      int64_t r = ResMul(p - val, inv2, p);
      if (r == 0 || kronecker_symbol<int64_t>(r, p) == 1) {
        return combine({{0, 1}, {1, y}, {2, SqrtModPrime(r, p)}});
      }
    }
  }
  return {};
}

/*
  A vector y of the even L, mod p, with L + Z y / p even: y G = 0 mod p,
  and y G y^T = 0 mod p^2 (odd p) or mod 8 (p = 2). Odd p: the form y G y^T / p mod p on the radical mod p is
  diagonalized. p = 2: the spans of the first j vectors of a basis of the
  radical are enumerated, j up to 20.
 */
template <typename T, typename Tint>
std::optional<std::vector<int64_t>>
IsotropicRadicalVector(MyMatrix<T> const &G, int64_t p) {
  std::vector<std::vector<int64_t>> R = KernelBasisModP<T, Tint>(G, p);
  int m = R.size();
  int n = G.rows();
  if (m == 0) {
    return {};
  }
  auto vec = [&](std::vector<int64_t> const &v) -> MyVector<T> {
    MyVector<T> w(n);
    for (int i = 0; i < n; i++) {
      w(i) = T(static_cast<long>(v[i]));
    }
    return w;
  };
  std::vector<MyVector<T>> RT;
  for (auto const &r : R) {
    RT.push_back(vec(r));
  }
  auto combine = [&](std::vector<int64_t> const &c) -> std::vector<int64_t> {
    std::vector<int64_t> y(n, 0);
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        y[j] = ResMod(y[j] + c[i] * R[i][j], p);
      }
    }
    return y;
  };
  if (p != 2) {
    T p_T(static_cast<long>(p));
    std::vector<std::vector<int64_t>> M(m, std::vector<int64_t>(m));
    for (int i = 0; i < m; i++) {
      for (int j = i; j < m; j++) {
        T val = RT[i].dot(G * RT[j]) / p_T;
        M[i][j] = M[j][i] = ResidueModInt<T, Tint>(val, p);
      }
    }
    std::optional<std::vector<int64_t>> c = IsotropicVectorFp(M, p);
    if (!c) {
      return {};
    }
    return combine(*c);
  }
  std::vector<std::vector<int64_t>> N(m, std::vector<int64_t>(m));
  for (int i = 0; i < m; i++) {
    for (int j = i; j < m; j++) {
      T val = RT[i].dot(G * RT[j]);
      N[i][j] = N[j][i] = ResidueModInt<T, Tint>(val, 8);
    }
  }
  int jmax = std::min(m, 20);
  for (int j = 0; j < jmax; j++) {
    // the vectors with c_j = 1, c_i = 0 for i > j
    for (int64_t mask = 0; mask < (int64_t(1) << j); mask++) {
      std::vector<int64_t> c(m, 0);
      c[j] = 1;
      for (int i = 0; i < j; i++) {
        c[i] = (mask >> i) & 1;
      }
      int64_t norm = 0;
      for (int a = 0; a <= j; a++) {
        if (c[a] == 0) {
          continue;
        }
        norm += N[a][a];
        for (int b = a + 1; b <= j; b++) {
          if (c[b] != 0) {
            norm += 2 * N[a][b];
          }
        }
      }
      if (ResMod(norm, 8) == 0) {
        return combine(c);
      }
    }
  }
  return {};
}

// L -> L + Z y / p, y_k made 1 and the k-th basis vector replaced by y / p.
template <typename T>
void GoUp(SymbolLattice<T> &L, std::vector<int64_t> y, int64_t p) {
  int n = L.G.rows();
  int k = 0;
  while (y[k] == 0) {
    k++;
  }
  int64_t inv = ResInv(y[k], p);
  for (auto &x : y) {
    x = ResMul(x, inv, p);
  }
  MyMatrix<T> E = IdentityMat<T>(n);
  T p_T(static_cast<long>(p));
  for (int j = 0; j < n; j++) {
    E(k, j) = T(static_cast<long>(y[j])) / p_T;
  }
  ApplyBasisChange(L, E);
}

/*
  A basis adapted to a Jordan splitting at p: rows P, p-integral and
  invertible mod p, with P S P^T block diagonal, S = G (odd p) or G / 2
  (p = 2), by the pivots of FindMinOrderEntry and Schur complements. Also
  the diagonal entries of P S P^T.
 */
template <typename T, typename Tint>
std::pair<MyMatrix<T>, std::vector<T>> PadicAdaptedBasis(MyMatrix<T> const &G,
                                                         int64_t p) {
  int n = G.rows();
  MyMatrix<T> S = (p == 2) ? MyMatrix<T>(G / T(2)) : G;
  MyMatrix<T> P = IdentityMat<T>(n);
  std::vector<int> active(n);
  for (int i = 0; i < n; i++) {
    active[i] = i;
  }
  while (!active.empty()) {
    int a = active.size();
    MyMatrix<T> Sact(a, a);
    for (int i = 0; i < a; i++) {
      for (int j = 0; j < a; j++) {
        Sact(i, j) = S(active[i], active[j]);
      }
    }
    auto [i0, j0] = FindMinOrderEntry<T, Tint>(Sact, p);
    int gi = active[i0], gj = active[j0];
    std::vector<int> I;
    if (i0 != j0 && p != 2) {
      // e_i <- e_i + e_j makes the diagonal of minimal order
      P.row(gi) += P.row(gj);
      S.row(gi) += S.row(gj);
      S.col(gi) += S.col(gj);
      I = {gi};
    } else if (i0 == j0) {
      I = {gi};
    } else {
      I = {gi, gj};
    }
    int b = I.size();
    MyMatrix<T> Bk(b, b);
    for (int u = 0; u < b; u++) {
      for (int v = 0; v < b; v++) {
        Bk(u, v) = S(I[u], I[v]);
      }
    }
    MyMatrix<T> Binv = Inverse(Bk);
    std::vector<int> rest;
    for (int r : active) {
      if (std::find(I.begin(), I.end(), r) == I.end()) {
        rest.push_back(r);
      }
    }
    std::vector<MyVector<T>> coef;
    for (int r : rest) {
      MyVector<T> sr(b);
      for (int u = 0; u < b; u++) {
        sr(u) = S(r, I[u]);
      }
      coef.push_back(Binv.transpose() * sr);
    }
    for (size_t x = 0; x < rest.size(); x++) {
      int r = rest[x];
      for (int u = 0; u < b; u++) {
        P.row(r) -= coef[x](u) * P.row(I[u]);
      }
      for (size_t y = 0; y < rest.size(); y++) {
        int s = rest[y];
        for (int u = 0; u < b; u++) {
          S(r, s) -= coef[x](u) * S(I[u], s);
        }
      }
    }
    for (int r : rest) {
      for (int u = 0; u < b; u++) {
        S(r, I[u]) = 0;
        S(I[u], r) = 0;
      }
    }
    active = rest;
  }
  std::vector<T> diag(n);
  for (int i = 0; i < n; i++) {
    diag[i] = S(i, i);
  }
  return {P, diag};
}

/*
  Linear forms f mod p on L whose kernels, the index p sublattices, are
  tried in the descent: in coordinates of a Jordan-adapted basis, the
  combinations of at most `weight` basis vectors, with all coefficients for
  p <= 13 and a few more those making a binary combination isotropic.
  Normalized with first non-zero coordinate 1, without repetition.
 */
template <typename T, typename Tint>
std::vector<std::vector<int64_t>> SublatticeForms(MyMatrix<T> const &G,
                                                  int64_t p, int weight) {
  int n = G.rows();
  auto [P, diag] = PadicAdaptedBasis<T, Tint>(G, p);
  MyMatrix<T> Q = Inverse(P);
  std::vector<std::vector<int64_t>> Qp(n, std::vector<int64_t>(n));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Qp[i][j] = ResidueModInt<T, Tint>(Q(i, j), p);
    }
  }
  std::vector<int64_t> lambdas;
  for (int64_t l = 1; l < p && l <= 12; l++) {
    lambdas.push_back(l);
  }
  std::vector<PadicClass> cls;
  for (int i = 0; i < n; i++) {
    cls.push_back(PadicFromRational<T, Tint>(diag[i] == 0 ? T(1) : diag[i], p));
  }
  std::set<std::vector<int64_t>> seen;
  std::vector<std::vector<int64_t>> res;
  auto push = [&](std::vector<std::pair<int, int64_t>> const &c) -> void {
    std::vector<int64_t> f(n, 0);
    for (auto const &[i, l] : c) {
      for (int r = 0; r < n; r++) {
        f[r] = ResMod(f[r] + ResMul(Qp[r][i], l, p), p);
      }
    }
    int k = 0;
    while (k < n && f[k] == 0) {
      k++;
    }
    if (k == n) {
      return;
    }
    int64_t inv = ResInv(f[k], p);
    for (auto &x : f) {
      x = ResMul(x, inv, p);
    }
    if (seen.insert(f).second) {
      res.push_back(f);
    }
  };
  for (int i = 0; i < n; i++) {
    push({{i, 1}});
  }
  if (weight >= 2) {
    for (int i = 0; i < n; i++) {
      for (int j = i + 1; j < n; j++) {
        std::vector<int64_t> ls = lambdas;
        if (p != 2 && diag[i] != 0 && diag[j] != 0 && cls[i].val == cls[j].val) {
          int64_t r = ResMul(p - cls[i].unit, ResInv(cls[j].unit, p), p);
          if (kronecker_symbol<int64_t>(r, p) == 1) {
            ls.push_back(SqrtModPrime(r, p));
          }
        }
        for (int64_t l : ls) {
          push({{i, 1}, {j, l}});
        }
      }
    }
  }
  if (weight >= 4) {
    // p = 2 only: the sums of four basis vectors
    for (int i = 0; i < n; i++) {
      for (int j = i + 1; j < n; j++) {
        for (int k = j + 1; k < n; k++) {
          for (int l = k + 1; l < n; l++) {
            push({{i, 1}, {j, 1}, {k, 1}, {l, 1}});
          }
        }
      }
    }
  }
  if (weight >= 3) {
    std::vector<int64_t> ls3;
    for (int64_t l = 1; l < p && l <= 4; l++) {
      ls3.push_back(l);
    }
    for (int i = 0; i < n; i++) {
      for (int j = i + 1; j < n; j++) {
        for (int k = j + 1; k < n; k++) {
          for (int64_t l1 : ls3) {
            for (int64_t l2 : ls3) {
              push({{i, 1}, {j, l1}, {k, l2}});
            }
          }
        }
      }
    }
  }
  return res;
}

// The kernel of f mod p as a basis change: rows e_i - f_i e_k (i != k) and
// p e_k, with f_k = 1.
template <typename T>
MyMatrix<T> IndexPSublatticeChange(std::vector<int64_t> const &f, int64_t p) {
  int n = f.size();
  int k = 0;
  while (f[k] == 0) {
    k++;
  }
  MyMatrix<T> E = IdentityMat<T>(n);
  for (int i = 0; i < n; i++) {
    if (i != k) {
      E(i, k) = -T(static_cast<long>(f[i]));
    }
  }
  E(k, k) = T(static_cast<long>(p));
  return E;
}

// The largest scale exponent of a Jordan component of G at p: minus the
// least valuation of an entry of G^{-1}.
template <typename T, typename Tint> int MaxScaleAt(MyMatrix<T> const &G, int64_t p) {
  MyMatrix<T> Ginv = Inverse(G);
  int e = 0;
  for (int i = 0; i < G.rows(); i++) {
    for (int j = 0; j < G.cols(); j++) {
      if (Ginv(i, j) != 0) {
        e = std::max(e, -PadicFromRational<T, Tint>(Ginv(i, j), p).val);
      }
    }
  }
  return e;
}

/*
  The local symbol key at p of the integral G; at p = 2 that of 2 G, which
  is even. e is the largest scale exponent of the target: the residues of
  the Gram matrix modulo p^{e+7} determine the Jordan components up to
  that scale and their units, which is much less than the p^{v_p(det)+5}
  of LocalSymbolKeyModular when there are many components (the 4 G of an
  odd genus at p = 2). A G with a component beyond that scale is not of
  the target class, and gets the empty key.
 */
// The precision of SymbolKeyAt for the largest scale exponent e, and the
// modulus p^{prec+5} of LocalSymbolKeyModular, when it fits.
inline std::optional<std::pair<int, int64_t>> SymbolKeyModulus(int64_t p,
                                                               int e) {
  int prec = (p == 2) ? e + 2 : e + 1;
  int64_t M = 1;
  for (int i = 0; i < prec + 5; i++) {
    if (M > (int64_t(1) << 30) / p) {
      return {};
    }
    M *= p;
  }
  return std::make_pair(prec, M);
}

template <typename T, typename Tint>
std::vector<int> SymbolKeyAt(MyMatrix<T> const &G, int64_t p, int e) {
  MyMatrix<T> G2 = (p == 2) ? MyMatrix<T>(T(2) * G) : G;
  std::optional<std::pair<int, int64_t>> mod = SymbolKeyModulus(p, e);
  if (!mod) {
    return LocalSymbolKey<T, Tint>(G2, p);
  }
  auto [prec, M] = *mod;
  int n = G.rows();
  MyMatrix<int64_t> Gm(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Gm(i, j) = ResidueModInt<T, Tint>(G2(i, j), M);
    }
  }
  std::optional<std::vector<int>> key = LocalSymbolKeyModular(Gm, p, prec);
  if (key) {
    return *key;
  }
  return {};
}

/*
  The residues modulo M of the Gram matrix of the kernel of f, from those of
  G, in the basis of IndexPSublatticeChange (f_k = 1): O(n^2) operations in
  place of the rational products.
 */
inline MyMatrix<int64_t> SublatticeGramModular(MyMatrix<int64_t> const &Gm,
                                               std::vector<int64_t> const &f,
                                               int64_t p, int64_t M) {
  int n = Gm.rows();
  int k = 0;
  while (f[k] == 0) {
    k++;
  }
  MyMatrix<int64_t> H(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      int64_t val;
      if (i == k && j == k) {
        val = ResMul(ResMul(p, p, M), Gm(k, k), M);
      } else if (i == k || j == k) {
        int o = (i == k) ? j : i;
        val = ResMul(p, Gm(o, k) - ResMul(f[o], Gm(k, k), M), M);
      } else {
        val = Gm(i, j) - ResMul(f[j], Gm(i, k), M) -
              ResMul(f[i], Gm(k, j), M) +
              ResMul(ResMul(f[i], f[j], M), Gm(k, k), M);
        val = ResMod(val, M);
      }
      H(i, j) = val;
      H(j, i) = val;
    }
  }
  return H;
}

/*
  The last resort of the descent: all the linear forms mod p, first
  non-zero coordinate 1, when there are at most 2^18 of them; otherwise
  2^18 random ones.
 */
inline std::vector<std::vector<int64_t>> AllLinearForms(int n, int64_t p) {
  std::vector<std::vector<int64_t>> res;
  int64_t total = 1;
  bool small = true;
  for (int i = 0; i < n; i++) {
    if (total > (int64_t(1) << 18) / p) {
      small = false;
      break;
    }
    total *= p;
  }
  if (small) {
    for (int64_t idx = 1; idx < total; idx++) {
      std::vector<int64_t> f(n);
      int64_t w = idx;
      for (int i = 0; i < n; i++) {
        f[i] = w % p;
        w /= p;
      }
      int k = 0;
      while (f[k] == 0) {
        k++;
      }
      if (f[k] == 1) {
        res.push_back(f);
      }
    }
    return res;
  }
  std::mt19937_64 gen(1);
  std::uniform_int_distribution<int64_t> dist(0, p - 1);
  for (int iter = 0; iter < (1 << 18); iter++) {
    std::vector<int64_t> f(n);
    int k = -1;
    for (int i = 0; i < n; i++) {
      f[i] = dist(gen);
      if (k == -1 && f[i] != 0) {
        k = i;
      }
    }
    if (k == -1) {
      continue;
    }
    int64_t inv = ResInv(f[k], p);
    for (auto &x : f) {
      x = ResMul(x, inv, p);
    }
    res.push_back(f);
  }
  return res;
}

/*
  The descent at p from L, maximal at p, to a sublattice of p-power index
  locally isometric to Gp. Going up from Gp by isotropic steps gives a
  chain Gp = T_0 < T_1 < ... < T_k of index p, T_k maximal, hence locally
  isometric to L. The descent follows the local symbol keys of that chain
  backwards: from a lattice locally isometric to T_{i+1}, some index p
  sublattice is locally isometric to T_i, and it is looked for among the
  kernels of SublatticeForms, of weight 2, 3 (and 4 at p = 2), and then of
  AllLinearForms. An unguided search over all the keys reachable at each
  level grows too fast (hundreds of keys after six steps at p = 2 in rank
  10).
 */
template <typename T, typename Tint>
SymbolLattice<T> DescendToLocalSymbol(SymbolLattice<T> const &L, int64_t p,
                                      MyMatrix<T> const &Gp,
                                      [[maybe_unused]] std::ostream &os) {
  int n = Gp.rows();
  int v_target = PadicFromRational<T, Tint>(DeterminantMat(Gp), p).val;
  // the chain of keys, from Gp up to a maximal lattice
  std::vector<std::vector<int>> chain;
  std::vector<int> scales;
  {
    SymbolLattice<T> M;
    M.B = IdentityMat<T>(n);
    M.G = Gp;
    int e = MaxScaleAt<T, Tint>(M.G, p);
    chain.push_back(SymbolKeyAt<T, Tint>(M.G, p, e));
    scales.push_back(e);
    while (true) {
      std::optional<std::vector<int64_t>> y =
          IsotropicRadicalVector<T, Tint>(M.G, p);
      if (!y) {
        break;
      }
      GoUp(M, *y, p);
      e = MaxScaleAt<T, Tint>(M.G, p);
      chain.push_back(SymbolKeyAt<T, Tint>(M.G, p, e));
      scales.push_back(e);
    }
  }
  int steps = static_cast<int>(chain.size()) - 1;
  int v0 = v_target - 2 * steps;
  if (PadicFromRational<T, Tint>(DeterminantMat(L.G), p).val != v0 ||
      SymbolKeyAt<T, Tint>(L.G, p, MaxScaleAt<T, Tint>(L.G, p)) !=
          chain.back()) {
    std::cerr << "GSYM: at p=" << p << " the maximal lattice of the rational "
              << "space is not locally isometric to that above the genus\n";
    throw TerminalException{1};
  }
#ifdef DEBUG_GENUS_SYMBOL
  os << "GSYM: p=" << p << " descent in " << steps << " steps\n";
#endif
  SymbolLattice<T> cur = L;
  for (int i = steps - 1; i >= 0; i--) {
    bool found = false;
    // weight 2, 3 (and 4 at p = 2) in an adapted basis, then all forms
    int max_weight = (p == 2) ? 4 : 3;
    for (int weight = 2; weight <= max_weight + 1 && !found; weight++) {
      std::vector<std::vector<int64_t>> forms =
          (weight <= max_weight) ? SublatticeForms<T, Tint>(cur.G, p, weight)
                                 : AllLinearForms(n, p);
#ifdef DEBUG_GENUS_SYMBOL
      os << "GSYM: p=" << p << " step " << i << " weight " << weight << ", "
         << forms.size() << " forms\n";
#endif
      // the keys of the candidates from the residues of the Gram matrix
      std::optional<std::pair<int, int64_t>> mod =
          SymbolKeyModulus(p, scales[i]);
      MyMatrix<int64_t> Gm;
      if (mod) {
        MyMatrix<T> G2 = (p == 2) ? MyMatrix<T>(T(2) * cur.G) : cur.G;
        Gm = MyMatrix<int64_t>(n, n);
        for (int a = 0; a < n; a++) {
          for (int b = 0; b < n; b++) {
            Gm(a, b) = ResidueModInt<T, Tint>(G2(a, b), mod->second);
          }
        }
      }
      for (auto const &f : forms) {
        bool match;
        if (mod) {
          std::optional<std::vector<int>> key = LocalSymbolKeyModular(
              SublatticeGramModular(Gm, f, p, mod->second), p, mod->first);
          match = key && *key == chain[i];
        } else {
          SymbolLattice<T> N = cur;
          ApplyBasisChange(N, IndexPSublatticeChange<T>(f, p));
          match = SymbolKeyAt<T, Tint>(N.G, p, scales[i]) == chain[i];
        }
        if (match) {
          ApplyBasisChange(cur, IndexPSublatticeChange<T>(f, p));
          found = true;
          break;
        }
      }
    }
    if (!found) {
      std::cerr << "GSYM: the descent at p=" << p << " found no sublattice "
                << "in the next local class, " << i << " steps from the end\n";
      throw TerminalException{1};
    }
  }
  return cur;
}

template <typename T, typename Tint>
void LLLReduceSymbolLattice(SymbolLattice<T> &L, std::ostream &os) {
  LLLreduction<T, Tint> red = LLLreducedBasis<T, Tint>(L.G, os);
  ApplyBasisChange(L, UniversalMatrixConversion<T, Tint>(red.Pmat));
}

/*
  A Gram matrix of a lattice in the genus of the symbol s, see the comment
  at the top. T is a field (mpq_class), Tint its ring of integers.
 */
template <typename T, typename Tint>
MyMatrix<T> GenusSymbolRepresentative(GenusSymbol const &s, std::ostream &os) {
  int n = s.rank;
  Tint d = GenusSymbolDeterminant<Tint>(s);
  T d_T = UniversalScalarConversion<T, Tint>(d);
  // 1. the local Gram matrices and their Hasse invariants. The scaling by
  // 4 of an odd genus is a square: it changes neither.
  std::map<int64_t, MyMatrix<T>> local;
  local.emplace(2, TwoAdicLocalGram<T>(s.sym2, n));
  for (auto const &[p, comps] : s.sym_odd) {
    local.emplace(p, OddPrimeLocalGram<T>(p, comps, n));
  }
  std::map<int64_t, int> hasse;
  int hasse_prod = 1;
  for (auto const &[p, Gp] : local) {
    std::vector<T> diag = RationalDiagonalForm(Gp);
    // det G_p / d must be a square of Q_p
    T ratio = DeterminantMat(Gp) / d_T;
    PadicClass c = PadicFromRational<T, Tint>(ratio, p);
    bool square = (c.val % 2 == 0) &&
                  ((p == 2) ? (c.unit == 1) : (kronecker_symbol<int64_t>(c.unit, p) == 1));
    if (!square) {
      std::cerr << "GSYM: the symbol at p=" << p
                << " does not match the determinant " << d << "\n";
      throw TerminalException{1};
    }
    hasse[p] = HasseInvariantDiagonal<T, Tint>(diag, p);
    hasse_prod *= hasse[p];
  }
  if (hasse_prod != 1) {
    std::cerr << "GSYM: the Hasse invariants of the symbol "
              << GenusSymbolString(s)
              << " violate the product formula: no positive definite "
                 "lattice has it (the oddity formula fails)\n";
    throw TerminalException{1};
  }
  // 2. the rational space
  std::vector<T> diag = RationalDiagonalWithInvariants<T, Tint>(n, d, hasse);
#ifdef DEBUG_GENUS_SYMBOL
  os << "GSYM: rational form diag(";
  for (auto const &x : diag) {
    os << " " << x;
  }
  os << " )\n";
#endif
  // 3. an even lattice of V, maximal at every prime of its determinant
  SymbolLattice<T> L;
  L.B = T(2) * IdentityMat<T>(n);
  L.G = ZeroMatrix<T>(n, n);
  for (int i = 0; i < n; i++) {
    L.G(i, i) = T(4) * diag[i];
  }
  {
    T det = DeterminantMat(L.G);
    Tint det_i = UniversalScalarConversion<Tint, T>(det);
    std::vector<int64_t> primes{2};
    for (int64_t p : OddPrimeDivisors<Tint>(det_i)) {
      primes.push_back(p);
    }
    for (int64_t p : primes) {
      while (true) {
        std::optional<std::vector<int64_t>> y =
            IsotropicRadicalVector<T, Tint>(L.G, p);
        if (!y) {
          break;
        }
        GoUp(L, *y, p);
      }
      LLLReduceSymbolLattice<T, Tint>(L, os);
      if (p != 2 && s.sym_odd.count(p) == 0) {
        T det_p = DeterminantMat(L.G);
        if (PadicFromRational<T, Tint>(det_p, p).val != 0) {
          std::cerr << "GSYM: the maximal lattice at p=" << p
                    << " is not unimodular\n";
          throw TerminalException{1};
        }
      }
    }
  }
  // 4. the descents at the primes of 2 d, to 2 L for an odd genus
  for (auto const &[p, Gp_L] : local) {
    MyMatrix<T> Gp = (p == 2 && !s.even) ? MyMatrix<T>(T(4) * Gp_L) : Gp_L;
    L = DescendToLocalSymbol<T, Tint>(L, p, Gp, os);
    LLLReduceSymbolLattice<T, Tint>(L, os);
  }
  if (!s.even) {
    L.G /= T(4);
    L.B /= T(2);
  }
  // 5. the check
  GenusSymbol s_L = GenusSymbolOfGram<T, Tint>(L.G);
  if (s_L != s) {
    std::cerr << "GSYM: the lattice built has the genus "
              << GenusSymbolString(s_L) << ", not " << GenusSymbolString(s)
              << "\n";
    throw TerminalException{1};
  }
  return L.G;
}

// clang-format off
#endif  // SRC_GENERA_GENUSSYMBOL_H_
// clang-format on
