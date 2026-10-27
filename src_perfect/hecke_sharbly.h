// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_HECKE_SHARBLY_H_
#define SRC_PERFECT_HECKE_SHARBLY_H_

// clang-format off
#include "hecke_symbols.h"
#include "hecke_dual.h"
// clang-format on

//
// The sharbly method: Ash-Rudolph reduction at the symbol level
//

/*
  The sharbly method computes the images at the symbol level: the
  well rounded cells whose cones have dimension n, which carry the top
  degree cohomology H^nu, with nu the virtual cohomological dimension.
  The image of a symbol cell F under a coset y is the modular symbol of
  the rows of EXT * y, rewritten as a sum of unimodular symbols by the
  continued fraction algorithm of Ash and Rudolph: for a symbol
  [v_1, ..., v_n] with |det| > 1 and a lattice point w = sum_i q_i v_i,
  0 <= q_i < 1, not all q_i integral,
    [v_1, ..., v_n] = sum_{q_i <> 0} [v_1, .., v_{i-1}, w, v_{i+1}, .., v_n]
  and each term has determinant |q_i| |det| < |det|. The identity is a
  relation of the Steinberg module, so the output is homologous to the
  input in the well rounded quotient complex: the operator induced on the
  homology of the symbol level is the Hecke operator. The level below
  is obtained by filling (see compute_sharbly_level_below) when the
  coverings are available (ComputePartialContractingHomotopy = T). The
  images of the other levels are left out (is_level_computed); the
  homology level sanity checks (the image of a cycle is a cycle and
  projects onto the harmonic space modulo boundaries) verify the result.

  The reducing points w are the lattice points of the box
  {sum_i q_i v_i : |q_i| < 1} minimizing max |q_i| (the largest new
  determinant), and the reduction is averaged over them. The box only
  depends on the symbol, so the reduction Phi is a G-equivariant function
  of the symbol: Phi(W g) = Phi(W) g. This makes the images of all the
  symbols consistent, which the level below relies on; the averaging over
  the coset stabilizer is kept as a safety net and the invariance is
  checked.

  A unimodular symbol is a cell of the unimodular frame orbit iOrb_0 (the
  orbit of the standard frame; GL_n(Z) is transitive on the frames). The
  orientation character of a frame cell on its stabilizer is the sign of
  the permutation of the rays, which is also the sign of the symbol
  relations, so inserting the symbol W as the cell (iOrb_0, E_0^{-1} W)
  is sign-consistent across all insertions.
 */

template <typename T, typename Tint, typename Tgroup>
struct SharblyReducer {
private:
  FullComplexEnumeration<T, Tint, Tgroup> const &fce;
  int symbol_level;
  int iOrb0;
  MyMatrix<Tint> E0_inv;
  // The cell iOrb0 is sign0 times the symbol of the rows of E0.
  int sign0;
  HeckeBudget budget;
  SymbolGeometry<T, Tint> geom;
  // The cells of the symbol level which are not unimodular (for example
  // the frames of determinant sqrt(-2) over Z[sqrt(-2)]): a symbol
  // equivalent to one of them is reduced, it is that cell. The cell is
  // sign times the symbol of the rows of E.
  struct SymbolCell {
    int iOrb;
    MyMatrix<Tint> E;
    MyMatrix<T> R_inv;
    Tint det;
    int sign;
  };
  std::vector<SymbolCell> l_symbol_cell;
  /*
    The reducing points of a symbol of a Hermitian T-space: the classes
    f of Z^{2m} modulo the rows of R(W), whose coordinates q_i are taken
    of the smallest norm modulo O_K (all of them if there are ties), with
    the classes ordered by max N(q_i), then the number of nonzero q_i,
    then sum N(q_i). The points w and u w for u a unit give the same
    symbols: one of them is kept. As for the classic case, the result only
    depends on the symbol, so the reduction is G-equivariant. Returned
    with the indices i of the nonzero q_i.
   */
  std::vector<std::pair<MyVector<Tint>, std::vector<int>>>
  get_reducing_points_hermitian(MyMatrix<Tint> const &W) const {
    using Elt = std::pair<T, T>;
    int m = geom.dim;
    MyMatrix<T> R_T = geom.realify(UniversalMatrixConversion<T, Tint>(W));
    MyMatrix<T> R_inv = Inverse(R_T);
    int n2 = 2 * m;
    std::vector<MyVector<T>> l_gen;
    for (int j = 0; j < n2; j++) {
      MyVector<T> g(n2);
      for (int i = 0; i < n2; i++) {
        g(i) = frac(R_inv(j, i));
      }
      l_gen.push_back(g);
    }
    std::set<MyVector<T>> set_class;
    std::vector<MyVector<T>> l_class;
    MyVector<T> zero = ZeroVector<T>(n2);
    set_class.insert(zero);
    l_class.push_back(zero);
    size_t pos = 0;
    while (pos < l_class.size()) {
      MyVector<T> f = l_class[pos];
      for (auto &g : l_gen) {
        MyVector<T> h(n2);
        for (int i = 0; i < n2; i++) {
          h(i) = frac(f(i) + g(i));
        }
        if (set_class.insert(h).second) {
          l_class.push_back(h);
        }
      }
      pos += 1;
    }
    std::optional<std::tuple<T, int, T>> best;
    std::vector<std::vector<std::vector<Elt>>> l_best;
    for (auto &f : l_class) {
      if (IsZeroVector(f)) {
        continue;
      }
      T score(0);
      T sum_norm(0);
      int n_nonzero = 0;
      std::vector<std::vector<Elt>> l_rep(m);
      for (int i = 0; i < m; i++) {
        Elt x{f(i), f(m + i)};
        if (ImagQuadRing<T>::is_zero(x)) {
          l_rep[i].push_back(x);
          continue;
        }
        n_nonzero += 1;
        for (auto &c : geom.ring->nearest_integers(x)) {
          l_rep[i].push_back({x.first - c.first, x.second - c.second});
        }
        T nrm = geom.ring->norm(l_rep[i][0]);
        sum_norm += nrm;
        if (nrm > score) {
          score = nrm;
        }
      }
      std::tuple<T, int, T> key{score, n_nonzero, sum_norm};
      if (!best || key < *best) {
        best = key;
        l_best.clear();
      }
      if (key == *best) {
        l_best.push_back(l_rep);
      }
    }
    if (!best || std::get<0>(*best) >= 1) {
      throw HeckeFailure{"sharbly: no reducing point of norm < 1, O_K is not norm-Euclidean"};
    }
    std::set<MyVector<T>> set_c;
    for (auto &l_rep : l_best) {
      std::vector<MyVector<T>> l_c{ZeroVector<T>(n2)};
      for (int i = 0; i < m; i++) {
        std::vector<MyVector<T>> l_c_new;
        for (auto &c : l_c) {
          for (auto &q : l_rep[i]) {
            MyVector<T> c2 = c;
            c2(i) = q.first;
            c2(m + i) = q.second;
            l_c_new.push_back(c2);
          }
        }
        l_c = std::move(l_c_new);
      }
      for (auto &c : l_c) {
        std::optional<MyVector<T>> c_can;
        for (auto &u : geom.l_unit) {
          MyVector<T> cu(n2);
          for (int i = 0; i < m; i++) {
            Elt q = geom.ring->mul(u, Elt{c(i), c(m + i)});
            cu(i) = q.first;
            cu(m + i) = q.second;
          }
          if (!c_can || cu < *c_can) {
            c_can = cu;
          }
        }
        set_c.insert(*c_can);
      }
    }
    std::vector<std::pair<MyVector<Tint>, std::vector<int>>> l_ret;
    for (auto &c : set_c) {
      MyVector<T> w_T = R_T.transpose() * c;
      std::vector<int> l_nz;
      for (int i = 0; i < m; i++) {
        if (c(i) != 0 || c(m + i) != 0) {
          l_nz.push_back(i);
        }
      }
      l_ret.push_back({geom.primitive_vector(w_T), l_nz});
    }
    return l_ret;
  }
  // The reducing points as primitive vectors, with the indices of the
  // rows they replace.
  std::vector<std::pair<MyVector<Tint>, std::vector<int>>>
  get_reducing_points_full(MyMatrix<Tint> const &W) const {
    if (geom.is_hermitian()) {
      return get_reducing_points_hermitian(W);
    }
    MyMatrix<T> W_T = UniversalMatrixConversion<T, Tint>(W);
    MyMatrix<T> W_inv = Inverse(W_T);
    std::vector<std::pair<MyVector<Tint>, std::vector<int>>> l_ret;
    for (auto &q : get_reducing_points(W_inv)) {
      MyVector<T> w_T = W_T.transpose() * q;
      std::vector<int> l_nz;
      for (int i = 0; i < q.size(); i++) {
        if (q(i) != 0) {
          l_nz.push_back(i);
        }
      }
      l_ret.push_back({UniversalVectorConversion<Tint, T>(RemoveFractionVector(w_T)), l_nz});
    }
    return l_ret;
  }
  std::ostream &os;
  // x - floor(x), in [0, 1).
  static T frac(T const &x) {
    Tint fl = UniversalFloorScalarInteger<Tint, T>(x);
    return x - T(fl);
  }
  /*
    The reducing points of the symbol of the rows of W: the nonzero lattice
    points w = sum_i q_i v_i with all |q_i| < 1 for which max_i |q_i| is the
    smallest, then the number of nonzero q_i, then sum_i q_i^2. They are returned as the
    vectors q, one for each pair {w, -w}. The lattice points of this box
    only depend on the symbol (the box is invariant under the permutations
    and the sign changes of the rows), so the reduction, which averages
    over all of them, is G-equivariant: Phi(W g) = Phi(W) g.
   */
  std::vector<MyVector<T>> get_reducing_points(MyMatrix<T> const &W_inv) const {
    int n = W_inv.rows();
    // The coefficients modulo 1 of the lattice points form the group
    // Z^n / L(W) of order |det W|, generated by those of the e_j, which are
    // the rows of W^{-1}.
    std::vector<MyVector<T>> l_gen;
    for (int j = 0; j < n; j++) {
      MyVector<T> g(n);
      for (int i = 0; i < n; i++) {
        g(i) = frac(W_inv(j, i));
      }
      l_gen.push_back(g);
    }
    std::set<MyVector<T>> set_class;
    std::vector<MyVector<T>> l_class;
    MyVector<T> zero = ZeroVector<T>(n);
    set_class.insert(zero);
    l_class.push_back(zero);
    size_t pos = 0;
    while (pos < l_class.size()) {
      MyVector<T> f = l_class[pos];
      for (auto &g : l_gen) {
        MyVector<T> h(n);
        for (int i = 0; i < n; i++) {
          h(i) = frac(f(i) + g(i));
        }
        if (set_class.insert(h).second) {
          l_class.push_back(h);
        }
      }
      pos += 1;
    }
    // In a class f, the smallest |q_i| is min(f_i, 1 - f_i). The classes
    // are ordered by max |q_i|, then the number of nonzero q_i, then
    // sum q_i^2 = w (W^T W)^{-1} w^T, all invariant under G.
    T half = T(1) / T(2);
    std::optional<std::tuple<T, int, T>> best;
    std::vector<MyVector<T>> l_best;
    for (auto &f : l_class) {
      if (IsZeroVector(f)) {
        continue;
      }
      T score(0);
      T sum_sqr(0);
      int n_nonzero = 0;
      for (int i = 0; i < n; i++) {
        if (f(i) != 0) {
          n_nonzero += 1;
          T a = (f(i) <= half) ? f(i) : T(1) - f(i);
          sum_sqr += a * a;
          if (a > score) {
            score = a;
          }
        }
      }
      std::tuple<T, int, T> key{score, n_nonzero, sum_sqr};
      if (!best || key < *best) {
        best = key;
        l_best.clear();
      }
      if (key == *best) {
        l_best.push_back(f);
      }
    }
    // The vectors q of the classes: f_i = 1/2 gives the two choices +-1/2.
    std::set<MyVector<T>> set_q;
    for (auto &f : l_best) {
      std::vector<MyVector<T>> l_q{ZeroVector<T>(n)};
      for (int i = 0; i < n; i++) {
        std::vector<MyVector<T>> l_q_new;
        for (auto &q : l_q) {
          auto f_push = [&](T const &val) -> void {
            MyVector<T> q2 = q;
            q2(i) = val;
            l_q_new.push_back(q2);
          };
          if (f(i) == half) {
            f_push(half);
            f_push(-half);
          } else if (f(i) < half) {
            f_push(f(i));
          } else {
            f_push(f(i) - T(1));
          }
        }
        l_q = std::move(l_q_new);
      }
      for (auto &q : l_q) {
        // w and -w give the same symbols: one of them is kept.
        MyVector<T> q_can = q;
        for (int i = 0; i < n; i++) {
          if (q_can(i) != 0) {
            if (q_can(i) < 0) {
              q_can = -q_can;
            }
            break;
          }
        }
        set_q.insert(q_can);
      }
    }
    return std::vector<MyVector<T>>(set_q.begin(), set_q.end());
  }
public:
  SharblyReducer(FullComplexEnumeration<T, Tint, Tgroup> const &_fce, int const &_symbol_level,
                 int const &_iOrb0, MyMatrix<Tint> const &_E0_inv,
                 int const &_sign0, SymbolGeometry<T, Tint> const &_geom, std::ostream &_os)
      : fce(_fce), symbol_level(_symbol_level), iOrb0(_iOrb0), E0_inv(_E0_inv), sign0(_sign0),
        geom(_geom), os(_os) {
    int n_orb = fce.levels[symbol_level].l_faces.size();
    for (int iOrb = 0; iOrb < n_orb; iOrb++) {
      FacePerfectComplex<T, Tint, Tgroup> const &face = fce.levels[symbol_level].l_faces[iOrb];
      if (!face.is_well_rounded) {
        continue;
      }
      MyMatrix<Tint> E = geom.symbol_rows(face.EXT);
      if (E.rows() != geom.dim) {
        continue;
      }
      Tint det = geom.det_norm(E);
      if (det > 1) {
        MyMatrix<T> R_inv = Inverse(geom.realify(UniversalMatrixConversion<T, Tint>(E)));
        int sign = symbol_orientation(face, E, fce.pctdi.LinSpa.ListMat);
        l_symbol_cell.push_back({iOrb, E, R_inv, det, sign});
      }
    }
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: sharbly |l_symbol_cell|=" << l_symbol_cell.size() << " (non unimodular symbol cells)\n";
#endif
  }
  SymbolGeometry<T, Tint> const &get_geometry() const { return geom; }
  int get_symbol_level() const { return symbol_level; }
  /*
    The cell of the symbol level equivalent to the symbol of the rows of W
    (|det| > 1), if any: (iOrb, sign, M) with [W] = sign F_iOrb M. The rows
    W' of W permuted by sigma and scaled by units are tried: E M = W' for
    the M = R(E)^{-1} R(W') which is integral (then of determinant +-1,
    so in G), and [W'] = sign(sigma) [W]. The scaling of all the rows by a
    unit is in G, so the first row is not scaled.
   */
  std::optional<std::tuple<int, int, MyMatrix<Tint>>> match_cell(MyMatrix<Tint> const &W) const {
    if (l_symbol_cell.empty()) {
      return {};
    }
    Tint det = geom.det_norm(W);
    int m = W.rows();
    MyMatrix<T> W_T = UniversalMatrixConversion<T, Tint>(W);
    std::vector<std::pair<T, T>> l_unit = geom.l_unit;
    if (!geom.is_hermitian()) {
      l_unit = {{T(1), T(0)}, {T(-1), T(0)}};
    }
    auto f_scale = [&](std::pair<T, T> const &u, MyVector<T> const &V) -> MyVector<T> {
      if (geom.is_hermitian()) {
        return geom.scale(u, V);
      }
      return u.first * V;
    };
    for (auto &sc : l_symbol_cell) {
      if (sc.det != det) {
        continue;
      }
      std::vector<int> sigma(m);
      for (int i = 0; i < m; i++) {
        sigma[i] = i;
      }
      do {
        int sign_sigma = 1;
        for (int i = 0; i < m; i++) {
          for (int j = i + 1; j < m; j++) {
            if (sigma[i] > sigma[j]) {
              sign_sigma = -sign_sigma;
            }
          }
        }
        // The unit choices for the rows 1, ..., m-1 as a counter.
        size_t n_unit = l_unit.size();
        std::vector<size_t> l_choice(m, 0);
        while (true) {
          MyMatrix<T> Wp(m, W.cols());
          for (int i = 0; i < m; i++) {
            MyVector<T> V = GetMatrixRow(W_T, sigma[i]);
            if (i > 0) {
              V = f_scale(l_unit[l_choice[i]], V);
            }
            AssignMatrixRow(Wp, i, V);
          }
          MyMatrix<T> M_T = sc.R_inv * geom.realify(Wp);
          if (IsIntegralMatrix(M_T)) {
            return std::tuple<int, int, MyMatrix<Tint>>{sc.iOrb, sign_sigma * sc.sign,
                                                         UniversalMatrixConversion<Tint, T>(M_T)};
          }
          int pos = 1;
          while (pos < m) {
            l_choice[pos] += 1;
            if (l_choice[pos] < n_unit) {
              break;
            }
            l_choice[pos] = 0;
            pos += 1;
          }
          if (pos >= m) {
            break;
          }
        }
      } while (std::next_permutation(sigma.begin(), sigma.end()));
    }
    return {};
  }
  // The size of a symbol for the reduction: 1 if it is reduced (a cell of
  // the symbol level), otherwise |det| (|N(det)| for a Hermitian T-space).
  Tint reduction_size(MyMatrix<Tint> const &W) const {
    Tint det = geom.det_norm(W);
    if (det > 1 && match_cell(W)) {
      return Tint(1);
    }
    return det;
  }
  // The reducing points of the symbol of the rows of W (|det| > 1), as
  // primitive integral vectors.
  std::vector<MyVector<Tint>> reducing_points(MyMatrix<Tint> const &W) const {
    std::vector<MyVector<Tint>> l_w;
    for (auto &pair : get_reducing_points_full(W)) {
      l_w.push_back(pair.first);
    }
    return l_w;
  }
  // Make the rows primitive integral vectors (generators of their lines
  // over O_K): the symbol is invariant under the scaling of its rows.
  MyMatrix<Tint> primitive_rows(MyMatrix<T> const &W_T) const {
    return geom.primitive_rows(W_T);
  }
  // The Ash-Rudolph reduction of the symbol of the rows of W, written
  // as a chain of cells of the symbol level. If l_sample is not null, the
  // symbols met and the 1-sharblies [w, v_1, ..., v_n] of the steps are
  // appended to it: they form a filling of the difference between the
  // symbol and its reduction in the sharbly complex.
  std::vector<PerfectFaceEntry<T, Tint>> reduce(MyMatrix<Tint> const &W_input,
                                                std::vector<MyMatrix<Tint>> *l_sample = nullptr) const {
    ChainBuilder<T, Tint, Tgroup> cb(symbol_level, fce, os);
    std::vector<std::pair<T, MyMatrix<Tint>>> l_work;
    l_work.push_back({T(1), W_input});
    size_t n_symbol = 0;
    while (!l_work.empty()) {
      n_symbol += 1;
      if (n_symbol > budget.max_sharbly_symbols) {
        throw HeckeFailure{"sharbly: more than " + std::to_string(budget.max_sharbly_symbols) +
                           " symbols in one reduction"};
      }
      std::pair<T, MyMatrix<Tint>> pair = l_work.back();
      l_work.pop_back();
      T const &coef = pair.first;
      MyMatrix<Tint> const &W = pair.second;
      Tint det = geom.det_norm(W);
      if (det == 0) {
        // A degenerate symbol vanishes in the quotient complex.
        continue;
      }
      if (det == 1) {
        // The symbol [W] = [E0] M is sign0 times the cell iOrb0 mapped by M.
        MyMatrix<Tint> M = E0_inv * geom.realify_int(W);
        cb.f_insert(coef * T(sign0), iOrb0, M);
        continue;
      }
      std::optional<std::tuple<int, int, MyMatrix<Tint>>> opt_cell = match_cell(W);
      if (opt_cell) {
        cb.f_insert(coef * T(std::get<1>(*opt_cell)), std::get<0>(*opt_cell), std::get<2>(*opt_cell));
        continue;
      }
      std::vector<std::pair<MyVector<Tint>, std::vector<int>>> l_q = get_reducing_points_full(W);
      if (l_q.empty()) {
        throw HeckeFailure{"sharbly: no reducing point for a symbol of determinant > 1"};
      }
      int n_q = l_q.size();
      T coef_q = coef / T(n_q);
      int n_row = W.rows();
      if (l_sample) {
        l_sample->push_back(W);
      }
      for (auto &pair : l_q) {
        MyVector<Tint> const &w = pair.first;
        if (l_sample) {
          MyMatrix<Tint> Wext(n_row + 1, W.cols());
          AssignMatrixRow(Wext, 0, w);
          for (int i = 0; i < n_row; i++) {
            AssignMatrixRow(Wext, i + 1, GetMatrixRow(W, i));
          }
          l_sample->push_back(std::move(Wext));
        }
        for (auto &i : pair.second) {
          MyMatrix<Tint> W_new = W;
          AssignMatrixRow(W_new, i, w);
          l_work.push_back({coef_q, std::move(W_new)});
        }
      }
    }
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: sharbly reduce n_symbol=" << n_symbol << " |chain|=" << cb.chain.size() << "\n";
#endif
    return cb.get_faces();
  }
};

// The orbit of the unimodular frames at the symbol level. GL_n(Z) is
// transitive on the frames so there is exactly one orbit; the standard
// frame is a face of the perfect cone of A_n, so the orbit is a cell.
template <typename T, typename Tint, typename Tgroup>
int get_unimodular_frame_orbit(FullComplexEnumeration<T, Tint, Tgroup> const &fce,
                               int const &symbol_level,
                               SymbolGeometry<T, Tint> const &geom,
                               [[maybe_unused]] std::ostream &os) {
  std::vector<int> l_found;
  int n_orb = fce.levels[symbol_level].l_faces.size();
  for (int iOrb = 0; iOrb < n_orb; iOrb++) {
    FacePerfectComplex<T, Tint, Tgroup> const &face = fce.levels[symbol_level].l_faces[iOrb];
    // The EXT of a face holds both v and -v (all the unit multiples for
    // a Hermitian T-space): the frame is one vector per line.
    MyMatrix<Tint> half = geom.symbol_rows(face.EXT);
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: sharbly symbol_level=" << symbol_level << " iOrb=" << iOrb
       << " |half|=" << half.rows() << " wr=" << face.is_well_rounded << "\n";
#endif
    if (face.is_well_rounded && half.rows() == geom.dim) {
      Tint det = geom.det_norm(half);
#ifdef DEBUG_HECKE_OPERATORS
      os << "HECKE: sharbly iOrb=" << iOrb << " det=" << det << "\n";
#endif
      if (T_abs(det) == 1) {
        l_found.push_back(iOrb);
      }
    }
  }
  if (l_found.size() != 1) {
    throw HeckeFailure{"sharbly: found " + std::to_string(l_found.size()) +
                       " orbits of unimodular frames at the symbol level instead of 1"};
  }
  return l_found[0];
}

/*
  The level below the symbol level (degree nu - 1 of the cohomology, the
  1-sharblies). The well rounded quotient complex (the cells modulo the
  cells which are not well rounded) is the compactly supported cochain
  complex of the well rounded retract W. Since Gamma is a duality group
  with dualizing module the Steinberg module, it is exact except at the
  symbol level, where its homology is St. For a cell F of the level below
  and a coset y, the chain R = Psi(d F, y) of the symbol level represents
  [d F] y = 0 in St, so R is a boundary in the quotient complex and the
  image of F is a filling X of it: d X = R modulo the cells which are not
  well rounded. Then d Psi = Psi d on the level below, which is all that
  is needed for the operator on its homology: the extension to the lower
  levels exists by exactness and does not change the induced map.

  The filling is searched by the covering method among the well rounded
  cells, which for the classic T-space are the cells whose vectors span
  R^n. The solution is averaged over the coset stabilizer as for the
  other methods.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<PerfectFaceEntry<T, Tint>>
chain_boundary_well_rounded(int const &index, std::vector<PerfectFaceEntry<T, Tint>> const &chain,
                            FullComplexEnumeration<T, Tint, Tgroup> const &fce, std::ostream &os) {
  std::vector<PerfectFaceEntry<T, Tint>> bnd = chain_boundary(index, chain, fce, os);
  std::vector<PerfectFaceEntry<T, Tint>> bnd_wr;
  for (auto &fe : bnd) {
    if (fce.levels[index + 1].l_faces[fe.iOrb].is_well_rounded) {
      bnd_wr.push_back(fe);
    }
  }
  return bnd_wr;
}

/*
  The perfect cones containing points of the simplex of the rank one forms
  of the rows of P: the points weighted towards a subset of its vertices
  (the vectors of the subset doubled, so their rank one forms have weight
  4), all the subsets for the small simplices, the singletons otherwise.
 */
template <typename T, typename Tint, typename Tgroup>
void locate_simplex_cones(MyMatrix<Tint> const &P, ConeWalker<T, Tint, Tgroup> const &walker,
                          std::vector<TopPerfectCone<Tint>> &l_top) {
  int max_row_all_subsets = 6;
  MyMatrix<T> P_T = UniversalMatrixConversion<T, Tint>(P);
  int n_row = P_T.rows();
  auto f_locate = [&](std::vector<int> const &l_row) -> void {
    MyMatrix<T> Q_T = P_T;
    for (auto &i_row : l_row) {
      for (int i_col = 0; i_col < Q_T.cols(); i_col++) {
        Q_T(i_row, i_col) *= 2;
      }
    }
    l_top.push_back(walker.find_containing_cone_T(Q_T));
  };
  f_locate({});
  if (n_row <= max_row_all_subsets) {
    for (int mask = 1; mask < (1 << n_row) - 1; mask++) {
      std::vector<int> l_row;
      for (int i_row = 0; i_row < n_row; i_row++) {
        if ((mask >> i_row) & 1) {
          l_row.push_back(i_row);
        }
      }
      f_locate(l_row);
    }
  } else {
    for (int i_row = 0; i_row < n_row; i_row++) {
      f_locate({i_row});
    }
  }
}

/*
  The seed of the covering for the filling of the image of F under y: the
  perfect cones containing the points of a filling in the sharbly complex,
  namely the cone F y with its facets and the symbols and the 1-sharblies
  met by the Ash-Rudolph reduction of the facets. The covering made of the
  cones containing the cells of the right hand side alone is too far from
  the filling.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<TopPerfectCone<Tint>>
get_sharbly_covering_seed(int const &index, int const &iOrb, MyMatrix<T> const &y,
                          FullComplexEnumeration<T, Tint, Tgroup> const &fce,
                          ConeWalker<T, Tint, Tgroup> const &walker,
                          SharblyReducer<T, Tint, Tgroup> const &reducer) {
  std::vector<MyMatrix<Tint>> l_sample;
  MyMatrix<Tint> const &EXT_F = fce.levels[index].l_faces[iOrb].EXT;
  SymbolGeometry<T, Tint> const &geom = reducer.get_geometry();
  l_sample.push_back(reducer.primitive_rows(UniversalMatrixConversion<T, Tint>(geom.symbol_rows(EXT_F)) * y));
  for (auto &ebnd : fce.boundaries[index].ll_bound[iOrb].l_bound) {
    FacePerfectComplex<T, Tint, Tgroup> const &face_G = fce.levels[index + 1].l_faces[ebnd.iOrb];
    if (!face_G.is_well_rounded) {
      continue;
    }
    MyMatrix<Tint> EXT_G = geom.symbol_rows(face_G.EXT) * ebnd.M;
    MyMatrix<Tint> W = reducer.primitive_rows(UniversalMatrixConversion<T, Tint>(EXT_G) * y);
    l_sample.push_back(W);
    (void)reducer.reduce(W, &l_sample);
  }
  std::vector<TopPerfectCone<Tint>> l_top;
  for (auto &P : l_sample) {
    locate_simplex_cones(P, walker, l_top);
  }
  return l_top;
}

/*
  The seed of the covering at the levels further below the symbol level
  (the k-sharblies, k >= 2): the cone F y and its facets. The right hand
  side is the image of the boundary of F by the chain map of the level
  above, itself a filling, and partial_contracting_homotopy_kernel adds
  the cones of its cells. There is no reduction of the k-sharblies to
  follow: the sharbly complex is not used below the level s - 1.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<TopPerfectCone<Tint>>
get_sharbly_covering_seed_lower(int const &index, int const &iOrb, MyMatrix<T> const &y,
                                FullComplexEnumeration<T, Tint, Tgroup> const &fce,
                                ConeWalker<T, Tint, Tgroup> const &walker,
                                SharblyReducer<T, Tint, Tgroup> const &reducer) {
  std::vector<MyMatrix<Tint>> l_sample;
  MyMatrix<Tint> const &EXT_F = fce.levels[index].l_faces[iOrb].EXT;
  SymbolGeometry<T, Tint> const &geom = reducer.get_geometry();
  l_sample.push_back(reducer.primitive_rows(UniversalMatrixConversion<T, Tint>(geom.symbol_rows(EXT_F)) * y));
  for (auto &ebnd : fce.boundaries[index].ll_bound[iOrb].l_bound) {
    FacePerfectComplex<T, Tint, Tgroup> const &face_G = fce.levels[index + 1].l_faces[ebnd.iOrb];
    if (!face_G.is_well_rounded) {
      continue;
    }
    MyMatrix<Tint> EXT_G = geom.symbol_rows(face_G.EXT) * ebnd.M;
    l_sample.push_back(reducer.primitive_rows(UniversalMatrixConversion<T, Tint>(EXT_G) * y));
  }
  std::vector<TopPerfectCone<Tint>> l_top;
  for (auto &P : l_sample) {
    locate_simplex_cones(P, walker, l_top);
  }
  return l_top;
}

/*
  The local filling of the image of a cell F of the level below the symbol
  level under a coset y (the 1-sharbly reduction).

  The 1-sharblies D = [u_0, ..., u_n] (n+1 vectors spanning R^n) are
  mapped to the chains Phi(dD) of the symbol level, with Phi the reduction
  of the symbols and dD = sum_k (-1)^k [u_0, .., u_k^, .., u_n]. Two
  identities are used:
  - The right hand side R of F y is sum_S c_S Phi(dD_S) for the
    1-sharblies D_S made of n+1 of the vectors of F y; the coefficients
    c_S are found by a small exact system, which also settles the
    orientations.
  - For any vector w, D = sum_k D[u_k <- w] + d[w, u_0, ..., u_n], so
    Phi(dD) = sum_k Phi(dD[u_k <- w]): a 1-sharbly can be split.
  A filling of R is then the sum of fillings of the pieces. A piece whose
  n-subsets have |det| <= 1 (its faces are reduced symbols, so Phi leaves
  them in place) is filled on a covering seeded by the cones met by its
  simplex. Otherwise it is split along its worst face sigma with the
  reducing points of sigma used by Phi: this follows the reduction of
  sigma by Phi one degree up, as the lifts of Gunnells do. The splitting
  is not equivariant, which does not matter since the result is averaged
  afterwards. Its termination is not proven (as for Gunnells), so the
  depth and the number of pieces are bounded, a piece at the maximal
  depth being filled directly.
 */
template <typename T, typename Tint, typename Tgroup>
struct SharblyLocalFiller {
private:
  int index;
  FullComplexEnumeration<T, Tint, Tgroup> const &fce;
  ConeWalker<T, Tint, Tgroup> const &walker;
  SharblyReducer<T, Tint, Tgroup> const &reducer;
  std::function<bool(MyMatrix<Tint> const &)> f_filter;
  CoveringOptions options;
  std::ostream &os;
  // The pieces are filled from the depth max_depth, and split further up to
  // max_depth_split when their covering is beyond the budget.
  int max_depth = 1;
  int max_depth_split = 3;
  size_t max_piece = 2000;
  size_t n_piece = 0;
  std::unordered_map<MyMatrix<Tint>, std::vector<PerfectFaceEntry<T, Tint>>> cache_phi;
  // The rows of D except the row k.
  static MyMatrix<Tint> remove_row(MyMatrix<Tint> const &D, int const &k) {
    MyMatrix<Tint> W(D.rows() - 1, D.cols());
    int pos = 0;
    for (int i = 0; i < D.rows(); i++) {
      if (i != k) {
        AssignMatrixRow(W, pos, GetMatrixRow(D, i));
        pos += 1;
      }
    }
    return W;
  }
  std::vector<PerfectFaceEntry<T, Tint>> const &phi(MyMatrix<Tint> const &W) {
    auto iter = cache_phi.find(W);
    if (iter != cache_phi.end()) {
      return iter->second;
    }
    std::vector<PerfectFaceEntry<T, Tint>> chain = reducer.reduce(W);
    return cache_phi.emplace(W, std::move(chain)).first->second;
  }
  // The largest |det| of the n-subsets of the rows of D, and their sum.
  std::pair<Tint, Tint> det_size(MyMatrix<Tint> const &D) const {
    Tint max_val(0), sum_val(0);
    for (int k = 0; k < D.rows(); k++) {
      Tint det = reducer.reduction_size(remove_row(D, k));
      sum_val += det;
      if (det > max_val) {
        max_val = det;
      }
    }
    return {max_val, sum_val};
  }
  std::vector<PerfectFaceEntry<T, Tint>> fill_piece(MyMatrix<Tint> const &D, int const &depth) {
    n_piece += 1;
    if (n_piece > max_piece) {
      throw HeckeFailure{"sharbly local filling: more than " + std::to_string(max_piece) + " pieces"};
    }
    std::vector<PerfectFaceEntry<T, Tint>> rhs = phi_boundary(D);
    if (rhs.empty()) {
      return {};
    }
    MyMatrix<T> D_T = UniversalMatrixConversion<T, Tint>(D);
    TopPerfectCone<Tint> tpc = walker.find_containing_cone_T(D_T);
    std::pair<Tint, Tint> size_D = det_size(D);
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: local piece depth=" << depth << " |D|=" << D.rows() << " max_det=" << size_D.first
       << " sum_det=" << size_D.second << " |rhs|=" << rhs.size() << "\n";
#endif
    if (size_D.first <= 1 || depth >= max_depth) {
      std::vector<TopPerfectCone<Tint>> l_top_seed{tpc};
      locate_simplex_cones(D, walker, l_top_seed);
      try {
        return partial_contracting_homotopy_kernel(index + 1, rhs, fce, f_filter, l_top_seed, options, os);
      } catch (CoveringFailure const &e) {
        // Too large for the budget of a piece: it is split further.
        if (size_D.first <= 1 || depth >= max_depth_split) {
          throw HeckeFailure{"sharbly local filling of a piece: " + e.message};
        }
      }
    }
    // The split along the worst face sigma = D minus u_j: with the reducing
    // points w of sigma used by Phi (averaged as Phi does),
    //   D = (1/|W|) sum_w sum_k D[u_k <- w]   modulo boundaries.
    // The pieces D[u_j <- w] have the faces sigma and its children in the
    // reduction of Phi, so their Phi(d .) sum to zero over w: they are
    // dropped (this is checked). The other pieces hold the children of
    // sigma, of smaller determinants, and are split further.
    int n_row = D.rows();
    int j_worst = -1;
    Tint det_worst(0);
    for (int k = 0; k < n_row; k++) {
      Tint det = reducer.reduction_size(remove_row(D, k));
      if (det > det_worst) {
        det_worst = det;
        j_worst = k;
      }
    }
    std::vector<MyVector<Tint>> l_w = reducer.reducing_points(remove_row(D, j_worst));
    if (l_w.empty()) {
      throw HeckeFailure{"sharbly local filling: no reducing point for a face"};
    }
    T coef_w = T(1) / T(static_cast<int>(l_w.size()));
    std::vector<MyMatrix<Tint>> l_piece_j;
    {
      ChainBuilder<T, Tint, Tgroup> cb_check(index + 1, fce, os);
      for (auto &w : l_w) {
        MyMatrix<Tint> Dj = D;
        AssignMatrixRow(Dj, j_worst, w);
        l_piece_j.push_back(Dj);
        for (auto &fe : phi_boundary(Dj)) {
          cb_check.f_insert(fe.value, fe.iOrb, fe.M);
        }
      }
      if (!cb_check.is_zero_chain()) {
        // Not expected; the pieces are then kept.
#ifdef DEBUG_HECKE_OPERATORS
        os << "HECKE: local piece, the pieces through the worst face do not cancel\n";
#endif
      } else {
        l_piece_j.clear();
      }
    }
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: local piece split along the face " << j_worst << " of det " << det_worst
       << " with " << l_w.size() << " reducing points\n";
#endif
    ChainBuilder<T, Tint, Tgroup> cb(index, fce, os);
    auto f_add = [&](MyMatrix<Tint> const &Dk) -> void {
      if (!reducer.get_geometry().is_spanning(Dk)) {
        return;
      }
      for (auto &fe : fill_piece(Dk, depth + 1)) {
        cb.f_insert(coef_w * fe.value, fe.iOrb, fe.M);
      }
    };
    for (auto &w : l_w) {
      for (int k = 0; k < n_row; k++) {
        if (k == j_worst) {
          continue;
        }
        MyMatrix<Tint> Dk = D;
        AssignMatrixRow(Dk, k, w);
        f_add(Dk);
      }
    }
    for (auto &Dj : l_piece_j) {
      f_add(Dj);
    }
    return cb.get_faces();
  }
public:
  SharblyLocalFiller(int const &_index, FullComplexEnumeration<T, Tint, Tgroup> const &_fce,
                     ConeWalker<T, Tint, Tgroup> const &_walker,
                     SharblyReducer<T, Tint, Tgroup> const &_reducer,
                     std::function<bool(MyMatrix<Tint> const &)> const &_f_filter,
                     CoveringOptions const &_options, std::ostream &_os)
      : index(_index), fce(_fce), walker(_walker), reducer(_reducer), f_filter(_f_filter),
        options(_options), os(_os) {
  }
  size_t get_n_piece() const { return n_piece; }
  // Phi(dD) for the (n+1) x n matrix D.
  std::vector<PerfectFaceEntry<T, Tint>> phi_boundary(MyMatrix<Tint> const &D) {
    ChainBuilder<T, Tint, Tgroup> cb(index + 1, fce, os);
    for (int k = 0; k < D.rows(); k++) {
      T sign = (k % 2 == 0) ? T(1) : T(-1);
      for (auto &fe : phi(remove_row(D, k))) {
        cb.f_insert(sign * fe.value, fe.iOrb, fe.M);
      }
    }
    return cb.get_faces();
  }
  // A filling of rhs, the image of the boundary of the cell whose vectors
  // (one per pair {v, -v}) mapped by y are the rows of V. Nothing if rhs
  // is not a combination of the Phi(dD_S).
  std::optional<std::vector<PerfectFaceEntry<T, Tint>>>
  fill(MyMatrix<Tint> const &V, std::vector<PerfectFaceEntry<T, Tint>> const &rhs) {
    n_piece = 0;
    int n = reducer.get_geometry().dim;
    int m = V.rows();
    // The (n+1)-subsets of the rows spanning R^n.
    std::vector<MyMatrix<Tint>> l_D;
    std::vector<int> subset(n + 1);
    std::function<void(int, int)> f_enum = [&](int start, int depth) -> void {
      if (depth == n + 1) {
        MyMatrix<Tint> D(n + 1, V.cols());
        for (int k = 0; k <= n; k++) {
          AssignMatrixRow(D, k, GetMatrixRow(V, subset[k]));
        }
        if (reducer.get_geometry().is_spanning(D)) {
          l_D.push_back(D);
        }
        return;
      }
      for (int i = start; i < m; i++) {
        subset[depth] = i;
        f_enum(i + 1, depth + 1);
      }
    };
    f_enum(0, 0);
    // The system sum_S c_S Phi(dD_S) = rhs, on the cells met.
    ComplexBuilder<T, Tint, Tgroup> cb2(index + 1, fce, os);
    std::vector<std::vector<PerfectFaceEntry<T, Tint>>> l_bnd;
    for (auto &D : l_D) {
      std::vector<PerfectFaceEntry<T, Tint>> bnd = phi_boundary(D);
      for (auto &fe : bnd) {
        cb2.f_insert(fe.iOrb, fe.M);
      }
      l_bnd.emplace_back(std::move(bnd));
    }
    for (auto &fe : rhs) {
      cb2.f_insert(fe.iOrb, fe.M);
    }
    int n_D = l_D.size();
    int n_cell = cb2.dimension();
    MyMatrix<T> A = ZeroMatrix<T>(n_D, n_cell);
    for (int i_D = 0; i_D < n_D; i_D++) {
      for (auto &fe : l_bnd[i_D]) {
        std::pair<size_t, int> p = cb2.get_position(fe.iOrb, fe.M);
        A(i_D, p.first) += fe.value * T(p.second);
      }
    }
    MyVector<T> b = ZeroVector<T>(n_cell);
    for (auto &fe : rhs) {
      std::pair<size_t, int> p = cb2.get_position(fe.iOrb, fe.M);
      b(p.first) += fe.value * T(p.second);
    }
    std::optional<MyVector<T>> opt = SolutionMat(A, b);
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: local fill |V|=" << m << " n_D=" << n_D << " n_cell=" << n_cell
       << " solvable=" << (opt ? 1 : 0) << "\n";
    if (!opt) {
      os << "HECKE: local fill V=\n";
      WriteMatrix(os, V);
      os << "HECKE: local fill b=" << StringVectorGAP(b) << "\n";
      for (int i_D = 0; i_D < n_D; i_D++) {
        os << "HECKE: local fill A(" << i_D << ")=" << StringVectorGAP(MyVector<T>(GetMatrixRow(A, i_D))) << "\n";
        for (int k = 0; k < l_D[i_D].rows(); k++) {
          MyMatrix<Tint> W = remove_row(l_D[i_D], k);
          os << "HECKE: local fill   k=" << k << " det_norm=" << reducer.get_geometry().det_norm(W)
             << " |phi|=" << phi(W).size() << "\n";
        }
      }
    }
#endif
    if (!opt) {
      return {};
    }
    ChainBuilder<T, Tint, Tgroup> cb(index, fce, os);
    for (int i_D = 0; i_D < n_D; i_D++) {
      T c = (*opt)(i_D);
      if (c == 0) {
        continue;
      }
      for (auto &fe : fill_piece(l_D[i_D], 0)) {
        cb.f_insert(c * fe.value, fe.iOrb, fe.M);
      }
    }
    return cb.get_faces();
  }
};

template <typename T, typename Tint, typename Tgroup>
void compute_sharbly_level_below(
    int const &index, FullComplexEnumeration<T, Tint, Tgroup> const &fce,
    HeckeCosets<T> const &cosets, ConeWalker<T, Tint, Tgroup> const &walker,
    SharblyReducer<T, Tint, Tgroup> const &reducer,
    std::vector<std::vector<std::vector<std::vector<PerfectFaceEntry<T, Tint>>>>> &images,
    std::ostream &os) {
  LinSpaceMatrix<T> const &LinSpa = fce.pctdi.LinSpa;
  int n = LinSpa.n;
  size_t n_coset = cosets.size();
  SymbolGeometry<T, Tint> const &geom = reducer.get_geometry();
  std::function<bool(MyMatrix<Tint> const &)> f_filter =
      [&geom](MyMatrix<Tint> const &EXT) -> bool {
    return geom.is_spanning(EXT);
  };
#ifdef TIMINGS_HECKE_OPERATORS
  MicrosecondTime time;
#endif
  HeckeBudget budget;
  // The level s - 1 (the 1-sharblies) has the seeds from the reduction of
  // the symbols and the local filling; the levels further below only the
  // fillings on the coverings.
  bool is_first_below = (index == reducer.get_symbol_level() - 1);
  auto f_seed = [&](int const &iOrb, MyMatrix<T> const &y) -> std::vector<TopPerfectCone<Tint>> {
    if (is_first_below) {
      return get_sharbly_covering_seed(index, iOrb, y, fce, walker, reducer);
    }
    return get_sharbly_covering_seed_lower(index, iOrb, y, fce, walker, reducer);
  };
  int n_orb = fce.levels[index].l_faces.size();
  for (int iOrb = 0; iOrb < n_orb; iOrb++) {
    FacePerfectComplex<T, Tint, Tgroup> const &face = fce.levels[index].l_faces[iOrb];
    if (!face.is_well_rounded) {
      continue;
    }
    std::vector<std::vector<std::pair<int, MyMatrix<Tint>>>> act =
        get_stab_coset_action(face.l_gens, cosets);
    std::vector<int> status(n_coset, 0);
    for (size_t j = 0; j < n_coset; j++) {
      if (status[j] == 1) {
        continue;
      }
      std::vector<PerfectFaceEntry<T, Tint>> rhs =
          get_hecke_boundary_image(index, iOrb, j, fce, cosets, images[index + 1], os);
      std::vector<PerfectFaceEntry<T, Tint>> sol;
      if (!rhs.empty()) {
#ifdef SANITY_CHECK_HECKE_OPERATORS
        for (auto &fe : rhs) {
          MyMatrix<Tint> EXT = fce.levels[index + 1].l_faces[fe.iOrb].EXT * fe.M;
          if (!f_filter(EXT) || !fce.levels[index + 1].l_faces[fe.iOrb].is_well_rounded) {
            std::cerr << "HECKE: sharbly, the right hand side should be made of well rounded cells\n";
            throw TerminalException{1};
          }
        }
#endif
        CoveringOptions options;
        options.guided_growth = true;
        options.clean_failure = true;
        options.max_entries_exact = budget.solve.max_entries_exact;
        options.l1_first = true;
        options.max_l1_iter = budget.solve.max_yall1_iter;
        // Three attempts: the filling of the whole right hand side on a
        // small covering, accepted if it is not bloated; the local filling
        // (the 1-sharbly reduction) with the same small coverings for the
        // pieces; the filling of the whole right hand side on a covering of
        // any size.
        CoveringOptions options_small = options;
        options_small.max_top = 100;
        options_small.max_unknowns = 50000;
        options_small.max_entries_exact = 500000;
        std::optional<std::vector<PerfectFaceEntry<T, Tint>>> opt_sol;
        {
          std::vector<TopPerfectCone<Tint>> l_top_seed = f_seed(iOrb, cosets.l_coset[j]);
          try {
            std::vector<PerfectFaceEntry<T, Tint>> sol_small =
                partial_contracting_homotopy_kernel(index + 1, rhs, fce, f_filter, l_top_seed, options_small, os);
            if (sol_small.size() <= 20 * rhs.size() + 50) {
              opt_sol = sol_small;
            }
          } catch (CoveringFailure const &e) {
#ifdef DEBUG_HECKE_OPERATORS
            os << "HECKE: sharbly small global filling failed: " << e.message << "\n";
#endif
          }
#ifdef TIMINGS_HECKE_OPERATORS
          os << "|HECKE: sharbly small global filling iOrb=" << iOrb << " j=" << j
             << " success=" << (opt_sol ? 1 : 0) << "|=" << time << "\n";
#endif
        }
        if (!opt_sol && is_first_below) {
          SharblyLocalFiller<T, Tint, Tgroup> filler(index, fce, walker, reducer, f_filter, options_small, os);
          MyMatrix<Tint> V = reducer.primitive_rows(
              UniversalMatrixConversion<T, Tint>(geom.symbol_rows(face.EXT)) * cosets.l_coset[j]);
          try {
            opt_sol = filler.fill(V, rhs);
          } catch (HeckeFailure const &e) {
#ifdef DEBUG_HECKE_OPERATORS
            os << "HECKE: sharbly local filling failed: " << e.message << "\n";
#endif
          }
#ifdef TIMINGS_HECKE_OPERATORS
          os << "|HECKE: sharbly local filling iOrb=" << iOrb << " j=" << j
             << " success=" << (opt_sol ? 1 : 0) << " n_piece=" << filler.get_n_piece() << "|=" << time << "\n";
#endif
        }
        if (opt_sol) {
          sol = *opt_sol;
        } else {
          std::vector<TopPerfectCone<Tint>> l_top_seed = f_seed(iOrb, cosets.l_coset[j]);
          try {
            sol = partial_contracting_homotopy_kernel(index + 1, rhs, fce, f_filter, l_top_seed, options, os);
          } catch (CoveringFailure const &e) {
            throw HeckeFailure{"sharbly, filling of the cell " + std::to_string(iOrb) + " of level " +
                               std::to_string(index) + " for the coset " + std::to_string(j) + ": " + e.message};
          }
        }
      }
#ifdef TIMINGS_HECKE_OPERATORS
      os << "|HECKE: sharbly level below iOrb=" << iOrb << "/" << n_orb << " j=" << j << "/" << n_coset
         << " |rhs|=" << rhs.size() << " |sol|=" << sol.size() << "|=" << time << "\n";
#endif
      CosetOrbit<Tint> co = get_coset_orbit(j, face.l_gens, act, n_coset, n);
      std::vector<std::pair<MyMatrix<Tint>, MyMatrix<Tint>>> l_schreier =
          get_coset_stabilizer_schreier(co, face.l_gens, act, n);
#ifdef SANITY_CHECK_HECKE_OPERATORS
      if (!is_chain_coset_stabilizer_invariant(index, iOrb, index + 1, rhs, l_schreier, fce, os)) {
        std::cerr << "HECKE: sharbly, the right hand side should be invariant under the coset stabilizer\n";
        throw TerminalException{1};
      }
#endif
      std::vector<PerfectFaceEntry<T, Tint>> sol_avg =
          average_chain_coset_stabilizer(index, iOrb, sol, l_schreier, fce, os);
#ifdef SANITY_CHECK_HECKE_OPERATORS
      {
        std::vector<PerfectFaceEntry<T, Tint>> bnd = chain_boundary_well_rounded(index, sol_avg, fce, os);
        if (!is_equal_chain(bnd, rhs, index + 1, fce, os)) {
          std::cerr << "HECKE: sharbly, the averaged chain is not a filling\n";
          throw TerminalException{1};
        }
        if (!is_chain_coset_stabilizer_invariant(index, iOrb, index, sol_avg, l_schreier, fce, os)) {
          std::cerr << "HECKE: sharbly, the averaged filling should be invariant\n";
          throw TerminalException{1};
        }
      }
#endif
      size_t n_orbit = co.l_orbit.size();
      for (size_t pos = 0; pos < n_orbit; pos++) {
        int j2 = co.l_orbit[pos];
        int sign = get_face_orientation(face.EXT, LinSpa.ListMat, face.or_info, co.transversal[pos]);
        MyMatrix<Tint> Minv = Inverse(co.transversal_conj[pos]);
        images[index][iOrb][j2] = chain_right_multiply(sol_avg, T(sign), Minv);
        status[j2] = 1;
      }
    }
#ifdef SANITY_CHECK_HECKE_OPERATORS
    // The chain map property on all the cosets, after the propagation.
    for (size_t j = 0; j < n_coset; j++) {
      std::vector<PerfectFaceEntry<T, Tint>> rhs =
          get_hecke_boundary_image(index, iOrb, j, fce, cosets, images[index + 1], os);
      std::vector<PerfectFaceEntry<T, Tint>> bnd =
          chain_boundary_well_rounded(index, images[index][iOrb][j], fce, os);
      if (!is_equal_chain(bnd, rhs, index + 1, fce, os)) {
        std::cerr << "HECKE: sharbly, the chain map property fails at index=" << index
                  << " iOrb=" << iOrb << " j=" << j << "\n";
        throw TerminalException{1};
      }
    }
#endif
  }
}

// The reducer of the symbols of the complex: its symbol level and its
// unimodular frame.
template <typename T, typename Tint, typename Tgroup>
SharblyReducer<T, Tint, Tgroup> get_sharbly_reducer(FullComplexEnumeration<T, Tint, Tgroup> const &fce,
                                                    SymbolGeometry<T, Tint> const &geom, std::ostream &os) {
  int n_levels = fce.levels.size();
  int symbol_level = get_symbol_level(fce, geom);
  if (symbol_level < 0 || symbol_level >= n_levels) {
    std::cerr << "HECKE: sharbly, the symbol level " << symbol_level
              << " is not in the complex (n_levels=" << n_levels << ")\n";
    throw TerminalException{1};
  }
  int iOrb0 = get_unimodular_frame_orbit(fce, symbol_level, geom, os);
  MyMatrix<Tint> E0 = geom.symbol_rows(fce.levels[symbol_level].l_faces[iOrb0].EXT);
  MyMatrix<Tint> E0_inv = Inverse(geom.realify_int(E0));
  int sign0 = symbol_orientation(fce.levels[symbol_level].l_faces[iOrb0], E0, fce.pctdi.LinSpa.ListMat);
  return SharblyReducer<T, Tint, Tgroup>(fce, symbol_level, iOrb0, E0_inv, sign0, geom, os);
}

template <typename T, typename Tint, typename Tgroup>
HeckeChainMap<T, Tint>
compute_hecke_chain_map_sharbly(FullComplexEnumeration<T, Tint, Tgroup> const &fce,
                                MyMatrix<T> const &x, std::ostream &os) {
  LinSpaceMatrix<T> const &LinSpa = fce.pctdi.LinSpa;
  int n = LinSpa.n;
  PerfectComplexOptions const &pco = fce.pctdi.pco;
  if (!pco.compute_boundary) {
    std::cerr << "HECKE: The sharbly method requires the boundary, please set ComputeBoundary = T\n";
    throw TerminalException{1};
  }
  if (!is_stab_space(x, LinSpa)) {
    std::cerr << "HECKE: The matrix x should preserve the T-space\n";
    throw TerminalException{1};
  }
#ifdef TIMINGS_HECKE_OPERATORS
  MicrosecondTime time;
#endif
  std::vector<MyMatrix<Tint>> l_gens_G = get_ambient_group_generators(fce, os);
  std::function<bool(MyMatrix<T> const &)> f_member =
      [&LinSpa](MyMatrix<T> const &g) -> bool {
    return is_in_ambient_group(LinSpa, g);
  };
  HeckeCosets<T> cosets = enumerate_hecke_cosets<T, Tint>(x, l_gens_G, f_member, os);
  size_t n_coset = cosets.size();
#ifdef TIMINGS_HECKE_OPERATORS
  os << "|HECKE: enumerate_hecke_cosets|=" << time << "\n";
#endif
  int n_levels = fce.levels.size();
  SymbolGeometry<T, Tint> geom(LinSpa);
#ifdef DEBUG_HECKE_OPERATORS
  os << "HECKE: sharbly is_hermitian=" << geom.is_hermitian() << " dim=" << geom.dim
     << " |units|=" << geom.l_unit.size() << "\n";
  // The structure of the complex: per level, the orbits, the well rounded
  // ones and the sizes of the half vector sets of the cells.
  for (int index = 0; index < n_levels; index++) {
    std::map<int, int> map_size;
    int n_wr = 0;
    for (auto &face : fce.levels[index].l_faces) {
      if (face.is_well_rounded) {
        n_wr += 1;
        map_size[geom.symbol_rows(face.EXT).rows()] += 1;
      }
    }
    os << "HECKE: sharbly level " << index << " n_orb=" << fce.levels[index].l_faces.size()
       << " n_wr=" << n_wr << " |half| of the well rounded:";
    for (auto &kv : map_size) {
      os << " " << kv.first << "^" << kv.second;
    }
    os << "\n";
  }
#endif
  SharblyReducer<T, Tint, Tgroup> reducer = get_sharbly_reducer(fce, geom, os);
  int symbol_level = reducer.get_symbol_level();
  ConeWalker<T, Tint, Tgroup> walker(fce, os);
  std::vector<std::vector<std::vector<std::vector<PerfectFaceEntry<T, Tint>>>>> images(n_levels);
  std::vector<uint8_t> l_level_computed(n_levels, 0);
  l_level_computed[symbol_level] = 1;
  for (int index = 0; index < n_levels; index++) {
    int n_orb = fce.levels[index].l_faces.size();
    images[index] = std::vector<std::vector<std::vector<PerfectFaceEntry<T, Tint>>>>(
        n_orb, std::vector<std::vector<PerfectFaceEntry<T, Tint>>>(n_coset));
  }
  int n_orb = fce.levels[symbol_level].l_faces.size();
  for (int iOrb = 0; iOrb < n_orb; iOrb++) {
    FacePerfectComplex<T, Tint, Tgroup> const &face = fce.levels[symbol_level].l_faces[iOrb];
    if (!face.is_well_rounded) {
      continue;
    }
    MyMatrix<Tint> half = geom.symbol_rows(face.EXT);
    if (half.rows() != geom.dim) {
      throw HeckeFailure{"sharbly: the cell " + std::to_string(iOrb) +
                         " of the symbol level has " + std::to_string(half.rows()) +
                         " vectors instead of " + std::to_string(geom.dim) +
                         ": not a symbol, use another method"};
    }
    MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tint>(half);
    // The cell is sign_F times the symbol of the rows of half.
    int sign_F = symbol_orientation(face, half, LinSpa.ListMat);
    std::vector<std::vector<std::pair<int, MyMatrix<Tint>>>> act =
        get_stab_coset_action(face.l_gens, cosets);
    std::vector<int> status(n_coset, 0);
    for (size_t j = 0; j < n_coset; j++) {
      if (status[j] == 1) {
        continue;
      }
#ifdef DEBUG_HECKE_OPERATORS
      os << "HECKE: sharbly iOrb=" << iOrb << "/" << n_orb << " j=" << j << "/" << n_coset << "\n";
#endif
      MyMatrix<T> W_T = EXT_T * cosets.l_coset[j];
      MyMatrix<Tint> W = reducer.primitive_rows(W_T);
      std::vector<PerfectFaceEntry<T, Tint>> sol = chain_right_multiply(reducer.reduce(W), T(sign_F), IdentityMat<Tint>(n));
      CosetOrbit<Tint> co = get_coset_orbit(j, face.l_gens, act, n_coset, n);
      std::vector<std::pair<MyMatrix<Tint>, MyMatrix<Tint>>> l_schreier =
          get_coset_stabilizer_schreier(co, face.l_gens, act, n);
      std::vector<PerfectFaceEntry<T, Tint>> sol_avg =
          average_chain_coset_stabilizer(symbol_level, iOrb, sol, l_schreier, fce, os);
#ifdef SANITY_CHECK_HECKE_OPERATORS
      if (!is_chain_coset_stabilizer_invariant(symbol_level, iOrb, symbol_level, sol_avg, l_schreier, fce, os)) {
        std::cerr << "HECKE: sharbly, the averaged chain should be invariant\n";
        throw TerminalException{1};
      }
#endif
      size_t n_orbit = co.l_orbit.size();
      for (size_t pos = 0; pos < n_orbit; pos++) {
        int j2 = co.l_orbit[pos];
        int sign = get_face_orientation(face.EXT, LinSpa.ListMat, face.or_info, co.transversal[pos]);
        MyMatrix<Tint> Minv = Inverse(co.transversal_conj[pos]);
        images[symbol_level][iOrb][j2] = chain_right_multiply(sol_avg, T(sign), Minv);
        status[j2] = 1;
      }
    }
  }
#ifdef TIMINGS_HECKE_OPERATORS
  os << "|HECKE: sharbly images|=" << time << "\n";
#endif
  // The level below needs the coverings by the perfect cones.
  if (symbol_level > 0 && pco.compute_partial_contracting_homotopy) {
    compute_sharbly_level_below(symbol_level - 1, fce, cosets, walker, reducer, images, os);
    l_level_computed[symbol_level - 1] = 1;
#ifdef TIMINGS_HECKE_OPERATORS
    os << "|HECKE: sharbly images level below|=" << time << "\n";
#endif
  }
  HeckeChainMap<T, Tint> hcm;
  hcm.x = x;
  hcm.cosets = std::move(cosets);
  hcm.images = std::move(images);
  hcm.l_level_computed = std::move(l_level_computed);
  hcm.method = "Sharbly";
  hcm.transposed = false;
  return hcm;
}

/*
  The extension of a chain map of the sharbly method to the levels s - 2,
  ..., s - depth below the symbol level s (the k-sharblies, k >= 2, which
  carry H^{nu - k}). The image of a cell F of the level k under a coset y
  is a filling X of the image R of its boundary by the chain map at the
  level k + 1: R is a cycle (d R = h(d d F) = 0), and the well rounded
  quotient complex is exact below the symbol level, so the filling exists.
  The levels already computed (for instance read from the cache) are kept.
  Returns whether a level was added.
 */
template <typename T, typename Tint, typename Tgroup>
bool extend_hecke_chain_map_sharbly(FullComplexEnumeration<T, Tint, Tgroup> const &fce,
                                    HeckeChainMap<T, Tint> &hcm, int const &depth, std::ostream &os) {
  if (hcm.method != "Sharbly") {
    std::cerr << "HECKE: only a chain map of the Sharbly method can be extended\n";
    throw TerminalException{1};
  }
  if (!fce.pctdi.pco.compute_partial_contracting_homotopy) {
    std::cerr << "HECKE: the levels below the symbol level need ComputePartialContractingHomotopy = T\n";
    throw TerminalException{1};
  }
  SymbolGeometry<T, Tint> geom(fce.pctdi.LinSpa);
  SharblyReducer<T, Tint, Tgroup> reducer = get_sharbly_reducer(fce, geom, os);
  int symbol_level = reducer.get_symbol_level();
  ConeWalker<T, Tint, Tgroup> walker(fce, os);
  bool added = false;
  for (int index = symbol_level - 2; index >= std::max(0, symbol_level - depth); index--) {
    if (hcm.is_level_computed(index)) {
      continue;
    }
    if (!hcm.is_level_computed(index + 1)) {
      throw HeckeFailure{"sharbly extension: the level " + std::to_string(index + 1) + " is not computed"};
    }
#ifdef TIMINGS_HECKE_OPERATORS
    MicrosecondTime time;
#endif
    compute_sharbly_level_below(index, fce, hcm.cosets, walker, reducer, hcm.images, os);
    hcm.l_level_computed[index] = 1;
    added = true;
#ifdef TIMINGS_HECKE_OPERATORS
    os << "|HECKE: sharbly images level " << index << "|=" << time << "\n";
#endif
  }
  return added;
}

// clang-format off
#endif  // SRC_PERFECT_HECKE_SHARBLY_H_
// clang-format on
