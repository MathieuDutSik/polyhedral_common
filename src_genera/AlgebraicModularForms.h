// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_ALGEBRAICMODULARFORMS_H_
#define SRC_GENERA_ALGEBRAICMODULARFORMS_H_

// clang-format off
#include "Group.h"
#include "Permutation.h"
#include "OrderLatticeNeighbor.h"
#include "LatticeStabEquiCan.h"
#include "LatticeAutomorphismVinberg.h"
#include "InvariantVectorFamily.h"
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_ALGEBRAIC_MODULAR_FORMS
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_ALGEBRAIC_MODULAR_FORMS
#endif

#ifdef TIMINGS
#define TIMINGS_ALGEBRAIC_MODULAR_FORMS
#endif

/*
  Algebraic modular forms (Gross) for the unitary group of a definite
  Hermitian lattice over O_K, or the orthogonal group of a positive definite
  even lattice over Z, of trivial weight, with a level structure at
  a prime N of O_K: the functions on the pairs (L_j, l) with L_j a class of
  the genus and l an orbit of Aut(L_j) on the lines of L_j / N L_j (all the
  lines for N split, the isotropic ones for N inert), and the Hecke operator
  at a prime P coprime to N: the sum over the neighbours of L_j at P.

  The classes are kept in canonical form (ComputeCanonicalFormMultiple on
  the family (G, G Omega^T)): the canonical basis of a neighbour is that of
  its class, so the class is found by a dictionary lookup, and the map
  sending the canonical basis of the neighbour to that of the class is the
  isometry that carries its level structure, without any isometry test.

  The genus is the closure of the seed under the neighbours at all the
  Hecke primes asked for. That is the genus when the neighbour graph is
  connected on it, which is not checked here; the mass sum_j 1/|Aut(L_j)|
  is returned for that.
 */

using AmfTidx = uint32_t;
using AmfTelt = permutalib::SingleSidedPerm<AmfTidx>;
using AmfTgroup = permutalib::Group<AmfTelt, mpz_class>;

struct AmfClass {
  // the canonical Z-basis, rows in the realified coordinates
  MyMatrix<AmfT> basis;
  MyMatrix<AmfT> basis_inv;
  MyMatrix<AmfTint> G;
  MyMatrix<AmfTint> Omega;
  std::vector<MyMatrix<AmfTint>> l_aut;
  mpz_class aut_order;
  // the level structure: the line space mod N, the orbit of each line
  // (-1 for a non isotropic one) and a representative line of each orbit
  LineSpace ls_level;
  std::vector<int64_t> orbit_of_line;
  std::vector<size_t> l_orbit_repr;
};

struct AmfHeckeMatrix {
  OrderPrime P;
  // the entries (row, column) -> multiplicity, rows the source states
  std::map<std::pair<size_t, size_t>, int64_t> entries;
  // the source states of each class (row sum = number of neighbours)
  size_t n_neighbor = 0;
};

struct AmfResult {
  std::vector<AmfClass> l_class;
  std::vector<size_t> l_offset;
  size_t n_state = 0;
  mpq_class mass;
  std::vector<AmfHeckeMatrix> l_hecke;
};

struct AmfCanonical {
  MyMatrix<AmfT> basis;
  MyMatrix<AmfTint> G;
  MyMatrix<AmfTint> Omega;
};

/*
  The canonical form of the lattice of basis A (realified coordinates).
 */
inline AmfCanonical amf_canonical(RealifiedRing<AmfT, AmfTint> const &ring, MyMatrix<AmfT> const &A,
                                  std::ostream &os) {
  // Over Z the family is {G} and Omega is empty.
  MyMatrix<AmfT> B_T = UniversalMatrixConversion<AmfT, AmfTint>(ring.B);
  MyMatrix<AmfT> G = A * B_T * A.transpose();
  std::vector<MyMatrix<AmfT>> l_mat{G};
  MyMatrix<AmfT> W_T;
  if (ring.d == 2) {
    W_T = UniversalMatrixConversion<AmfT, AmfTint>(ring.W);
    MyMatrix<AmfT> Om = A * W_T * Inverse(A);
    l_mat.push_back(G * Om.transpose());
  }
  MyMatrix<AmfTint> Pc = ComputeCanonicalFormMultiple<AmfT, AmfTint, AmfTgroup>(l_mat, os);
  MyMatrix<AmfT> Pc_T = UniversalMatrixConversion<AmfT, AmfTint>(Pc);
  MyMatrix<AmfT> A_can = Pc_T * A;
  MyMatrix<AmfT> G_can = A_can * B_T * A_can.transpose();
  MyMatrix<AmfT> Om_can = ring.d == 2 ? MyMatrix<AmfT>(A_can * W_T * Inverse(A_can)) : MyMatrix<AmfT>(0, 0);
#ifdef SANITY_CHECK_ALGEBRAIC_MODULAR_FORMS
  if (!IsIntegralMatrix(G_can) || !IsIntegralMatrix(Om_can)) {
    std::cerr << "AMF: the lattice is not integral or not an O_K-module\n";
    throw TerminalException{1};
  }
#endif
  return {A_can, UniversalMatrixConversion<AmfTint, AmfT>(G_can),
          UniversalMatrixConversion<AmfTint, AmfT>(Om_can)};
}

inline bool amf_is_isotropic(MyMatrix<AmfTint> const &G, std::vector<int64_t> const &z,
                             int64_t const &ell) {
  MyVector<AmfTint> x(z.size());
  for (size_t k = 0; k < z.size(); k++) {
    x(k) = amf_int(z[k]);
  }
  return order_residue(order_hermitian_norm(G, x), ell) == 0;
}

// The residue vector of the row vector z under the matrix g, mod ell.
inline std::vector<int64_t> amf_apply(std::vector<int64_t> const &z, MyMatrix<AmfTint> const &g,
                                      int64_t const &ell) {
  int n = z.size();
  std::vector<int64_t> w(n, 0);
  for (int j = 0; j < n; j++) {
    int64_t sum = 0;
    for (int k = 0; k < n; k++) {
      sum = (sum + z[k] * order_residue(g(k, j), ell)) % ell;
    }
    w[j] = sum;
  }
  return w;
}

/*
  A new class: its automorphism group and the orbits of the level lines.
 */
inline AmfClass amf_new_class(RealifiedRing<AmfT, AmfTint> const &ring, AmfCanonical const &can,
                              std::optional<OrderPrime> const &N, std::ostream &os) {
  AmfClass cl;
  cl.basis = can.basis;
  cl.basis_inv = Inverse(can.basis);
  cl.G = can.G;
  cl.Omega = can.Omega;
  MyMatrix<AmfT> G_T = UniversalMatrixConversion<AmfT, AmfTint>(cl.G);
  MyMatrix<AmfT> Om_T = UniversalMatrixConversion<AmfT, AmfTint>(cl.Omega);
  std::vector<MyMatrix<AmfT>> l_mat{G_T};
  if (ring.d == 2) {
    l_mat.push_back(G_T * Om_T.transpose());
  }
  cl.l_aut = ArithmeticAutomorphismGroupMultiple<AmfT, AmfTint, AmfTgroup>(l_mat, os);
#ifdef SANITY_CHECK_ALGEBRAIC_MODULAR_FORMS
  for (auto &g : cl.l_aut) {
    if (g * cl.G * g.transpose() != cl.G || (ring.d == 2 && g * cl.Omega != cl.Omega * g)) {
      std::cerr << "AMF: an automorphism does not preserve the Hermitian structure\n";
      throw TerminalException{1};
    }
  }
#endif
  MyMatrix<AmfTint> SHV = GetCanonicVectorFamily<AmfT, AmfTint>(G_T, os).get_full();
  cl.aut_order = OrderFromGens<AmfTint, AmfTgroup>(SHV, cl.l_aut);
  if (!N) {
    cl.orbit_of_line = {0};
    cl.l_orbit_repr = {0};
    return cl;
  }
  cl.ls_level = get_line_space(ring, cl.Omega, *N);
  size_t n_line = cl.ls_level.n_line();
  cl.orbit_of_line.assign(n_line, -1);
  std::vector<bool> is_valid(n_line, true);
  if (N->f == 2) {
    for (size_t idx = 0; idx < n_line; idx++) {
      is_valid[idx] = amf_is_isotropic(cl.G, cl.ls_level.line_vector(idx), N->ell);
    }
  }
  int64_t n_orbit = 0;
  for (size_t idx = 0; idx < n_line; idx++) {
    if (!is_valid[idx] || cl.orbit_of_line[idx] != -1) {
      continue;
    }
    cl.l_orbit_repr.push_back(idx);
    cl.orbit_of_line[idx] = n_orbit;
    std::vector<size_t> l_active{idx};
    while (!l_active.empty()) {
      size_t pos = l_active.back();
      l_active.pop_back();
      std::vector<int64_t> z = cl.ls_level.line_vector(pos);
      for (auto &g : cl.l_aut) {
        std::optional<size_t> opt = cl.ls_level.line_index(amf_apply(z, g, N->ell));
        size_t img = *opt;
        if (cl.orbit_of_line[img] == -1) {
          cl.orbit_of_line[img] = n_orbit;
          l_active.push_back(img);
        }
      }
    }
    n_orbit++;
  }
  return cl;
}

inline AmfResult compute_algebraic_modular_forms(RealifiedRing<AmfT, AmfTint> const &ring,
                                                 std::optional<OrderPrime> const &N,
                                                 std::vector<OrderPrime> const &l_prime,
                                                 std::ostream &os) {
  AmfResult res;
  std::map<std::pair<MyMatrix<AmfTint>, MyMatrix<AmfTint>>, size_t> map_class;
  auto f_insert = [&](AmfCanonical const &can) -> size_t {
    std::pair<MyMatrix<AmfTint>, MyMatrix<AmfTint>> key{can.G, can.Omega};
    auto iter = map_class.find(key);
    if (iter != map_class.end()) {
      return iter->second;
    }
    size_t idx = res.l_class.size();
    res.l_class.push_back(amf_new_class(ring, can, N, os));
    map_class[key] = idx;
#ifdef DEBUG_ALGEBRAIC_MODULAR_FORMS
    os << "AMF: class " << idx << " |Aut|=" << res.l_class[idx].aut_order << " level orbits "
       << res.l_class[idx].l_orbit_repr.size() << "\n";
#endif
    return idx;
  };
  int n = ring.n();
  (void)f_insert(amf_canonical(ring, IdentityMat<AmfT>(n), os));
  // The transitions (class, orbit) -> (class, orbit), per prime.
  struct Transition {
    size_t i;
    size_t o;
    size_t j;
    size_t o2;
  };
  std::vector<std::vector<Transition>> l_trans(l_prime.size());
  std::vector<size_t> l_n_neighbor(l_prime.size(), 0);
  for (size_t i = 0; i < res.l_class.size(); i++) {
    for (size_t ip = 0; ip < l_prime.size(); ip++) {
      OrderPrime const &P = l_prime[ip];
      // copies: f_insert may reallocate l_class
      MyMatrix<AmfTint> G_i = res.l_class[i].G;
      MyMatrix<AmfTint> Om_i = res.l_class[i].Omega;
      MyMatrix<AmfT> basis_i = res.l_class[i].basis;
      std::vector<size_t> repr_i = res.l_class[i].l_orbit_repr;
      LineSpace ls_level_i = res.l_class[i].ls_level;
      LineSpace ls = get_line_space(ring, Om_i, P);
      size_t n_line = ls.n_line();
      size_t n_nb = 0;
      for (size_t idx = 0; idx < n_line; idx++) {
        std::optional<MyMatrix<AmfT>> opt = order_neighbor(ring, G_i, Om_i, P, ls.line_vector(idx));
        if (!opt) {
          continue;
        }
        n_nb++;
        MyMatrix<AmfT> A_nb = (*opt) * basis_i;
        AmfCanonical can = amf_canonical(ring, A_nb, os);
        size_t j = f_insert(can);
        MyMatrix<AmfT> can_inv_tr = Inverse(can.basis).transpose();
        for (size_t o = 0; o < repr_i.size(); o++) {
          size_t o2 = 0;
          if (N) {
            // the level line of L_i, carried to L_j: its vector in the
            // coordinates of the canonical basis of the neighbour, read in
            // those of L_j
            std::vector<int64_t> z = ls_level_i.line_vector(repr_i[o]);
            MyVector<AmfT> v(n);
            for (int k = 0; k < n; k++) {
              v(k) = amf_rat(z[k]);
            }
            MyVector<AmfT> v_amb = basis_i.transpose() * v;
            MyVector<AmfT> v_nb = can_inv_tr * v_amb;
            std::vector<int64_t> w = order_residue_vector(v_nb, N->ell);
            std::optional<size_t> opt_idx = res.l_class[j].ls_level.line_index(w);
            if (!opt_idx) {
              std::cerr << "AMF: a level line is sent to zero\n";
              throw TerminalException{1};
            }
            int64_t orb = res.l_class[j].orbit_of_line[*opt_idx];
            if (orb == -1) {
              std::cerr << "AMF: a level line is sent to a non isotropic line\n";
              throw TerminalException{1};
            }
            o2 = static_cast<size_t>(orb);
          }
          l_trans[ip].push_back({i, o, j, o2});
        }
      }
      if (i == 0) {
        l_n_neighbor[ip] = n_nb;
      } else if (l_n_neighbor[ip] != n_nb) {
        std::cerr << "AMF: the classes do not have the same number of neighbours\n";
        throw TerminalException{1};
      }
#ifdef DEBUG_ALGEBRAIC_MODULAR_FORMS
      os << "AMF: class " << i << " prime " << P.name << " done, " << res.l_class.size()
         << " classes known\n";
#endif
    }
  }
  size_t offset = 0;
  res.mass = 0;
  for (auto &cl : res.l_class) {
    res.l_offset.push_back(offset);
    offset += cl.l_orbit_repr.size();
    res.mass += mpq_class(1) / mpq_class(cl.aut_order);
  }
  res.n_state = offset;
  for (size_t ip = 0; ip < l_prime.size(); ip++) {
    AmfHeckeMatrix hm;
    hm.P = l_prime[ip];
    hm.n_neighbor = l_n_neighbor[ip];
    for (auto &t : l_trans[ip]) {
      hm.entries[{res.l_offset[t.i] + t.o, res.l_offset[t.j] + t.o2}] += 1;
    }
    res.l_hecke.push_back(hm);
  }
  return res;
}

// clang-format off
#endif  // SRC_GENERA_ALGEBRAICMODULARFORMS_H_
// clang-format on
