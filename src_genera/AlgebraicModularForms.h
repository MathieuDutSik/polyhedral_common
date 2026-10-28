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
#include "ClassicLLL.h"
#include "LatticeRootDecomposition.h"
#include "LocalDensity.h"
#include "GenusMass.h"
#include "HermitianMass.h"
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
  std::vector<size_t> l_orbit_size;
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
  // the mass of the genus by the mass formula, when it is covered
  std::optional<mpq_class> mass_formula;
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
inline AmfCanonical amf_canonical(RealifiedRing<AmfT, AmfTint> const &ring, MyMatrix<AmfT> const &A_in,
                                  std::ostream &os) {
  // Over Z the family is {G} and Omega is empty.
  MyMatrix<AmfT> B_T = UniversalMatrixConversion<AmfT, AmfTint>(ring.B);
  /*
    The basis of a neighbour (from GetZbasis) and the canonical bases are
    skewed, and the canonical form starts by enumerating the shells of short
    vectors, which a skewed Gram matrix makes very slow (99% of the time of a
    quinary genus of 389 classes was in that enumeration). The canonical form
    does not depend on the basis, so it is reduced first.
   */
  LLLreduction<AmfT, AmfTint> rec = LLLreducedBasis<AmfT, AmfTint>(MyMatrix<AmfT>(A_in * B_T * A_in.transpose()), os);
  MyMatrix<AmfT> A = UniversalMatrixConversion<AmfT, AmfTint>(rec.Pmat) * A_in;
  MyMatrix<AmfT> G = rec.GramMatRed;
  std::vector<MyMatrix<AmfT>> l_mat{G};
  MyMatrix<AmfT> W_T;
  if (ring.d == 2) {
    W_T = UniversalMatrixConversion<AmfT, AmfTint>(ring.W);
    MyMatrix<AmfT> Om = A * W_T * Inverse(A);
    l_mat.push_back(G * Om.transpose());
  }
  // Over Z the Weyl-chamber canonical form of GenusEnumeration; over O_K the
  // canonical form of the family.
  MyMatrix<AmfTint> Pc = ring.d == 1 ? ComputeCanonicalFormWeylChamber<AmfT, AmfTint, AmfTgroup>(G, os)
                                     : ComputeCanonicalFormMultiple<AmfT, AmfTint, AmfTgroup>(l_mat, os);
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

/*
  The orbits of the group generated by l_gen on the lines of ls that pass
  the filter (the filter has to be invariant), by their indices, with
  their sizes. If orbit_of_line is given, it receives the orbit of each
  line (-1 for those failing the filter); otherwise only one bit per line
  is kept.
 */
struct LineOrbits {
  std::vector<size_t> l_repr;
  std::vector<size_t> l_size;
};

template <typename F>
LineOrbits amf_line_orbits(LineSpace const &ls, std::vector<MyMatrix<AmfTint>> const &l_gen, F f_filter,
                           std::vector<int64_t> *orbit_of_line) {
  int64_t ell = ls.F.ell;
  int n = ls.n;
  // the generators mod ell, those acting as scalars on the lines dropped
  std::vector<std::vector<int64_t>> l_gen_red;
  for (auto &g : l_gen) {
    std::vector<int64_t> g_red(n * n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        g_red[i * n + j] = order_residue(g(i, j), ell);
      }
    }
    l_gen_red.push_back(g_red);
  }
  auto f_apply = [&](std::vector<int64_t> const &z, std::vector<int64_t> const &g_red) -> size_t {
    std::vector<int64_t> w(n, 0);
    for (int j = 0; j < n; j++) {
      int64_t sum = 0;
      for (int k = 0; k < n; k++) {
        sum = (sum + z[k] * g_red[k * n + j]) % ell;
      }
      w[j] = sum;
    }
    return *ls.line_index(w);
  };
  size_t n_line = ls.n_line();
  Face status(n_line);
  if (orbit_of_line) {
    orbit_of_line->assign(n_line, -1);
  }
  LineOrbits orb;
  for (size_t idx = 0; idx < n_line; idx++) {
    if (status[idx] || !f_filter(idx)) {
      continue;
    }
    int64_t i_orb = orb.l_repr.size();
    orb.l_repr.push_back(idx);
    size_t size = 1;
    status[idx] = 1;
    if (orbit_of_line) {
      (*orbit_of_line)[idx] = i_orb;
    }
    std::vector<size_t> l_active{idx};
    while (!l_active.empty()) {
      size_t pos = l_active.back();
      l_active.pop_back();
      std::vector<int64_t> z = ls.line_vector(pos);
      for (auto &g_red : l_gen_red) {
        size_t img = f_apply(z, g_red);
        if (!status[img]) {
          status[img] = 1;
          size++;
          if (orbit_of_line) {
            (*orbit_of_line)[img] = i_orb;
          }
          l_active.push_back(img);
        }
      }
    }
    orb.l_size.push_back(size);
  }
  return orb;
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
  // The automorphisms are computed on a reduced basis Q (rows in the
  // canonical basis) and brought back: g = Q^{-1} g_red Q.
  MyMatrix<AmfT> G_can = UniversalMatrixConversion<AmfT, AmfTint>(cl.G);
  LLLreduction<AmfT, AmfTint> rec = LLLreducedBasis<AmfT, AmfTint>(G_can, os);
  MyMatrix<AmfT> Q = UniversalMatrixConversion<AmfT, AmfTint>(rec.Pmat);
  MyMatrix<AmfT> Q_inv = Inverse(Q);
  MyMatrix<AmfT> G_T = rec.GramMatRed;
  std::vector<MyMatrix<AmfT>> l_mat{G_T};
  if (ring.d == 2) {
    MyMatrix<AmfT> Om_T = Q * UniversalMatrixConversion<AmfT, AmfTint>(cl.Omega) * Q_inv;
    l_mat.push_back(G_T * Om_T.transpose());
  }
  std::vector<MyMatrix<AmfTint>> l_aut_red;
  if (ring.d == 1) {
    // Over Z the Vinberg decomposition of GenusEnumeration gives the order
    // and the generators without enumerating the group.
    LatticeAutInfo<AmfT, AmfTint, AmfTgroup> info =
        GetLatticeAutInfo<AmfT, AmfTint, AmfTgroup>(G_T, GenusScheme::CanonicalForm, os);
    l_aut_red = info.ListGenMat;
    cl.aut_order = info.order;
  } else {
    l_aut_red = ArithmeticAutomorphismGroupMultiple<AmfT, AmfTint, AmfTgroup>(l_mat, os);
    MyMatrix<AmfTint> SHV = GetCanonicVectorFamily<AmfT, AmfTint>(G_T, os).get_full();
    cl.aut_order = OrderFromGens<AmfTint, AmfTgroup>(SHV, l_aut_red);
  }
  for (auto &g_red : l_aut_red) {
    MyMatrix<AmfT> g = Q_inv * UniversalMatrixConversion<AmfT, AmfTint>(g_red) * Q;
    cl.l_aut.push_back(UniversalMatrixConversion<AmfTint, AmfT>(g));
  }
#ifdef SANITY_CHECK_ALGEBRAIC_MODULAR_FORMS
  for (auto &g : cl.l_aut) {
    if (g * cl.G * g.transpose() != cl.G || (ring.d == 2 && g * cl.Omega != cl.Omega * g)) {
      std::cerr << "AMF: an automorphism does not preserve the Hermitian structure\n";
      throw TerminalException{1};
    }
  }
#endif
  if (!N) {
    cl.orbit_of_line = {0};
    cl.l_orbit_repr = {0};
    cl.l_orbit_size = {1};
    return cl;
  }
  cl.ls_level = get_line_space(ring, cl.Omega, *N);
  auto f_filter = [&](size_t const &idx) -> bool {
    return N->f == 1 || amf_is_isotropic(cl.G, cl.ls_level.line_vector(idx), N->ell);
  };
  LineOrbits orb = amf_line_orbits(cl.ls_level, cl.l_aut, f_filter, &cl.orbit_of_line);
  cl.l_orbit_repr = orb.l_repr;
  cl.l_orbit_size = orb.l_size;
  return cl;
}

/*
  The mass of the genus of the seed lattice by the mass formula: GenusMass
  over Z (even lattices), HermitianGenusMass over O_K (the lattices that are
  unimodular at every prime). When it is known, the classes found account
  for the whole genus iff the sum of the 1 / |Aut(L_j)| is equal to it.
 */
inline std::optional<mpq_class> amf_genus_mass(RealifiedRing<AmfT, AmfTint> const &ring, std::ostream &os) {
  if (ring.d == 1) {
    return GenusMass<AmfT, AmfTint>(UniversalMatrixConversion<AmfT, AmfTint>(ring.B));
  }
  return HermitianGenusMass(ring, os);
}

/*
  The Hecke operator at P on the class L_i: one neighbour per orbit X of
  Aut(L_i) on the lines of L_i / P L_i. For a level orbit O of L_i, the
  entry (O, O') counts the lines x with the image of a fixed l in O by the
  neighbour of x in O'; counting the pairs (x, l) of X x O with that
  property in two ways, it is |X| / |O| times the number of l in O sent to
  O' by the neighbour of the representative of X. So the level lines are
  all carried once per orbit representative, which is a product mod N,
  while the neighbours, which cost a canonical form each, are one per
  orbit.
 */
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
  // The weighted transitions (class, orbit) -> (class, orbit), per prime.
  struct Transition {
    size_t i;
    size_t o;
    size_t j;
    size_t o2;
    mpq_class weight;
  };
  std::vector<std::vector<Transition>> l_trans(l_prime.size());
  std::vector<size_t> l_n_neighbor(l_prime.size(), 0);
  for (size_t i = 0; i < res.l_class.size(); i++) {
    for (size_t ip = 0; ip < l_prime.size(); ip++) {
      OrderPrime const &P = l_prime[ip];
#ifdef TIMINGS_ALGEBRAIC_MODULAR_FORMS
      MicrosecondTime time;
#endif
      // copies: f_insert may reallocate l_class
      MyMatrix<AmfTint> G_i = res.l_class[i].G;
      MyMatrix<AmfTint> Om_i = res.l_class[i].Omega;
      MyMatrix<AmfT> basis_i = res.l_class[i].basis;
      std::vector<int64_t> orbit_of_line_i = res.l_class[i].orbit_of_line;
      std::vector<size_t> size_i = res.l_class[i].l_orbit_size;
      LineSpace ls_level_i = res.l_class[i].ls_level;
      std::vector<MyMatrix<AmfTint>> l_aut_i = res.l_class[i].l_aut;
      LineSpace ls = get_line_space(ring, Om_i, P);
      LineOrbits orb = amf_line_orbits(ls, l_aut_i, [](size_t const &) -> bool { return true; }, nullptr);
#ifdef TIMINGS_ALGEBRAIC_MODULAR_FORMS
      os << "|AMF: line orbits|=" << time << "\n";
#endif
      size_t n_nb = 0;
      for (size_t i_orb = 0; i_orb < orb.l_repr.size(); i_orb++) {
        size_t size_X = orb.l_size[i_orb];
        std::optional<MyMatrix<AmfT>> opt = order_neighbor(ring, G_i, Om_i, P, ls.line_vector(orb.l_repr[i_orb]));
        if (!opt) {
          continue;
        }
        n_nb += size_X;
        MyMatrix<AmfT> A_nb = (*opt) * basis_i;
        AmfCanonical can = amf_canonical(ring, A_nb, os);
        size_t j = f_insert(can);
        if (!N) {
          l_trans[ip].push_back({i, 0, j, 0, mpq_class(static_cast<long>(size_X))});
          continue;
        }
        // The level lines of L_i carried to L_j: the coordinates in the
        // basis of L_i, then in the canonical basis of the neighbour, which
        // are those of L_j. The matrix has denominators dividing ell(P),
        // which is prime to N.
        MyMatrix<AmfT> R = basis_i * Inverse(can.basis);
        std::vector<int64_t> R_red(n * n);
        for (int a = 0; a < n; a++) {
          for (int b = 0; b < n; b++) {
            R_red[a * n + b] = order_residue(R(a, b), N->ell);
          }
        }
        AmfClass const &cl_j = res.l_class[j];
        std::map<std::pair<size_t, size_t>, size_t> count;
        for (size_t idx = 0; idx < orbit_of_line_i.size(); idx++) {
          if (orbit_of_line_i[idx] == -1) {
            continue;
          }
          std::vector<int64_t> z = ls_level_i.line_vector(idx);
          std::vector<int64_t> w(n, 0);
          for (int b = 0; b < n; b++) {
            int64_t sum = 0;
            for (int a = 0; a < n; a++) {
              sum = (sum + z[a] * R_red[a * n + b]) % N->ell;
            }
            w[b] = sum;
          }
          std::optional<size_t> opt_idx = cl_j.ls_level.line_index(w);
          if (!opt_idx) {
            std::cerr << "AMF: a level line is sent to zero\n";
            throw TerminalException{1};
          }
          int64_t orb2 = cl_j.orbit_of_line[*opt_idx];
          if (orb2 == -1) {
            std::cerr << "AMF: a level line is sent to a non isotropic line\n";
            throw TerminalException{1};
          }
          count[{static_cast<size_t>(orbit_of_line_i[idx]), static_cast<size_t>(orb2)}]++;
        }
        for (auto &kv : count) {
          size_t o = kv.first.first;
          mpq_class weight = mpq_class(static_cast<long>(size_X * kv.second)) /
                             mpq_class(static_cast<long>(size_i[o]));
          l_trans[ip].push_back({i, o, j, kv.first.second, weight});
        }
      }
      if (i == 0) {
        l_n_neighbor[ip] = n_nb;
      } else if (l_n_neighbor[ip] != n_nb) {
        std::cerr << "AMF: the classes do not have the same number of neighbours\n";
        throw TerminalException{1};
      }
#ifdef TIMINGS_ALGEBRAIC_MODULAR_FORMS
      os << "|AMF: class " << i << " prime " << P.name << ", " << orb.l_repr.size() << " orbits|=" << time << "\n";
#endif
#ifdef DEBUG_ALGEBRAIC_MODULAR_FORMS
      os << "AMF: class " << i << " prime " << P.name << " done, " << orb.l_repr.size() << " orbits of lines, "
         << res.l_class.size() << " classes known\n";
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
  res.mass_formula = amf_genus_mass(ring, os);
  for (size_t ip = 0; ip < l_prime.size(); ip++) {
    AmfHeckeMatrix hm;
    hm.P = l_prime[ip];
    hm.n_neighbor = l_n_neighbor[ip];
    std::map<std::pair<size_t, size_t>, mpq_class> entries;
    for (auto &t : l_trans[ip]) {
      entries[{res.l_offset[t.i] + t.o, res.l_offset[t.j] + t.o2}] += t.weight;
    }
    for (auto &kv : entries) {
      if (kv.second.get_den() != 1) {
        std::cerr << "AMF: a Hecke matrix entry is not an integer\n";
        throw TerminalException{1};
      }
      hm.entries[kv.first] = kv.second.get_num().get_si();
    }
    res.l_hecke.push_back(hm);
  }
  return res;
}

// clang-format off
#endif  // SRC_GENERA_ALGEBRAICMODULARFORMS_H_
// clang-format on
