// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_MATRIXCANONICALFORM_H_
#define SRC_LATT_MATRIXCANONICALFORM_H_

// clang-format off
#include "InvariantVectorFamily.h"
#include "MAT_MatrixInt.h"
#include "MatrixGroup.h"
#include "Shvec_exact.h"
#include "PolytopeEquiStabInt.h"
#include <utility>
#include <string>
#include <vector>
// clang-format on

/*
  This is code for computing equivalence, canonical form and automorphism
  group of positive definite quadratic forms.
  This is not the T-space relevant code.
 */

/*
  This is the set of algorithms for computing:
  * Arithmetic equivalence under GL_n(Z) of positive definite quadratic forms.
  * Stabilizer under GL_n(Z) action of positive definite quadratic forms.
  * Canonical form under GL_n(Z) action of positive definite form.

  The classic algorithm for iterm 1 and 2 is the Plesken-Souvignier algorithm:
  W. Plesken, B. Souvignier, Computing isometries of lattices, Journal of
     Symbolic Computation (1997) 24, 327-344.
  It works by using the standard basis (e_1, ...., e_n) with e_i the vector
  having
  (e_i)_j = | 1 is i=j
            | 0 if i<>j
  That standard basis is then looked for possible images in the list of short
  vectors. This is very efficient. The algorithm of this section are different.
  Instead we compute the set of pairwise scalar products and use partition
  backtrack for testing isomorphism and stabilizer:
  * This approach is not scalable for large lattice such as Leech lattice.
    Here PS algorithm works better.
  * But it may be more efficient for small dimensional lattices.

  For the canonical form computation, we have
  Mathieu Dutour Sikirić, Anna Haensch, John Voight, Wessel Van Woerden,
  A canonical form for positive definite matrices, Proceedings of the
  Fourteenth Algorithmic Number Theory Symposium (ANTS-XIV), edited by
  Steven Galbraith, Open Book Series 4, Mathematical Sciences Publishers,
  Berkeley, 2020, prepring at https://arxiv.org/abs/2004.14022

  The approach of Plesken-Souvignier does not appear to be feasible for
  computing canonical form.

  So, the approach of this file is to compute configuration of short
  vectors and then compute automorphism, isomorphism and canonical forms.

  Computing those configurations of short vectors is a difficult
  business. Two tricks are available to accelerate it:
  * Computing also for the dual. The inverse matrix A^(-1) might be
    easier to compute with.
  * We do not necessarily need a spanning configuration. Having a full
    rank one (but not neceessarily spanning) that can then be used for
    rational stabilizer and equivalence is good enough. Then we can
    use the finite index algorithm for concluding.

 */


#ifdef DEBUG
#define DEBUG_LATTICE_STAB_EQUI_CAN
#endif

#ifdef DISABLE_DEBUG_LATTICE_STAB_EQUI_CAN
#undef DEBUG_LATTICE_STAB_EQUI_CAN
#endif

#ifdef TIMINGS
#define TIMINGS_LATTICE_STAB_EQUI_CAN
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_LATTICE_STAB_EQUI_CAN
#endif

template <typename T, typename Tint> struct Canonic_PosDef {
  MyMatrix<Tint> Basis;
  MyMatrix<Tint> SHV;
  MyMatrix<T> Mat;
};

//
// The canonic form
//

/*
  Whether every matrix of the configuration has integral entries, and the
  conversion that follows. The weight matrices below are built out of the
  scalar products v A w only, so an integral configuration lets them be built
  over the ring instead of over the field.
 */
template <typename T>
bool IsIntegralListMat(std::vector<MyMatrix<T>> const &ListMat) {
  for (auto &eMat : ListMat) {
    if (!IsIntegralMatrix(eMat)) {
      return false;
    }
  }
  return true;
}

template <typename T, typename Tring>
std::vector<MyMatrix<Tring>>
ConvertListMatToRing(std::vector<MyMatrix<T>> const &ListMat) {
  std::vector<MyMatrix<Tring>> ListMatRet;
  ListMatRet.reserve(ListMat.size());
  for (auto &eMat : ListMat) {
    ListMatRet.push_back(UniversalMatrixConversion<Tring, T>(eMat));
  }
  return ListMatRet;
}

template <typename Tval, typename Tint>
MyMatrix<Tint>
CanonicallyReorder_SHV_kernel(std::vector<MyMatrix<Tval>> const &ListMat,
                              MyMatrix<Tint> const &SHV, std::ostream &os) {
  using T = Tval;
  using Tgr = GraphListAdj;
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
  os << "LSEC: Begining of ComputeCanonicalForm\n";
#endif
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  MicrosecondTime time;
#endif
  int nbRow = SHV.rows();
  int n = SHV.cols();
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
  os << "LSEC: nbRow=" << nbRow << " n=" << n << "\n";
  if (!CheckCentralSymmetry(SHV)) {
    std::cerr
        << "LSEC: The set of vector does not respect the central symmetry "
           "condition\n";
    throw TerminalException{1};
  }
#endif
  //
  // Computing the scalar product matrix
  //
  using Tidx_value = int16_t;
  // The matrices are Gram matrices, so the weight matrix is symmetric. That
  // halves its entries and halves the number of vertices of the graph, which
  // is a quarter of the edges, and the edges are what nauty is bounded by.
  const bool is_symm = true;
  WeightMatrix<is_symm, std::vector<T>, Tidx_value> WMat =
      T_TranslateToMatrix_ListMat_SHV<is_symm, T, Tint, Tidx_value>(ListMat,
                                                                    SHV, os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: WMat|=" << time << "\n";
#endif
  WMat.ReorderingSetWeight();
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: ReorderingSetWeight|=" << time << "\n";
#endif
  //
  // Computing the canonicalization of the scalar product matrix
  //
  std::vector<int> CanonicOrd =
    GetCanonicalizationVector_Kernel<std::vector<T>, Tgr, int>(WMat, os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: GetCanonicalizationVector_Kernel|=" << time << "\n";
#endif
  //
  // Building the canonical basis
  //
  MyMatrix<Tint> SHVcan(n, nbRow);
  for (int iRowCan = 0; iRowCan < nbRow; iRowCan++) {
    int iRowNative = CanonicOrd[iRowCan];
    MyVector<Tint> eRow_Tint = GetMatrixRow(SHV, iRowNative);
    AssignMatrixCol(SHVcan, iRowCan, eRow_Tint);
  }
  return SHVcan;
}

// Over the ring rather than over the field when the configuration allows it,
// see CanonicallyReorder_SHV_AbsTrick for why that leaves the answer alone.
template <typename T, typename Tint>
MyMatrix<Tint> CanonicallyReorder_SHV(std::vector<MyMatrix<T>> const &ListMat,
                                      MyMatrix<Tint> const &SHV,
                                      std::ostream &os) {
  using Tring = typename underlying_ring<T>::ring_type;
  if (IsIntegralListMat(ListMat)) {
    std::vector<MyMatrix<Tring>> ListMat_ring =
        ConvertListMatToRing<T, Tring>(ListMat);
    return CanonicallyReorder_SHV_kernel<Tring, Tint>(ListMat_ring, SHV, os);
  }
  return CanonicallyReorder_SHV_kernel<T, Tint>(ListMat, SHV, os);
}

template<typename Tint>
MyMatrix<Tint> get_canonicalization_matrix(MyMatrix<Tint> const& SHVcan, [[maybe_unused]] std::ostream &os) {
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  MicrosecondTime time;
#endif
  MyMatrix<Tint> BasisCan_pre = ComputeRowHermiteNormalForm_first(SHVcan);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: get_canonicalization_matrix, BasisCan_pre|=" << time << "\n";
#endif
  MyMatrix<Tint> BasisCan = TransposedMat(Inverse(BasisCan_pre));
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: get_canonicalization_matrix, BasisCan|=" << time << "\n";
#endif
  return BasisCan;
}



template <typename T, typename Tint>
MyMatrix<Tint> ComputeCanonicalForm_inner(std::vector<MyMatrix<T>> const &ListMat,
                                          MyMatrix<Tint> const &SHV,
                                          std::ostream &os) {
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  MicrosecondTime time;
#endif
  MyMatrix<Tint> SHVcan = CanonicallyReorder_SHV<T,Tint>(ListMat, SHV, os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: SHVcan|=" << time << "\n";
#endif
  MyMatrix<Tint> BasisCan = get_canonicalization_matrix(SHVcan, os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: BasisCan|=" << time << "\n";
#endif
  return BasisCan;
}

/*
  The canonical reordering obtained from the absolute trick: the graph is
  built on the antipodal pairs rather than on the vectors, which is a quarter
  of the vertices, and the sign of each vector is recovered afterwards.

  The lifting test and the sign propagation are the ones of
  src_group/PolytopeEquiStab.h, shared with the polytope side. Both can
  decline, in which case this returns nothing and the caller falls back to
  the reordering on the whole family. Whether they decline depends only on
  the configuration, so two isometric lattices always take the same branch
  and their canonical forms remain comparable, which is what the
  deduplication of the genus enumeration needs.

  SHVhalf holds one vector per antipodal pair. The returned matrix has the
  chosen representative of each pair as a column, in canonical order and with
  its canonical sign. It spans the same lattice as the whole family, so the
  Hermite normal form built from it is the same.
 */
template <typename Tval, typename Tint>
std::optional<MyMatrix<Tint>> CanonicallyReorder_SHV_AbsTrick_kernel(
    std::vector<MyMatrix<Tval>> const &ListMat, MyMatrix<Tint> const &SHVhalf,
    std::ostream &os) {
  using T = Tval;
  using Tgr = GraphListAdj;
  using Tidx = uint32_t;
  using Tidx_value = int16_t;
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  MicrosecondTime time;
#endif
  size_t nbPair = SHVhalf.rows();
  int n = SHVhalf.cols();
  WeightMatrixAbs<std::vector<T>, Tidx_value> WMatAbs =
      T_TranslateToMatrixAntipodal_AbsTrick_ListMat_SHV<T, Tint, Tidx_value>(
          ListMat, SHVhalf, os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: WMatAbs|=" << time << "\n";
#endif
  std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> ePair =
      GetGroupCanonicalizationVector_Kernel<std::vector<T>, Tgr, Tidx,
                                            Tidx_value, true>(WMatAbs.WMat, os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: GetGroupCanonicalizationVector_Kernel|=" << time << "\n";
#endif
  if (!AbsTrick_TestLiftGenerators<std::vector<T>, Tidx, Tidx_value>(
          WMatAbs, ePair.second, nbPair)) {
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
    os << "LSEC: a generator does not lift, falling back\n";
#endif
    return {};
  }
  std::vector<Tidx> const &CanonicOrd = ePair.first;
  std::optional<std::vector<int>> opt_signs =
      AbsTrick_GetSigns<std::vector<T>, Tidx, Tidx_value>(WMatAbs, CanonicOrd,
                                                          nbPair);
  if (!opt_signs) {
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
    os << "LSEC: the signs are not determined, falling back\n";
#endif
    return {};
  }
  std::vector<int> const &ListSigns = *opt_signs;
  MyMatrix<Tint> SHVcan(n, nbPair);
  for (size_t iPair = 0; iPair < nbPair; iPair++) {
    size_t jPair = CanonicOrd[iPair];
    Tint eSign = ListSigns[iPair];
    for (int i = 0; i < n; i++) {
      SHVcan(i, iPair) = eSign * SHVhalf(jPair, i);
    }
  }
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: CanonicallyReorder_SHV_AbsTrick|=" << time << "\n";
#endif
  return SHVcan;
}

/*
  The same, over the ring rather than over the field whenever the matrices
  allow it.

  The entries of the weight matrix are the scalar products v A w, which are
  integers as soon as A is one, and the weights are only ever compared with
  one another, an order that is the same in the ring as in the field. So the
  graph, its canonical ordering and the reordered family are unchanged.

  What changes is the price. Over Q every one of the O(|V|^2 n^2) products
  and sums normalizes its result by a gcd, and a profile of a rank-14
  determinant-351 genus has __gmpq_mul, __gmpq_set_z, __gmpz_gcd and
  __gmpn_gcd_11 at the top, with the weight matrix construction taking half
  of the canonicalization. T is a field here only because ComputeCanonicalForm
  needs one for the basis it returns, not because the scalar products need
  one.
 */
template <typename T, typename Tint>
std::optional<MyMatrix<Tint>> CanonicallyReorder_SHV_AbsTrick(
    std::vector<MyMatrix<T>> const &ListMat, MyMatrix<Tint> const &SHVhalf,
    std::ostream &os) {
  using Tring = typename underlying_ring<T>::ring_type;
  if (IsIntegralListMat(ListMat)) {
    std::vector<MyMatrix<Tring>> ListMat_ring =
        ConvertListMatToRing<T, Tring>(ListMat);
    return CanonicallyReorder_SHV_AbsTrick_kernel<Tring, Tint>(ListMat_ring,
                                                               SHVhalf, os);
  }
  return CanonicallyReorder_SHV_AbsTrick_kernel<T, Tint>(ListMat, SHVhalf, os);
}

/*
  Permutation generators of the automorphism group of an antipodal family,
  through the absolute trick: the graph is built on the pairs, which is a
  quarter of the vertices, and each of its generators is lifted to a signed
  permutation of the whole family.

  The generators returned act on the family ordered as
  CanonicVectorFamily::get_full gives it, the representatives first and their
  negatives after. The map v -> -v is among them.

  Nothing when a generator of the graph does not lift, or when the signs are
  not connected enough to determine the lift; the caller then computes the
  automorphisms from the whole family as before.
 */
template <typename Tval, typename Tint, typename Tgroup>
std::optional<std::vector<std::vector<typename Tgroup::Telt::Tidx>>>
GetListGenAutomorphism_AbsTrick_kernel(
    std::vector<MyMatrix<Tval>> const &ListMat, MyMatrix<Tint> const &SHVhalf,
    std::ostream &os) {
  using T = Tval;
  using Tgr = GraphListAdj;
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tidx_value = int16_t;
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  MicrosecondTime time;
#endif
  size_t nbPair = SHVhalf.rows();
  WeightMatrixAbs<std::vector<T>, Tidx_value> WMatAbs =
      T_TranslateToMatrixAntipodal_AbsTrick_ListMat_SHV<T, Tint, Tidx_value>(
          ListMat, SHVhalf, os);
  std::vector<std::vector<Tidx>> ListGen =
      GetStabilizerWeightMatrix_Kernel<std::vector<T>, Tgr, Tidx, Tidx_value,
                                       true>(WMatAbs.WMat, os);
  std::optional<std::vector<std::vector<Tidx>>> opt =
      AbsTrick_LiftGenerators<std::vector<T>, Tidx, Tidx_value>(WMatAbs,
                                                                ListGen,
                                                                nbPair);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: GetListGenAutomorphism_AbsTrick|=" << time << "\n";
#endif
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
  if (!opt) {
    os << "LSEC: the absolute trick does not lift the automorphisms, falling "
       << "back\n";
  }
#endif
  return opt;
}

// Over the ring rather than over the field when the configuration allows it,
// see CanonicallyReorder_SHV_AbsTrick. The generators are permutations, so
// they do not depend on which of the two the weights were computed in.
template <typename T, typename Tint, typename Tgroup>
std::optional<std::vector<std::vector<typename Tgroup::Telt::Tidx>>>
GetListGenAutomorphism_AbsTrick(std::vector<MyMatrix<T>> const &ListMat,
                                MyMatrix<Tint> const &SHVhalf,
                                std::ostream &os) {
  using Tring = typename underlying_ring<T>::ring_type;
  if (IsIntegralListMat(ListMat)) {
    std::vector<MyMatrix<Tring>> ListMat_ring =
        ConvertListMatToRing<T, Tring>(ListMat);
    return GetListGenAutomorphism_AbsTrick_kernel<Tring, Tint, Tgroup>(
        ListMat_ring, SHVhalf, os);
  }
  return GetListGenAutomorphism_AbsTrick_kernel<T, Tint, Tgroup>(ListMat,
                                                                 SHVhalf, os);
}

/*
  Canonical basis from a family that generates Z^n. The Hermite normal form
  of the canonically reordered family is then already a canonical basis.
 */
template <typename T, typename Tint>
MyMatrix<Tint>
ComputeCanonicalFormSpanning_family(std::vector<MyMatrix<T>> const &ListMat,
                                    CanonicVectorFamily<Tint> const &fam,
                                    std::ostream &os) {
  // The graph is built on the antipodal pairs, which is a quarter of the
  // vertices. When the signs cannot be recovered from it we fall back to the
  // reordering on the whole family.
  std::optional<MyMatrix<Tint>> opt =
      CanonicallyReorder_SHV_AbsTrick<T, Tint>(ListMat, fam.SHVhalf, os);
  if (opt) {
    return get_canonicalization_matrix(*opt, os);
  }
  return ComputeCanonicalForm_inner<T, Tint>(ListMat, fam.get_full(), os);
}

template <typename T, typename Tint, typename Tgroup>
MyMatrix<Tint> ComputeCanonicalFormMultiple(std::vector<MyMatrix<T>> const &ListMat,
                                            std::ostream &os) {
  //
  // Computing the Z-basis on which the computation relies.
  //
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
  os << "LSEC: Begining of ComputeCanonicalForm\n";
#endif
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  MicrosecondTime time;
#endif
  MyMatrix<T> const &inpMat = ListMat[0];
  CanonicVectorFamily<Tint> fam = GetCanonicVectorFamily<T, Tint>(inpMat, os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: GetCanonicVectorFamily|=" << time << "\n";
#endif
  return ComputeCanonicalForm_family<T, Tint, Tgroup>(ListMat, fam, os);
}

/*
  Canonical form computed from a full rank invariant vector family that does
  not necessarily span the full lattice Z^n. The vectors are canonically
  reordered from the weight matrix exactly as in ComputeCanonicalForm; the
  ambiguity left by the reordering is the full group of rational
  transformations preserving the family, and the position of the ambient
  lattice Z^n relative to the lattice spanned by the family is canonicalized
  under that finite group by LinPolytopeIntegral_Canonicalization_Subspaces.

  The returned matrix B belongs to GL_n(Z) and B * inpMat * B^T is the
  canonical form. The canonical forms of this function and of
  ComputeCanonicalForm are both canonical but differ in general: reductions
  computed with one method can only be compared with reductions computed by
  the same method.
 */
template <typename T, typename Tint, typename Tgroup>
MyMatrix<Tint>
ComputeCanonicalFormFullRank_family(std::vector<MyMatrix<T>> const &ListMat,
                                    MyMatrix<Tint> const &SHV,
                                    std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  MicrosecondTime time;
#endif
  [[maybe_unused]] MyMatrix<T> const &inpMat = ListMat[0];
  MyMatrix<Tint> SHVcan = CanonicallyReorder_SHV<T, Tint>(ListMat, SHV, os);
  MyMatrix<Tint> SHVord = TransposedMat(SHVcan);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: CanonicallyReorder_SHV|=" << time << "\n";
#endif
  MyMatrix<T> SHVord_T = UniversalMatrixConversion<T, Tint>(SHVord);
  int n_row = SHVord_T.rows();
  std::vector<T> Vdiag(n_row, T(0));
  std::vector<std::vector<Tidx>> ListGen =
      GetListGenAutomorphism_ListMat_Vdiag<T, T, Tgroup>(SHVord_T, ListMat,
                                                         Vdiag, os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: GetListGenAutomorphism_ListMat_Vdiag|=" << time << "\n";
#endif
  std::vector<MyMatrix<T>> ListMatrGens;
  for (auto &eList : ListGen) {
    Telt ePerm(eList);
    std::optional<MyMatrix<T>> opt =
        FindTransformationGeneral(SHVord_T, SHVord_T, ePerm);
    MyMatrix<T> eMatrGen =
        unfold_opt(opt, "the transformation of the family should exist");
    ListMatrGens.emplace_back(std::move(eMatrGen));
  }
  MyMatrix<T> B_T = LinPolytopeIntegral_Canonicalization_Subspaces<T, Tgroup>(
      ListMatrGens, SHVord_T, os);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: LinPolytopeIntegral_Canonicalization_Subspaces|=" << time
     << "\n";
#endif
#ifdef SANITY_CHECK_LATTICE_STAB_EQUI_CAN
  MyMatrix<T> eProd = B_T * inpMat * B_T.transpose();
  if (!IsSymmetricMatrix(eProd)) {
    std::cerr << "LSEC: the canonical form should be symmetric\n";
    throw TerminalException{1};
  }
#endif
  return UniversalMatrixConversion<Tint, T>(B_T);
}

/*
  Canonical basis from a family, routed on whether the family generates Z^n:
  the Hermite normal form when it does, the subspace canonicalization when
  it only has full rank.
 */
template <typename T, typename Tint, typename Tgroup>
MyMatrix<Tint>
ComputeCanonicalForm_family(std::vector<MyMatrix<T>> const &ListMat,
                            CanonicVectorFamily<Tint> const &fam,
                            std::ostream &os) {
  if (fam.spans_lattice) {
    return ComputeCanonicalFormSpanning_family<T, Tint>(ListMat, fam, os);
  }
  // The subspace route needs the automorphisms of the configuration as a
  // permutation group, so it needs both signs.
  return ComputeCanonicalFormFullRank_family<T, Tint, Tgroup>(
      ListMat, fam.get_full(), os);
}

/*
  The canonical form. Works from the smaller of the two full rank families
  and takes whichever of the two canonicalizations that family allows.

  T has to be a field: a family that does not span Z^n is canonicalized
  through LinPolytopeIntegral_Canonicalization_Subspaces, which divides.
  Instantiating with a ring such as mpz_class fails on a static assertion.

  Different families give different canonical forms, so a reduction computed
  here may only be compared with another computed here, never with one from
  ComputeCanonicalFormSymplectic.
 */
template <typename T, typename Tint, typename Tgroup>
MyMatrix<Tint> ComputeCanonicalForm(MyMatrix<T> const &inpMat,
                                    std::ostream &os) {
  CanonicVectorFamily<Tint> fam = GetCanonicVectorFamily<T, Tint>(inpMat, os);
  std::vector<MyMatrix<T>> ListMat{inpMat};
  return ComputeCanonicalForm_family<T, Tint, Tgroup>(ListMat, fam, os);
}

template <typename T, typename Tint>
MyMatrix<Tint> ComputeCanonicalFormSymplectic(MyMatrix<T> const &inpMat, std::ostream &os) {
  int n_tot = inpMat.rows();
  if (n_tot % 2 == 1) {
    std::cerr << "LSEC: The dimension is odd\n";
    throw TerminalException{1};
  }
  int n = n_tot / 2;
  MyMatrix<T> SympFormMat = ZeroMatrix<T>(2 * n, 2 * n);
  for (int i = 0; i < n; i++) {
    SympFormMat(i, n + i) = 1;
    SympFormMat(n + i, i) = -1;
  }
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(inpMat, os);
  std::vector<MyMatrix<T>> ListMat{inpMat, SympFormMat};
  MyMatrix<Tint> SHVvan = CanonicallyReorder_SHV<T,Tint>(ListMat, SHV, os);
  MyMatrix<Tint> BasisSymp = SYMPL_ComputeSymplecticBasis(SHVvan);
  return BasisSymp;
}

//
// Automorphism code
//

/*
  The integral automorphisms of the configuration, read off the graph on the
  antipodal pairs when the absolute trick concludes and off the whole family
  otherwise. SHV_T has to be the family in the order CanonicVectorFamily
  gives it, the representatives followed by their negatives, since that is
  what the lifted generators index.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<T>>
GetIntAutomorphism_Family(std::vector<MyMatrix<T>> const &ListMat,
                          CanonicVectorFamily<Tint> const &fam,
                          MyMatrix<T> const &SHV_T, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  std::optional<std::vector<std::vector<Tidx>>> opt =
      GetListGenAutomorphism_AbsTrick<T, Tint, Tgroup>(ListMat, fam.SHVhalf,
                                                       os);
  if (opt) {
    return GetIntAutomorphism_FromPermGens<T, Tgroup>(SHV_T, ListMat, *opt, os);
  }
  std::vector<T> Vdiag(SHV_T.rows(), T(0));
  return GetIntAutomorphism_ListMat_Vdiag<T, Tgroup>(SHV_T, ListMat, Vdiag, os);
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>> ArithmeticAutomorphismGroupMultiple_inner(
    std::vector<MyMatrix<T>> const &ListMat, MyMatrix<Tint> const &SHV,
    std::ostream &os) {
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T, Tint>(SHV);
  int n_row = SHV_T.rows();
  std::vector<T> Vdiag(n_row, T(0));
  std::vector<MyMatrix<T>> LGen =
      GetIntAutomorphism_ListMat_Vdiag<T, Tgroup>(SHV_T, ListMat, Vdiag, os);
  std::vector<MyMatrix<Tint>> ListGenRet;
  for (auto &M_T : LGen) {
    MyMatrix<Tint> M = UniversalMatrixConversion<Tint, T>(M_T);
    ListGenRet.push_back(M);
  }
  return ListGenRet;
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>>
ArithmeticAutomorphismGroup_inner(MyMatrix<T> const &inpMat,
                                  MyMatrix<Tint> const &SHV, std::ostream &os) {
  std::vector<MyMatrix<T>> ListMat{inpMat};
  return ArithmeticAutomorphismGroupMultiple_inner<T, Tint, Tgroup>(ListMat,
                                                                    SHV, os);
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>>
ArithmeticAutomorphismGroupMultiple(std::vector<MyMatrix<T>> const &ListMat,
                                    std::ostream &os) {
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  MicrosecondTime time;
#endif
  // The smaller of the two families rather than the shells, and the
  // automorphisms off the antipodal pairs when the trick concludes.
  CanonicVectorFamily<Tint> fam =
      GetCanonicVectorFamily<T, Tint>(ListMat[0], os);
  MyMatrix<Tint> SHV = fam.get_full();
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: GetCanonicVectorFamily|=" << time << "\n";
#endif
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T, Tint>(SHV);
  std::vector<MyMatrix<T>> LGen =
      GetIntAutomorphism_Family<T, Tint, Tgroup>(ListMat, fam, SHV_T, os);
  std::vector<MyMatrix<Tint>> ListGenRet;
  for (auto &M_T : LGen) {
    ListGenRet.push_back(UniversalMatrixConversion<Tint, T>(M_T));
  }
  return ListGenRet;
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>>
ArithmeticAutomorphismGroup(MyMatrix<T> const &inpMat, std::ostream &os) {
  std::vector<MyMatrix<T>> ListMat{inpMat};
  return ArithmeticAutomorphismGroupMultiple<T, Tint, Tgroup>(ListMat, os);
}

//
// Equivalence code
//

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>> ArithmeticEquivalenceMultiple_inner(
    std::vector<MyMatrix<T>> const &ListMat1, MyMatrix<Tint> const &SHV1,
    std::vector<MyMatrix<T>> const &ListMat2, MyMatrix<Tint> const &SHV2,
    std::ostream &os) {
#ifdef DEBUG_LATTICE_STAB_EQUI_CAN
  os << "LSEC: |SHV1|=" << SHV1.rows() << " |SHV2|=" << SHV2.rows() << "\n";
#endif
  if (SHV1.rows() != SHV2.rows())
    return {};
  MyMatrix<T> SHV1_T = UniversalMatrixConversion<T, Tint>(SHV1);
  MyMatrix<T> SHV2_T = UniversalMatrixConversion<T, Tint>(SHV2);
  int n_rows = SHV1_T.rows();
  std::vector<T> Vdiag(n_rows, T(0));
  std::optional<MyMatrix<T>> opt = TestIntEquivalence_ListMat_Vdiag<T, Tgroup>(
      SHV1_T, ListMat1, Vdiag, SHV2_T, ListMat2, Vdiag, os);
  if (!opt) {
    return {};
  }
  MyMatrix<Tint> M = UniversalMatrixConversion<Tint, T>(*opt);
  return M;
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>> ArithmeticEquivalence_inner(
    MyMatrix<T> const &inpMat1, MyMatrix<Tint> const &SHV1,
    MyMatrix<T> const &inpMat2, MyMatrix<Tint> const &SHV2, std::ostream &os) {
  std::vector<MyMatrix<T>> ListMat1{inpMat1};
  std::vector<MyMatrix<T>> ListMat2{inpMat2};
  return ArithmeticEquivalenceMultiple_inner<T,Tint,Tgroup>(ListMat1, SHV1, ListMat2, SHV2, os);
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
ArithmeticEquivalenceMultiple(std::vector<MyMatrix<T>> const &ListMat1,
                              std::vector<MyMatrix<T>> const &ListMat2,
                              std::ostream &os) {
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  MicrosecondTime time;
#endif
  /*
    The smaller of the two families on each side. The rule picking it depends
    only on the isometry class, both candidate sizes being invariants, so two
    isometric lattices pick corresponding families and the comparison below
    of the two sizes cannot reject an equivalence that exists.
   */
  MyMatrix<Tint> SHV1 =
      GetCanonicVectorFamily<T, Tint>(ListMat1[0], os).get_full();
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: GetCanonicVectorFamily1|=" << time << "\n";
#endif
  MyMatrix<Tint> SHV2 =
      GetCanonicVectorFamily<T, Tint>(ListMat2[0], os).get_full();
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "|LSEC: GetCanonicVectorFamily2|=" << time << "\n";
#endif
  return ArithmeticEquivalenceMultiple_inner<T,Tint,Tgroup>(ListMat1, SHV1, ListMat2, SHV2, os);
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>> ArithmeticEquivalence(MyMatrix<T> const &inpMat1,
                                                    MyMatrix<T> const &inpMat2,
                                                    std::ostream &os) {
  std::vector<MyMatrix<T>> ListMat1{inpMat1};
  std::vector<MyMatrix<T>> ListMat2{inpMat2};
  return ArithmeticEquivalenceMultiple<T, Tint, Tgroup>(ListMat1, ListMat2, os);
}

// clang-format off
#endif  // SRC_LATT_MATRIXCANONICALFORM_H_
// clang-format on
