// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_PERFECTFORM_H_
#define SRC_PERFECT_PERFECTFORM_H_

// clang-format off>
#include "MatrixGroup.h"
#include "PolytopeEquiStab.h"
#include "EquiStabMemoization.h"
#include "shortest_flipping.h"
#include "Positivity.h"
#include "Tspace_General.h"
#include "fractions.h"
#include "POLY_RecursiveDualDesc.h"
#include <map>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_PERFECT_FORM
#define DEBUG_INITIAL_PERFECT
#define DEBUG_PERFECT_REPR
#define DEBUG_BOUNDED_FACE
#endif

#ifdef DISABLE_DEBUG_PERFECT_FORM
#undef DEBUG_PERFECT_FORM
#endif

#ifdef DISABLE_DEBUG_INITIAL_PERFECT
#undef DEBUG_INITIAL_PERFECT
#endif

#ifdef DISABLE_DEBUG_PERFECT_REPR
#undef DEBUG_PERFECT_REPR
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_INITIAL_PERFECT
#define SANITY_CHECK_PERFECT_REPR
#define SANITY_CHECK_BOUNDED_FACE
#endif

#ifdef TIMINGS
#define TIMINGS_PERFECT_FORM
#endif

template <typename T, typename Tint> struct NakedPerfect {
  MyMatrix<T> eGram;
  MyMatrix<Tint> SHV;
  MyMatrix<Tint> SHVred;
  MyMatrix<T> PerfDomEXT;
  std::vector<std::vector<int>> ListBlock;
  std::vector<int> ListPos;
};

template <typename T, typename Tint>
MyMatrix<T> GetNakedPerfectConeClassical(MyMatrix<Tint> const &M) {
  int nbRow = M.rows();
  int n = M.cols();
  int dimSymm = n * (n + 1) / 2;
  MyMatrix<T> RetMat(nbRow, dimSymm);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    MyVector<Tint> V = GetMatrixRow(M, iRow);
    MyMatrix<T> M(n, n);
    for (int u = 0; u < n; u++)
      for (int v = 0; v < n; v++)
        M(u, v) = V(u) * V(v);
    MyVector<T> Vm = SymmetricMatrixToVector(M);
    AssignMatrixRow(RetMat, iRow, Vm);
  }
  return RetMat;
}

template <typename T, typename Tint>
NakedPerfect<T, Tint> GetNakedPerfectCone(LinSpaceMatrix<T> const &LinSpa,
                                          MyMatrix<T> const &eGram,
                                          Tshortest<T, Tint> const &rec_shv,
                                          [[maybe_unused]] std::ostream &os) {
  int nbSHV = rec_shv.SHV.rows();
  std::vector<int> ListPos(nbSHV);
  int nbMat = LinSpa.ListMat.size();
  int n = eGram.rows();
  MyMatrix<T> RyshkovLoc(nbSHV, nbMat);
  for (int iSHV = 0; iSHV < nbSHV; iSHV++) {
    MyVector<Tint> eVect = GetMatrixRow(rec_shv.SHV, iSHV);
    for (int iMat = 0; iMat < nbMat; iMat++) {
      T eSum = EvaluationQuadForm<T, Tint>(LinSpa.ListMat[iMat], eVect);
      RyshkovLoc(iSHV, iMat) = eSum;
    }
  }
  //  os << "RyshkovLoc=\n";
  //  WriteMatrix(os, RyshkovLoc);
  std::map<MyVector<T>, std::vector<int>> map;
  for (int iSHV = 0; iSHV < nbSHV; iSHV++) {
    MyVector<T> V = GetMatrixRow(RyshkovLoc, iSHV);
    map[V].push_back(iSHV);
  }
  std::vector<std::vector<int>> ListBlock;
  int i_block = 0;
  for (auto &kv : map) {
    ListBlock.push_back(kv.second);
    for (auto &iSHV : kv.second) {
      ListPos[iSHV] = i_block;
    }
    i_block += 1;
  }
  int nbBlock = i_block;
#ifdef DEBUG_PERFECT_FORM
  os << "m=" << n << " nbBlock=" << nbBlock << " nbSHV=" << nbSHV
     << " nbMat=" << nbMat << "\n";
#endif
  MyMatrix<T> PerfDomEXT(nbBlock, nbMat);
  MyMatrix<Tint> SHVred(nbBlock, n);
  for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
    int iSHV = ListBlock[iBlock][0];
    MyVector<Tint> eVect = GetMatrixRow(rec_shv.SHV, iSHV);
    AssignMatrixRow(SHVred, iBlock, eVect);
    for (int iMat = 0; iMat < nbMat; iMat++) {
      PerfDomEXT(iBlock, iMat) = RyshkovLoc(iSHV, iMat);
    }
  }
  return {eGram, rec_shv.SHV, SHVred, PerfDomEXT, ListBlock, ListPos};
}

template <typename T, typename Tgroup> struct RyshkovGRP {
  MyMatrix<T> PerfDomEXT;
  Tgroup GRPsub;
  vectface ListIncd;
  std::vector<int> ListPos;
};

template <typename T, typename Tgroup>
RyshkovGRP<T, Tgroup>
GetNakedPerfectCone_GRP(LinSpaceMatrix<T> const &LinSpa,
                        MyMatrix<T> const &SHV_T,
                        Tgroup const &GRP, [[maybe_unused]] std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int nbSHV = SHV_T.rows();
  std::vector<int> ListPos(nbSHV);
  int nbMat = LinSpa.ListMat.size();
  int n = LinSpa.n;
  MyMatrix<T> RyshkovLoc(nbSHV, nbMat);
  for (int iSHV = 0; iSHV < nbSHV; iSHV++) {
    MyVector<T> eVect = GetMatrixRow(SHV_T, iSHV);
    for (int iMat = 0; iMat < nbMat; iMat++) {
      T eSum = EvaluationQuadForm<T, T>(LinSpa.ListMat[iMat], eVect);
      RyshkovLoc(iSHV, iMat) = eSum;
    }
  }
  //  os << "RyshkovLoc=\n";
  //  WriteMatrix(os, RyshkovLoc);
  std::map<MyVector<T>, std::vector<int>> map;
  for (int iSHV = 0; iSHV < nbSHV; iSHV++) {
    MyVector<T> V = GetMatrixRow(RyshkovLoc, iSHV);
    map[V].push_back(iSHV);
  }
  std::vector<std::vector<int>> ListBlock;
  int i_block = 0;
  for (auto &kv : map) {
    ListBlock.push_back(kv.second);
    for (auto &iSHV : kv.second) {
      ListPos[iSHV] = i_block;
    }
    i_block += 1;
  }
  int nbBlock = i_block;
#ifdef DEBUG_PERFECT_FORM
  os << "m=" << n << " nbBlock=" << nbBlock << " nbSHV=" << nbSHV
     << " nbMat=" << nbMat << "\n";
#endif
  MyMatrix<T> PerfDomEXT(nbBlock, nbMat);
  MyMatrix<T> SHVred(nbBlock, n);
  for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
    int iSHV = ListBlock[iBlock][0];
    MyVector<T> eVect = GetMatrixRow(SHV_T, iSHV);
    AssignMatrixRow(SHVred, iBlock, eVect);
    for (int iMat = 0; iMat < nbMat; iMat++) {
      PerfDomEXT(iBlock, iMat) = RyshkovLoc(iSHV, iMat);
    }
  }
  std::vector<Telt> l_gens;
  for (auto &eGen : GRP.GeneratorsOfGroup()) {
    std::vector<Tidx> eList(nbBlock);
    for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
      Tidx iSHV = ListBlock[iBlock][0];
      Tidx jSHV = eGen.at(iSHV);
      Tidx jBlock = ListPos[jSHV];
      eList[iBlock] = jBlock;
    }
    Telt eElt(eList);
    l_gens.emplace_back(std::move(eElt));
  }
  Tgroup GRPsub(l_gens, nbBlock);
  vectface ListIncd =
    DualDescriptionStandard<T, Tgroup>(PerfDomEXT, GRPsub, os);
  return {std::move(PerfDomEXT), std::move(GRPsub), std::move(ListIncd), std::move(ListPos)};
}

template <typename T, typename Tgroup>
Face get_big_incd(RyshkovGRP<T,Tgroup> const& ryshk, Face const& f_sma) {
  int n_shv = ryshk.ListPos.size();
  Face f_big(n_shv);
  for (int i=0; i<n_shv; i++) {
    int i_block = ryshk.ListPos[i];
    if (f_sma[i_block] == 1) {
      f_big[i] = 1;
    }
  }
  return f_big;
}

template <typename T, typename Tint, typename Tgroup>
Tgroup MapLatticeGroupToConeGroup(NakedPerfect<T, Tint> const &eNaked,
                                  Tgroup const &GRPshv) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int nbBlock = eNaked.ListBlock.size();
  std::vector<Telt> ListGen;
  std::vector<Telt> LGen = GRPshv.GeneratorsOfGroup();
  for (auto &eGen : LGen) {
    std::vector<Tidx> v(nbBlock);
    for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
      int iSHV = eNaked.ListBlock[iBlock][0];
      int jSHV = OnPoints(iSHV, eGen);
      int jBlock = eNaked.ListPos[jSHV];
      v[iBlock] = jBlock;
    }
    ListGen.push_back(Telt(v));
  }
  return Tgroup(ListGen, nbBlock);
}




template <typename T, typename Tint>
std::pair<MyMatrix<T>, Tshortest<T, Tint>>
Flipping_Perfect(MyMatrix<T> const &eMatIn, MyMatrix<T> const &eMatDir,
                 std::ostream &os) {
  auto f_admissible = [&os](MyMatrix<T> const &eMat) -> bool {
    return IsPositiveDefinite<T>(eMat, os);
  };
  auto f_shortest = [&os](MyMatrix<T> const &eMat) -> Tshortest<T, Tint> {
    return T_ShortestVectorHalf<T, Tint>(eMat, os);
  };
  return Kernel_Flipping_Perfect<T, Tint, decltype(f_admissible),
                                 decltype(f_shortest)>(f_admissible, f_shortest,
                                                       eMatIn, eMatDir, os);
}

template <typename T, typename Tint>
MyMatrix<T> get_scal_mat(LinSpaceMatrix<T> const &LinSpa,
                         Tshortest<T, Tint> const &rec_shv) {
  int nbMat = LinSpa.ListMat.size();
  int nbShort = rec_shv.SHV.rows();
  MyMatrix<T> ScalMat(nbShort, nbMat);
  for (int iShort = 0; iShort < nbShort; iShort++) {
    MyVector<Tint> eVectShort = rec_shv.SHV.row(iShort);
    for (int iMat = 0; iMat < nbMat; iMat++) {
      T eNorm = EvaluationQuadForm<T, Tint>(LinSpa.ListMat[iMat], eVectShort);
      ScalMat(iShort, iMat) = eNorm;
    }
  }
  return ScalMat;
}

template <typename T, typename Tint>
bool is_perfect_in_space(LinSpaceMatrix<T> const &LinSpa,
                         Tshortest<T, Tint> const &rec_shv) {
  int nbMat = LinSpa.ListMat.size();
  MyMatrix<T> ScalMat = get_scal_mat<T, Tint>(LinSpa, rec_shv);
  return RankMat(ScalMat) == nbMat;
}

/*
  Look for a positive definite matrix in the T-space such that
  lambda = A[v] for v in SHV
 */
template<typename T, typename Tint>
bool is_bounded_face_iterative(LinSpaceMatrix<T> const &LinSpa, MyMatrix<Tint> const& SHV, std::ostream& os) {
  int n = LinSpa.n;
  int n_vect = SHV.rows();
  int n_mat = LinSpa.ListMat.size();
  MyMatrix<T> MatScal(n_mat, n_vect);
  for (int i_mat=0; i_mat<n_mat; i_mat++) {
    MyVector<Tint> V0 = GetMatrixRow(SHV, 0);
    T val0 = EvaluationQuadForm<T, Tint>(LinSpa.ListMat[i_mat], V0);
    for (int i_vect=0; i_vect<n_vect-1; i_vect++) {
      MyVector<Tint> V = GetMatrixRow(SHV, i_vect+1);
      T val = EvaluationQuadForm<T, Tint>(LinSpa.ListMat[i_mat], V);
      MatScal(i_mat, i_vect) = val - val0;
    }
  }
  MyMatrix<T> NSP = NullspaceMat(MatScal);
  std::vector<MyMatrix<T>> BasisSpace;
  for (int i_nsp=0; i_nsp<NSP.rows(); i_nsp++) {
    MyMatrix<T> mat = ZeroMatrix<T>(n, n);
    for (int i_mat=0; i_mat<n_mat; i_mat++) {
      mat += NSP(i_nsp, i_mat) * LinSpa.ListMat[i_mat];
    }
    BasisSpace.push_back(mat);
  }
  std::optional<MyMatrix<T>> opt = GetOnePositiveDefiniteMatrix<T,Tint>(BasisSpace, os);
  if (opt) {
    return false;
  } else {
    return true;
  }
}

/*
  Look for a positive semi-definite matrix in the T-space such that
  0 = A[v] for v in SHV
 */
template<typename T, typename Tint>
bool is_bounded_face_iterative_bis(LinSpaceMatrix<T> const &LinSpa, MyMatrix<Tint> const& SHV, std::ostream& os) {
#ifdef DEBUG_BOUNDED_FACE
  os << "PERF: is_bounded_face_iterative_bis, step 1\n";
#endif
  int n = LinSpa.n;
  int n_vect = SHV.rows();
  int n_mat = LinSpa.ListMat.size();
#ifdef DEBUG_BOUNDED_FACE
  os << "PERF: is_bounded_face_iterative_bis, step 2\n";
#endif
  MyMatrix<T> MatScal(n_mat, n_vect);
  for (int i_mat=0; i_mat<n_mat; i_mat++) {
    for (int i_vect=0; i_vect<n_vect; i_vect++) {
      MyVector<Tint> V = GetMatrixRow(SHV, i_vect);
      T val = EvaluationQuadForm<T, Tint>(LinSpa.ListMat[i_mat], V);
      MatScal(i_mat, i_vect) = val;
    }
  }
#ifdef DEBUG_BOUNDED_FACE
  os << "PERF: is_bounded_face_iterative_bis, step 3\n";
#endif
  MyMatrix<T> NSP = NullspaceMat(MatScal);
  std::vector<MyMatrix<T>> BasisSpace;
#ifdef DEBUG_BOUNDED_FACE
  os << "PERF: is_bounded_face_iterative_bis, step 4\n";
#endif
  for (int i_nsp=0; i_nsp<NSP.rows(); i_nsp++) {
    MyMatrix<T> mat = ZeroMatrix<T>(n, n);
    for (int i_mat=0; i_mat<n_mat; i_mat++) {
      mat += NSP(i_nsp, i_mat) * LinSpa.ListMat[i_mat];
    }
#ifdef SANITY_CHECK_BOUNDED_FACE
    for (int i_vect=0; i_vect<n_vect; i_vect++) {
      MyVector<Tint> V = GetMatrixRow(SHV, i_vect);
      T val = EvaluationQuadForm<T, Tint>(mat, V);
      if (val != 0) {
        std::cerr << "PERF: The vector should have value 0\n";
        throw TerminalException{1};
      }
    }
#endif
    BasisSpace.push_back(mat);
  }
#ifdef DEBUG_BOUNDED_FACE
  os << "PERF: is_bounded_face_iterative_bis, step 5\n";
#endif
  std::optional<MyMatrix<T>> opt = GetOnePositiveSemiDefiniteMatrix<T,Tint>(BasisSpace, os);
  if (opt) {
    return false;
  } else {
    return true;
  }
}

/*
  Look for a positive semi-definite matrix in the T-space such that
  Min(A) = SHV
 */
template<typename T, typename Tint>
bool is_bounded_face_iterative_thi(LinSpaceMatrix<T> const &LinSpa, MyMatrix<Tint> const& SHV, std::ostream& os) {
  // We look for the full stabilizer of the subspace and the quadratic forms.
  MyMatrix<T> TheGramMat = SHORT_GetGram<T,Tint>(SHV, os);
  std::vector<MyMatrix<T>> ListMat = LinSpa.Basis;
  ListMat.push_back(TheGramMat);

  

}






/*
  A set of vectors define a bounded face if the conditions
  A in LinSpa, Min(A) = SHV
  defines a bounded face of the complex.
  ---
  What do we mean by bounded face?
  * The condition A >= 0,
    A = t_1 A_1 + .... + t_m A_m
    and A[v] = 1 for v in SHV
    implies that Tr(A) <= C for some constant C.
  * An equivalent condition is that
    A in LinSpa, A[v] = 0 for v in SHV
    A >= 0 implies that A = 0.

  For the cone of positive definite matrices, the necessary
  and sufficient condition is that SHV is of full rank.

  For self-dual, we have a fast algorithm.

  For non self-dual, we have a more expensive iterative
  algorithm.
 */
template<typename T, typename Tint>
bool is_bounded_face(LinSpaceMatrix<T> const &LinSpa, MyMatrix<Tint> const& SHV, std::ostream& os) {
  int n = LinSpa.n;
  int rnk = RankMat(SHV);
  if (rnk == n) {
    // For the case of classic Voronoi, that would be a necessary and
    // sufficient condition.
    return true;
  }
  int n_vect = SHV.rows();
  int n_mat = LinSpa.ListMat.size();
  if (LinSpa.self_dual_info) {
#ifdef SANITY_CHECK_BOUNDED_FACE
    bool is_bounded_iterative = is_bounded_face_iterative(LinSpa, SHV, os);
# ifdef DEBUG_BOUNDED_FACE
    os << "TSPACE: is_bounded_iterative=" << is_bounded_iterative << "\n";
# endif
#endif
    // For self-dual, we can apply a more direct algorithm.
    SelfDualInfo<T> const& self_dual_info = *LinSpa.self_dual_info;
    MyVector<T> SumRay = ZeroVector<T>(n_mat);
    for (int i_vect=0; i_vect<n_vect; i_vect++) {
      MyVector<Tint> V = GetMatrixRow(SHV, i_vect);
      for (int i_mat=0; i_mat<n_mat; i_mat++) {
        T val = EvaluationQuadForm<T, Tint>(LinSpa.ListMat[i_mat], V);
        SumRay(i_mat) += val;
      }
    }
    MyVector<T> ExprBasis = self_dual_info.PairwiseScalarInv * SumRay;
    MyMatrix<T> SumMat = ZeroMatrix<T>(n,n);
    for (int i_mat=0; i_mat<n_mat; i_mat++) {
      SumMat += ExprBasis(i_mat) * LinSpa.ListMat[i_mat];
    }
#ifdef SANITY_CHECK_BOUNDED_FACE
    for (int i_mat=0; i_mat<n_mat; i_mat++) {
      MyMatrix<T> const& eMat = LinSpa.ListMat[i_mat];
      T val1(0);
      for (int i_vect=0; i_vect<n_vect; i_vect++) {
        MyVector<Tint> V = GetMatrixRow(SHV, i_vect);
        val1 += EvaluationQuadForm<T, Tint>(LinSpa.ListMat[i_mat], V);
      }
      T val2 = frobenius_inner(eMat, SumMat);
      if (val1 != val2) {
        std::cerr << "PERF: inconsistency in the matrix values\n";
        //        throw TerminalException{1};
      }
    }
    if (!IsPositiveSemiDefinite(SumMat, os)) {
      std::cerr << "The matrix should be positive semi-definite\n";
      std::cerr << "So, maybe the cone is not self-dual\n";
      throw TerminalException{1};
    }
#endif
    int rnk = RankMat(SumMat);
    bool is_bounded_self_dual = rnk == n;
#ifdef SANITY_CHECK_BOUNDED_FACE
    if (is_bounded_self_dual != is_bounded_iterative) {
      std::cerr << "TSPACE: Incoherent results\n";
      std::cerr << "TSPACE: rnk=" << rnk << " n=" << n << "\n";
      std::cerr << "TSPACE: is_bounded_self_dual=" << is_bounded_self_dual << "\n";
      std::cerr << "TSPACE: is_bounded_iterative=" << is_bounded_iterative << "\n";
      throw TerminalException{1};
    }
#endif
    return is_bounded_self_dual;
  }
  // Not self-dual. So use instead the iterative
  return is_bounded_face_iterative(LinSpa, SHV, os);
}



/*
  Finds a perfect form by starting from the super matrix of the space
  and then doing some moves if that is not perfect.
  ---
  The direction of the movements are obtained by computing the
  shortest vectors, getting the kernel.
  The vector in the kernel gives direction of change.
  However, if that direction is positive definite then we will not obtain
  a new form. This scenario is explained in "Enumerating Perfect Forms"
  arXiv:0901.1587
  ---
  During the running of this test, the minimal norm of the matrix remains the
  same, but the number of vectors attaining it is increasing.
 */
template <typename T, typename Tint>
std::pair<MyMatrix<T>, Tshortest<T, Tint>>
GetOnePerfectForm(LinSpaceMatrix<T> const &LinSpa, std::ostream &os) {
  int nbMat = LinSpa.ListMat.size();
  MyMatrix<T> ThePerfMat = LinSpa.SuperMat;
  Tshortest<T, Tint> rec_shv = T_ShortestVectorHalf<T, Tint>(ThePerfMat, os);
#ifdef SANITY_CHECK_INITIAL_PERFECT
  T TheMin = rec_shv.min;
#endif
#ifdef DEBUG_INITIAL_PERFECT
  int iter = 0;
#endif
  while (true) {
    MyMatrix<T> ScalMat = get_scal_mat<T, Tint>(LinSpa, rec_shv);
    SelectionRowCol<T> eSelect = TMat_SelectRowCol(ScalMat);
    int TheRank = eSelect.TheRank;
#ifdef DEBUG_INITIAL_PERFECT
    os << "PERF: GetOnePerfectForm, iter=" << iter << " min=" << rec_shv.min
       << " |SHV|=" << rec_shv.SHV.rows() << "\n";
#endif
    if (TheRank == nbMat) {
#ifdef DEBUG_INITIAL_PERFECT
      os << "PERF: GetOnePerfectForm, returning at iter=" << iter << "\n";
#endif
      return {std::move(ThePerfMat), std::move(rec_shv)};
    }
    MyVector<T> V = eSelect.NSP.row(0);
    auto iife_get_dir = [&]() -> MyMatrix<T> {
      MyMatrix<T> M = LINSPA_GetMatrixInTspace(LinSpa, V);
      if (IsPositiveDefinite<T>(M, os)) {
        // For a positive definite matrix, we need to take the opposite
        // because we need a direction outside of the cone.
        return -M;
      }
      return M;
    };
    MyMatrix<T> DirMat = iife_get_dir();
#ifdef DEBUG_INITIAL_PERFECT
    os << "PERF: GetOnePerfectForm, iter=" << iter << " DirMat=\n";
    WriteMatrix(os, DirMat);
#endif
    auto pair = Flipping_Perfect<T, Tint>(ThePerfMat, DirMat, os);
    ThePerfMat = pair.first;
    rec_shv = pair.second;
#ifdef SANITY_CHECK_INITIAL_PERFECT
    if (rec_shv.min != TheMin) {
      std::cerr << "PERF: The rec_shv minimum should remain invariant\n";
      throw TerminalException{1};
    }
#endif
#ifdef DEBUG_INITIAL_PERFECT
    iter += 1;
#endif
  }
}

template <typename T, typename Tint> struct SimplePerfect {
  MyMatrix<T> Gram;
  Tshortest<T, Tint> rec_shv;
};

template <typename T, typename Tint>
std::istream &operator>>(std::istream &is, SimplePerfect<T, Tint> &obj) {
  MyMatrix<T> eG = ReadMatrix<T>(is);
  MyMatrix<Tint> SHV = ReadMatrix<Tint>(is);
  MyVector<Tint> V = GetMatrixRow(SHV, 0);
  T min = EvaluationQuadForm<T, Tint>(eG, V);
  Tshortest<T, Tint> rec_shv{min, SHV};
  obj = {eG, rec_shv};
  return is;
}

template <typename T, typename Tint>
std::ostream &operator<<(std::ostream &os, SimplePerfect<T, Tint> const &obj) {
  WriteMatrix(os, obj.Gram);
  WriteMatrix(os, obj.rec_shv.SHV);
  return os;
}

template <typename T, typename Tint>
MyMatrix<T> conversion_and_duplication(MyMatrix<Tint> const &SHV) {
  int dim = SHV.cols();
  int nbSHV = SHV.rows();
  MyMatrix<T> SHV_T(2 * nbSHV, dim);
  for (int i_row = 0; i_row < nbSHV; i_row++) {
    for (int i = 0; i < dim; i++) {
      T val = UniversalScalarConversion<T, Tint>(SHV(i_row, i));
      SHV_T(2 * i_row, i) = val;
      SHV_T(2 * i_row + 1, i) = -val;
    }
  }
  return SHV_T;
}

template <typename T, typename Tint>
MyMatrix<T> get_shv_t(MyMatrix<T> const &eMat, Tshortest<T, Tint> const &rec_shv,
                      std::ostream &os) {
  MyMatrix<T> SHVorig_T = UniversalMatrixConversion<T, Tint>(rec_shv.SHV);
  if (IsFullDimZbasis(rec_shv.SHV)) {
    return conversion_and_duplication<T, Tint>(rec_shv.SHV);
  }
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(eMat, os);
  return UniversalMatrixConversion<T, Tint>(SHV);
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>> SimplePerfect_TestEquivalence(
    LinSpaceMatrix<T> const &LinSpa, MyMatrix<T> const &eMat1,
    MyMatrix<T> const &eMat2, Tshortest<T, Tint> const &rec_shv1,
    Tshortest<T, Tint> const &rec_shv2, std::ostream &os) {
  MyMatrix<T> SHV1_T = get_shv_t(eMat1, rec_shv1, os);
  MyMatrix<T> SHV2_T = get_shv_t(eMat2, rec_shv2, os);
#ifdef SANITY_CHECK_PERFECT_REPR
  if (has_duplication(SHV1_T)) {
    std::cerr << "PERF: SHV1 has duplication\n";
    throw TerminalException{1};
  }
  if (has_duplication(SHV2_T)) {
    std::cerr << "PERF: SHV2 has duplication\n";
    throw TerminalException{1};
  }
  if (!is_antipodal(SHV1_T)) {
    std::cerr << "PERF: SHV1 is not antipodal\n";
    throw TerminalException{1};
  }
  if (!is_antipodal(SHV2_T)) {
    std::cerr << "PERF: SHV2 is not antipodal\n";
    throw TerminalException{1};
  }
#endif
#ifdef DEBUG_PERFECT_REPR
  os << "PERF: TestEquivalence, before LINSPA_TestEquivalenceGramMatrix_SHV\n";
  os << "PERF: TestEquivalence, eMat1=\n";
  WriteMatrix(os, RemoveFractionMatrix(eMat1));
  os << "PERF: TestEquivalence, eMat2=\n";
  WriteMatrix(os, RemoveFractionMatrix(eMat2));
#endif
  std::optional<MyMatrix<T>> opt =
      LINSPA_TestEquivalenceGramMatrix_SHV<T, Tgroup>(LinSpa, eMat1, eMat2,
                                                      SHV1_T, SHV2_T, os);
#ifdef DEBUG_PERFECT_REPR
  os << "PERF: TestEquivalence, det1=" << DeterminantMat(eMat1)
     << " det2=" << DeterminantMat(eMat2)
     << " opt.has_value()=" << opt.has_value() << "\n";
#endif
  if (!opt) {
    return {};
  }
  MyMatrix<T> const &M_T = *opt;
  MyMatrix<Tint> M = UniversalMatrixConversion<Tint, T>(M_T);
  return M;
}

template <typename T, typename Tint>
size_t
SimplePerfect_Invariant(size_t const &seed, LinSpaceMatrix<T> const &LinSpa,
                        MyMatrix<T> const &eMat,
                        Tshortest<T, Tint> const &rec_shv, std::ostream &os) {
  MyMatrix<T> SHV_T = get_shv_t(eMat, rec_shv, os);
  return LINSPA_Invariant_SHV<T>(seed, LinSpa, eMat, SHV_T, os);
}

template <typename T, typename Tint, typename Tgroup>
std::pair<Tgroup, std::vector<MyMatrix<Tint>>>
SimplePerfect_Stabilizer(LinSpaceMatrix<T> const &LinSpa,
                         MyMatrix<T> const &eMat,
                         Tshortest<T, Tint> const &rec_shv, std::ostream &os) {
  //
  // Functionality for checking quality of equivalences
  //
  MyMatrix<T> SHVorig_T = conversion_and_duplication<T, Tint>(rec_shv.SHV);
  MyMatrix<T> SHV_T = get_shv_t(eMat, rec_shv, os);
  Result_ComputeStabilizer_SHV<T, Tgroup> result =
      LINSPA_ComputeStabilizer_SHV<T, Tgroup>(LinSpa, eMat, SHV_T, os);
#ifdef DEBUG_PERFECT_FORM
  os << "PERFECT: SimplePerfect_Stabilizer, we have result\n";
  os << "PERFECT: SimplePerfect_Stabilizer |SHVorig_T|=" << SHVorig_T.rows()
     << " |SHV_T|=" << SHV_T.rows() << "\n";
#endif
  if (TestEqualityMatrix(SHVorig_T, SHV_T)) {
    // This is the most likely scenario: The original
    // SHVorig is adequate
#ifdef DEBUG_PERFECT_FORM
    os << "PERFECT: SimplePerfect_Stabilizer(A), Case SHVorig_T == SHV_T\n";
#endif
    std::vector<MyMatrix<T>> l_matr_t =
        result.get_list_matrix(SHV_T, eMat, LinSpa, os);
#ifdef DEBUG_PERFECT_FORM
    os << "PERFECT: SimplePerfect_Stabilizer(A), We have l_matr_t\n";
#endif
    std::vector<MyMatrix<Tint>> l_matr =
        UniversalStdVectorMatrixConversion<Tint, T>(l_matr_t);
#ifdef DEBUG_PERFECT_FORM
    os << "PERFECT: SimplePerfect_Stabilizer(A), We have l_matr\n";
#endif
    Tgroup GRP = result.get_perm_group(SHV_T, os);
#ifdef DEBUG_PERFECT_FORM
    os << "PERFECT: SimplePerfect_Stabilizer(A), We have GRP\n";
#endif
    return {std::move(GRP), std::move(l_matr)};
  } else {
#ifdef DEBUG_PERFECT_FORM
    os << "PERFECT: SimplePerfect_Stabilizer(B), Case SHVorig_T != SHV_T\n";
#endif
    std::vector<MyMatrix<T>> l_matr_t =
        result.get_list_matrix(SHV_T, eMat, LinSpa, os);
#ifdef DEBUG_PERFECT_FORM
    os << "PERFECT: SimplePerfect_Stabilizer(B), We have l_matr_t\n";
#endif
    std::vector<MyMatrix<Tint>> l_matr =
        UniversalStdVectorMatrixConversion<Tint, T>(l_matr_t);
#ifdef DEBUG_PERFECT_FORM
    os << "PERFECT: SimplePerfect_Stabilizer(B), We have l_matr\n";
#endif
    Tgroup GRP =
        get_perm_group_from_list_matrices<T, Tgroup>(l_matr_t, SHVorig_T, os);
#ifdef DEBUG_PERFECT_FORM
    os << "PERFECT: SimplePerfect_Stabilizer(B), We have GRP\n";
#endif
    return {std::move(GRP), std::move(l_matr)};
  }
}

// clang-format off
#endif  // SRC_PERFECT_PERFECTFORM_H_
// clang-format on
