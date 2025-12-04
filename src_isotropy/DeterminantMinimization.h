// Copyright (C) 2023 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_DETERMINANTMINIMIZATION_H_
#define SRC_INDEFINITE_DETERMINANTMINIMIZATION_H_

// clang-format off
#include "factorizations.h"
#include "MAT_MatrixMod.h"
#include "ClassicLLL.h"
#include <algorithm>
#include <limits>
#include <utility>
#include <vector>
#include <map>
// clang-format on

template <typename T> struct ResultDetMin {
  MyMatrix<T> P;
  MyMatrix<T> Mred;
};

#ifdef DEBUG
#define DEBUG_DETERMINANT_MINIMIZATION
#endif

#ifdef DISABLE_DEBUG_DETERMINANT_MINIMIZATION
#undef DEBUG_DETERMINANT_MINIMIZATION
#endif

template <typename T> struct ResultNullspaceMod {
  int dimNSP;
  MyMatrix<T> BasisTot;
};

template <typename T>
ResultNullspaceMod<T> GetAdjustedBasis(MyMatrix<T> const &M, T const &TheMod,
                                       [[maybe_unused]] std::ostream &os) {
  int n_row = M.rows();
  // We avoid a copy by having NullspaceTrMat instead of NullspaceMat, but that
  // does not matter since M is symmetric here
#ifdef DEBUG_DETERMINANT_MINIMIZATION
  os << "DETMIN: GetAdjustedBasis, before NullspaceTrMatMod p=" << TheMod
     << "\n";
  os << "DETMIN: M=\n";
  WriteMatrix(os, M);
#endif
  MyMatrix<T> NSP = NullspaceTrMatMod(M, TheMod);
  int dimNSP = NSP.rows();
#ifdef DEBUG_DETERMINANT_MINIMIZATION
  os << "DETMIN: GetAdjustedBasis, NSP=\n";
  WriteMatrix(os, NSP);
#endif
  MyMatrix<T> Orth = NullspaceTrMat(NSP);
#ifdef DEBUG_DETERMINANT_MINIMIZATION
  os << "DETMIN: GetAdjustedBasis, Orth=\n";
  WriteMatrix(os, Orth);
#endif
  MyMatrix<T> OrthTr = Orth.transpose();
#ifdef DEBUG_DETERMINANT_MINIMIZATION
  os << "DETMIN: GetAdjustedBasis, OrthTr=\n";
  WriteMatrix(os, OrthTr);
#endif
  MyMatrix<T> NSP_b = SublatticeBasisReduction(NullspaceIntMat(OrthTr), os);
#ifdef DEBUG_DETERMINANT_MINIMIZATION
  os << "DETMIN: GetAdjustedBasis, Now NSP_b=\n";
  WriteMatrix(os, NSP_b);
#endif
  MyMatrix<T> BasisComp = SubspaceCompletionInt(NSP_b, n_row);
#ifdef DEBUG_DETERMINANT_MINIMIZATION
  os << "DETMIN: GetAdjustedBasis, After SubspaceCompletionInt\n";
#endif
  MyMatrix<T> BasisTot = Concatenate(NSP_b, BasisComp);
  return {dimNSP, BasisTot};
}

template <typename T>
bool IsMatrixNonZeroMultiple(MyMatrix<T> const &M1, MyMatrix<T> const &M2) {
  int n_row = M1.rows();
  int n_col = M1.cols();
  for (int i = 0; i < n_row; i++) {
    for (int j = 0; j < n_col; j++) {
      if (M1(i, j) != 0) {
        T scal = M2(i, j) / M1(i, j);
        if (scal == 0) {
          return false;
        }
        for (int i1 = 0; i1 < n_row; i1++) {
          for (int j1 = 0; j1 < n_col; j1++) {
            T val1 = scal * M1(i1, j1);
            if (val1 != M2(i1, j1)) {
              return false;
            }
          }
        }
        return true;
      }
    }
  }
  // M1 is zero only possibility is if M2 is zero as well
  return IsZeroMatrix(M2);
}

/*
  We apply a numbr of ideas from the preprint
  "Quadratic equations in dimensions 4, 5 and more" (P1)
  Denis Simon
  ---
  The following notions are used
  * v_p(Q) to be the multiplicity of p as a prime factor of the determinant of
  Q.
  * d = dim Ker_{F_p}(Q).
  We have the basic result d <= v.
  ---
  The ideas are clear from the manuscript except for Lemma 9, which I do not
  understand yet.
  Things are not done strictly in the same way, but the same ideas
  are implemented.
  ----------------------------
  When having "only_equivariant = true", only the operations that do not depend
  on a choice are applied.
 */
template <typename T>
ResultDetMin<T> DeterminantMinimization(MyMatrix<T> const &Q,
                                        bool only_equivariant,
                                        std::ostream &os) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  using Tring = typename underlying_ring<T>::ring_type;
  if (!IsIntegralMatrix(Q)) {
    std::cerr << "The matrix M should be integral\n";
    throw TerminalException{1};
  }
  int n = Q.rows();
#ifdef DEBUG_DETERMINANT_MINIMIZATION
  bool is_n_even = true;
  if (ResInt(n, 2) == 1) {
    is_n_even = false;
  }
#endif
  T det = DeterminantMat(Q);
  T det_abs = T_abs(det);
#ifdef DEBUG_DETERMINANT_MINIMIZATION
  os << "DETMIN: is_n_even=" << is_n_even << "\n";
  os << "DETMIN: n=" << n << " det_abs=" << det_abs << "\n";
  os << "DETMIN: Q=\n";
  WriteMatrix(os, Q);
  int iter = 0;
#endif
  Tring det_ai = UniversalScalarConversion<Tring, T>(det_abs);
  std::map<Tring, size_t> map = FactorsIntMap(det_ai);
#ifdef DEBUG_DETERMINANT_MINIMIZATION
  for (auto &kv : map) {
    os << "DETMIN: p=" << kv.first << " mult=" << kv.second << "\n";
  }
#endif
  MyMatrix<T> Qw = Q;
  MyMatrix<T> Pw = IdentityMat<T>(n);
#ifdef DEBUG_DETERMINANT_MINIMIZATION
  auto check_state = [&]() -> void {
    MyMatrix<T> ResQ = Pw * Q * Pw.transpose();
    if (!IsMatrixNonZeroMultiple(ResQ, Qw)) {
      std::cerr << "The matrix are not scalar multiple\n";
      throw TerminalException{1};
    }
  };
#endif
  std::map<T, size_t> map_lemma6;
  while (true) {
    std::vector<Tring> list_P_erase;
    bool DoSomethingGlobal = false;
#ifdef DEBUG_DETERMINANT_MINIMIZATION
    T TheDet = DeterminantMat(Qw);
    Tring TheDet_abs = UniversalScalarConversion<Tring, T>(T_abs(TheDet));
    os << "DETMIN: TheDet=" << TheDet << " TheDet_abs=" << TheDet_abs << "\n";
    Tring eProd(1);
    for (auto &kv : map) {
      for (size_t u = 0; u < kv.second; u++) {
        eProd *= kv.first;
      }
    }
    os << "DETMIN: eProd=" << eProd << "\n";
    if (eProd != TheDet_abs) {
      std::cerr << "DETMIN: Incoherence of determinant\n";
      throw TerminalException{1};
    }
#endif
    for (auto &kv : map) {
      Tring p_ring = kv.first;
      T p = UniversalScalarConversion<T, Tring>(p_ring);
      size_t &v_mult_s = kv.second;
      int v_mult_i = v_mult_s;
#ifdef DEBUG_DETERMINANT_MINIMIZATION
      T p_sqr = p * p;
      os << "DETMIN: iter=" << iter << " p=" << p << " mult=" << v_mult_i
         << "\n";
#endif
      ResultNullspaceMod<T> res = GetAdjustedBasis(Qw, p, os);
      int d_mult_i = res.dimNSP;
#ifdef DEBUG_DETERMINANT_MINIMIZATION
      os << "DETMIN: d_mult_i=" << d_mult_i << " v_mult_i=" << v_mult_i << "\n";
      os << "DETMIN: BasisTot=\n";
      WriteMatrix(os, res.BasisTot);
#endif
      auto change_basis = [&]() -> void {
        Pw = res.BasisTot * Pw;
        Qw = res.BasisTot * Qw * res.BasisTot.transpose();
#ifdef DEBUG_DETERMINANT_MINIMIZATION
        check_state();
        for (int i = 0; i < d_mult_i; i++) {
          for (int j = 0; j < n; j++) {
            T res = ResInt(Qw(i, j), p);
            if (res != 0) {
              std::cerr << "The coefficient M(i,j) is not reduced modulo p\n";
              std::cerr << "i=" << i << " j=" << j << " p=" << p << "\n";
              throw TerminalException{1};
            }
          }
        }
        int p_dim = n - d_mult_i;
        MyMatrix<T> U(p_dim, p_dim);
        for (int i = 0; i < p_dim; i++) {
          for (int j = 0; j < p_dim; j++) {
            U(i, j) = Qw(i + d_mult_i, j + d_mult_i);
          }
        }
        T det = DeterminantMat(U);
        T res_det = ResInt(det, p);
        if (res_det == 0) {
          std::cerr << "The matrix U is not invertible\n";
          throw TerminalException{1};
        }
#endif
      };
      auto get_qtilde = [&]() -> MyMatrix<T> {
        MyMatrix<T> Qtilde(d_mult_i, d_mult_i);
        for (int i = 0; i < d_mult_i; i++) {
          for (int j = 0; j < d_mult_i; j++) {
            Qtilde(i, j) = Qw(i, j) / p;
          }
        }
        return Qtilde;
      };
      bool DoSomething = false;
      //
      // Lemma 4 is eminently equivariant and do not depend a choice. Apply it if possible.
      //
      bool apply_lemma4 = true;
      if (apply_lemma4) {
        if (d_mult_i == n && !DoSomething) {
#ifdef DEBUG_DETERMINANT_MINIMIZATION
          os << "DETMIN: Apply Lemma 4\n";
#endif
          Qw = Qw / p;
          DoSomething = true;
          DoSomethingGlobal = true;
          v_mult_i -= 2 * n;
          v_mult_s -= 2 * n;
#ifdef DEBUG_DETERMINANT_MINIMIZATION
          for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
              T res = ResInt(Qw(i, j), p);
              if (res != 0) {
                std::cerr << "Qw is not divisible by p as expected\n";
                os << "Qw=\n";
                WriteMatrix(os, Qw);
                throw TerminalException{1};
              }
            }
          }
#endif
        }
      }
      //
      // Lemma 5 reduces the matrix by the kernel of the nullspace mod p.
      // So, yes, it is equivariant. More precisely, any way to apply it
      // lead us to an equivalent result.
      //
      bool apply_lemma5 = true;
      if (apply_lemma5) {
#ifdef DEBUG_DETERMINANT_MINIMIZATION
        os << "DETMIN: d_mult_i=" << d_mult_i << " v_mult_i=" << v_mult_i << "\n";
#endif
        if (d_mult_i < v_mult_i && !DoSomething) {
#ifdef DEBUG_DETERMINANT_MINIMIZATION
          os << "DETMIN: Apply Lemma 5\n";
#endif
          change_basis();
          MyMatrix<T> Qtilde = get_qtilde();
          ResultNullspaceMod<T> resB = GetAdjustedBasis(Qtilde, p, os);
          MyMatrix<T> Hmat = IdentityMat<T>(n);
          for (int i = 0; i < d_mult_i; i++) {
            for (int j = 0; j < d_mult_i; j++) {
              Hmat(i, j) = resB.BasisTot(i, j);
            }
          }
          Pw = Hmat * Pw;
          Qw = Hmat * Qw * Hmat.transpose();
          int dimNSPB = resB.dimNSP;
#ifdef DEBUG_DETERMINANT_MINIMIZATION
          check_state();
          for (int i = 0; i < dimNSPB; i++) {
            for (int j = 0; j < dimNSPB; j++) {
              T res = ResInt(Qw(i, j), p_sqr);
              if (res != 0) {
                std::cerr << "Qw[1..dimNSPB][1..dimNSPB] is not divisible by p * "
                  "p as expected\n";
                throw TerminalException{1};
              }
            }
          }
#endif
          MyMatrix<T> U = IdentityMat<T>(n);
          for (int i = 0; i < dimNSPB; i++)
            U(i, i) = 1 / p;
          Pw = U * Pw;
          Qw = U * Qw * U.transpose();
#ifdef DEBUG_DETERMINANT_MINIMIZATION
          check_state();
          if (!IsIntegralMatrix(Qw)) {
            std::cerr << "The matrix Qw is not integral\n";
            throw TerminalException{1};
          }
#endif
          DoSomething = true;
          DoSomethingGlobal = true;
          v_mult_i -= 2 * dimNSPB;
          v_mult_s -= 2 * dimNSPB;
        }
      }
      //
      // Lemma 6 also applies a nullspace reduction. So, same as Lemma 5.
      // 
      bool apply_lemma6 = true;
      if (apply_lemma6) {
        // Apply Lemma 6
        size_t &ref = map_lemma6[p];
        if (d_mult_i == v_mult_i && d_mult_i >= 2 && n > 2 * d_mult_i &&
            !DoSomething && ref == 0) {
          ref = 1;
#ifdef DEBUG_DETERMINANT_MINIMIZATION
          os << "DETMIN: Apply Lemma 6\n";
          os << "DETMIN: Before det=" << DeterminantMat(Qw)
             << " d_mult_i=" << d_mult_i << "\n";
#endif
          change_basis();
          MyMatrix<T> U = IdentityMat<T>(n);
          for (int u = d_mult_i; u < n; u++) {
            U(u, u) = p;
          }
          Pw = U * Pw;
          Qw = U * Qw * U.transpose() / p;
#ifdef DEBUG_DETERMINANT_MINIMIZATION
          check_state();
#endif
          DoSomething = true;
          DoSomethingGlobal = true;
          int dec = 2 * (n - d_mult_i) - n;
#ifdef DEBUG_DETERMINANT_MINIMIZATION
          os << "DETMIN: Before v_mult_i=" << v_mult_i << " dec=" << dec << "\n";
          os << "DETMIN: After  det(Qw)=" << DeterminantMat(Qw)
             << " det(U)=" << DeterminantMat(U) << "\n";
#endif
          v_mult_i += dec;
          v_mult_s += dec;
#ifdef DEBUG_DETERMINANT_MINIMIZATION
          os << "DETMIN: After  v_mult_i=" << v_mult_i << "\n";
#endif
        }
      }
      //
      // This Lemma relies on finding an isotropic vector mod p.
      // It very much depends on a random process and is not equivariant
      // at all.
      //
      bool apply_lemma78 = !only_equivariant;
      if (apply_lemma78) {
        // Apply Lemma 7 (joined with Lemma 8)
        if (d_mult_i == v_mult_i && !DoSomething) {
#ifdef DEBUG_DETERMINANT_MINIMIZATION
          os << "DETMIN: Apply Lemma 7 (or 8)\n";
#endif
          change_basis();
#ifdef DEBUG_DETERMINANT_MINIMIZATION
          os << "DETMIN: After change_basis\n";
#endif
          MyMatrix<T> Qtilde = get_qtilde();
#ifdef DEBUG_DETERMINANT_MINIMIZATION
          os << "DETMIN: After get_qtilde\n";
          os << "DETMIN: Qtilde=\n";
          WriteMatrix(os, Qtilde);
#endif
          std::optional<MyVector<T>> opt = FindIsotropicVectorMod(Qtilde, p);
#ifdef DEBUG_DETERMINANT_MINIMIZATION
          os << "DETMIN: After FindIsotropicVectorMod\n";
#endif
          if (opt) {
            MyVector<T> const &eV_pre = *opt;
            MyVector<T> eV = RemoveFractionVectorPlusCoeff(eV_pre).TheVect;
#ifdef DEBUG_DETERMINANT_MINIMIZATION
            os << "DETMIN: We have |eV|=" << eV.size() << "\n";
            os << "DETMIN: eV=" << StringVectorGAP(eV) << "\n";
            T val = EvaluationQuadForm(Qtilde, eV);
            T val_mod = ResInt(val, p);
            if (val_mod != 0) {
              std::cerr << "DETMIN: We have val=" << val
                        << " but val_mod=" << val_mod << "\n";
              throw TerminalException{1};
            }
#endif
            MyMatrix<T> M = ZeroMatrix<T>(1, n);
            for (int i = 0; i < d_mult_i; i++)
              M(0, i) = eV(i);
#ifdef DEBUG_DETERMINANT_MINIMIZATION
            os << "DETMIN: We have M, n=" << n << "\n";
            os << "DETMIN: M=\n";
            WriteMatrix(os, M);
#endif
            MyMatrix<T> BasisCompl = SubspaceCompletionInt(M, n);
#ifdef DEBUG_DETERMINANT_MINIMIZATION
            os << "DETMIN: We have BasisCompl\n";
#endif
            MyMatrix<T> U = Concatenate(M, BasisCompl);
#ifdef DEBUG_DETERMINANT_MINIMIZATION
            os << "DETMIN: We have U\n";
#endif
            Pw = U * Pw;
            Qw = U * Qw * U.transpose();
#ifdef DEBUG_DETERMINANT_MINIMIZATION
            check_state();
            T res = ResInt(Qw(0, 0), p_sqr);
            if (res != 0) {
              std::cerr << "We do not have Qtilde(0,0) divisible by p^2\n";
              std::cerr << "Qw=\n";
              WriteMatrix(std::cerr, Qw);
              throw TerminalException{1};
            }
#endif
            MyMatrix<T> V = IdentityMat<T>(n);
            V(0, 0) = 1 / p;
            Pw = V * Pw;
            Qw = V * Qw * V.transpose();
            DoSomething = true;
            DoSomethingGlobal = true;
            v_mult_i -= 2;
            v_mult_s -= 2;
          }
        }
      }
#ifdef DEBUG_DETERMINANT_MINIMIZATION
      os << "DETMIN: --------------------------------------\n";
#endif
      if (v_mult_s == 0) {
        list_P_erase.push_back(p_ring);
      }
    }
    for (auto &p : list_P_erase) {
      map.erase(p);
    }
    if (!DoSomethingGlobal) {
      break;
    }
#ifdef DEBUG_DETERMINANT_MINIMIZATION
    iter += 1;
#endif
  }
#ifdef DEBUG_DETERMINANT_MINIMIZATION
  os << "DETMIN: EXITING DETERMINANT MINIMIZATION\n";
  check_state();
#endif
  return {Pw, Qw};
}

// clang-format off
#endif  //  SRC_INDEFINITE_DETERMINANTMINIMIZATION_H_
// clang-format on
