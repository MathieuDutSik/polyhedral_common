// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_SINGLE_DELAUNAY_ADJACENCY_SCHEME_H_
#define SRC_SINGLE_DELAUNAY_ADJACENCY_SCHEME_H_

// clang-format off
#include "GRAPH_BitsetType.h"
#include "GRAPH_GraphicalBasic.h"
#include "POLY_LinearProgramming.h"
#include "POLY_DirectDualDesc.h"
#include "POLY_Fundamental.h"
// clang-format on

#ifdef SANITY_CHECK
#define SANITY_CHECK_SINGLE_DELAUNAY
#endif

template<typename T>
struct SingleDelaunaySpace {
  std::vector<MyMatrix<T>> l_matrices;
};

template<typename T, typename Tint>
struct SingleDelaunay {
  MyMatrix<T> QuadFunc; // Of dimension n+1
  MyMatrix<Tint> EXT; // Of dimension n+1
  MyMatrix<Tint> NSP; // Of dimension n (not n+1)
};


template<typename T>
MyMatrix<T> sd_get_gram_matrix(MyMatrix<T> const& QuadFunc, std::ostream& os) {
  int n = QuadFunc.rows() - 1;
  MyMatrix<T> G(n, n);
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      G(i, j) = QuadFunc(i+1, j+1);
    }
  }
  return G;
}

template<typename T>
MyMatrix<T> sd_get_linear_term(MyMatrix<T> const& QuadFunc, std::ostream& os) {
  int n = sd.QuadFunc.rows() - 1;
  MyVector<T> V(n);
  for (int i=0; i<n; i++) {
    V(i) = QuadFunc(0, i+1);
  }
  return V;
}

template<typename T, typename Tint>
void check_coherency_single_delaunay(SingleDelaunay<T,Tint> const& sd, std::ostream& os) {
  int n_ext = sd.EXT.rows();
  for (int i_ext=0; i_ext<n_ext; i_ext++) {
    MyVector<T> eEXT = GetMatrixRow(sd.EXT, i_ext);
    T val = EvaluationQuadForm(sd.QuadFunc, eEXT);
    if (val != 0) {
      std::cerr << "SD: The vector is not a vertex\n";
      throw TerminalException{1};
    }
  }
  int dim_nsp = sd.NSP.rows();
  for (int i_nsp=0; i_nsp<dim_nsp; i_nsp++) {
    MyVector<T> eNSP = GetMatrixRow(sd.NSP, i_nsp);
    T val = EvaluationQuadForm(sd.QuadFunc, eNSP);
    if (val != 0) {
      std::cerr << "SD: The direction is not a vertex\n";
      throw TerminalException{1};
    }
  }
}

template<typename Tint>
struct SingleDelaunayError {
  std::string reason;
  MyVector<Tint> V;
};

template<typename T, typename Tint>
struct ResultSingleDelaunay {
  std::variant<SingleDelaunay<T,Tint>,SingleDelaunayDefect<Tint>> res;
};

template<typename T>
MyVector<T> ConcatenateScalVec(T const& val, MyVector<T> const& v) {
  int len = v.size();
  MyVector<T> w(len + 1);
  w(0) = val;
  for (int i=0; i<len; i++) {
    w(i + 1) = v(i);
  }
  return w;
}

/*
  Find a vector x such that
  a + 2 v.x < 0
 */
template<typename T, typename Tint>
MyVector<Tint> sd_get_negative_vector(T const& a, MyVector<T> const& v_lin, std::ostream& os) {
  int len = v.size();
  if (a < 0) {
    return ZeroVector<Tint>(len);
  }
  // We have a >= 0.
  Myvector<Tint> Vret = ZeroVector<Tint>(len);
  Tint norm(0);
  auto get_val=[](T const& val) -> Tint {
    Tint val_i(1);
    T val_T(1);
    while(true) {
      if (val_T > val) {
        return val_i;
      }
      val_i += 1;
      val_T += 1;
    }
  };
  for (int i=0; i<len; i++) {
    if (v_lin(i) != 0) {
      T v_abs = T_abs(v_lin(i));
      int sign = T_sign(v_lin(i));
      T val = a / (2 * v_abs);
      // Suppose sign=1.
      // a + 2 v_abs x(i) < 0
      // a < 2 v_abs ( - x(i) )
      // a / (2 v_abs) < - x(i)
      //
      Tint val_i = get_val(val);
      Myvector<Tint> Vtest = ZeroVector<Tint>(len);
      Vtest(i) = - sign * val_i;
#ifdef SANITY_CHECK_SINGLE_DELAUNAY
      MyVector<T> Vtest_T = UniversalVectorConversion<T,Tint>(Vtest);
      T scal = Vtest_T.dot(v_lin);
      T obj_val = a + 2 * scal;
      if (obj_val >= 0) {
        std::cerr << "SD: The value should be strictly negative\n";
        throw TerminalException{1};
      }
#endif
      Tint test_norm = L1_norm(Vtest);
      if (norm == 0) {
        Vret = Vtest;
        norm = test_norm;
      } else {
        if (test_norm < norm) {
          Vret = Vtest;
          norm = test_norm;
        }
      }
    }
  }
  return Vret;
}

template<typename T, typename Tint>
ResultQuadFuncDelaunay<T,Tint> get_single_delaunay(MyMatrix<T> const& QuadFunc, std::ostream& os) {
  MyMatrix<T> GramMat = sd_get_gram_matrix(QuadFunc, os);
  int n = GramMat.rows();
  // Testing for positive semidefiniteness
  if (!IsPositiveSemidefinite(sd, os)) {
    T MaxNorm(0);
    MyVector<Tint> V = GetShortVector<T,Tint>(GramMat, MaxNorm, os);
    MyVector<T> Vret = ConcatenateScalVec(T(0), V);
    std::string reason = "not semi-definite matrix";
    SingleDelaunayError<Tint> sde{reason, Vret};
    return sde;
  }
  // Computing the subspace
  MyMatrix<T> NSP_T = NullspaceIntMat(GramMat);
  MyMatrix<Tint> NSP = UniversalMatrixConversion<Tint,T>(NSP_T);
  MyVector<T> v_lin = sd_get_linear_term(QuadFunc, os);
  MyVector<T> v_l_nsp = NSP_T * v_lin;
  T cst = QuadFunc(0,0);
  if (!IsZeroVector(v_l_nsp) || cst < 0) {
    MyVector<Tint> V1 = sd_get_negative_vector(cst, v_l_nsp, os);
    MyVector<Tint> V2 = NSP.transpose() * V1;
    MyVector<T> V3 = ConcatenateScalVec(T(1), V2);
    std::string reason = "semi-definite, but kernel leads to a non-zero entry";
    SingleDelaunayError<Tint> sde{reason, V3};
    return sde;
  }
  // Forming the complement
  MyMatrix<Tint> Compl_nsp = SubspaceCompletionInt(NSP, n);
  MyMatrix<T> Compl_nsp_T = UniversalMatrixConversion<T,Tint>(Compl_nsp);
  MyMatrix<T> GramMat_c = Compl_nsp_T * GramMat * Compl_nsp_T.transpose();
  MyMatrix<T> GramMat_c_inv = Inverse(GramMat_c);
  MyMatrix<T> v_lc = Compl_nsp_T * v_lin;
  // The function is now cst + 2 * v_lc * x + Q[x].
  MyVector<T> v_lc_i = - GramMat_c_inv * v_lc;
  T delta = v_lc_i.dot(v_lc);
  T cst_b = cst + delta;
  resultCVP<T, Tint> result = NearestVectors<T,Tint>(GramMat_c, v_lc_i, os);
  T cst_c = cst_b + result.TheNorm;
  if (cst_c < 0) {
    MyVector<Tint> V1 = GetMatrixRow(result.ListVect, 0);
    MyVector<Tint> V2 = NSP.transpose() * V1;
    MyVector<T> V3 = ConcatenateScalVec(T(1), V2);
    std::string reason = "found a shortest vectors of negative norm";
    SingleDelaunayError<Tint> sde{reason, V3};
    return sde;
  }
  // Now we are correct.
  if (cst_c > 0) {
    // No vector. Odd, but it can technicall happen.
    MyMatrix<Tint> EXT(0, n+1);
    SingleDelaunay sd{QuadFunc, EXT, NSP};
    return sd;
  }
  // The vectors to be returned
  int n_ext = result.ListVect.rows();
  MyMatrix<Tint> EXTred = result.ListVect * Compl_nsp;
  MyMatrix<Tint> EXT(n_ext, n + 1);
  for (int i_ext=0; i_ext<n_ext; i_ext++) {
    EXT(i_ext, 0) = 1;
    for (int i=0; i<n; i++) {
      EXT(i_ext, i+1) = EXTred(i_ext, i);
    }
  }
  SingleDelaunay sd{QuadFunc, EXT, NSP};
  return sd;
}



template<typename T, typename Tint>
SingleDeleunay<T,Tint> flip_evolution(SingleDelaunay<T,Tint> const& sd, MyMatrix<T> const& quad_evol) {
  int n_ext = sd.EXT.rows();
  for (int i_ext=0; i_ext<n_ext; i_ext++) {
    MyVector<T> eEXT = GetMatrixRow(sd.EXT, i_ext);
  }

}



/*
  This is a scheme for given a polytope
  Thi

  
 */
template<typename T>
std::vector<





// clang-format off
#endif  // SRC_SINGLE_DELAUNAY_ADJACENCY_SCHEME_H_
// clang-format on

