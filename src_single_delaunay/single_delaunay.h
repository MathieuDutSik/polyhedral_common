// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_SINGLE_DELAUNAY_ADJACENCY_SCHEME_H_
#define SRC_SINGLE_DELAUNAY_ADJACENCY_SCHEME_H_

// clang-format off
#include "GRAPH_BitsetType.h"
#include "GRAPH_GraphicalBasic.h"
#include "POLY_LinearProgramming.h"
#include "POLY_DirectDualDesc.h"
#include "POLY_Fundamental.h"
#include "InvariantVectorFamily.h"
// clang-format on

#ifdef SANITY_CHECK
#define SANITY_CHECK_SINGLE_DELAUNAY
#endif

#ifdef DEBUG
#define DEBUG_SINGLE_DELAUNAY
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
MyMatrix<T> sd_get_gram_matrix(MyMatrix<T> const& QuadFunc) {
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
MyMatrix<T> sd_get_linear_term(MyMatrix<T> const& QuadFunc) {
  int n = sd.QuadFunc.rows() - 1;
  MyVector<T> V(n);
  for (int i=0; i<n; i++) {
    V(i) = QuadFunc(0, i+1);
  }
  return V;
}

template<typename T, typename Tint>
struct SingleDelaunayReduction {
  MyMatrix<T> NSP_T;
  MyMatrix<Tint> NSP;
  MyMatrix<Tint> Compl_nsp;
  MyMatrix<T> Compl_nsp_T;
  MyMatrix<T> GramMat_c;
  MyMatrix<T> v_lc;
  MyMatrix<T> QuadFuncRed;
};

template<typename T, typename Tint>
SingleDelaunayReduction<T,Tint> sd_get_single_delaunay_reduction(MyMatrix<T> const& QuadFunc) {
  int n = QuadFunc.rows() - 1;
  MyMatrix<T> GramMat = sd_get_gram_matrix(QuadFunc);
  MyMatrix<T> NSP_T = NullspaceIntMat(GramMat);
  MyMatrix<Tint> NSP = UniversalMatrixConversion<Tint,T>(NSP_T);
  MyVector<T> v_lin = sd_get_linear_term(QuadFunc);
  MyMatrix<Tint> Compl_nsp = SubspaceCompletionInt(NSP, n);
  MyMatrix<T> Compl_nsp_T = UniversalMatrixConversion<T,Tint>(Compl_nsp);
  MyMatrix<T> GramMat_c = Compl_nsp_T * GramMat * Compl_nsp_T.transpose();
  MyMatrix<T> v_lc = Compl_nsp_T * v_lin;
  int dim_compl = NSP.rows();
  MyMatrix<T> QuadFuncRed(dim_comp+1, dim_compl+1);
  QuadFuncRed(0, 0) = QuadFunc(0, 0);
  for (int i_compl=0; i_compl<dim_compl; i_compl++) {
    QuadFuncRed(0, i_compl+1) = v_lc(i);
    QuadFuncRed(i_compl+1, 0) = v_lc(i);
  }
  for (int i_compl=0; i_compl<dim_compl; i_compl++) {
    for (int j_compl=0; j_compl<dim_compl; j_compl++) {
      QuadFuncRed(i_compl+1, j_compl+1) = GramMat_c(i_compl, j_compl);
    }
  }
  return {NSP_T, NSP, Compl_nsp, Compl_nsp_T, GramMat_c, v_lc, QuadFuncRed};
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
  std::vector<MyVector<Tint>> ListV;
};

template<typename T, typename Tint>
struct ResultSingleDelaunay {
  std::variant<SingleDelaunay<T,Tint>,SingleDelaunayError<Tint>> res;
  std::optional<SingleDelaunayError<Tint>> get_error() {
    try {
      return res.get<SingleDelaunayError<Tint>>();
    }
    catch (const std::bad_variant_access& ex) {
      return {};
    }
  }
  std::optional<SingleDelaunay<T,Tint>> get_result() {
    try {
      return res.get<SingleDelaunay<T,Tint>>();
    }
    catch (const std::bad_variant_access& ex) {
      return {};
    }
  }
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



/*
  In the GAP code, we compute the facets of the polytope, and then the adjacent Delaunay
  polytopes. 
 */
template<typename T, typename Tint>
std::vector<MyVector<Tint>> get_adjacent_vertices(SingleDelaunay<T,Tint> const& sd, std::ostream& os) {
  int n = sd.QuadFunc.rows() - 1;
  SingleDelaunayReduction<T,Tint> sdr = sd_get_single_delaunay_reduction<T,Tint>(QuadFunc);
  int dim_compl = sdr.Compl_nsp.rows();
  MyMatrix<Tint> InvBasis = ExtractInvariantVectorFamilyZbasis<T,Tint>(sdr.GramMat_c, os);
  MyMatrix<Tint> Compl_p_NSP = Concatenate(sdr.Compl_nsp, sdr.NSP);
#ifdef SANITY_CHECK_SINGLE_DELAUNAY
  Tint det = DeterminantMat(Compl_p_NSP);
  if (T_abs(det) == 1) {
    std::cerr << "SD: The determinant should be 1\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<Tint> Compl_p_NSP_inv = Inverse(Compl_p_NSP);
  int n_ext = sd.EXT.rows();
  std::unordered_set<MyVector<Tint>> set_EXT;
  for (int i_ext=0; i_ext<n_ext; i_ext++) {
    MyVector<Tint> V1(n);
    for (int i=0; i<n; i++) {
      V1(i) = sd.EXT(i_ext, i+1);
    }
    MyVector<Tint> V2 = Compl_p_NSP_inv.transpose() * V1;
    MyVector<Tint> V3(dim_compl);
    for (int i_compl=0; i_compl<dim_compl; i_compl++) {
      V3(i_compl) = V2(i_compl);
    }
    set_EXT.insert(V3);
  }
  std::unordered_set<MyVector<Tint>> adj_set;
  int n_basis = InvBasis.rows();
  for (int i_basis=0; i_basis<n_basis; i_basis++) {
    MyVector<Tint> v_basis = GetMatrixRow(InvBasis, i_basis);
    for (auto & eV: set_EXT) {
      MyVector<Tint> eW = eV + v_basis;
      if (set_vert.count(eW) == 0) {
        adj_set.insert(eW);
      }
    }
  }
  std::vector<MyVector<Tint>> adj_list;
  for (auto & V1: adj_set) {
    MyVector<Tint> V2 = sdr.Compl_nsp.transpose() * V1;
    MyVector<Tint> V3 = ConcatenateScalVec(Tint(1), V2);
    adj_list.push_back(V3);
  }
  return adj_list;
}



template<typename T, typename Tint>
ResultQuadFuncDelaunay<T,Tint> get_single_delaunay(MyMatrix<T> const& QuadFunc, std::ostream& os) {
  MyMatrix<T> GramMat = sd_get_gram_matrix(QuadFunc);
  int n = GramMat.rows();
  // Testing for positive semidefiniteness
  if (!IsPositiveSemidefinite(sd, os)) {
    T MaxNorm(0);
    MyVector<Tint> V = GetShortVector<T,Tint>(GramMat, MaxNorm, os);
    MyVector<Tint> Vret = ConcatenateScalVec(Tint(0), V);
    std::vector<MyVector<Tint>> ListV{Vret};
    std::string reason = "not semi-definite matrix";
    SingleDelaunayError<Tint> sde{reason, ListV};
    return {sde};
  }
  // Computing the subspace
  SingleDelaunayReduction<T,Tint> sdr = sd_get_single_delaunay_reduction<T,Tint>(QuadFunc);
  MyVector<T> v_lin = sd_get_linear_term(QuadFunc);
  MyVector<T> v_l_nsp = sdr.NSP_T * v_lin;
  T cst = QuadFunc(0,0);
  if (!IsZeroVector(v_l_nsp) || cst < 0) {
    MyVector<Tint> V1 = sd_get_negative_vector(cst, v_l_nsp, os);
    MyVector<Tint> V2 = sdr.Compl_nsp.transpose() * V1;
    MyVector<Tint> V3 = ConcatenateScalVec(T(1), V2);
    std::vector<MyVector<Tint>> ListV{V3};
    std::string reason = "semi-definite, but kernel leads to a non-zero entry";
    SingleDelaunayError<Tint> sde{reason, ListV};
    return {sde};
  }
  // Forming the complement
  MyMatrix<T> GramMat_c_inv = Inverse(sdr.GramMat_c);
  // The function is now cst + 2 * v_lc * x + Q[x].
  MyVector<T> v_lc_i = - GramMat_c_inv * sdr.v_lc;
  T delta = v_lc_i.dot(sdr.v_lc);
  T cst_b = cst + delta;
  resultCVP<T, Tint> result = NearestVectors<T,Tint>(sdr.GramMat_c, v_lc_i, os);
  T cst_c = cst_b + result.TheNorm;
  if (cst_c < 0) {
    std::vector<MyVector<Tint>> ListV;
    int n_vect = result.ListVect.rows();
    for (int i_vect=0; i_vect<n_evct; i_vect++) {
      MyVector<Tint> V1 = GetMatrixRow(result.ListVect, i_vect);
      MyVector<Tint> V2 = sdr.Compl_nsp.transpose() * V1;
      MyVector<Tint> V3 = ConcatenateScalVec(Tint(1), V2);
      ListV.push_back(V3);
    }
    std::string reason = "found a shortest vectors of negative norm";
    SingleDelaunayError<Tint> sde{reason, ListV};
    return {sde};
  }
  // Now we are correct.
  if (cst_c > 0) {
    // No vector. Odd, but it can technicall happen.
    MyMatrix<Tint> EXT(0, n+1);
    SingleDelaunay sd{QuadFunc, EXT, sdr.NSP};
    return {sd};
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
  SingleDelaunay sd{QuadFunc, EXT, sdr.NSP};
  return {sd};
}





/*
  GAP source code for the flip: InfDel_LiftingDelaunay
  
  
 */
template<typename T, typename Tint>
SingleDeleunay<T,Tint> flip_evolution(SingleDelaunay<T,Tint> const& sd, MyMatrix<T> const& dir_change, std::ostream& os) {
  int n_ext = sd.EXT.rows();
  std::vector<MyVector<Tint>> ListOutVertices;
  // Finding the vertices that are not in the intersection
  for (int i_ext=0; i_ext<n_ext; i_ext++) {
    MyVector<T> eEXT = GetMatrixRow(sd.EXT, i_ext);
    T val = EvaluationQuadForm(sd.QuadFunc, eEXT);
    if (val < 0) {
      std::cerr << "SD: The direction is incorrect\n";
      throw TerminalException{1};
    }
    if (val > 0) {
      ListOutVertices.push_back(eEXT);
    }
  }
  // Finding the vertices that are not
  for (auto &eEXT: get_adjacent_vertices(sd, os)) {
    ListOutVertices.push_back(eEXT);
  }
  auto get_det=[](MyVector<T> const& EXT1, MyVector<T> const& EXT2) -> T {
    return EXT1(0) * EXT2(1) - EXT2(0) * EXT1(1);
  };
  auto get_matrix=[&](T const& a1, T const& a2) -> MyMatrix<T> {
    MyMatrix<T> M = a1 * sd.QuadFunc + a2 * dir_change;
    return M;
  };
  auto get_quad_funcs=[&]() -> std::vector<MyMatrix<T>> {
    std::unordered_set<MyVector<T>> set_ineq;
    for (auto &eEXT: ListOutVertices) {
      T scal1 = EvaluationQuadForm(sd.QuadFunc, eEXT);
      T scal2 = EvaluationQuadForm(dir_change, eEXT);
      MyVector<T> V1(2);
      V1(0) = scal1;
      V1(1) = scal2;
      MyVector<T> V2 = ScalarCanonicalizationVector(V1);
      set_ineq.insert(V2);
    }
    std::vector<MyVector<T>> l_ineq.insert(l_ineq.begin(), set_ineq.begin(), set_ineq.end());
    std::vect<T> EXT1 = l_ineq[0];
    MyVector<T> EXT2 = l_ineq[1];
    T det12 = get_det(EXT1, EXT2);
    if (det12 == 0) {
      std::cerr << "SD: The determinant (det12) should be non-zero\n";
      throw TerminalException{1};
    }
    for (size_t i_elt=2; i_elt<l_ineq.size(); i_elt++) {
      MyVector<T> EXTN = l_ineq[i_elt];
      T det1N = get_det(EXT1, EXTN);
      T det2N = get_det(EXT2, EXTN);
      if (det1N == 0 || det2N == 0) {
        std::cerr << "SD: The determinants (det1N, det2N) should be non-zero\n";
        throw TerminalException{1};
      }
      if (det1N * det2N > 0) {
        if (det12 * det1N > 0) {
          EXT2 = EXTN;
          det12 = det1N;
        } else {
          EXT1 = EXTN;
          det12 = - det2N;
        }
      }
    }
    if (det12 > 0) {
      MyMatrix<T> sd1 = get_matrix(-EXT1[1], EXT1[0]);
      MyMatrix<T> sd2 = get_matrix( EXT2[1],-EXT2[0]);
      return {sd1, sd2};
    } else {
      MyMatrix<T> sd1 = get_matrix(-EXT2[1], EXT2[0]);
      MyMatrix<T> sd2 = get_matrix( EXT1[1],-EXT1[0]);
      return {sd1, sd2};
    }
  };
  auto get_new_quad_func=[&]() -> MyMatrix<T> {
    std::vector<MyMatrix<T>> quad_funcs = get_quad_funcs();
    MyMatrix<T> m_func = ScalarCanonicalizationMatrix(sd.QuadFunc);
    for (auto &eM: quad_funcs) {
      MyMatrix<T> fM = ScalarCanonicalizationMatrix(eM);
      if (fM != m_func) {
        return fM;
      }
    }
    std::cerr << "SD: Failed to find a matrix\n";
    throw TerminalException{1};
  };
  while(true) {
    MyMatrix<T> new_quad_func = get_new_quad_func();
    ResultQuadFuncDelaunay<T,Tint> res = get_single_delaunay(new_quad_func, os);
    std::optional<SingleDelaunayError<Tint>> opt = res.get_error();
    if (opt) {
      SingleDelaunayError<Tint> const& error = *opt;
#ifdef DEBUG_SINGLE_DELAUNAY
      os << "SD: Failing for reason=" << error.reason << "\n";
#endif
      for (auto &eEXT: error.ListV) {
        ListOutVertices.push_back(eEXT);
      }
    } else {
      std::optional<SingleDelaunay<T,Tint>> opt = res.get_result();
      if (!opt) {
        std::cerr << "SD: We should have a result\n";
        throw TerminalException{1};
      }
      return *opt;
    }
  }
}


template<typename T, typename Tint, typename Tgroup>
Tgroup sd_get_group(




/*
  This is a scheme for given a polytope
  Thi
 */
template<typename T, typename Tint>
std::vector<MyMatrix<T>> 





// clang-format off
#endif  // SRC_SINGLE_DELAUNAY_ADJACENCY_SCHEME_H_
// clang-format on

