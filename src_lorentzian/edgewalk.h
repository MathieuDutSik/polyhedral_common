#ifndef INCLUDE_EDGEWALK_H
#define INCLUDE_EDGEWALK_H


#include "two_dim_lorentzian.h"
#include "coxeter_dynkin.h"
#include "vinberg_code.h"



template<typename T>
MyMatrix<T> ComputeLattice_LN(MyMatrix<T> const& G, T const& N)
{
  int n = G.rows();
  MyMatrix<T> M1 = IdentityMat<T>(n);
  MyMatrix<T> M2 = (N / 2) * Inverse(G);
  return IntersectionLattice(M1, M2);
}


template<typename T, typename Tint>
struct FundDomainVertex {
  MyVector<T> gen;
  std::vector<MyVector<Tint>> l_roots;
};



template<typename T, typename Tint>
struct PossibleExtension {
  int sign; // 0 for 0, 1 for positive, -1 for negative
  T quant1; // this is (k.alpha_{N,\Delta'})^2 / R_{N,\Delta'}
  T quant2; // this is (k.alpha_{N,\Delta'})^2 / N
  T e_norm;
  MyVector<Tint> alpha;
};

template<typename T>
int get_sign_sing(T const& val)
{
  if (val > 0)
    return 1;
  if (val < 0)
    return -1;
  return 0;
}


template<typename T>
int get_sign_pair_t(T const& p1, T const& p2)
{
  if (p1 < p2)
    return 1;
  if (p1 > p2)
    return -1;
  return 0;
}



// return 0 is p1 == p2 :  1 if p1 < p2 : -1 if p1 > p2
template<typename T>
int get_sign_pair_stdpair(std::pair<int,T> const& p1, std::pair<int,T> const& p2)
{
  if (p1.first == p2.first)
    return get_sign_pair_t(p1.second, p2.second);
  return get_sign_pair_t(p1.first, p2.first);
}


template<typename T>
PossibleExtension<T> gen_possible_extension(MyMatrix<T> const& G, MyVector<T> const& k, MyVector<Tint> const& alpha, T const& res_norm, T const& e_norm)
{
  MyVector<T> alpha_T = UniversalVectorConversion<T,Tint>(alpha);
  T scal = - k.dot(G * alpha_T);
  T quant1 = (scal * scal) / res_norm;
  T quant2 = (scal * scal) / e_norm;
  return {get_sign_sing(scal), quant1, quant2, e_norm, alpha};
}


// sign as before + means preferable according to page 27.
// return 1 if poss1 is preferable to poss2
template<typename T, typename Tint>
int get_sign_poss_ext(PossibleExtension<T,Tint> const& poss1, PossibleExtension<T,Tint> const& poss2)
{
  int sign1 = get_sign_pair_stdpair({poss1.sign, poss1.quant1}, {poss2.sign, poss2.quant1});
  if (sign1 != 0)
    return sign1; // because -k.alpha1 / sqrt(R1)    <     -k.alpha2 / sqrt(R2)   correspond to 1 in the above.
  int sign2 = get_sign_pair_stdpair({poss1.sign, poss1.quant2}, {poss2.sign, poss2.quant2});
  if (sign2 != 0)
    return sign2; // because -k.alpha1 / sqrt(N1)    <     -k.alpha2 / sqrt(N2)   correspond to 1 in the above.
  int sign3 = get_sign_pair_t(poss1.e_norm, poss2.e_norm);
  return sign3; // because N1 < N2 corresponds to 1
}




/*
  We take the notations as in EDGEWALK paper.
  ---The dimension is n+1
  ---We have an edge e between two rays k and k'.
  We have u1, .... u(n-1) roots that are orthogonal and U the real span of those vectors.
  ---P is the real plane orthogonal to U.
  ---pi_U and pi_P are the corresponding projectors
  ---How is (1/2) P defined and correspond to k (typo correction)
  ---


 */
template<typename T, typename Tint>
std::option<FundDomainVertex<T,Tint>> EdgewalkProcedure(MyMatrix<T> const& G, MyVector<T> const& k, std::vector<MyVector<Tint>> l_ui, std::vector<T> const& l_norms)
{
  int dim = G.size();
  size_t n_root = l_ui.size();
  MyMatrix<T> Space(n_root,dim);
  MyMatrix<T> EquaB(n_root+1,dim);
  for (size_t i_root=0; i_root<n_root; i_root++) {
    MyVector<T> eV = UniversalVectorConversion<T,Tint>(l_ui[i_root]);
    MyVector<T> eP = G * eV;
    T eScal = k.dot(eP);
    if (eScal != 0) {
      std::cerr << "The scalar product should be 0\n";
      throw TerminalException{1};
    }
    for (int i=0; i<dim; i++) {
      Space(i_root,i) = eV(i);
      EquaB(i_root,i) = eP(i);
    }
  }
  MyVector<T> eP = G * k;
  for (int i=0; i<dim; i++)
    EquaB(n_root,i) = eP(i);
  MyMatrix<T> NSP = NullspaceMat(EquaB);
  if (NSP.rows() != 1) {
    std::cerr << "The dimension should be exactly 2\n";
    throw TerminalException{1};
  }
  MyVector<T> r0 = GetMatrixRow(NSP,0);
  std::vector<FundDomainVertex<T,Tint>> l_candidates;
  bool allow_euclidean = false;
  std::vector<Possible_Extension<T>> l_extension = ComputePossibleExtensions(G, l_ui, l_norm, allow_euclidean);
  for (auto & e_extension : l_extension) {
    T e_norm = e_extension.e_norm;
    MyMatrix<T> Latt = ComputeLattice_LN(G, e_norm);
    // Now getting into the LN space
    MyMatrix<T> Space_LN = Space * Inverse(Latt);
    MyMatrix<T> G_LN = Latt * G * Latt.transpose();
    MyMatrix<T> Equas = Space_LN * G_LN;
    MyMatrix<T> NSP = NullspaceIntMat(TransposedMat(Equas));
    MyMatrix<T> GP_LN = NSP * G_LN * NSP.transpose();
    MyVector<T> r0_LN = Inverse(Latt).transpose() * r0;
    MyVector<T> r0_NSP = SolutionMat(NSP, r0_LN);
    MyVector<Tint> r0_work = UniversalVectorConversion<Tint,T>(RemoveVectorFraction(r0_NSP));
    std::optional<MyVector<Tint>> opt_v = get_first_next_vector(GP_LN, r0_work, e_extension.res_norm);
    if (opt_v) {
      MyVector<T> v = UniversalVectorConversion<T,Tint>(*opt_v);
      
      MyVector<T> root_T = e_extension.u + v * NSP;
      MyVector<Tint> root = UniversalVectorConversion<Tint,T>(root_T);
      
      
     
    }
  }

}







template<typename T, typename Tint>
struct PairVertices {
  FundDomainVertex<T,Tint> vert1;
  FundDomainVertex<T,Tint> vert2;
};


template<typename T, typename Tint>
struct ResultEdgewalk {
  std::vector<MyMatrix<Tint>> l_gen_isom_cox;
  std::vector<PairVertices<T,Tint>> l_orbit_pair_vertices;
};





template<typename T, typename Tint>
  ResultEdgeWalk<T,Tint> LORENTZ_RunEdgewalkAlgorithm(MyMatrix<T> const& G, std::vector<T> const& Nlist_norms, FundDomainVertex<T,Tint> const& eVert)
{
  std::vector<MyMatrix<Tint>> l_gen_isom_cox;
  std::vector<PairVertices<T,Tint>> l_orbit_pair_vertices;
  auto func_insert_pair_vertices=[&](PairVertices<T,Tint> const& e_pair) -> void {
  };
  
  while(true) {
  }
}









#endif
