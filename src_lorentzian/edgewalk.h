#ifndef INCLUDE_EDGEWALK_H
#define INCLUDE_EDGEWALK_H


#include "two_dim_lorentzian.h"
#include "coxeter_dynkin.h"
#include "vinberg_code.h"



template<typename T>
struct DataLattice {
  // Basis of L \cap N/2 L*
  MyMatrix<T> Basis;
  T eNorm;
};


template<typename T>
DataLattice<T> ComputeDataLattice(MyMatrix<T> const& G, T const& N)
{
  int n = G.rows();
  MyMatrix<T> M1 = IdentityMat<T>(n);
  MyMatrix<T> M2 = (N / 2) * Inverse(G);
  MyMatrix<T> M1_inter_M2 = IntersectionLattice(M1, M2);
  return {M1_inter_M2, N};
}


template<typename T, typename Tint>
struct FundDomainVertex {
  MyVector<T> gen;
  std::vector<MyVector<Tint>> l_roots;
};



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
  MyMatrix<T> Equa(n_root+1,dim);
  for (size_t i_root=0; i_root<n_root; i_root++) {
    MyVector<T> eP = G * UniversalVectorConversion<T,Tint>(l_ui[i_root]);
    T eScal = k.dot(eP);
    if (eScal != 0) {
      std::cerr << "The scalar product should be 0\n";
      throw TerminalException{1};
    }
    for (int i=0; i<dim; i++)
      Equa(i_root,i) = eP(i);
  }
  MyVector<T> eP = G * k;
  for (int i=0; i<dim; i++)
    Equa(n_root,i) = eP(i);
  MyMatrix<T> NSP = NullspaceMat(Equa);
  if (NSP.rows() != 1) {
    std::cerr << "The dimension should be exactly 2\n";
    throw TerminalException{1};
  }
  MyVector<T> r0 = GetMatrixRow(NSP,0);
  std::vector<FundDomainVertex<T,Tint>> l_candidates;
  bool allow_euclidean = false;
  std::vector<Possible_Extension<T>> l_extension = ComputePossibleExtensions(G, l_ui, l_norm, allow_euclidean);
  for (auto & e_extension : l_extension) {
    
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
