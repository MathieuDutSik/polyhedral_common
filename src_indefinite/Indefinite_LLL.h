#ifndef INCLUDE_INDEFINITE_LLL_H
#define INCLUDE_INDEFINITE_LLL_H

template<typename T>
struct ResultGramSchmidt_Indefinite {
  bool success; // true means we have a basis. False that we have an isotropic vector
  MyMatrix<T> Bstar;
  MyMatrix<T> mu;
  MyVector<T> Xisotrop;
};



template<typename T>
ResultGramSchmidt_Indefinite<T> GramSchmidtOrthonormalization(MyMatrix<T> const& M)
{
  int n = M.rows();
  MyMatrix<T> Bstar = IdentityMat<T>(n);
  MyMatrix<T> mu(n,n);
  struct PairInf {
    MyVector<T> Bistar_M;
    T Bistar_norm;
  };
  std::vector<PairInf> l_inf;
  for (int i=0; i<n; i++) {
    MyVector<T> Bistar = GetMatrixRow(Bstar, i);
    for (int j=0; j<i; j++) {
      T muij = (Bi.dot(l_inf[j].Bistar_M)) / (l_inf[j].Bistar_norm);
      mu(i,j) = muij;
      BiStar -= muij * Bistar;
    }
    MyVector<T> Bistar_M = M * BiStar;
    T scal = Bistar_M.dot(Bistar);
    if (scal == 0) {
      return {false, {}, {}, BiStar};
    }
    AssignMatrixRow(Bstar, Bistar, i);
    PairInf epair{std::move(Bistar_M), std::move(scal)};
    l_inf.push_back(epair);
  }
  return {true, Bstar, mu, {}};
}




template<typename T, typename Tint>
struct IndefiniteLLL {
  MyMatrix<T> Mred;
  MyMatrix<Tint> Pmat;
};


template<typename T, typename Tint>
struct ResultIndefiniteLLL {
  bool success; // true if we obtained the reduced matrix. false if we found an isotropic vector
  IndefiniteLLL<T,Tint> res;
  MyVector<T> Xisotrop;
};



// Adapted from Denis Simon, Solving Quadratic Equations Using Reduced
// Unimodular Quadratic Forms, Math. Comp. 74(251) 1531--1543
template<typename T, typename Tint>
ResultIndefiniteLLL<T,Tint> Indefinite_LLL(MyMatrix<T> const& M)
{
  T c = T(7) / T(8); // The c constant of the LLL algorithm
  int k = 1;
  while(true) {
    
  }
}





#endif
