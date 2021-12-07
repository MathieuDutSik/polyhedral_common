#ifndef INCLUDE_INDEFINITE_LLL_H
#define INCLUDE_INDEFINITE_LLL_H

template<typename T>
struct ResultGramSchmidt_Indefinite {
  bool success; // true means we have a basis. False that we have an isotropic vector
  std::vector<MyVector<T>> Bstar;
  std::vector<T> Bstar_norms;
  MyMatrix<T> mu;
  MyVector<T> Xisotrop;
};



template<typename T, typename Tint>
ResultGramSchmidt_Indefinite<T> GramSchmidtOrthonormalization(MyMatrix<T> const& M, std::vector<MyVector<Tint>> const& B)
{
  int n = M.rows();
  MyMatrix<T> mu(n,n);
  struct PairInf {
    MyVector<T> Bistar_M;
    T Bistar_norm;
  };
  std::vector<PairInf> l_inf;
  std::vector<MyVector<T>> Bstar;
  std::vector<T> Bstar_norms;
  for (int i=0; i<n; i++) {
    MyVector<T> Bistar = UniversalVectorConversion<T,Tint>(B[i]);
    for (int j=0; j<i; j++) {
      T muij = (Bistar.dot(l_inf[j].Bistar_M)) / (l_inf[j].Bistar_norm);
      mu(i,j) = muij;
      Bistar -= muij * Bistar;
    }
    MyVector<T> Bistar_M = M * Bistar;
    T scal = Bistar_M.dot(Bistar);
    if (scal == 0) {
      return {false, {}, {}, {}, Bistar};
    }
    PairInf epair{std::move(Bistar_M), scal};
    Bstar.push_back(Bistar);
    Bstar_norms.push_back(scal);
    l_inf.push_back(epair);
  }
  return {true, Bstar, mu, Bstar_norms, {}};
}






template<typename T, typename Tint>
struct ResultIndefiniteLLL {
  bool success; // true if we obtained the reduced matrix. false if we found an isotropic vector
  std::vector<Myector<Tint>> B;
  MyVector<T> Xisotrop;
};



// Adapted from Denis Simon, Solving Quadratic Equations Using Reduced
// Unimodular Quadratic Forms, Math. Comp. 74(251) 1531--1543
template<typename T, typename Tint>
ResultIndefiniteLLL<T,Tint> Indefinite_LLL(MyMatrix<T> const& M)
{
  int n = M.rows();
  T c = T(7) / T(8); // The c constant of the LLL algorithm
  std::vector<MyVector<Tint>> B;
  for (int i=0; i<n; i++) {
    MyVector<Tint> V = ZeroVector<Tint>(n);
    V(i) = 1;
    B.emplace_back(std::move(V));
  }
  int k = 1;
  while(true) {
    ResultGramSchmidt_Indefinite<T> ResGS = GramSchmidtOrthonormalization(M, B);
    if (!ResGS.success) {
      return {false, {}, ResGS.Xisotrop};
    }
    for (int i=n-1; i >=0; i--) {
      for (int j=0; j<i; j++) {
        Tint q = UniversalNearestScalarInteger<Tint,T>(ResGS.mu(i,j));
        B[i] -= q * B[j];
      }
    }
    T mu = ResGS.mu(k, k-1);
    T sum1 = T_abs(ResGS.Bstar_norms[k] + mu * mu * ResGS.Bstar_norms[k-1]);
    T sum2 = c * T_abs(ResGS.Bstar_norms[k-1]);
    if (sum1 < sum2) {
      std::swap(B[k], B[k-1]);
      k = std::max(k-1, 1);
    } else {
      k++;
    }
    if (k >= n)
      break;
  }
  return {true, B, {}};
}





#endif
