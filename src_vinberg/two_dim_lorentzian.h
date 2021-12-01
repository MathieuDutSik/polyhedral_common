#ifndef INCLUDE_TWO_DIM_LORENTZIAN_H
#define INCLUDE_TWO_DIM_LORENTZIAN_H

/*
  We implement here algorithms from Section 8 in
  "AN ALTERNATIVE TO VINBERG’S ALGORITHM"

  We work with a form of signature (n,1) (also in Allcock paper).
  In this article:
  ---"timelike" vectors are vector of negative norm
  (which are contained in the cone).
  ---"lightlike" vectors have zero norm
  ---"spacelike" have positive norm

  R k is the boundary of the half plane.
  k has N(k) <= 0

  S is a set of vectors of posiive norm.
  What means "l must not lie in the interior of the
  convex hull of the norm M hyperbola in S"

 */

template<typename T, typename Tint>
T eval(const MyMatrix<T>& G, const MyVector<Tint>& v1, const MyVector<Tint>& v2)
{
  return v1(0) * v2(0) * G(0,0) + G(0,1) * (v1(0) * v2(1) + v1(1) * v2(0)) + v1(1) * v2(1) * G(1,1);
}



template<typename T, typename Tint>
std::pair<MyVector<Tint>, MyVector<Tint>> Promised(const MyMatrix<T>& G, const T&M, const MyVector<Tint>& r, const MyVector<Tint>& l)
{
  MyVector<Tint> m = r + l;
  T norm_mm = eval(G, m, m);
  T scal_ml = eval(G, m, l);
  if (norm_mm <= M || scal_ml < 0)
    return Promised(G, M, r, m);
  if (eval(G, l, l) >= 0 && eval(G, r, l) > 0)
    return {l, -r};
  return Promised(G, M, m, r);
}


template<typename T, typename Tint>
Tint LowerSquareRoot(const T& val)
{
  auto is_lower=[&](const Tint& x) -> bool {
    T x_T = x;
    return val - x_T * x_T >= 0;
  };
  double val_d = UniversalScalarConversion<double,T>(val);
  long val_sqrt = long( floor(sqrt(val_d)) );
  Tint x_guess = UniversalScalarConversion<Tint,long>(val_sqrt);
  while(true) {
    bool test1 = is_lower(x_guess);
    bool test2 = is_lower(x_guess+1);
    if (test1 && !test2) {
      return x_guess;
    } else {
      if (test2) {
        x_guess++;
      } else {
        if (!test1)
          x_guess--;
      }
    }
  }
}


template<typename T, typename Tint>
MyVector<Tint> Canonical(const MyMatrix<T>& G, const T&M, const MyVector<Tint>& r, const MyVector<Tint>& l)
{
  T scal1 = eval(G, r, l);
  T scal2 = eval(G, r, r);
  T scal3 = eval(G, l, l);
  T e_ent = scal1 * scal1 - scal2 * (scal3 - M);
  if (e_ent < 0) {
    std::cerr << "We cannot compute square root of a negative number\n";
    throw TerminalException{1};
  }
  Tint low_sqrt_tint =  LowerSquareRoot<T,Tint>(e_ent);
  T quot = (-scal1 + low_sqrt_tint) / scal2;
  Tint K = UniversalFloorScalarInteger<Tint,T>(quot);
  MyVector<Tint> l_new = l + K * r;
  return l_new;
}


template<typename T, typename Tint>
std::optional<std::pair<MyVector<Tint>, MyVector<Tint>>> Shorter(const MyMatrix<T>& G, MyVector<Tint> r, MyVector<Tint> l)
{
  auto get_char_mat=[&](const MyVector<Tint>& v1, const MyVector<Tint>& v2) -> std::vector<T> {
    return {eval(G, v1, v1), eval(G, v1, v2), eval(G, v2, v2)};
  };
  T M = eval(G, r, r);
  std::vector<T> char_ent1 = get_char_mat(r, l);
  T scal_rr;
  while(true) {
    // Step 2
    scal_rr = eval(G, r, r);
    std::pair<MyVector<Tint>, MyVector<Tint>> e_pair = Promised(G, scal_rr, r, l);
    r = e_pair.first;
    l = e_pair.second;
    // Step 3
    scal_rr = eval(G, r, r);
    l = Canonical(G, scal_rr, r, l);
    // Step 4
    if (scal_rr < M) {
      std::pair<MyVector<Tint>, MyVector<Tint>> pair{r, l};
      return pair;
    }
    // Step 5
    std::vector<T> char_ent2 = get_char_mat(r, l);
    if (char_ent1 == char_ent2)
      return {};
  }
}



template<typename T, typename Tint>
std::optional<std::pair<MyVector<Tint>, MyVector<Tint>>> NotPromised(const MyMatrix<T>& G, const T& M, MyVector<Tint> r, MyVector<Tint> l)
{
  T scal_rr = eval(G, r, r);
  if (scal_rr <= M) {
    std::pair<MyVector<Tint>, MyVector<Tint>> e_pair = Promised(G, M, r, l);
    r = e_pair.first;
    l = e_pair.second;
    l = Canonical(G, M, r, l);
    std::pair<MyVector<Tint>, MyVector<Tint>> pair{r, l};
    return pair;
  }
  while(true) {
    std::optional<std::pair<MyVector<Tint>, MyVector<Tint>>> opt_pair = Shorter(G, r, l);
    if (!opt_pair)
      return {};
    r = opt_pair->first;
    l = opt_pair->second;
    if (eval(G, r, r) <= M) {
      l = Canonical(G, M, r, l);
      std::pair<MyVector<Tint>, MyVector<Tint>> pair{r, l};
      return pair;
    }
  }
}


template<typename T, typename Tint>
std::optional<std::pair<MyMatrix<Tint>,std::vector<MyVector<Tint>>>> Anisotropic(const MyMatrix<T>& G, const T& M, MyVector<Tint> r0, MyVector<Tint> l0)
{
  std::optional<std::pair<MyVector<Tint>, MyVector<Tint>>> pair_opt = NotPromised(G, M, r0, l0);
  if (!pair_opt)
    return {};
  auto get_char_mat=[&](const MyVector<Tint>& v1, const MyVector<Tint>& v2) -> std::vector<T> {
    return {eval(G, v1, v1), eval(G, v1, v2), eval(G, v2, v2)};
  };
  MyVector<Tint> r1 = pair_opt->first;
  MyVector<Tint> l1 = pair_opt->second;
  std::vector<MyVector<Tint>> list_r{r1};
  MyVector<Tint> r = r1;
  MyVector<Tint> l = l1;
  std::vector<T> A_vect = get_char_mat(r, l);
  while(true) {
    std::pair<MyVector<Tint>,MyVector<Tint>> pair = Promised(G, M, r, l);
    r = pair.first;
    l = pair.second;
    l = Canonical(G, M, r, l);
    if (A_vect == get_char_mat(r,l)) { // Concluding step
      MyMatrix<Tint> M1_Tint = MatrixFromVectorFamily<Tint>({r1, l1});
      MyMatrix<T> M1_T = UniversalMatrixConversion<T,Tint>(M1_Tint);
      MyMatrix<Tint> M_Tint = MatrixFromVectorFamily<Tint>({r, l});
      MyMatrix<T> M_T = UniversalMatrixConversion<T,Tint>(M_Tint);
      MyMatrix<T> gMat_T = Inverse(M1_T) * M_T;
      MyMatrix<Tint> gMat_Tint = UniversalMatrixConversion<Tint,T>(gMat_T);
      std::pair<MyMatrix<Tint>,std::vector<MyVector<Tint>>> pair{gMat_Tint, list_r};
      return pair;
    }
    list_r.push_back(r);
  }

}



#endif
