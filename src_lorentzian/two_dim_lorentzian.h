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

#define DEBUG_TWO_DIM_LORENTZIAN
//#undef DEBUG_TWO_DIM_LORENTZIAN


template<typename T, typename Tint>
T eval_quad(const MyMatrix<T>& G, const MyVector<Tint>& v)
{
  return v(0) * v(0) * G(0,0) + 2 * G(0,1) * v(0) * v(1) + v(1) * v(1) * G(1,1);
}

template<typename T, typename Tint>
T eval_sval(const MyMatrix<T>& G, const MyVector<Tint>& v1, const MyVector<Tint>& v2)
{
  return v1(0) * v2(0) * G(0,0) + G(0,1) * (v1(0) * v2(1) + v1(1) * v2(0)) + v1(1) * v2(1) * G(1,1);
}


template<typename Tint>
Tint det_two(const MyVector<Tint>& r, const MyVector<Tint>& l)
{
  return r(0) * l(1) - r(1) * l(0);
}


template<typename T, typename Tint>
std::pair<MyVector<Tint>, MyVector<Tint>> Promised(const MyMatrix<T>& G, const T&M, const MyVector<Tint>& r, const MyVector<Tint>& l)
{
  MyVector<Tint> m = r + l;
#ifdef DEBUG_TWO_DIM_LORENTZIAN
  T norm_rr = eval_quad(G, r);
  if (norm_rr <= 0) {
    std::cerr << "norm_rr=" << norm_rr << "\n";
    std::cerr << "The vector r must be of positive norm\n";
    throw TerminalException{1};
  }
  Tint det = det_two(r,l);
  if (det != 1) {
    std::cerr << "det=" << det << "\n";
    std::cerr << "The configuration of vectors should satisfy det(r,l) = 1\n";
    throw TerminalException{1};
  }
  T norm_ll = eval_quad(G, l);
  if (norm_ll > M) {
    std::cerr << "norm_ll=" << norm_ll << " M=" << M << "\n";
    std::cerr << "l must not lie in the interior of the norm M hyperbola in S, that is its norm must not go above M\n";
    throw TerminalException{1};
  }
#endif
  std::cerr << "Promised : r=" << r << " l=" << l << " |l|=" << eval_quad(G, l) << " det=" << det_two(r,l) << "\n";
  T norm_mm = eval_quad(G, m);
  T scal_ml = eval_scal(G, m, l);
  if (norm_mm <= M || scal_ml < 0) {
    std::cerr << "Promised : Branching at Go Right norm_mm=" << norm_mm << " scal_ml=" << scal_ml << "\n";
    return Promised(G, M, r, m);
  }
  if (eval_quad(G, l) >= 0 && eval_scal(G, r, l) > 0) {
    std::cerr << "Promised : Exiting at the Done\n";
    return {l, -r};
  }
  std::cerr << "Promised : Exiting at the Go Left\n";
  return Promised(G, M, m, l);
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
  T scal1 = eval_scal(G, r, l);
  T scal2 = eval_quad(G, r);
  T scal3 = eval_quad(G, l);
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
    return {eval_quad(G, v1), eval_scal(G, v1, v2), eval_quad(G, v2)};
  };
  T M = eval_quad(G, r);
  std::vector<T> char_ent1 = get_char_mat(r, l);
  T scal_rr;
  while(true) {
    // Step 2
    scal_rr = eval_quad(G, r);
    std::pair<MyVector<Tint>, MyVector<Tint>> e_pair = Promised(G, scal_rr, r, l);
    r = e_pair.first;
    l = e_pair.second;
    // Step 3
    scal_rr = eval_quad(G, r);
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
  T scal_rr = eval_quad(G, r);
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
    if (eval_quad(G, r) <= M) {
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
    return {eval_quad(G, v1), eval_scal(G, v1, v2), eval_quad(G, v2)};
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
    std::cerr << "Finding |r|=" << eval_quad(G, r) << "\n";
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



template<typename T, typename Tint>
MyVector<Tint> GetPositiveVector(const MyMatrix<T>& G)
{
  MyVector<Tint> r(2);
  if (G(0,0) > 0) {
    r(0) = 1;
    r(1) = 0;
    return r;
  }
  if (G(1,1) > 0) {
    r(0) = 0;
    r(1) = 1;
    return r;
  }
  double a_d = UniversalScalarConversion<double,T>(G(0,0));
  double b_d = UniversalScalarConversion<double,T>(G(0,1));
  double c_d = UniversalScalarConversion<double,T>(G(1,1));
  // Equation for eigenvalue is
  // (a - lambda) (c - lambda) = b*b
  // l^2 - (a + c) l + [a c - b b] = 0
  // discriminant is D = (a+c)^2 - 4 (ac - b^2) = (a - c)^2 + 4 b^2
  // positive root is l_1 = [ (a + c) + sqrt(D) ] / 2
  // eigenvector equation is (a - l1) x + b y = 0
  // Thus (x,y) = (-b , a - l1)
  double disc = sqrt( (a_d - c_d) * (a_d - c_d) + 4 * b_d * b_d);
  double l1 = ( (a_d + c_d) + sqrt(disc)) / 2;
  double x = -b_d;
  double y = a_d - l1;
  size_t mult = 1;
  while(true) {
    Tint x_i = UniversalNearestScalarInteger<Tint,double>(mult * x);
    Tint y_i = UniversalNearestScalarInteger<Tint,double>(mult * y);
    T norm = G(0,0) * x_i * y_i + 2 * G(0,1) * x_i * y_i + G(1,1) * y_i * y_i;
    if (norm > 0) {
      Tint eGcd = GcdPair(x_i, y_i);
      r(0) = x_i / eGcd;
      r(1) = y_i / eGcd;
      return r;
    }
    mult++;
  }
}


template<typename Tint>
MyVector<Tint> GetTwoComplement(const MyVector<Tint>& r)
{
  GCD_int<Tint> RecGCD = ComputePairGcd(r(0), r(1));
  MyVector<Tint> l(2);
  l(0) = -RecGCD.Pmat(1,0);
  l(1) =  RecGCD.Pmat(0,0);
  return l;
}




template<typename T, typename Tint>
std::optional<std::pair<MyMatrix<Tint>,std::vector<MyVector<Tint>>>> AnisotropicComplete(const MyMatrix<T>& G, const T& M)
{
  MyVector<Tint> r = GetPositiveVector<T,Tint>(G);
  MyVector<Tint> l_A = GetTwoComplement(r);
  MyVector<Tint> l_B = Canonical(G, M, r, l_A);
  std::cerr << "r=" << r << "  l_A=" << l_A << " l_B=" << l_B << "\n";
  return Anisotropic(G, M, r, l_B);
}




#endif
