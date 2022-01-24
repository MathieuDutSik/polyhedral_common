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

//#define CHECK_TWO_DIM_LORENTZIAN
//#define DEBUG_TWO_DIM_LORENTZIAN


template<typename T, typename Tint>
T eval_quad(const MyMatrix<T>& G, const MyVector<Tint>& v)
{
  return v(0) * v(0) * G(0,0) + 2 * G(0,1) * v(0) * v(1) + v(1) * v(1) * G(1,1);
}

template<typename T, typename Tint>
T eval_scal(const MyMatrix<T>& G, const MyVector<Tint>& v1, const MyVector<Tint>& v2)
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
  std::cerr << "Promised G=" << G(0,0) << "," << G(1,0) << "," << G(1,1) << " M=" << M << "\n";
#endif
#ifdef CHECK_TWO_DIM_LORENTZIAN
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
#ifdef DEBUG_TWO_DIM_LORENTZIAN
  std::cerr << "Promised : r=" << r << " l=" << l << " |l|=" << eval_quad(G, l) << " det=" << det_two(r,l) << "\n";
#endif
  T norm_mm = eval_quad(G, m);
  T scal_mr = eval_scal(G, m, r);
  if (norm_mm <= M || scal_mr < 0) {
#ifdef DEBUG_TWO_DIM_LORENTZIAN
    std::cerr << "Promised : Branching at Go Right norm_mm=" << norm_mm << " scal_mr=" << scal_mr << "\n";
#endif
    return Promised(G, M, r, m);
  }
  if (eval_quad(G, l) >= 0 && eval_scal(G, r, l) > 0) {
#ifdef DEBUG_TWO_DIM_LORENTZIAN
    std::cerr << "Promised : Exiting at the Done\n";
#endif
    return {l, -r};
  }
#ifdef DEBUG_TWO_DIM_LORENTZIAN
  std::cerr << "Promised : Exiting at the Go Left\n";
#endif
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
  //  l = Canonical(G, M, r, l);
  std::vector<T> char_ent1 = get_char_mat(r, l);
#ifdef DEBUG_TWO_DIM_LORENTZIAN
  std::cerr << "Shorter M=" << M << "\n";
#endif
  T scal_rr;
  while(true) {
    // Step 2
    scal_rr = eval_quad(G, r);
#ifdef DEBUG_TWO_DIM_LORENTZIAN
    std::cerr << "Shorter : Step 2, scal_rr=" << scal_rr << "\n";
#endif
    std::pair<MyVector<Tint>, MyVector<Tint>> e_pair = Promised(G, scal_rr, r, l);
    r = e_pair.first;
    l = e_pair.second;
    // Step 3
#ifdef DEBUG_TWO_DIM_LORENTZIAN
    std::cerr << "Shorter : Step 3\n";
#endif
    scal_rr = eval_quad(G, r);
    l = Canonical(G, scal_rr, r, l);
    // Step 4
#ifdef DEBUG_TWO_DIM_LORENTZIAN
    std::cerr << "Shorter : Step 4\n";
#endif
    if (scal_rr < M) {
      std::pair<MyVector<Tint>, MyVector<Tint>> pair{r, l};
      return pair;
    }
    // Step 5
#ifdef DEBUG_TWO_DIM_LORENTZIAN
    std::cerr << "Shorter : Step 5\n";
#endif
    std::vector<T> char_ent2 = get_char_mat(r, l);
#ifdef DEBUG_TWO_DIM_LORENTZIAN
    std::cerr << "Shorter : char_ent2=" << char_ent2 << " char_ent1=" << char_ent1 << "\n";
#endif
    if (char_ent1 == char_ent2)
      return {};
  }
}



template<typename T, typename Tint>
std::optional<std::pair<MyVector<Tint>, MyVector<Tint>>> NotPromised(const MyMatrix<T>& G, const T& M, MyVector<Tint> r, MyVector<Tint> l)
{
  T scal_rr = eval_quad(G, r);
#ifdef DEBUG_TWO_DIM_LORENTZIAN
  std::cerr << "NotPromised : scal_rr=" << scal_rr << " M=" << M << "\n";
#endif
  if (scal_rr <= M) {
    std::pair<MyVector<Tint>, MyVector<Tint>> e_pair = Promised(G, M, r, l);
    r = e_pair.first;
    l = e_pair.second;
    l = Canonical(G, M, r, l);
    std::pair<MyVector<Tint>, MyVector<Tint>> pair{r, l};
    return pair;
  }
#ifdef DEBUG_TWO_DIM_LORENTZIAN
  std::cerr << "NotPromised : Before infinite loop\n";
#endif
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
#ifdef DEBUG_TWO_DIM_LORENTZIAN
  std::cerr << "We pass the NotPromised stage\n";
#endif
  auto get_char_mat=[&](const MyVector<Tint>& v1, const MyVector<Tint>& v2) -> std::vector<T> {
    return {eval_quad(G, v1), eval_scal(G, v1, v2), eval_quad(G, v2)};
  };
  MyVector<Tint> r1 = pair_opt->first;
  MyVector<Tint> l1 = pair_opt->second;
#ifdef DEBUG_TWO_DIM_LORENTZIAN
  std::cerr << "Anisotropic r1=" << StringVectorGAP(r1) << " / " << StringVectorGAP(l1) << "\n";
#endif
  std::vector<MyVector<Tint>> list_r{r1};
  MyVector<Tint> r = r1;
  MyVector<Tint> l = l1;
#ifdef DEBUG_TWO_DIM_LORENTZIAN
  std::cerr << "First : Anisotropic r=" << StringVectorGAP(r) << " / " << StringVectorGAP(l) << "\n";
#endif
  std::vector<T> A_vect = get_char_mat(r, l);
  while(true) {
    std::pair<MyVector<Tint>,MyVector<Tint>> pair = Promised(G, M, r, l);
    r = pair.first;
#ifdef DEBUG_TWO_DIM_LORENTZIAN
    std::cerr << "Finding |r|=" << eval_quad(G, r) << "\n";
#endif
    l = pair.second;
    l = Canonical(G, M, r, l);
#ifdef DEBUG_TWO_DIM_LORENTZIAN
    std::cerr << "Now : Anisotropic r=" << StringVectorGAP(r) << " / " << StringVectorGAP(l) << "\n";
    std::cerr << "After canonical\n";
#endif
    if (A_vect == get_char_mat(r,l)) { // Concluding step
#ifdef DEBUG_TWO_DIM_LORENTZIAN
      std::cerr << "Exiting case\n";
#endif
      MyMatrix<Tint> M1_Tint = MatrixFromVectorFamily<Tint>({r1, l1});
      MyMatrix<T> M1_T = UniversalMatrixConversion<T,Tint>(M1_Tint);
      MyMatrix<Tint> M_Tint = MatrixFromVectorFamily<Tint>({r, l});
      MyMatrix<T> M_T = UniversalMatrixConversion<T,Tint>(M_Tint);
      MyMatrix<T> gMat_T = Inverse(M1_T) * M_T;
      MyMatrix<Tint> gMat_Tint = UniversalMatrixConversion<Tint,T>(gMat_T);
      std::pair<MyMatrix<Tint>,std::vector<MyVector<Tint>>> pair{gMat_Tint, list_r};
#ifdef DEBUG_TWO_DIM_LORENTZIAN
      std::cerr << "Before returning pair\n";
#endif
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




/*
  The matrix is written as
  | a b |
  | b c |
  The associated quadratic form is
  Q(x,y) = a x^2 + 2 b xy + c y^2
  The associated polynomial is
  P(x)   = a x^2 + 2 b x + c
  The disciminant is Delta =  4 b^2 - 4 a c
  If that disciminant is a square then no isotropic decomposition is possible
  ---
  The roots of the polynomial are
  x1 = (-2b + sqrt(4 Delta)) / 2a = (-b + sqrt(Delta)) / a
  x2 = (-2b = sqrt(4 Delta)) / 2a = (-b - sqrt(Delta)) / a
 */
template<typename T>
std::optional<MyMatrix<T>> GetIsotropicFactorization(MyMatrix<T> const& G)
{
  if (G.rows() != 2 || G.cols() != 2) {
    std::cerr << "The matrix is not square of order 2\n";
    throw TerminalException{1};
  }
  T a = G(0,0);
  T b = G(1,0);
  T c = G(1,1);
  T Delta = b * b - a * c;
  std::optional<T> opt_root = UniversalSquareRoot(Delta);
  if (!opt_root)
    return {};
  T Delta_root = *opt_root;
  MyMatrix<T> F(2,2);
  if (a != 0) {
    T x1 = (-b + Delta_root) / a;
    T x2 = (-b - Delta_root) / a;
    // So we have a factorization as a (x - x1) (x - x2)
    // and so of the quadratic form as
    // (a x - a x1 y) (x - x2 y)
    F(0,0) = a;
    F(0,1) = -a * x1;
    F(1,0) = 1;
    F(1,1) = -x2;
    return F;
  }
  if (c != 0) {
    // Polynomial becomes c y^2 + 2 b y + a
    // roots are
    // y1 = (-b + sqrt(Delta)) / c
    // y2 = (-b - sqrt(Delta)) / c
    T y1 = (-b + Delta_root) / c;
    T y2 = (-b - Delta_root) / c;
    // So now we have a factorization of the form
    // c (y - y1) ( y - y2)
    // and so the quadratic form as
    // Q(x,y) = (c y - cy1 x) (y - y2 x)
    F(0,0) = -c * y1;
    F(0,1) = c;
    F(1,0) = -y2;
    F(1,1) = 1;
    return F;
  }
  // We are now in the case a = 0, c=0
  // Factorization is obvious: (2b x) ( y )
  F(0,0) = 2 * b;
  F(0,1) = 0;
  F(1,0) = 0;
  F(1,1) = 1;
  return F;
}



/*
  F is the factorization with each row representing one term of the factorization.
  Q(x,y) = (a1 x + b1 y) (a2 x + b2 y)

 */
template<typename T, typename Tint>
std::vector<MyVector<Tint>> EnumerateVectorFixedNorm_Factorization(MyMatrix<T> const& F, T const& M)
{
#ifdef DEBUG_TWO_DIM_LORENTZIAN
  std::cerr << "EnumerateVectorFixedNorm_Factorization M=" << M << " F=\n";
  WriteMatrix(std::cerr, F);
#endif
  MyVector<T> v1(2);
  v1(0) = F(0,0);
  v1(1) = F(0,1);
  FractionVector<T> Fr1 = RemoveFractionVectorPlusCoeff(v1);
  //
  MyVector<T> v2(2);
  v2(0) = F(1,0);
  v2(1) = F(1,1);
  FractionVector<T> Fr2 = RemoveFractionVectorPlusCoeff(v2);
  //
  T M_scal = M * Fr1.TheMult * Fr2.TheMult;
#ifdef DEBUG_TWO_DIM_LORENTZIAN
  std::cerr << "M_scal=" << M_scal << "\n";
#endif
  if (!IsInteger(M_scal))
    return {};
  MyMatrix<T> A(2,2);
  A(0,0) = Fr1.TheVect(0);
  A(0,1) = Fr1.TheVect(1);
  A(1,0) = Fr2.TheVect(0);
  A(1,1) = Fr2.TheVect(1);
#ifdef DEBUG_TWO_DIM_LORENTZIAN
  std::cerr << "EnumerateVectorFixedNorm_Factorization A=\n";
  WriteMatrix(std::cerr, A);
#endif
  MyMatrix<T> Ainv = Inverse(A);
#ifdef DEBUG_TWO_DIM_LORENTZIAN
  std::cerr << "EnumerateVectorFixedNorm_Factorization Ainv=\n";
  WriteMatrix(std::cerr, Ainv);
#endif
  MyVector<T> v(2);
  std::vector<MyVector<Tint>> l_sol;
  std::vector<T> list_div = GetAllFactors(T_abs(M_scal));
#ifdef DEBUG_TWO_DIM_LORENTZIAN
  std::cerr << "|list_div|=" << list_div.size() << "\n";
#endif
  for (auto & e_div : list_div) {
    v(0) = e_div;
    v(1) = M_scal / e_div;
    MyVector<T> e_sol = Ainv * v;
#ifdef DEBUG_TWO_DIM_LORENTZIAN
    std::cerr << "v=" << StringVectorGAP(v) << "  e_sol=" << StringVectorGAP(e_sol) << "\n";
#endif
    if (IsIntegralVector(e_sol)) {
      MyVector<Tint> e_sol_i = UniversalVectorConversion<Tint,T>(e_sol);
#ifdef DEBUG_TWO_DIM_LORENTZIAN
      std::cerr << "  e_sol_i=" << StringVectorGAP(e_sol_i) << "\n";
#endif
      l_sol.push_back( e_sol_i);
      l_sol.push_back(-e_sol_i);
    }
  }
#ifdef DEBUG_TWO_DIM_LORENTZIAN
  std::cerr << "|l_sol|=" << l_sol.size() << "\n";
#endif
  return l_sol;
}







/*
  We are happy for the purpose of this code to have only an upper bound of the
  minimum
  Finding the exact minimum appears relatively complex. See 
 */
template<typename T, typename Tint>
T GetUpperBoundMinimum(const MyMatrix<T>& G)
{
  MyVector<Tint> eV = GetPositiveVector(G);
  T eNorm = eval_quad(G, eV);
  auto get_spann=[&](MyVector<Tint> const& v) -> std::vector<MyVector<Tint>> {
    MyVector<Tint> v1{v(0), v(1)+1};
    MyVector<Tint> v2{v(0), v(1)-1};
    MyVector<Tint> v3{v(0)+1, v(1)};
    MyVector<Tint> v4{v(0)-1, v(1)};
    return {v1,v2,v3,v4};
  };
  auto try_better=[&]() -> bool {
    for (auto & newV : get_spann(eV)) {
      T newNorm = eval_quad(G, newV);
      if (newNorm < eNorm) {
        eV = newV;
        eNorm = newNorm;
        return true;
      }
    }
    return false;
  };
  while(try_better()) {}
  return eNorm;
}

template<typename T>
T GetLowerBoundMinimum(const MyMatrix<T>& G)
{
  T val1 = G(0,0);
  T val2 = 2 * G(0,1);
  T val3 = G(1,1);
  T eGcd1 = GcdPair(val1, val2);
  return GcdPair(eGcd1, val3);
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
#ifdef DEBUG_TWO_DIM_LORENTZIAN
  std::cerr << "r=" << r << "  l_A=" << l_A << " l_B=" << l_B << "\n";
#endif
  return Anisotropic(G, M, r, l_B);
}


template<typename T, typename Tint>
std::optional<MyVector<Tint>> get_first_next_vector_isotropic(MyMatrix<T> const& G, MyVector<Tint> const& r0, T const& SearchNorm, MyMatrix<T> const& F)
{
  std::vector<MyVector<Tint>> l_sol = EnumerateVectorFixedNorm_Factorization<T,Tint>(F, SearchNorm);
  auto is_corr=[&](MyVector<Tint> const& x) -> bool {
    T norm = eval_quad(G, x);
    if (norm != SearchNorm)
      return false;
    T scal = eval_scal(G, r0, x);
    if (scal <= 0)
      return false;
    return det_two(r0, x) > 0;
  };
  std::vector<MyVector<Tint>> l_sol_red;
  std::optional<MyVector<Tint>> e_sol;
  size_t n_match=0;
  for (auto & e_v : l_sol) {
    if (is_corr(e_v)) {
      n_match++;
      if (e_sol) {
        if (det_two(e_v, *e_sol) > 0)
          e_sol = e_v;
      }
    } else {
      e_sol = e_v;
    }
  }
#ifdef DEBUG_TWO_DIM_LORENTZIAN
  std::cerr << "n_match=" << n_match << "\n";
#endif
  return e_sol;
}

template<typename T, typename Tint>
std::optional<MyVector<Tint>> get_first_next_vector_anisotropic(MyMatrix<T> const& G, MyVector<Tint> const& r0, T const& SearchNorm)
{
  T M = eval_quad(G, r0);
#ifdef DEBUG_TWO_DIM_LORENTZIAN
  std::cerr << "M=" << M << " SearchNorm=" << SearchNorm << " r0=" << StringVectorGAP(r0) << " G=\n";
  WriteMatrix(std::cerr, G);
#endif
  MyVector<Tint> l_A = GetTwoComplement(r0);
  MyVector<Tint> l_B = Canonical(G, M, r0, l_A);
  std::optional<std::pair<MyMatrix<Tint>,std::vector<MyVector<Tint>>>> opt = Anisotropic(G, M, r0, l_B);
#ifdef DEBUG_TWO_DIM_LORENTZIAN
  std::cerr << "We have opt\n";
#endif
  if (!opt)
    return {};
  std::vector<MyVector<Tint>> const& l_vect = opt->second;
#ifdef DEBUG_TWO_DIM_LORENTZIAN
  std::cerr << "|l_vect|=" << l_vect.size() << "\n";
#endif
  for (auto & e_v : l_vect) {
    T norm = eval_quad(G, e_v);
#ifdef DEBUG_TWO_DIM_LORENTZIAN
    std::cerr << "norm=" << norm << " SearchNorm=" << SearchNorm << " e_v=" << StringVectorGAP(e_v) << "\n";
#endif
    if (norm == SearchNorm)
      return e_v;
  }
#ifdef DEBUG_TWO_DIM_LORENTZIAN
  std::cerr << "returning missing\n";
#endif
  return {};
}





template<typename T, typename Tint>
std::optional<MyVector<Tint>> get_first_next_vector(MyMatrix<T> const& G, MyVector<Tint> const& r0, T const& SearchNorm)
{
#ifdef DEBUG_TWO_DIM_LORENTZIAN
  std::cerr << "SearchNorm=" << SearchNorm << "\n";
#endif
  if (SearchNorm <= 0)
    return {};
  T lower_bnd = GetLowerBoundMinimum(G);
#ifdef DEBUG_TWO_DIM_LORENTZIAN
  std::cerr << "lower_bnd=" << lower_bnd << "\n";
#endif
  if (SearchNorm < lower_bnd) // no solution possible
    return {};
  std::optional<MyMatrix<T>> opt = GetIsotropicFactorization(G);
  if (opt) {
#ifdef DEBUG_TWO_DIM_LORENTZIAN
    std::cerr << "get_first_next_vector : isotropic\n";
#endif
    return get_first_next_vector_isotropic(G, r0, SearchNorm, *opt);
  }
#ifdef DEBUG_TWO_DIM_LORENTZIAN
  std::cerr << "get_first_next_vector : anisotropic\n";
#endif
  return get_first_next_vector_anisotropic(G, r0, SearchNorm);
}







#endif
