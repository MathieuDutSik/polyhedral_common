#ifndef TEMP_SPARSE_MATRIX_SOLVER_GAMPMATLAB
#define TEMP_SPARSE_MATRIX_SOLVER_GAMPMATLAB

#include "MAT_Matrix.h"
#include <vector>
#include <string>
#include <algorithm>

//#define DEBUG
#undef DEBUG

template <typename T> T L1_Norm(MyVector<T> const &b) {
  T eNorm = 0;
  int siz = b.size();
  for (int i = 0; i < siz; i++) {
    T eAbs = std::abs(b(i));
    eNorm += eAbs;
  }
  return eNorm;
}

template <typename T> T Linf_Norm(MyVector<T> const &b) {
  T eNorm = 0;
  int siz = b.size();
  for (int i = 0; i < siz; i++) {
    T eAbs = std::abs(b(i));
    if (eAbs > eNorm) {
      eNorm = eAbs;
    }
  }
  return eNorm;
}

template <typename T> T L2_Norm(MyVector<T> const &b) {
  T eNorm = 0;
  int siz = b.size();
  for (int i = 0; i < siz; i++) {
    T eVal = b(i);
    eNorm += eVal * eVal;
  }
  return sqrt(eNorm);
}

template <typename T> struct RecSparse {
  int n;
  int m;
  std::function<MyVector<T>(MyVector<T> const &)> A;
  std::function<MyVector<T>(MyVector<T> const &)> At;
  MyVector<T> weights;
};

template <typename T> struct RecOptSparse {
  T delta;
  T tol;
  T mu;
  int maxit;
  bool print;
  T nu;
  T eps;
  T rho;
  bool nonneg;
  T gamma;
  T xs;
  bool nonorth;
  int stepfreq;
  bool UseWeight;
  bool basis;
};

template <typename T>
RecSparse<T> AMP_linear_operators(MySparseMatrix<T> const &preA0) {
  const int n = preA0.rows();
  const int m = preA0.cols();
  //  std::cerr << "(rows/cols)preA0=" << preA0.rows() << " / " << preA0.cols()
  //  << "\n";
  MySparseMatrix<T> A0 = preA0.transpose();
  //  std::cerr << "(rows/cols)A0=" << A0.rows() << " / " << A0.cols() << "\n";
  std::function<MyVector<T>(MyVector<T> const &)> f1 =
      [=](MyVector<T> const &x) -> MyVector<T> {
    //    std::cerr << "Before product A0*x\n";
    //    std::cerr << "(rows/cols)x = " << x.rows() << " / " << x.cols() <<
    //    "\n"; std::cerr << "(rows/cols)A0 = " << A0.rows() << " / " <<
    //    A0.cols() << "\n";
    MyVector<T> yRet = A0 * x;
    //    std::cerr << "L1(x)=" << L1_norm(x) << " L1(yRet)=" << L1_norm(yRet)
    //    << "\n";
    return yRet;
  };
  std::function<MyVector<T>(MyVector<T> const &)> f2 =
      [=](MyVector<T> const &x) -> MyVector<T> {
    //    std::cerr << "f2: n=" << n << " m=" << m << "\n";
    //    std::cerr << "(rows/cols)x = " << x.rows() << " / " << x.cols() <<
    //    "\n"; std::cerr << "(rows/cols)A0 = " << A0.rows() << " / " <<
    //    A0.cols() << "\n"; std::cerr << "Before product x*A0\n";
    MyMatrix<T> xTr = x.transpose();
    //    std::cerr << "We have xTr\n";
    //    std::cerr << "(rows/cols)xTr = " << xTr.rows() << " / " << xTr.cols()
    //    << "\n";
    MyMatrix<T> eProduct = xTr * A0;
    //    std::cerr << "We have eProduct\n";
    MyVector<T> yRet = eProduct.transpose();
    //    std::cerr << "L1(x)=" << L1_norm(x) << " L1(yRet)=" << L1_norm(yRet)
    //    << "\n";
    return yRet;
  };
  MyVector<T> weight(n);
  size_t siz = n;
  for (size_t i = 0; i < siz; i++)
    weight(i) = 1;
  return RecSparse<T>{n, m, f1, f2, weight};
}

template <typename T> struct OutSolver {
  MyVector<T> x;
  MyVector<T> z;
  int iter;
  int cntAt;
  int cntA;
  std::string exit;
  int cputime;
  std::vector<T> error;
  std::vector<T> optim;
};

template <typename T>
OutSolver<T> AMP_yall1(RecSparse<T> const &eRecSparse, MyVector<T> const &b,
                       RecOptSparse<T> const &eRecOpt) {
  /*
A solver for L1-minimization models:

min ||Wx||_{w,1}, st Ax = b
min ||Wx||_{w,1} + (1/nu)||Ax - b||_1
min ||Wx||_{w,1} + (1/2*rho)||Ax - b||_2^2
min ||x||_{w,1}, st Ax = b                and x > = 0
min ||x||_{w,1} + (1/nu)||Ax - b||_1,      st x > = 0
min ||x||_{w,1} + (1/2*rho)||Ax - b||_2^2, st x > = 0

where (A,b,x) can be complex or real
(but x must be real in the last 3 models)

Copyright(c) 2009 Yin Zhang
Test Version: please do NOT distribute
--------------------------------------

--- Input:
    A --- either an m x n matrix or
          a structure with 2 fields:
          1) A.times: a function handle for A*x
          2) A.trans: a function handle for A'*y
    b --- an m-vector, real or complex
 opts --- a structure with fields:
          opts.tol   -- tolerance *** required ***
          opts.nu    -- values > 0 for L1/L1 model
          opts.rho   -- values > 0 for L1/L2 model
          opts.basis -- sparsifying unitary basis W (W*W = I)
                       a struct with 2 fields:
                       1) times: a function handle for W*x
                       2) trans: a function handle for W'*y
          opts.nonneg  -- 1 for nonnegativity constraints
          opts.nonorth -- 1 for A with non-orthonormal rows
          see the User's Guide for all other options

--- Output:
    x --- last iterate (hopefully an approximate solution)
  Out --- a structure with fields:
          Out.status  --- exit information
          Out.iter    --- #iterations taken
          Out.cputime --- solver CPU time
          Out.z       --- final dual slack value
  */

  int m = eRecSparse.m;
  bool L1L1 = false;
  if (eRecOpt.nu > 0) {
    L1L1 = true;
  }
  if (L1L1 && eRecOpt.UseWeight) {
    // opts.weights = [opts.weights(:); ones(m,1)];
  }

  bool posrho = eRecOpt.rho > 0;
  bool posdel = eRecOpt.delta > 0;
  bool posnu = eRecOpt.nu > 0;
  bool nonneg = eRecOpt.nonneg;
  MyVector<T> x0;
  MyVector<T> z0;
  RecOptSparse<T> NewRecOpt = eRecOpt;
  if ((posdel && posrho) || (posdel && posnu) || (posrho && posnu)) {
    std::cerr << "Model parameter conflict! YALL1: set delta = 0 && nu = 0;\n";
    NewRecOpt.delta = 0;
    posdel = false;
    NewRecOpt.nu = 0;
    posnu = false;
  }
  std::string prob = "the basis pursuit problem";
  if (posrho)
    prob = "the unconstrained L1L2 problem";
  if (posdel)
    prob = "the constrained L1L2 problem";
  if (posnu)
    prob = "the unconstrained L1L1 problem";
  MyVector<T> Atb = eRecSparse.At(b);
  int n = Atb.size();
  T bmax = Linf_Norm(b);
  T tol = eRecOpt.tol;
  bool L2Unc_zsol = posrho && Linf_Norm(Atb) <= eRecOpt.rho;
  bool L2Con_zsol = posdel && L2_Norm(b) <= eRecOpt.delta;
  bool L1L1_zsol = posnu && bmax < tol;
  bool BP_zsol = (!posrho) && (!posdel) && (!posnu) && bmax < tol;
  bool zsol = L2Unc_zsol || L2Con_zsol || BP_zsol || L1L1_zsol;
  if (zsol) {
    MyVector<T> x(n);
    for (int i = 0; i < n; i++)
      x(i) = 0;
    OutSolver<T> Out;
    Out.x = x;
    Out.iter = 0;
    Out.cntAt = 1;
    Out.cntA = 0;
    Out.exit = "Data b = 0";
    //    std::cerr << "Return AMP_yall1, case 1\n";
    return Out;
  };
  MyVector<T> b1 = b / bmax;
  if (posrho)
    NewRecOpt.rho = NewRecOpt.rho / bmax;
  if (posdel)
    NewRecOpt.delta = NewRecOpt.delta / bmax;
  std::time_t t0 = std::time(0);

  OutSolver<T> Out = yall1_solve(eRecSparse, b1, x0, z0, NewRecOpt);

  Out.cputime = std::time(0) - t0;

  MyVector<T> x = Out.x * bmax;
  if (L1L1) {
    int xSize = Out.x.size();
    MyVector<T> xShort(xSize - m);
    for (int i = 0; i < xSize - m; i++) {
      xShort(i) = Out.x(i);
    }
    Out.x = xShort;
  }
#ifdef DEBUG
  if (eRecOpt.basis) {
    // x = opts.basis.trans(x);
    std::cerr << "Missing code here\n";
  }
#endif
  if (nonneg) {
    int xSize = Out.x.size();
    for (int i = 0; i < xSize; i++) {
      if (Out.x(i) <= 0)
        Out.x(i) = 0;
    }
  }
  //  std::cerr << "Return AMP_yall1, case 2\n";
  return Out;
}

/* code for generating random double
  const double lower_bound = 0;
  const double upper_bound = 1;
  std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
  std::random_device rand_dev;          // Use random_device to get a random
  seed. std::mt19937 rand_engine(rand_dev()); // mt19937 is a good pseudo-random
  number std::cout << x << std::endl; for (int i=0; i<n; i++) { double x =
  unif(rand_engine);
*/

template <typename T> bool check_orth(RecSparse<T> const &eRecSparse) {
  int n = eRecSparse.n;
  MyVector<T> s1(n);
  for (int i = 0; i < n; i++) {
    int eVal = rand() % 100;
    T eValT = eVal;
    s1(i) = eValT;
  }
  MyVector<T> s2 = eRecSparse.A(eRecSparse.At(s1));
  T err = L2_Norm(s1 - s2) / L2_Norm(s1);
  T tol = 1 / 100;
  if (err > tol)
    return false;
  return true;
}

template <typename T>
MyVector<T> proj2box(MyVector<T> const &z, MyVector<T> const &w,
                     bool const &nonneg, T const &nu, int const &m) {
  size_t siz = z.size();
  MyVector<T> zRet = MyVector<T>(siz);
  //  std::cerr << "nonneg=" << nonneg << " siz=" << siz << "\n";
  //  std::cerr << "(rows/cols)w=" << w.rows() << " / " << w.cols() << "\n";
  //  std::cerr << "(rows/cols)z=" << z.rows() << " / " << z.cols() << "\n";
  if (nonneg) {
    for (size_t i = 0; i < siz; i++)
      zRet(i) = std::min(w(i), z(i));
    if (nu > 0) {
      for (size_t i = siz - m; i < siz; i++) {
        zRet(i) = std::max(T(-1), zRet(i));
      }
    }
  } else {
    for (size_t i = 0; i < siz; i++) {
      zRet(i) = z(i) * w(i) / std::max(w(i), std::abs(z(i)));
    }
  }
  //  std::cerr << "Returning from proj2box\n";
  return zRet;
}

template <typename T> T T_sumAbs(MyVector<T> const &a, MyVector<T> const &b) {
  T eSum;
  eSum = 0;
  size_t siz = a.size();
  for (size_t i = 0; i < siz; i++) {
    T eProd = a(i) * b(i);
    T eAbs = std::abs(eProd);
    eSum += eAbs;
  }
  return eSum;
}

template <typename T>
T T_ScalarProduct(MyVector<T> const &a, MyVector<T> const &b) {
  T eSum;
  eSum = 0;
  size_t siz = a.size();
  for (size_t i = 0; i < siz; i++) {
    T eProd = a(i) * b(i);
    eSum += eProd;
  }
  return eSum;
}

template <typename T> T MeanAbs(MyVector<T> const &V) {
  T eSum = 0;
  size_t siz = V.size();
  for (size_t i = 0; i < siz; i++)
    eSum += std::abs(V(i));
  return eSum / T(siz);
}

template <typename T>
OutSolver<T> yall1_solve(RecSparse<T> const &eRecSparse, MyVector<T> const &b,
                         MyVector<T> const &x0, MyVector<T> const &z0,
                         RecOptSparse<T> const &eRecOpt) {
  /*
    yall1_solve version beta-6 (Aug. 2, 2009)
    Copyright(c) 2009 Yin Zhang
  */
  int m = eRecSparse.m;
  int n = eRecSparse.n;
  //  std::cerr << "m=" << m << " n=" << n << "\n";
  T sqrtm = sqrt(m);
  T bnrm = L2_Norm(b);
  T tol = eRecOpt.tol;
  T mu = eRecOpt.mu;
  if (mu == 0)
    mu = MeanAbs(b);
  T gamma = eRecOpt.gamma;
  if (gamma == 0)
    gamma = 1.618;
  T eps = eRecOpt.eps;
  int maxit = eRecOpt.maxit;
#ifdef DEBUG
  bool print = eRecOpt.print;
#endif
  T nu = eRecOpt.nu;
  T rho = eRecOpt.rho;
  T delta = eRecOpt.delta;
  bool nonneg = eRecOpt.nonneg;
  bool nonorth = eRecOpt.nonorth;
  int stepfreq = eRecOpt.stepfreq;
  int iter;
  OutSolver<T> Out;
  MyVector<T> rd;
  MyVector<T> xp;
  MyVector<T> w = eRecSparse.weights;
  MyVector<T> x = x0;
  //  std::cerr << "(rows/cols)x0=" << x0.rows() << " / " << x0.cols() << "\n";
  if (x0.rows() == 0)
    x = eRecSparse.At(b);
  //  std::cerr << "(rows/cols)x=" << x.rows() << " / " << x.cols() << "\n";
  MyVector<T> z = z0;
  if (z0.rows() == 0)
    z = ZeroVector<T>(n);
  MyVector<T> y(m);
  MyVector<T> Aty(n);
  /*
  auto iprint1_0=[&]() -> void {
    MyVector<T> rp = eRecSparse.A(x) - b;
    T rpnrm = L2_Norm(rp);
    std::cerr << " norm( A*x0 - b ) = " << rpnrm << "\n";
    };*/
#ifdef DEBUG
  auto iprint1_1 = [&]() -> void {
    MyVector<T> rp = eRecSparse.A(x) - b;
    T objp = T_sumAbs(w, x);
    T objd = T_ScalarProduct(b, y);
    if (rho > 0) {
      T eNorm1 = T_ScalarProduct(rp, rp);
      objp = objp + (1 / (2 * rho)) * eNorm1;
      T eNorm2 = T_ScalarProduct(y, y);
      objd = objd - (rho / 2) * eNorm2;
    }
    T dgap = std::abs(objd - objp);
    T rel_gap = dgap / std::abs(objp);
    T rdnrm = L2_Norm(rd);
    T rel_rd = rdnrm / sqrtm;
    T rpnrm = L2_Norm(rp);
    T rel_rp = rpnrm / bnrm;
    std::cerr << " Rel_Dgap  Rel_ResD  Rel_ResP\n";
    std::cerr << " " << rel_gap << "  " << rel_rd << "  " << rel_rp << "\n";
  };
#endif
#ifdef DEBUG
  auto iprint2 = [&]() -> void {
    //    std::cerr << "Beginning of iprint2\n";
    T rdnrm = L2_Norm(rd);
    //    std::cerr << "iprint2, step 1\n";
    MyVector<T> rp = eRecSparse.A(x) - b;
    //    std::cerr << "iprint2, step 2\n";
    T rpnrm = L2_Norm(rp);
    //    std::cerr << "iprint2, step 3\n";
    T objp = T_sumAbs(w, x);
    //    std::cerr << "iprint2, step 4\n";
    T objd = T_ScalarProduct(b, y);
    //    std::cerr << "iprint2, step 5\n";
    if (rho > 0) {
      T eNorm1 = T_ScalarProduct(rp, rp);
      objp = objp + (1 / (2 * rho)) * eNorm1;
      T eNorm2 = T_ScalarProduct(y, y);
      objd = objd - (rho / 2) * eNorm2;
    }
    //    std::cerr << "iprint2, step 6\n";
    T dgap = std::abs(objd - objp);
    //    std::cerr << "iprint2, step 7 iter=" << iter << "\n";
#ifdef DEBUG
    if (iter % 50 == 0) {
      std::cerr << "  Iter " << iter << ":\n";
      std::cerr << "  Dgap " << dgap << "\n";
      std::cerr << "  ResD " << rdnrm << "\n";
      std::cerr << "  ResP " << rpnrm << "\n";
      std::cerr << "\n";
    }
#endif
    //    std::cerr << "xs=" << eRecOpt.xs << "\n";
    //    std::cerr << "nu=" << eRecOpt.nu << "\n";
    if (eRecOpt.xs > 0 && eRecOpt.nu == 0) {
      T optim = std::max(dgap / std::abs(objp), rdnrm / sqrtm);
      if (rho == 0) {
        optim = std::max(optim, rpnrm / bnrm);
      };
      Out.optim.push_back(optim);
      size_t siz = x.size();
      MyVector<T> xShift(siz);
      for (size_t i = 0; i < siz; i++)
        xShift(i) = x(i) - eRecOpt.xs;
      Out.error.push_back(L2_Norm(xShift));
    }
  };
#endif
  auto check_stopping = [&]() -> bool {
    bool stop = false;
    T q = 1 / 10; // q in [0,1)
    if (delta > 0)
      q = 0;
    // check relative change
    //    std::cerr << "(rows/cols) x=" << x.rows() << " / " << x.cols() <<
    //    "\n"; std::cerr << "(rows/cols)xp=" << xp.rows() << " / " << xp.cols()
    //    << "\n";
    MyVector<T> TheDiff = x - xp;
    T xrel_chg = L2_Norm(TheDiff) / L2_Norm(x);
    if (xrel_chg < tol * (1 - q)) {
      Out.exit = "Exit: Stablized";
      stop = true;
      //      std::cerr << "L2_Norm(TheDiff)=" << L2_Norm(TheDiff) << "\n";
      //      std::cerr << "L2_Norm(x)=" << L2_Norm(x) << "\n";
      //      std::cerr << "xrel_chg=" << xrel_chg << "\n";
      //      std::cerr << "tol=" << tol << " q=" << q << "\n";
      //      std::cerr << "check_stopping, exit, case 1\n";
      return stop;
    }
    if (xrel_chg >= tol * (1 + q)) {
      //      std::cerr << "check_stopping, exit, case 2\n";
      return stop;
    }
    // check dual residual
    T rdnrm = L2_Norm(rd);
    bool d_feasible = rdnrm < tol * sqrtm;
    if (!d_feasible) {
      //      std::cerr << "check_stopping, exit, case 3\n";
      return stop;
    }
    // check duality gap
    T objp = T_sumAbs(w, x);
    T objd = T_ScalarProduct(b, y);
    MyVector<T> rp;
    if (rho > 0) {
      rp = eRecSparse.A(x) - b;
      Out.cntA = Out.cntA + 1;
      T eNorm1 = T_ScalarProduct(rp, rp);
      objp = objp + (1 / (2 * rho)) * eNorm1;
      T eNorm2 = T_ScalarProduct(y, y);
      objd = objd - (rho / 2) * eNorm2;
    }
    bool gap_small = std::abs(objd - objp) < tol * std::abs(objp);
    if (!gap_small) {
      //      std::cerr << "check_stopping, exit, case 4\n";
      return stop;
    }
    // check primal residual
    if (rho == 0) {
      rp = eRecSparse.A(x) - b;
      Out.cntA++;
    }
    T rpnrm = L2_Norm(rp);
    bool p_feasible;
    if (rho > 0) {
      p_feasible = true;
    } else {
      p_feasible = rpnrm < tol * bnrm;
    }
    if (p_feasible) {
      stop = true;
      Out.exit = "Exit: Converged";
    }
    //    std::cerr << "check_stopping, exit, case 5\n";
    return stop;
  };
  if (nonorth) {
    for (size_t i = 0; i < size_t(m); i++)
      y(i) = 0;
    for (size_t i = 0; i < size_t(n); i++)
      Aty(i) = 0;
  }
#ifdef DEBUG
  if (print)
    std::cerr << "--- YALL1 vb6 ---\n";
#endif
  T rdmu = rho / mu;
  T stp = 0;
  T rdmu1 = rdmu + 1;
  MyVector<T> bdmu = b / mu;
  T ddmu = delta / mu;
  Out.cntA = 0;
  Out.cntAt = 0;
  //  std::cerr << "Before central loop\n";
  for (iter = 1; iter <= maxit; iter++) {
    //    std::cerr << "iter=" << iter << " / " << maxit << "\n";
    //    std::cerr << "nonorth=" << nonorth << "\n";
    MyVector<T> xdmu = x / mu;
    //    std::cerr << "mu=" << mu << "\n";
    //    std::cerr << "MeanAbs(x)=" << MeanAbs(x) << "\n";
    //    std::cerr << "MeanAbs(xdmu)=" << MeanAbs(xdmu) << "\n";
    //    std::cerr << "stepfreq=" << stepfreq << "\n";
    //    std::cerr << "rho=" << rho << "\n";
    //    std::cerr << "(rows/cols)x=" << x.rows() << " / " << x.cols() << "\n";
    //    std::cerr << "(rows/cols)xdmu=" << xdmu.rows() << " / " << xdmu.cols()
    //    << "\n";
    if (!nonorth) {
      //      std::cerr << "Passing in orthonormal case\n";
      y = eRecSparse.A(z - xdmu) + bdmu;
      if (rho > 0) {
        y = y / rdmu1;
      } else {
        if (delta > 0) {
          T alpha = std::max(T(0), 1 - ddmu / L2_Norm(y));
          y = alpha * y;
        }
      }
    } else {
      //      std::cerr << "Passing in non-orthonormal case\n";
      MyVector<T> ry = eRecSparse.A(Aty - z + xdmu) - bdmu;
      //      std::cerr << "MeanAbs(ry)=" << MeanAbs(ry) << "\n";
      if (rho > 0)
        ry = ry + rdmu * y;
      if (iter <= 1 || iter % stepfreq == 0) {
        MyVector<T> Atry = eRecSparse.At(ry);
        T eNorm = L2_Norm(Atry);
        T denom = eNorm * eNorm;
        T eScal = T_ScalarProduct(ry, ry);
        if (rho > 0) {
          denom += rdmu * eScal;
        }
        //	std::cerr << "denom=" << denom << "\n";
        //	std::cerr << "eps=" << eps << "\n";
        stp = eScal / (denom + eps);
        Out.cntAt++;
      }
      //      std::cerr << "stp=" << stp << "\n";
      y = y - stp * ry;
    }
    //    std::cerr << "L1_Norm(y)=" << L1_Norm(y) << "\n";
    Aty = eRecSparse.At(y);
    //    std::cerr << "1: MeanAbs(Aty)=" << MeanAbs(Aty) << "\n";

    z = Aty + xdmu;
    //    std::cerr << "2: MeanAbs(z)=" << MeanAbs(z) << "\n";
    z = proj2box(z, w, nonneg, nu, m);
    //    std::cerr << "3: MeanAbs(z)=" << MeanAbs(z) << "\n";

    Out.cntA++;
    Out.cntAt++;

    rd = Aty - z;
    //    std::cerr << "4: MeanAbs(rd)=" << MeanAbs(rd) << "\n";
    xp = x;
    //    std::cerr << "5: MeanAbs(xp)=" << MeanAbs(xp) << "\n";
    //    std::cerr << "gamma=" << gamma << "\n";
    //    std::cerr << "mu=" << mu << "\n";
    x = x + (gamma * mu) * rd;
    //    std::cerr << "MeanAbs(x)=" << MeanAbs(x) << "\n";

    bool stop = check_stopping();
    //    if (iter == 2)
    //      exit(1);

    //    std::cerr << "print=" << print << "\n";
#ifdef DEBUG
    if (print) {
      //      std::cerr << "Before calls to iprint2\n";
      iprint2();
    }
#endif
    //    std::cerr << "stop=" << stop << "\n";
    if (stop)
      break;
  }
  Out.z = z;
  Out.iter = iter + 1;
  if (iter == maxit) {
    Out.exit = "Exit: maxiter";
  }
#ifdef DEBUG
  if (print)
    iprint1_1();
#endif
  Out.x = x;
  //  std::cerr << "Exiting yall1_solve, final\n";
  return Out;
}

template <typename T>
MyVector<T> AMP_SolutionSparseSystem(MySparseMatrix<T> const &SpMat,
                                     MyVector<T> const &eVect) {
  RecSparse<double> eRecSparse = AMP_linear_operators(SpMat);
  //
  RecOptSparse<double> eRecOpt;
  eRecOpt.delta = 0;
  eRecOpt.mu = 0;
  eRecOpt.nu = 0;
  eRecOpt.eps = 2.2204e-16; // The value of eps in matlab which is a default
                            // value (CRAZY Matlab!)
  eRecOpt.rho = 0;
  eRecOpt.tol = 1e-10;
  eRecOpt.gamma = 0;
  eRecOpt.nonneg = false;
  eRecOpt.nonorth = true;
  eRecOpt.print = true;
  eRecOpt.UseWeight = false;
  eRecOpt.basis = false;
  eRecOpt.stepfreq = 1;
  eRecOpt.maxit = 99999;
  eRecOpt.xs = -1;
  OutSolver<double> eRecOut = AMP_yall1(eRecSparse, eVect, eRecOpt);
  //  std::cerr << "After AMP_yall1\n";
  //
  return eRecOut.x;
}

#endif
