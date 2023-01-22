// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "generic_polarization.h"

// The single-variable function for which derivatives are needed

int main() {
  int n = 2;
  int TotDim = n * (n + 1) / 2;
  const double expo = 2;
  std::function<dual2nd(dual2nd const &)> f =
      [&expo](dual2nd const &Aval) -> dual2nd { return exp(-expo * Aval); };
  double onethird = 1.0 / 3.0;
  //  double onethird=0;
  std::vector<dual2nd> cA{dual2nd(onethird), dual2nd(onethird), dual2nd(2),
                          dual2nd(1), dual2nd(2)};
  std::vector<int> Ranges{1, 1};

  //  dual2nd u = Compute_Polarization(2, Ranges, c, A, f);

  dual2nd u = Compute_Polarization(n, Ranges, cA, f);
  double u_d = u.val.val;
  std::cout << "u=" << u << " u_d=" << u_d << "\n";

  int nPtot = n + TotDim;
  MyVector<double> gradC(n);
  for (int i = 0; i < nPtot; i++) {
    autodiff::dual ux =
        derivative(Compute_Polarization, wrt(cA[i]), at(2, Ranges, cA, f));
    double ux_val = ux.val;
    double ux_grad = ux.grad;
    std::cout << "grad_cA(" << i << "), ux=" << ux << " ux_val=" << ux_val
              << " ux_grad=" << ux_grad << "\n";
  }
  MyMatrix<double> Hess(nPtot, nPtot);
  for (int i = 0; i < nPtot; i++)
    for (int j = 0; j < nPtot; j++) {
      double eVal;
      //      std::cout << "Before i=" << i << " j=" << j << "\n";
      if (i == j) {
        eVal = derivative(Compute_Polarization, autodiff::wrt<2>(cA[i]),
                          at(n, Ranges, cA, f));
      } else {
        eVal = derivative(Compute_Polarization, autodiff::wrt(cA[i], cA[j]),
                          at(n, Ranges, cA, f));
      }
      //      std::cout << "After i=" << i << " j=" << j << "\n";
      Hess(i, j) = eVal;
    }
  while (true) {
  }
}
