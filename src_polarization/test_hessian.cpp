#include "generic_polarization.h"

autodiff::var SumFprod(int const &n, double const &a,
                       const autodiff::VectorXvar &v,
                       std::function<dual2nd(dual2nd const &)> const &f) {
  autodiff::var RetVal = 0;
  for (int i = 0; i < n; i++) {
    RetVal += f(a * v[i] * v[j]);
  }
  return RetVal;
}

int main() {
  std::function<dual2nd(dual2nd const &)> f =
      [](dual2nd const &Aval) -> dual2nd { return exp(Aval); };
  int n = 4;
  double a = 4.5;
  //
  autodiff::VectorXvar v(n);
  for (int i = 0; i < n; i++)
    v[i] = i + 1;
  //
  dual u;
  autodiff::VectorXd gx = gradient(SumFprod, wrt(v), at(n, a, v, f), u);
}
