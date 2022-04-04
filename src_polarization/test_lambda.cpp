#include "generic_polarization.h"

dual2nd ShiftValues(int const &n, dual2nd const &x,
                    std::function<dual2nd(dual2nd const &)> const &f) {
  dual2nd RetVal = 0;
  for (int i = -n; i <= n; i++) {
    RetVal += f(x - i);
  }
  return RetVal;
}

int main() {
  std::function<dual2nd(dual2nd const &)> f =
      [](dual2nd const &Aval) -> dual2nd { return exp(Aval); };
  dual2nd x = 2;
  autodiff::dual ux = derivative(ShiftValues, wrt(x), at(2, x, f));
  std::cout << "grad_x=" << ux << "\n";
}
