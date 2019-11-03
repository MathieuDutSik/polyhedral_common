#include "generic_polarization.h"





// The single-variable function for which derivatives are needed

int main()
{
  int n=2;
  int TotDim=n*(n+1)/2;
  const double expo = 2;
  std::function<dual2nd(dual2nd const&)> f = [&expo](dual2nd const& Aval) -> dual2nd {
    return Aval;
    //    return exp(-expo * Aval);
  };
  double onethird=double(1) / double(3);
  //  double onethird=0;
  std::vector<dual2nd> cA{dual2nd(onethird), dual2nd(onethird), dual2nd(2), dual2nd(1), dual2nd(2)};
  std::vector<int> Ranges{1, 1};

  //  dual2nd u = Compute_Polarization(2, Ranges, c, A, f);


  for (int i=0; i<n; i++) {
    autodiff::dual ux = derivative(Compute_Polarization, wrt(cA[i]), at(2, Ranges, cA, f));
    std::cout << "grad_c(" << i << ")=" << ux << "\n";
  }
  for (int idx=0; idx<TotDim; idx++) {
    autodiff::dual ux = derivative(Compute_Polarization, wrt(cA[n+idx]), at(2, Ranges, cA, f));
    std::cout << "grad_A(" << idx << ")=" << ux << "\n";
  }
}
