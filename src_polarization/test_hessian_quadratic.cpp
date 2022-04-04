#include "generic_polarization.h"

// The single-variable function for which derivatives are needed

dual2nd f(MyMatrix<double> const &A, std::vector<dual2nd> const &v) {
  int n = A.rows();
  dual2nd RetVal = 0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      RetVal += A(i, j) * v[i] * v[j];
    }
  return RetVal;
}

int main() {
  int n = 3;
  MyMatrix<double> A(n, n);
  //
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++) {
      int val = 1 + i * j + i * i + j;
      A(i, j) = double(val);
      A(j, i) = double(val);
    }
  //
  std::vector<dual2nd> v;
  for (int i = 0; i < n; i++) {
    int val = 10 + i;
    dual2nd val_d = double(val);
    v.push_back(val_d);
    double val_dd = val_d.val.val;
    std::cerr << "i=" << i << " val=" << val << " val_dd=" << val_dd << "\n";
  }
  //
  double val = 0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      //      val += A(i,j) * v(i).val.val * v(j).val.val;
    }
  std::cout << "val=" << val << "\n";
  //
  dual2nd u = f(A, v);
  double u_d = u.val.val;
  std::cout << "u=" << u << " u_d=" << u_d << "\n";
  //
  MyVector<double> gradC(n);
  for (int i = 0; i < n; i++) {
    autodiff::dual ux = derivative(f, wrt(v[i]), at(A, v));
    double ux_val = ux.val;
    double ux_grad = ux.grad;
    std::cout << "grad_cA(" << i << "), ux=" << ux << " ux_val=" << ux_val
              << " ux_grad=" << ux_grad << "\n";
  }
  MyMatrix<double> Hess(n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      double eVal;
      if (i == j) {
        eVal = derivative(f, autodiff::wrt<2>(v[i]), at(A, v));
      } else {
        eVal = derivative(f, autodiff::wrt(v[i], v[j]), at(A, v));
      }
      Hess(i, j) = eVal;
    }
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      std::cout << "i=" << i << " j=" << j << " A(i,j)=" << A(i, j)
                << " Hess(i,j)=" << Hess(i, j) << "\n";
    }
}
