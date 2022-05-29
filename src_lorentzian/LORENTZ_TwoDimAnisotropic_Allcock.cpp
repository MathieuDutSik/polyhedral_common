// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "MAT_Matrix.h"
#include "MatrixCanonicalForm.h"
#include "NumberTheory.h"
#include "Temp_PolytopeEquiStab.h"
#include "two_dim_lorentzian.h"

int main(int argc, char *argv[]) {
  try {
    if (argc != 5) {
      std::cerr << "TEST_Anisotropic a b c M\n";
      throw TerminalException{1};
    }
    using T = mpq_class;
    using Tint = mpz_class;
    //
    int a, b, c, M;
    sscanf(argv[1], "%d", &a);
    sscanf(argv[2], "%d", &b);
    sscanf(argv[3], "%d", &c);
    sscanf(argv[4], "%d", &M);
    T a_T = a;
    T b_T = b;
    T c_T = c;
    T M_T = M;
    std::cerr << "Working with the quadratic form\n";
    std::cerr << "q = " << a << " x^2 + 2 * " << b << " x y + " << c
              << " y^2\n";
    std::cerr << "M = " << M << "\n";
    std::cerr << "\n";
    T discriminant = b * b - a * c;
    if (discriminant <= 0) {
      std::cerr << "discriminant = " << discriminant << "\n";
      throw TerminalException{1};
    }
    //
    MyVector<Tint> r(2);
    auto set_r = [&]() -> void {
      if (a_T > 0) { // Then (1,0) is of positive norm
        r(0) = 1;
        r(1) = 0;
        return;
      }
      if (b_T > 0) {
        r(0) = 0;
        r(1) = 1;
        return;
      }
      std::cerr
          << "More cases need to be covered if test were to be exhaustive\n";
      throw TerminalException{1};
    };
    set_r();
    std::cerr << "r = " << r(0) << "," << r(1) << "\n";
    //
    // A very basic algorithm for computing
    MyVector<Tint> l(2);
    auto set_l = [&]() -> void {
      int shift = 1;
      while (true) {
        for (int i = -shift; i <= shift; i++)
          for (int j = -shift; j <= shift; j++) {
            l(0) = i;
            l(1) = j;
            MyMatrix<Tint> eMat = MatrixFromVectorFamily<Tint>({r, l});
            if (DeterminantMat(eMat) > 0) {
              T eNorm = a_T * T(i * i) + 2 * b_T * T(i * j) + c_T * T(j * j);
              if (eNorm < M_T)
                return;
            }
          }
        shift++;
      }
    };
    set_l();
    std::cerr << "l = " << l(0) << "," << l(1) << "\n";
    //
    MyMatrix<T> G(2, 2);
    G(0, 0) = a_T;
    G(0, 1) = b_T;
    G(1, 0) = b_T;
    G(1, 1) = c_T;
    //    std::optional<std::pair<MyMatrix<Tint>,std::vector<MyVector<Tint>>>>
    //    pair_opt = Anisotropic(G, M_T, r, l);
    std::optional<std::pair<MyMatrix<Tint>, std::vector<MyVector<Tint>>>>
        pair_opt = AnisotropicComplete<T, Tint>(G, M_T);
    if (pair_opt) {
      std::cerr << "g =\n";
      WriteMatrix(std::cerr, pair_opt->first);
      std::cerr << "|R| = " << pair_opt->second.size() << "\n";
      std::cerr << "R =";
      for (auto &eVect : pair_opt->second) {
        Tint x = eVect(0);
        Tint y = eVect(1);
        T norm = a_T * x * x + 2 * b_T * x * y + c_T * y * y;
        std::cerr << " [" << x << "," << y << "] : " << norm;
      }
      std::cerr << "\n";
    } else {
      std::cerr << "No solution were found\n";
    }
    std::cerr << "Normal end of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something went wrong\n";
    exit(e.eVal);
  }
}
