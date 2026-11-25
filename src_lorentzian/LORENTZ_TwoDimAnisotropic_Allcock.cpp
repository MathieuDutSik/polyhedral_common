// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MAT_Matrix.h"
#include "LatticeStabEquiCan.h"
#include "PolytopeEquiStab.h"
#include "two_dim_lorentzian.h"
// clang-format on

// The algorithms of Section 8 of Allcock, "An alternative to Vinberg's
// algorithm", on the 2-dimensional Lorentzian forms a x^2 + 2 b x y + c y^2
// listed in an input file, the results being written in GAP format.

template <typename T> MyMatrix<T> GetGram(int a, int b, int c) {
  MyMatrix<T> G(2, 2);
  G(0, 0) = a;
  G(0, 1) = b;
  G(1, 0) = b;
  G(1, 1) = c;
  return G;
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 4) {
      std::cerr << "LORENTZ_TwoDimAnisotropic_Allcock [Mode] [FileIn] "
                   "[FileOut]\n";
      std::cerr << "\n";
      std::cerr << "Mode anisotropic: FileIn has the number of cases, then "
                   "one line a b c M per\n";
      std::cerr << "  anisotropic form a x^2 + 2 b x y + c y^2. For each, the "
                   "positive vector r0 of\n";
      std::cerr << "  GetPositiveVector, and the generator g and the list R "
                   "of AnisotropicComplete\n";
      std::cerr << "  (Algorithm 8.7) for the norm bound M.\n";
      std::cerr << "Mode first_next: FileIn has the number of cases, then "
                   "one line a b c x y N per\n";
      std::cerr << "  form. For each, get_first_next_vector from r0 = (x,y) "
                   "for the norm N.\n";
      std::cerr << "FileOut: the GAP readable list of the results\n";
      throw TerminalException{1};
    }
    using T = mpq_class;
    using Tint = mpz_class;
    std::string Mode = argv[1];
    std::string FileIn = argv[2];
    std::string FileOut = argv[3];
    std::ifstream is(FileIn);
    int n_case;
    is >> n_case;
    std::ofstream os(FileOut);
    os << "return [";
    for (int i_case = 0; i_case < n_case; i_case++) {
      if (i_case > 0)
        os << ",\n";
      if (Mode == "anisotropic") {
        int a, b, c, M;
        is >> a >> b >> c >> M;
        MyMatrix<T> G = GetGram<T>(a, b, c);
        MyVector<Tint> r0 = GetPositiveVector<T, Tint>(G);
        std::optional<std::pair<MyMatrix<Tint>, std::vector<MyVector<Tint>>>>
            opt = AnisotropicComplete<T, Tint>(G, T(M));
        os << "rec(r0:=" << StringVectorGAP(r0);
        if (opt) {
          os << ", g:=" << StringMatrixGAP(opt->first)
             << ", R:=" << StringMatrixGAP(MatrixFromVectorFamily(opt->second));
        }
        os << ")";
        continue;
      }
      if (Mode == "first_next") {
        int a, b, c, x, y, N;
        is >> a >> b >> c >> x >> y >> N;
        MyMatrix<T> G = GetGram<T>(a, b, c);
        MyVector<Tint> r0(2);
        r0(0) = x;
        r0(1) = y;
        std::optional<MyVector<Tint>> opt =
            get_first_next_vector<T, Tint>(G, r0, T(N));
        if (opt) {
          os << "rec(v:=" << StringVectorGAP(*opt) << ")";
        } else {
          os << "rec()";
        }
        continue;
      }
      std::cerr << "Mode=" << Mode << " but allowed modes are anisotropic and "
                   "first_next\n";
      throw TerminalException{1};
    }
    os << "];\n";
    std::cerr << "Normal end of LORENTZ_TwoDimAnisotropic_Allcock\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LORENTZ_TwoDimAnisotropic_Allcock\n";
    exit(e.eVal);
  }
  runtime(time);
}
