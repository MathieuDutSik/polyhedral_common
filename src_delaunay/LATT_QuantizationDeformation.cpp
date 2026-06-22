// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
#include "QuantizationDeformation.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

template <typename T, typename Tint, typename Tgroup>
void process(std::string const &Qfile, std::string const &Hfile,
             std::ostream &os) {
  MyMatrix<T> Q = ReadMatrixFile<T>(Qfile);
  MyMatrix<T> H = ReadMatrixFile<T>(Hfile);
  int n = Q.rows();
  os << "QDEF: n=" << n << "\n";
  //
  // 1. The common symmetry group of (Q, H).
  //
  std::vector<MyMatrix<T>> gens_T =
      compute_qh_symmetry_gens<T, Tint, Tgroup>(Q, H, os);
  os << "QDEF: |symmetry generators|=" << gens_T.size() << "\n";
  //
  // 2. The T-space of the invariant forms.
  //
  LinSpaceMatrix<T> LinSpa = build_qh_tspace<T, Tint, Tgroup>(Q, gens_T, os);
  os << "QDEF: T-space dimension=" << LinSpa.ListMat.size() << "\n";
  // Check that Q and H are in the T-space.
  MyVector<T> cQ = LINSPA_GetVectorOfMatrixExpression(LinSpa, Q);
  MyVector<T> cH = LINSPA_GetVectorOfMatrixExpression(LinSpa, H);
  os << "QDEF: Q and H are in the T-space\n";
  //
  // 3. The iso-Delaunay segment.
  //
  T t_init(1);
  IsoDelaunaySegment<T, Tgroup> seg =
      find_iso_delaunay_segment<T, Tint, Tgroup>(LinSpa, Q, H, t_init, os);
  os << "QDEF: number of Delaunay orbits=" << seg.DT.l_dels.size() << "\n";
  os << "QDEF: number of defining inequalities=" << seg.ListIneq.size() << "\n";
  os << "QDEF: bounded segment=" << seg.bounded << "\n";
  os << "QDEF: tmax=" << seg.tmax << "\n";
}

template <typename T, typename Tint> void process_B(std::string const &Qfile,
                                                    std::string const &Hfile) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tint_grp = mpz_class;
  using Tgroup = permutalib::Group<Telt, Tint_grp>;
  return process<T, Tint, Tgroup>(Qfile, Hfile, std::cerr);
}

void process_C(std::string const &arith, std::string const &Qfile,
               std::string const &Hfile) {
  if (arith == "gmp") {
    using T = mpq_class;
    using Tint = mpz_class;
    return process_B<T, Tint>(Qfile, Hfile);
  }
  if (arith == "safe") {
    using T = Rational<SafeInt64>;
    using Tint = SafeInt64;
    return process_B<T, Tint>(Qfile, Hfile);
  }
  std::cerr << "LATT_QuantizationDeformation: no match for arith=" << arith
            << "\n";
  std::cerr << "Available types: gmp, safe\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_QuantizationDeformation [arith] [Qfile] [Hfile]\n";
      std::cerr << "     arith   gmp, safe\n";
      std::cerr << "     Qfile   the positive definite form Q\n";
      std::cerr << "     Hfile   the symmetric perturbation H\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string Qfile = argv[2];
    std::string Hfile = argv[3];
    process_C(arith, Qfile, Hfile);
    std::cerr << "Normal termination of LATT_QuantizationDeformation\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_QuantizationDeformation\n";
    exit(e.eVal);
  }
  runtime(time);
}
