// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "erdahl_polyhedron.h"
// clang-format on

/*
  The zero set of the degree 2 function f(x) = F[(1, x)] on Z^n, with F the
  (n+1) x (n+1) symmetric matrix of the file, written in GAP format:
  --- nonnegative:=true when f >= 0 on Z^n, with min the minimum of f over
      Z^n and the zero set EXT + L: EXT the representatives as rows (1, x)
      and L a basis of the lattice. EXT has no rows when min > 0.
  --- nonnegative:=false when f takes negative values, with the reason,
      lattice points (1, x) where f is negative and, when f is unbounded
      below, the directions along which it decreases. min is given when
      f is bounded below.
 */

template <typename Tint>
std::string StringListVectorGAP(std::vector<MyVector<Tint>> const &l_v) {
  std::string str = "[";
  bool IsFirst = true;
  for (auto &v : l_v) {
    if (!IsFirst) {
      str += ", ";
    }
    IsFirst = false;
    str += StringVectorGAP(v);
  }
  return str + "]";
}

template <typename T, typename Tint>
void process(std::string const &FileF, std::ostream &os_out) {
  MyMatrix<T> F = ReadMatrixFile<T>(FileF);
  if (F.rows() != F.cols() || !IsSymmetricMatrix(F)) {
    std::cerr << "ERDAHL_zero_set: the matrix F should be square and "
                 "symmetric\n";
    throw TerminalException{1};
  }
  auto res_min = erdahl_lattice_minimum<T, Tint>(F, std::cerr);
  auto res = erdahl_zero_set<T, Tint>(F, std::cerr);
  if (res.is_ok()) {
    ErdahlLatticeSet<Tint> const &zero_set = res.get_ok();
    os_out << "return rec(nonnegative:=true, min:="
           << res_min.get_ok().min_value
           << ", EXT:=" << StringMatrixGAP(zero_set.EXT)
           << ", L:=" << StringMatrixGAP(zero_set.L) << ");\n";
    return;
  }
  ErdahlNegativeWitness<Tint> const &witness = res.get_err();
  os_out << "return rec(nonnegative:=false";
  if (res_min.is_ok()) {
    os_out << ", min:=" << res_min.get_ok().min_value;
  }
  os_out << ", reason:=\"" << witness.reason << "\""
         << ", points:=" << StringListVectorGAP(witness.points)
         << ", directions:=" << StringListVectorGAP(witness.directions)
         << ");\n";
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 2 && argc != 3) {
      std::cerr << "ERDAHL_zero_set [FileF] [OutFile]\n";
      std::cerr << "or\n";
      std::cerr << "ERDAHL_zero_set [FileF]\n";
      std::cerr << "\n";
      std::cerr << "FileF   : the (n+1) x (n+1) symmetric matrix F of the "
                   "function f(x) = F[(1, x)]\n";
      std::cerr << "OutFile : stdout, stderr or the name of a file "
                   "(default stderr)\n";
      std::cerr << "\n";
      std::cerr << "The output is the zero set of f on Z^n in GAP format\n";
      return -1;
    }
    std::string FileF = argv[1];
    std::string OutFile = "stderr";
    if (argc == 3) {
      OutFile = argv[2];
    }
    // The FLINT rationals and integers when available, being the faster.
#ifdef ENABLE_FLINT_SUPPORT
    using T = fmpq_class;
    using Tint = fmpz_class;
#else
    using T = mpq_class;
    using Tint = mpz_class;
#endif
    auto f_print = [&](std::ostream &os_out) -> void {
      process<T, Tint>(FileF, os_out);
    };
    FILE_PrintStderrStdoutFile(OutFile, f_print);
    std::cerr << "Normal termination of ERDAHL_zero_set\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in ERDAHL_zero_set\n";
    exit(e.eVal);
  }
  runtime(time);
}
