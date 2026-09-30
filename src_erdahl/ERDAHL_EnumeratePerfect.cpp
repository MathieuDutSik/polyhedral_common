// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "erdahl_enumeration.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

/*
  The perfect Delaunay polyhedra of Z^n relative to a space W of degree 2
  functions, up to the affine integral transformations preserving W.
  * space = "full": all the functions a + b.x + C[x].
  * space = "centered": the functions a + C[x - c], whose zero sets are the
    Delaunay polyhedra symmetric around c. The center c is read from a file
    in the format of ReadVectorFile (n followed by the n entries, rational
    entries written as p/q).
 */

template <typename T, typename Tint, typename Tgroup>
void process(int n, std::string const &space, std::string const &FileCenter,
             std::string const &OutFile, std::string const &FileDualDesc) {
  ErdahlFunctionSpace<T> W = [&]() -> ErdahlFunctionSpace<T> {
    if (space == "full") {
      return erdahl_full_space<T>(n);
    }
    if (space == "centered") {
      MyVector<T> c = ReadVectorFile<T>(FileCenter);
      if (c.size() != n) {
        std::cerr << "ERDAHL_EnumeratePerfect: the center should be of "
                     "length n\n";
        throw TerminalException{1};
      }
      return erdahl_centered_space(c);
    }
    std::cerr << "ERDAHL_EnumeratePerfect: space should be full or centered, "
                 "not "
              << space << "\n";
    throw TerminalException{1};
  }();
  std::vector<DelaunayPolyhedron<T, Tint>> l_perf =
      erdahl_enumerate_perfect<T, Tint, Tgroup>(W, FileDualDesc, std::cerr);
  auto f = [&](std::ostream &os) -> void {
    os << "return [";
    for (size_t i = 0; i < l_perf.size(); i++) {
      if (i > 0) {
        os << ",\n";
      }
      DelaunayPolyhedron<T, Tint> const &D = l_perf[i];
      os << "rec(EXT:=" << StringMatrixGAP(D.EXT)
         << ", L:=" << StringMatrixGAP(D.L)
         << ", F:=" << StringMatrixGAP(D.F) << ")";
    }
    os << "];\n";
  };
  if (OutFile == "stderr") {
    f(std::cerr);
  } else {
    if (OutFile == "stdout") {
      f(std::cout);
    } else {
      std::ofstream os(OutFile);
      f(os);
    }
  }
  std::cerr << "ERDAHL_EnumeratePerfect: n=" << n << " space=" << space
            << " |l_perf|=" << l_perf.size() << "\n";
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    auto print_usage = [&]() -> void {
      std::cerr << "ERDAHL_EnumeratePerfect [n] full [OutFile] "
                   "[FileDualDesc]\n";
      std::cerr << "or\n";
      std::cerr << "ERDAHL_EnumeratePerfect [n] centered [FileCenter] "
                   "[OutFile] [FileDualDesc]\n";
      std::cerr << "OutFile can be stderr or stdout\n";
      std::cerr << "FileDualDesc (optional) is the namelist file of the "
                   "heuristics of the dual descriptions\n";
    };
    if (argc < 4) {
      print_usage();
      return -1;
    }
    int n = ParseScalar<int>(argv[1]);
    std::string space = argv[2];
    std::string FileCenter = "unset";
    std::string OutFile;
    std::string FileDualDesc = "unset";
    int pos = 3;
    if (space == "centered") {
      FileCenter = argv[pos];
      pos++;
    }
    if (argc != pos + 1 && argc != pos + 2) {
      print_usage();
      return -1;
    }
    OutFile = argv[pos];
    if (argc == pos + 2) {
      FileDualDesc = argv[pos + 1];
    }
    using T = mpq_class;
    using Tint = mpz_class;
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint_grp = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint_grp>;
    process<T, Tint, Tgroup>(n, space, FileCenter, OutFile, FileDualDesc);
    std::cerr << "Normal termination of ERDAHL_EnumeratePerfect\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in ERDAHL_EnumeratePerfect\n";
    exit(e.eVal);
  }
  runtime(time);
}
