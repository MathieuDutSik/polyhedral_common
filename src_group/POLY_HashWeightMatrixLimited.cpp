// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "GRP_GroupFct.h"
#include "GRP_GroupFile.h"
#include "Group.h"
#include "Permutation.h"
#include "PolytopeEquiStab.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

/*
  The hash of faces given by the limited weight matrix of WeightMatrixLimited.h,
  built either from the polytope (GetWeightMatrixLimited) or from a
  permutation group acting on the vertices. The hash is invariant under the
  symmetry group of the polytope, resp. under the group, which is what the
  CI checks.
 */

template <typename T, typename Tgroup>
WeightMatrixLimited<true, T> GetWMat(std::string const &input,
                                     std::string const &FileIn,
                                     size_t max_offdiag) {
  if (input == "EXT") {
    MyMatrix<T> EXT = ReadMatrixFile<T>(FileIn);
    return GetWeightMatrixLimited(EXT, max_offdiag, std::cerr);
  }
  if (input == "GRP") {
    Tgroup GRP = ReadGroupFile<Tgroup>(FileIn);
    return WeightMatrixLimited<true, T>(GRP, max_offdiag, std::cerr);
  }
  std::cerr << "POLY_HashWeightMatrixLimited: input=" << input
            << " but allowed values are EXT and GRP\n";
  throw TerminalException{1};
}

template <typename T, typename Tgroup>
void process(std::string const &input, std::string const &FileIn,
             std::string const &FileFaces, size_t max_offdiag,
             std::string const &OutFormat, std::ostream &os) {
  WeightMatrixLimited<true, T> WMatLimited =
      GetWMat<T, Tgroup>(input, FileIn, max_offdiag);
  MyMatrix<int> ListFace = ReadMatrixFile<int>(FileFaces);
  size_t n_face = ListFace.rows();
  size_t len = ListFace.cols();
  std::vector<size_t> ListHash;
  for (size_t i_face = 0; i_face < n_face; i_face++) {
    Face eFace(len);
    for (size_t i = 0; i < len; i++) {
      eFace[i] = ListFace(i_face, i);
    }
    ListHash.push_back(WMatLimited.get_hash(eFace));
  }
  if (OutFormat == "GAP") {
    os << "return [";
    for (size_t i_face = 0; i_face < n_face; i_face++) {
      if (i_face > 0) {
        os << ",";
      }
      os << ListHash[i_face];
    }
    os << "];\n";
    return;
  }
  if (OutFormat == "CPP") {
    os << n_face << "\n";
    for (auto &eHash : ListHash) {
      os << eHash << "\n";
    }
    return;
  }
  std::cerr << "POLY_HashWeightMatrixLimited: No matching entry for "
               "OutFormat="
            << OutFormat << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    std::string input;
    std::string FileIn;
    std::string FileFaces;
    size_t max_offdiag = 0;
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    add_rational_arithmetics(dispatch, [&]<typename T>() -> void {
      process<T, Tgroup>(input, FileIn, FileFaces, max_offdiag, OutFormat,
                         *os_ptr);
    });
    if (argc != 6 && argc != 8) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_HashWeightMatrixLimited [arith] [input] [FileIn] "
                   "[FileFaces] [max_offdiag] [OutFormat] [OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "POLY_HashWeightMatrixLimited [arith] [input] [FileIn] "
                   "[FileFaces] [max_offdiag]\n";
      std::cerr << "\n";
      std::cerr << "input       : EXT (FileIn is the list of vertices) or GRP "
                   "(FileIn is\n";
      std::cerr << "              a permutation group acting on the "
                   "vertices)\n";
      std::cerr << "FileFaces   : the faces, one 0/1 row per face\n";
      std::cerr << "max_offdiag : the bound on the number of off diagonal "
                   "weights\n";
      std::cerr << "OutFormat   : GAP (default) or CPP\n";
      std::cerr << "OutFile     : the output file, or stderr (default) / "
                   "stdout\n";
      std::cerr << "\n";
      std::cerr << "     ------- arith -------\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      return -1;
    }
    std::string arith = argv[1];
    input = argv[2];
    FileIn = argv[3];
    FileFaces = argv[4];
    max_offdiag = ParseScalar<size_t>(argv[5]);
    if (argc == 8) {
      OutFormat = argv[6];
      OutFile = argv[7];
    }
    auto f = [&](std::ostream &os) -> void {
      os_ptr = &os;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    std::cerr << "Normal termination of POLY_HashWeightMatrixLimited\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_HashWeightMatrixLimited\n";
    exit(e.eVal);
  }
  runtime(time1);
}
