// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Tspace_General.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    using T = mpq_class;
    using Tint = mpz_class;
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint_grp = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint_grp>;
    FullNamelist eFull = NAMELIST_GetOneTSPACE();
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "TSPACE_Equivalence [FileLinSpa] [FileMat1] [FileMat2]\n";
      std::cerr << "\n";
      return -1;
    }
    std::string FileTspace = argv[1];
    std::string FileGram1 = argv[2];
    std::string FileGram2 = argv[3];
    //
    LinSpaceMatrix<T> LinSpa = ReadLinSpaceFile<T>(FileTspace);
    MyMatrix<T> eMat1 = ReadMatrixFile<T>(FileGram1);
    MyMatrix<T> eMat2 = ReadMatrixFile<T>(FileGram2);
    std::optional<MyMatrix<T>> opt = LINSPA_TestEquivalenceGramMatrix<T,Tint,Tgroup>(LinSpa, eMat1, eMat2, std::cerr);
    if (opt) {
      std::cerr << "It is equivalence\n";
    } else {
      std::cerr << "It is not equivalence\n";
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TSPACE_Equivalence\n";
    exit(e.eVal);
  }
  runtime(time);
}
