// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Tspace_General.h"
#include "Permutation.h"
#include "Group.h"
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
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "TSPACE_Stabilizer [FileTspace] [FileGram]\n";
      return -1;
    }
    std::string FileTspace = argv[1];
    std::string FileGram = argv[2];
    //
    LinSpaceMatrix<T> LinSpa = ReadLinSpaceFile<T>(FileTspace);
    MyMatrix<T> eMat = ReadMatrixFile<T>(FileGram);
    std::vector<MyMatrix<T>> ListGen = LINSPA_ComputeStabilizer<T, Tint, Tgroup>(LinSpa, eMat, std::cerr);
    std::cerr << "|ListGen|=" << ListGen.size() << "\n";
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TSPACE_Stabilizer\n";
    exit(e.eVal);
  }
  runtime(time);
}
