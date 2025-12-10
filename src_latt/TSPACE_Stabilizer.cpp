// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Tspace_General.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

template <typename T>
void write_group(std::vector<MyMatrix<T>> const &LGen,
                 std::string const &OutFormat, std::ostream &os) {
  if (OutFormat == "count") {
    os << "number of generators=" << LGen.size() << "\n";
    return;
  }
  if (OutFormat == "GAP") {
    os << "return Group([";
    bool IsFirst = true;
    for (auto &eGen : LGen) {
      if (!IsFirst) {
        os << ",\n";
      }
      os << StringMatrixGAP_line(eGen);
      IsFirst = false;
    }
    os << "]);\n";
    return;
  }
  std::cerr << "Failed to find a matching format\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    using T = mpq_class;
    using Tint = mpz_class;
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint_grp = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint_grp>;
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "TSPACE_Stabilizer [FileTspace] [FileGram]\n";
      std::cerr << "or\n";
      std::cerr << "TSPACE_Stabilizer [FileTspace] [FileGram] [OutFormat] "
                   "[FileOut]\n";
      return -1;
    }
    std::string FileTspace = argv[1];
    std::string FileGram = argv[2];
    std::string OutFormat = "count";
    std::string FileOut = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      FileOut = argv[4];
    }
    //
    LinSpaceMatrix<T> LinSpa = ReadLinSpaceFile<T>(FileTspace, std::cerr);
    MyMatrix<T> eMat = ReadMatrixFile<T>(FileGram);
    std::vector<MyMatrix<T>> ListGen =
        LINSPA_ComputeStabilizer<T, Tint, Tgroup>(LinSpa, eMat, std::cerr);
    auto f = [&](std::ostream &os_out) -> void {
      write_group(ListGen, OutFormat, os_out);
    };
    print_stderr_stdout_file(FileOut, f);
    std::cerr << "Normal termination of TSPACE_Stabilizer\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TSPACE_Stabilizer\n";
    exit(e.eVal);
  }
  runtime(time);
}
