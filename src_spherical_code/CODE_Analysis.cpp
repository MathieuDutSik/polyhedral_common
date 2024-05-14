// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheoryQuadField.h"
#include "Group.h"
#include "Permutation.h"
// clang-format on

template<typename T, typename Tgroup>
void process_entry_type(std::string const& FileCode) {
  MyMatrix<T> CODE = ReadMatrixFile<T>(FileCode);
  int nbEnt = CODE.rows();
  int dim = CODE.cols();
  for (int iEnt=0; iEnt<nbEnt; iEnt++) {
    std::map<T, size_t> map;
    for (int jEnt=0; jEnt<nbEnt; jEnt++) {
      if (iEnt != jEnt) {
        T scal(0);
        for (int i=0; i<dim; i++) {
          scal += CODE(iEnt, i) * CODE(jEnt, i);
        }
        map[scal] += 1;
      }
    }
    std::cerr << "iEnt=" << iEnt;
    for (auto & kv : map) {
      std::cerr << " (" << kv.first << "/" << kv.second << ")";
    }
    std::cerr << "\n";
  }

}





template <typename Tgroup>
void process_entry(std::string const& arith, std::string const& FileCode) {
  if (arith == "rational") {
    using T = mpq_class;
    return process_entry_type<T, Tgroup>(FileCode);
  }
  if (arith == "Qsqrt3") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 3>;
    return process_entry_type<T, Tgroup>(FileCode);
  }
  std::cerr << "Failed to find a matching arithmetic\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3) {
      std::cerr << "CODE_Analysis [arith] [FileCODE]\n";
      throw TerminalException{1};
    }
    using Tidx = int32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    std::string arith = argv[1];
    std::string FileCode = argv[2];
    process_entry<Tgroup>(arith, FileCode);
    std::cerr << "Normal termination of the program time=" << time << "\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in CODE_Analysis time=" << time << "\n";
    exit(e.eVal);
  }
}
