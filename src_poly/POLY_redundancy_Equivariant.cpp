// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Group.h"
#include "NumberTheory.h"
#include "POLY_RedundancyElimination.h"
#include "Permutation.h"

int main(int argc, char *argv[]) {
  SingletonTime time1;
  try {
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_redundancy_Equivariant [DATAIN] [GRPIN] [DATAOUT]\n";
      std::cerr << "\n";
      std::cerr << "DATAIN  : The polyhedral cone inequalities\n";
      std::cerr << "GRPIN   : The permutation group on the inequalities\n";
      std::cerr << "DATAOUT : The list of irredundant facets\n";
      return -1;
    }
    //
    using T = mpq_class;
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    //
    std::string FileIneq = argv[1];
    std::ifstream is1(FileIneq);
    MyMatrix<T> EXT = ReadMatrix<T>(is1);
    std::cerr << "We have EXT\n";
    std::string FileGRP = argv[2];
    std::ifstream is2(FileGRP);
    Tgroup GRP = ReadGroup<Tgroup>(is2);
    std::cerr << "We have GRP\n";
    //
    Face status = GetNonRedundant_Equivariant(EXT, GRP);
    //
    std::ofstream os(argv[3]);
    os << "return [";
    size_t nbIrred = status.count();
    boost::dynamic_bitset<>::size_type pos = status.find_first();
    for (size_t i = 0; i < nbIrred; i++) {
      if (i > 0)
        os << ",";
      int eVal = pos + 1;
      os << eVal;
      pos = status.find_next(pos);
    }
    os << "];\n";
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_redundancy_Equivariant\n";
    exit(e.eVal);
  }
  runtime(time1);
}
