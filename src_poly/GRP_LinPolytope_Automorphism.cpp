// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "GRP_GroupFct.h"
#include "Group.h"
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "QuadField.h"
#include "Permutation.h"
#include "Temp_PolytopeEquiStab.h"
// clang-format on


template<typename T, typename Tgroup>
void full_process_A(std::string const& eFile, std::ostream & os) {
  MyMatrix<T> EXT = ReadMatrixFile<T>(eFile);
  int nbCol = EXT.cols();
  int nbRow = EXT.rows();
  std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
  //
  const bool use_scheme = true;
  Tgroup GRP = LinPolytope_Automorphism<T, use_scheme, Tgroup>(EXT);
  std::cerr << "|GRP|=" << GRP.size() << "\n";
  WriteGroup(os, GRP);
}

template<typename Tgroup>
void full_process_B(std::string const& arith, std::string const& eFile, std::ostream & os) {
  if (arith == "rational") {
    using T = mpq_class;
    return full_process_A<T,Tgroup>(eFile, os);
  }
  if (arith == "Qsqrt5") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 5>;
    return full_process_A<T,Tgroup>(eFile, os);
  }
  if (arith == "Qsqrt2") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 2>;
    return full_process_A<T,Tgroup>(eFile, os);
  }
  std::optional<std::string> opt_realalgebraic =
      get_postfix(arith, "RealAlgebraic=");
  if (opt_realalgebraic) {
    std::string const &FileAlgebraicField = *opt_realalgebraic;
    if (!IsExistingFile(FileAlgebraicField)) {
      std::cerr << "FileAlgebraicField=" << FileAlgebraicField
                << " is missing\n";
      throw TerminalException{1};
    }
    using T_rat = mpq_class;
    HelperClassRealField<T_rat> hcrf(FileAlgebraicField);
    int const idx_real_algebraic_field = 1;
    insert_helper_real_algebraic_field(idx_real_algebraic_field, hcrf);
    using T = RealField<idx_real_algebraic_field>;
    return full_process_A<T,Tgroup>(eFile, os);
  }
  std::cerr << "Failed to find a matching arithmetic\n";
  throw TerminalException{1};
}



int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3 && argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_LinPolytope_Automorphism Arith [EXTIN] [OutGroup]\n";
      std::cerr << "or\n";
      std::cerr << "POLY_LinPolytope_Automorphism Arith [EXTIN]\n";
      std::cerr << "\n";
      std::cerr
          << "EXTIN : The list of vertices (or inequalities for that matter)\n";
      std::cerr << "OutGroup : The automorphism group file\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      std::cerr << "integer  : integer arithmetic on input\n";
      std::cerr << "rational : rational arithmetic on input\n";
      std::cerr << "Qsqrt2   : arithmetic over the field Q(sqrt(2))\n";
      std::cerr << "Qsqrt5   : arithmetic over the field Q(sqrt(5))\n";
      std::cerr << "RealAlgebraic=FileDesc  : For the real algebraic case of a ";
      std::cerr << "field whose description is in FileDesc\n";
      return -1;
    }
    //
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    std::string arith = argv[1];
    std::string eFile = argv[2];
    //
    if (argc == 3) {
      full_process_B<Tgroup>(arith, eFile, std::cerr);
    }
    if (argc == 4) {
      std::ofstream os(argv[3]);
      full_process_B<Tgroup>(arith, eFile, os);
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_LinPolytope_Automorphism\n";
    exit(e.eVal);
  }
  runtime(time1);
}
