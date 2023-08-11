// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
# include "NumberTheoryBoostGmpInt.h"
#else
# include "NumberTheory.h"
#endif
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "POLY_PolytopeFct.h"
#include "POLY_DirectDualDesc.h"
#include "GRP_DoubleCoset.h"
#include "GRP_GroupFct.h"
#include "Group.h"
#include "Permutation.h"
// clang-format on

template <typename T, typename Tgroup>
void process(std::string const &eFileI, std::string const& eFileG, std::string const& ansProg, std::string const& OutFormat, std::ostream &os) {
  MyMatrix<T> EXT = ReadMatrixFile<T>(eFileI);
  vectface vf = DirectFacetComputationIncidence(EXT, ansProg, os);
  Tgroup GRP = ReadGroupFile<Tgroup>(eFileG);
  std::unordered_set<Face> SetFace;
  for (auto & eFace : vf) {
    SetFace.insert(eFace);
  }
  vectface vfo = OrbitSplittingSet_T(SetFace, GRP);
  if (OutFormat == "GAP") {
    os << "return ";
    VectVectInt_Gap_Print(os, vfo);
    os << ";\n";
    return;
  }
  if (OutFormat == "Oscar") {
    MyMatrix<int> M = VectfaceAsMatrix(vfo);
    WriteMatrix(os, M);
    return;
  }
  std::cerr << "POLY_dual_description_group : Failed to find a matching entry for OutFormat=" << OutFormat << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 5 && argc != 7) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_dual_description_group arith command [FileEXT] [FileGRP] [OutFormat] [FileOUT]\n";
      std::cerr << "or\n";
      std::cerr << "POLY_dual_description_group arith command [FileEXT] [FileGRP]\n";
      std::cerr << "\n";
      std::cerr << "with:\n";
      std::cerr << "arith     : the chosen arithmetic\n";
      std::cerr << "command   : the program used for computing the dual description\n";
      std::cerr << "FileEXT   : The polyhedral cone inequalities\n";
      std::cerr << "FileGRP   : The permutation group used for reduction\n";
      std::cerr << "OutFormat : The file format in output, GAP or Oscar\n";
      std::cerr << "FileOUT   : The file of output (if present, otherwise std::cerr)\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      std::cerr << "rational : rational arithmetic on input\n";
      std::cerr << "Qsqrt2   : arithmetic over the field Q(sqrt(2))\n";
      std::cerr << "Qsqrt5   : arithmetic over the field Q(sqrt(5))\n";
      std::cerr << "RealAlgebraic=FileDesc  : For the real algebraic case of a "
                   "field whose description is in FileDesc\n";
      std::cerr << "\n";
      std::cerr << "        --- command ---\n";
      std::cerr << "\n";
      std::cerr << "cdd      : the cdd program\n";
      std::cerr << "lrs      : the lrs program\n";
      std::cerr << "lrs_ring : the lrs program but reduced to ring computation (remove denominators)\n";
      std::cerr << "glrs     : the external program glrs\n";
      std::cerr << "ppl_ext  : the external program ppl_ext\n";
      std::cerr << "cdd_ext  : the external program cdd_ext\n";
      std::cerr << "normaliz : the external program normaliz\n";
      return -1;
    }
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
    using Tint = boost::multiprecision::mpz_int;
    using Trat = boost::multiprecision::mpq_rational;
#else
    using Tint = mpz_class;
    using Trat = mpq_class;
#endif
    //
    using Tidx = int32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt, Tint>;
    std::string arith = argv[1];
    std::string command = argv[2];
    std::string eFileI = argv[3];
    std::string eFileG = argv[4];
    std::string OutFormat = "GAP";
    std::string eFileO = "stderr";
    if (argc == 7) {
      OutFormat = argv[5];
      eFileO = argv[6];
    }
    auto call_dualdesc = [&](std::ostream &os) -> void {
      if (arith == "safe_rational") {
        using T = Rational<SafeInt64>;
        return process<T,Tgroup>(eFileI, eFileG, command, OutFormat, os);
      }
      if (arith == "rational") {
        using T = Trat;
        return process<T,Tgroup>(eFileI, eFileG, command, OutFormat, os);
      }
      if (arith == "Qsqrt5") {
        using T = QuadField<Trat, 5>;
        return process<T,Tgroup>(eFileI, eFileG, command, OutFormat, os);
      }
      if (arith == "Qsqrt2") {
        using T = QuadField<Trat, 2>;
        return process<T,Tgroup>(eFileI, eFileG, command, OutFormat, os);
      }
      std::optional<std::string> opt_realalgebraic =
          get_postfix(arith, "RealAlgebraic=");
      if (opt_realalgebraic) {
        std::string FileAlgebraicField = *opt_realalgebraic;
        if (!IsExistingFile(FileAlgebraicField)) {
          std::cerr << "FileAlgebraicField=" << FileAlgebraicField
                    << " is missing\n";
          throw TerminalException{1};
        }
        HelperClassRealField<Trat> hcrf(FileAlgebraicField);
        int const idx_real_algebraic_field = 1;
        insert_helper_real_algebraic_field(idx_real_algebraic_field, hcrf);
        using T = RealField<idx_real_algebraic_field>;
        return process<T,Tgroup>(eFileI, eFileG, command, OutFormat, os);
      }
      std::cerr << "Failed to find a matching field for arith=" << arith
                << "\n";
      std::cerr << "Available possibilities: rational, Qsqrt5, Qsqrt2, "
                   "RealAlgebraic\n";
      throw TerminalException{1};
    };
    if (eFileO == "stderr") {
      call_dualdesc(std::cerr);
    } else {
      if (eFileO == "stdout") {
        call_dualdesc(std::cout);
      } else {
        std::ofstream os(eFileO);
        call_dualdesc(os);
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_dual_description\n";
    exit(e.eVal);
  }
  runtime(time1);
}
