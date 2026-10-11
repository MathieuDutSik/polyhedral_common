// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#ifdef QUADRATIC_FIELDS
#include "NumberTheoryQuadField.h"
#endif
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "POLY_DirectTriangulation.h"
// clang-format on

// Each simplex comes with the signed determinant of its vertices, taken in
// the order of the simplex (sorted).
template <typename T>
void process(std::string const &eFileI, std::string const &method,
             std::string const &OutFormat, std::ostream &os_out) {
  MyMatrix<T> EXT = ReadMatrixFile<T>(eFileI);
  TriangulationProgram prog = triangulation_program_from_string(method);
  std::vector<std::pair<std::vector<int>, T>> l_trig =
      DirectTriangulationDet<T>(EXT, prog, std::cerr);
  if (OutFormat == "NbTrig") {
    os_out << "NbTrig=" << l_trig.size() << "\n";
    return;
  }
  if (OutFormat == "GAP") {
    os_out << "return [";
    bool is_first = true;
    for (auto &[trig, det] : l_trig) {
      if (!is_first) {
        os_out << ",";
      }
      is_first = false;
      os_out << "rec(LV:=[";
      for (size_t i = 0; i < trig.size(); i++) {
        if (i > 0)
          os_out << ",";
        os_out << (trig[i] + 1);
      }
      os_out << "], det:=" << det << ")";
    }
    os_out << "];\n";
    return;
  }
  std::cerr << "Failed to find a matching output\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    if (argc != 4 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_triangulation_det arith method [DATAIN] [OutFormat] "
                   "[FileOut]\n";
      std::cerr << "\n";
      std::cerr << "POLY_triangulation_det arith method [DATAIN]\n";
      std::cerr << "\n";
      std::cerr << "with:\n";
      std::cerr << "arith     : the chosen arithmetic\n";
      std::cerr << "method    : lrs (reverse search) or beneath_beyond\n";
      std::cerr << "DATAIN    : The polytope vertices\n";
      std::cerr << "OutFormat : NbTrig (default) or GAP, the list of\n";
      std::cerr << "            rec(LV, det) with LV the sorted vertices and\n";
      std::cerr << "            det the signed determinant of the simplex\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      std::cerr << "mpz_class : integer arithmetic over GMP mpz_class\n";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << "mpz_int : integer arithmetic over "
                   "boost::multiprecision::mpz_int\n";
#endif
#ifdef ENABLE_BOOST_TYPES
      std::cerr << "cpp_int : integer arithmetic over "
                   "boost::multiprecision::cpp_int\n";
#endif
      std::cerr << "mpq_class : rational arithmetic over GMP mpq_class\n";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << "mpq_rational : rational arithmetic over "
                   "boost::multiprecision::mpq_rational\n";
#endif
#ifdef ENABLE_BOOST_TYPES
      std::cerr << "cpp_rational : rational arithmetic over "
                   "boost::multiprecision::cpp_rational\n";
#endif
#ifdef QUADRATIC_FIELDS
      std::cerr << "Qsqrt2   : arithmetic over the field Q(sqrt(2))\n";
      std::cerr << "Qsqrt5   : arithmetic over the field Q(sqrt(5))\n";
#endif
      std::cerr
          << "RealAlgebraic=FileDesc  : For the real algebraic case of a ";
      std::cerr << "    field whose description is in FileDesc\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string method = argv[2];
    std::string eFileI = argv[3];
    std::string OutFormat = "NbTrig";
    std::string FileO = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      FileO = argv[5];
    }
    auto call_trig = [&](std::ostream &os) -> void {
      if (arith == "mpz_class") {
        using T = mpz_class;
        return process<T>(eFileI, method, OutFormat, os);
      }
#ifdef ENABLE_BOOST_TYPES
      if (arith == "mpz_int") {
        using T = boost::multiprecision::mpz_int;
        return process<T>(eFileI, method, OutFormat, os);
      }
      if (arith == "cpp_int") {
        using T = boost::multiprecision::cpp_int;
        return process<T>(eFileI, method, OutFormat, os);
      }
#endif
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "fmpz_class") {
        using T = fmpz_class;
        return process<T>(eFileI, method, OutFormat, os);
      }
#endif
      if (arith == "mpq_class") {
        using T = mpq_class;
        return process<T>(eFileI, method, OutFormat, os);
      }
#ifdef ENABLE_BOOST_TYPES
      if (arith == "mpq_rational") {
        using T = boost::multiprecision::mpq_rational;
        return process<T>(eFileI, method, OutFormat, os);
      }
      if (arith == "cpp_rational") {
        using T = boost::multiprecision::cpp_rational;
        return process<T>(eFileI, method, OutFormat, os);
      }
#endif
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "fmpq_class") {
        using T = fmpq_class;
        return process<T>(eFileI, method, OutFormat, os);
      }
#endif
      if (arith == "double") {
        using T = double;
        return process<T>(eFileI, method, OutFormat, os);
      }
#ifdef QUADRATIC_FIELDS
      if (arith == "Qsqrt5") {
        using Trat = Trat_quad_field;
        using T = QuadField<Trat, 5>;
        return process<T>(eFileI, method, OutFormat, os);
      }
      if (arith == "Qsqrt2") {
        using Trat = Trat_quad_field;
        using T = QuadField<Trat, 2>;
        return process<T>(eFileI, method, OutFormat, os);
      }
#endif
      std::optional<std::string> opt_realalgebraic =
          get_postfix(arith, "RealAlgebraic=");
      if (opt_realalgebraic) {
        std::string FileAlgebraicField = *opt_realalgebraic;
        if (!FILE_IsExistingFile(FileAlgebraicField)) {
          std::cerr << "FileAlgebraicField=" << FileAlgebraicField
                    << " is missing\n";
          throw TerminalException{1};
        }
        HelperClassRealField<Trat_real_field> hcrf(FileAlgebraicField);
        int const idx_real_algebraic_field = 1;
        insert_helper_real_algebraic_field(idx_real_algebraic_field, hcrf);
        using T = RealField<idx_real_algebraic_field>;
        return process<T>(eFileI, method, OutFormat, os);
      }
      std::cerr << "Failed to find a matching field for arith=" << arith
                << "\n";
#ifdef QUADRATIC_FIELDS
      std::cerr << "Available possibilities: mpq_class, mpz_class, Qsqrt5, "
                   "Qsqrt2, RealAlgebraic\n";
#else
      std::cerr << "Available possibilities: mpq_class, mpz_class, "
                   "RealAlgebraic\n";
#endif
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(FileO, call_trig);
    std::cerr << "Normal termination of POLY_triangulation_det\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_triangulation_det\n";
    exit(e.eVal);
  }
  runtime(time1);
}
