// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheoryQuadField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "POLY_lrslib.h"
// clang-format on

template <typename T>
void process(std::string const &eFileI, std::string const &OutFormat,
             std::ostream &os_out) {
  MyMatrix<T> EXT = ReadMatrixFile<T>(eFileI);
  vectface vf = lrs::GetTriangulation(EXT);
  if (OutFormat == "Volume") {
    T sum_det(0);
    for (auto &trig : vf) {
      MyMatrix<T> EXTtrig = SelectRow(EXT, trig);
      T det = T_abs(DeterminantMat(EXTtrig));
      std::cerr << "det=" << det << "\n";
      sum_det += T_abs(det);
    }
    int dim = EXT.cols() - 1;
    T fact(1);
    for (int u = 1; u <= dim; u++) {
      fact *= u;
    }
    os_out << "sum_det=" << sum_det << " fact=" << fact << "\n";
    T volume = sum_det / fact;
    os_out << "volume=" << volume << "\n";
    return;
  }
  if (OutFormat == "NbTrig") {
    os_out << "NbTrig=" << vf.size() << "\n";
    return;
  }
  if (OutFormat == "Trigs") {
    os_out << "return " << StringVectfaceGAP(vf) << ";\n";
    return;
  }
  if (OutFormat == "GAP") {
    os_out << "return ";
    VectVectInt_Gap_Print(os_out, vf);
    os_out << ";\n";
    return;
  }
  if (OutFormat == "GAPtrig_det") {
    os_out << "return [";
    bool is_first = true;
    for (auto &trig : vf) {
      MyMatrix<T> EXTtrig = SelectRow(EXT, trig);
      T det = T_abs(DeterminantMat(EXTtrig));
      if (!is_first) {
        os_out << ",";
      }
      is_first = false;
      size_t siz = trig.count();
      os_out << "rec(LV:=[";
      boost::dynamic_bitset<>::size_type eVal = trig.find_first();
      for (size_t i = 0; i < siz; i++) {
        if (i > 0)
          os_out << ",";
        os_out << int(eVal + 1);
        eVal = trig.find_next(eVal);
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
  HumanTime time1;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr
          << "POLY_lrs_triangulation arith [DATAIN] [OutFormat] [FileOut]\n";
      std::cerr << "\n";
      std::cerr << "POLY_lrs_triangulation arith [DATAIN]\n";
      std::cerr << "\n";
      std::cerr << "with:\n";
      std::cerr << "arith   : the chosen arithmetic\n";
      std::cerr << "DATAIN  : The polytope vertices\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      std::cerr
          << "safe_integer  : integer arithmetic based on int64_t that fails\n";
      std::cerr << "    gracefully on overflow\n";
      std::cerr << "safe_rational : rational arithmetic based on int64_t that "
                   "fails\n";
      std::cerr << "    gracefully on overflow\n";
      std::cerr << "mpz_class : integer arithmetic over GMP mpz_class\n";
      std::cerr << "mpz_int : integer arithmetic over "
                   "boost::multiprecision::mpz_int\n";
      std::cerr << "cpp_int : integer arithmetic over "
                   "boost::multiprecision::cpp_int\n";
      std::cerr << "mpq_class : rational arithmetic over GMP mpq_class\n";
      std::cerr << "mpq_rational : rational arithmetic over "
                   "boost::multiprecision::mpq_rational\n";
      std::cerr << "cpp_rational : rational arithmetic over "
                   "boost::multiprecision::cpp_rational\n";
      std::cerr << "Qsqrt2   : arithmetic over the field Q(sqrt(2))\n";
      std::cerr << "Qsqrt5   : arithmetic over the field Q(sqrt(5))\n";
      std::cerr
          << "RealAlgebraic=FileDesc  : For the real algebraic case of a ";
      std::cerr << "    field whose description is in FileDesc\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string eFileI = argv[2];
    std::string OutFormat = "NbTrig";
    std::string FileO = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      FileO = argv[4];
    }
    auto call_lrs = [&](std::ostream &os) -> void {
      if (arith == "safe_integer") {
        using T = SafeInt64;
        return process<T>(eFileI, OutFormat, os);
      }
      /*
      if (arith == "safe_rational") {
        using T = Rational<SafeInt64>;
        return process<T>(eFileI, OutFormat, os);
      }
      */
      if (arith == "mpz_class") {
        using T = mpz_class;
        return process<T>(eFileI, OutFormat, os);
      }
      if (arith == "mpz_int") {
        using T = boost::multiprecision::mpz_int;
        return process<T>(eFileI, OutFormat, os);
      }
      if (arith == "cpp_int") {
        using T = boost::multiprecision::cpp_int;
        return process<T>(eFileI, OutFormat, os);
      }
      if (arith == "mpq_class") {
        using T = mpq_class;
        return process<T>(eFileI, OutFormat, os);
      }
      if (arith == "mpq_rational") {
        using T = boost::multiprecision::mpq_rational;
        return process<T>(eFileI, OutFormat, os);
      }
      if (arith == "cpp_rational") {
        using T = boost::multiprecision::cpp_rational;
        return process<T>(eFileI, OutFormat, os);
      }
      if (arith == "double") {
        using T = double;
        return process<T>(eFileI, OutFormat, os);
      }
      if (arith == "Qsqrt5") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 5>;
        return process<T>(eFileI, OutFormat, os);
      }
      if (arith == "Qsqrt2") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 2>;
        return process<T>(eFileI, OutFormat, os);
      }
      std::optional<std::string> opt_realalgebraic =
          get_postfix(arith, "RealAlgebraic=");
      if (opt_realalgebraic) {
        using T_rat = mpq_class;
        std::string FileAlgebraicField = *opt_realalgebraic;
        if (!IsExistingFile(FileAlgebraicField)) {
          std::cerr << "FileAlgebraicField=" << FileAlgebraicField
                    << " is missing\n";
          throw TerminalException{1};
        }
        HelperClassRealField<T_rat> hcrf(FileAlgebraicField);
        int const idx_real_algebraic_field = 1;
        insert_helper_real_algebraic_field(idx_real_algebraic_field, hcrf);
        using T = RealField<idx_real_algebraic_field>;
        return process<T>(eFileI, OutFormat, os);
      }
      std::cerr << "Failed to find a matching field for arith=" << arith
                << "\n";
      std::cerr << "Available possibilities: mpq_class, mpz_class, Qsqrt5, "
                   "Qsqrt2, RealAlgebraic\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(FileO, call_lrs);
    std::cerr << "Normal termination of POLY_lrs_triangulation\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_lrs_triangulation\n";
    exit(e.eVal);
  }
  runtime(time1);
}
