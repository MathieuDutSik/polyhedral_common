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
#include "Arithmetic_dispatch_types.h"
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
    std::string method;
    std::string eFileI;
    std::string OutFormat = "NbTrig";
    std::string FileO = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T>() -> void {
      process<T>(eFileI, method, OutFormat, *os_ptr);
    };
    add_integer_arithmetics(dispatch, f_arith);
    add_rational_arithmetics(dispatch, f_arith);
    add_algebraic_field_arithmetics(dispatch, f_arith);
    dispatch.add("double", "floating point arithmetic (inexact)",
                 [f_arith]() -> void { f_arith.template operator()<double>(); });
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
      dispatch.print_descriptions(std::cerr);
      return -1;
    }
    //
    std::string arith = argv[1];
    method = argv[2];
    eFileI = argv[3];
    if (argc == 6) {
      OutFormat = argv[4];
      FileO = argv[5];
    }
    auto call_trig = [&](std::ostream &os) -> void {
      os_ptr = &os;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(FileO, call_trig);
    std::cerr << "Normal termination of POLY_triangulation_det\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_triangulation_det\n";
    exit(e.eVal);
  }
  runtime(time1);
}
