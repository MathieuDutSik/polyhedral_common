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
#include "POLY_DualDesc_reverse_search.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void process(std::string const &eFileI, std::string const &choice,
             std::ostream &os) {
  auto read_file = [&]() -> MyMatrix<T> {
    std::optional<std::string> opt = get_prefix(eFileI, "lrs");
    std::ifstream is(eFileI);
    if (opt) {
      return ReadMatrixLrsCdd<T>(is);
    } else {
      return ReadMatrix<T>(is);
    }
  };
  MyMatrix<T> EXT = read_file();
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  //
  if (choice == "lrs") {
    os << "V-representation\n";
    os << "begin\n";
    // The facets are rays: lrs writes them with a leading zero.
    os << "****** " << (nbCol + 1) << " rational\n";
    size_t nVertices = 0;
    auto fPrint = [&]([[maybe_unused]] rev_search::Dictionary<T> &dict,
                      [[maybe_unused]] rev_search::Problem &prob,
                      [[maybe_unused]] int const &col, T *out) -> void {
      os << " 0";
      for (int iCol = 0; iCol < nbCol; iCol++)
        os << " " << out[iCol];
      os << "\n";
      nVertices++;
    };
    rev_search::Kernel_DualDescription(EXT, fPrint);
    os << "end\n";
    os << "*Total: nvertices=" << nVertices << "\n";
    return;
  }
  if (choice == "GAP") {
    os << "return [";
    bool IsFirst = true;
    auto fPrint = [&]([[maybe_unused]] rev_search::Dictionary<T> &dict,
                      [[maybe_unused]] rev_search::Problem &prob,
                      [[maybe_unused]] int const &col, T *out) -> void {
      if (!IsFirst)
        os << ",\n";
      os << "[";
      for (int iCol = 0; iCol < nbCol; iCol++) {
        if (iCol > 0)
          os << ",";
        os << out[iCol];
      }
      os << "]";
      IsFirst = false;
    };
    rev_search::Kernel_DualDescription(EXT, fPrint);
    os << "];\n";
  }
  if (choice == "vertex_incidence") {
    std::vector<size_t> VertexIncd(nbRow, 0);
    T eScal;
    auto fUpdateIncd = [&]([[maybe_unused]] rev_search::Dictionary<T> &dict,
                           [[maybe_unused]] rev_search::Problem &prob,
                           [[maybe_unused]] int const &col, T *out) -> void {
      for (int iRow = 0; iRow < nbRow; iRow++) {
        eScal = 0;
        for (int iCol = 0; iCol < nbCol; iCol++)
          eScal += out[iCol] * EXT(iRow, iCol);
        if (eScal == 0)
          VertexIncd[iRow] += 1;
      }
    };
    rev_search::Kernel_DualDescription(EXT, fUpdateIncd);
    os << "VertexIncd=[";
    for (int iRow = 0; iRow < nbRow; iRow++) {
      if (iRow > 0)
        os << ",";
      os << VertexIncd[iRow];
    }
    os << "]\n";
    return;
  }
  if (choice == "number_facet") {
    size_t nFacets = 0;
    auto fIncrement = [&]([[maybe_unused]] rev_search::Dictionary<T> &dict,
                          [[maybe_unused]] rev_search::Problem &prob,
                          [[maybe_unused]] int const &col,
                          [[maybe_unused]] T *out) -> void { nFacets++; };
    rev_search::Kernel_DualDescription(EXT, fIncrement);
    os << "nFacets=" << nFacets << "\n";
    return;
  }
  if (choice == "qhull_incidence") {
    T eScal;
    auto fPrintIncd = [&]([[maybe_unused]] rev_search::Dictionary<T> &dict,
                          [[maybe_unused]] rev_search::Problem &prob,
                          [[maybe_unused]] int const &col, T *out) -> void {
      bool IsFirst = true;
      for (int iRow = 0; iRow < nbRow; iRow++) {
        eScal = 0;
        for (int iCol = 0; iCol < nbCol; iCol++)
          eScal += out[iCol] * EXT(iRow, iCol);
        if (eScal == 0) {
          if (!IsFirst)
            os << " ";
          IsFirst = false;
          os << iRow;
        }
      }
      os << "\n";
    };
    rev_search::Kernel_DualDescription(EXT, fPrintIncd);
    return;
  }
  if (choice == "structure_vertex_facets") {
    std::vector<std::vector<size_t>> VertexIncd(nbRow);
    std::vector<Face> ListFace;
    size_t idx_facet = 0;
    auto f_insert = [&]([[maybe_unused]] rev_search::Dictionary<T> &dict,
                        [[maybe_unused]] rev_search::Problem &prob,
                        [[maybe_unused]] int const &col, T *out) -> void {
      std::cerr << "idx_facet=" << idx_facet << "\n";
      std::vector<size_t> eIncd;
      Face f(nbRow);
      for (int iRow = 0; iRow < nbRow; iRow++) {
        T eScal(0);
        for (int iCol = 0; iCol < nbCol; iCol++)
          eScal += out[iCol] * EXT(iRow, iCol);
        if (eScal == 0) {
          f[iRow] = 1;
          VertexIncd[iRow].push_back(idx_facet);
        }
      }
      ListFace.push_back(f);
      idx_facet++;
    };
    rev_search::Kernel_DualDescription(EXT, f_insert);
    for (int iRow = 0; iRow < nbRow; iRow++) {
      os << "iRow=" << iRow << " |Contained Facet|=" << VertexIncd[iRow].size()
         << "\n";
      std::map<size_t, size_t> map;
      for (auto &idx_facet : VertexIncd[iRow]) {
        size_t len = ListFace[idx_facet].count();
        map[len]++;
      }
      os << " |map(len(face)) =";
      for (auto &kv : map)
        os << " (" << kv.first << "," << kv.second << ")";
      os << "\n";
    }
    return;
  }
  std::cerr << "Failed to find a matching entry in POLY_lrs\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    std::string choice;
    std::string eFileI;
    std::string eFileO = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T>() -> void {
      process<T>(eFileI, choice, *os_ptr);
    };
    add_integer_arithmetics(dispatch, f_arith);
    add_rational_arithmetics(dispatch, f_arith);
    add_algebraic_field_arithmetics(dispatch, f_arith);
    dispatch.add("double", "floating point arithmetic (inexact)",
                 [f_arith]() -> void { f_arith.template operator()<double>(); });
    if (argc != 4 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_lrs choice arith [DATAIN] [DATAOUT]\n";
      std::cerr << "or\n";
      std::cerr << "POLY_lrs choice arith [DATAIN]\n";
      std::cerr << "\n";
      std::cerr << "with:\n";
      std::cerr << "choice  : the chosen processing option\n";
      std::cerr << "arith   : the chosen arithmetic\n";
      std::cerr << "DATAIN  : The polyhedral cone inequalities\n";
      std::cerr
          << "DATAOUT : The file of output (if present, otherwise std::cerr)\n";
      std::cerr << "\n";
      std::cerr << "        --- choice ---\n";
      std::cerr << "\n";
      std::cerr << "lrs: exact behavior as in lrs\n";
      std::cerr << "vertex_incidence: the incidence of the vertices\n";
      std::cerr << "number_facet: total number of facets\n";
      std::cerr << "qhull_incidence: print the list of incidences\n";
      std::cerr << "structure_vertex_facets: number of facets contained by "
                   "vertices and the number of vertices of such facets\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      return -1;
    }
    //
    choice = argv[1];
    std::string arith = argv[2];
    eFileI = argv[3];
    if (argc == 5)
      eFileO = argv[4];
    auto call_lrs = [&](std::ostream &os) -> void {
      os_ptr = &os;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(eFileO, call_lrs);
    std::cerr << "Normal termination of POLY_lrs\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_lrs\n";
    exit(e.eVal);
  }
  runtime(time1);
}
