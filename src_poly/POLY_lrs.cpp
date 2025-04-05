// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheoryQuadField.h"
#include "NumberTheorySafeInt.h"
#include "POLY_lrslib.h"
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
  MyMatrix<T> EXT_pre = read_file();
  std::pair<MyMatrix<T>, int> pair = lrs::FirstColumnZeroCond(EXT_pre);
  MyMatrix<T> const &EXT = pair.first;
  int shift = pair.second;
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  //
  if (choice == "lrs") {
    os << "V-representation\n";
    os << "begin\n";
    os << "****** " << nbCol << " rational\n";
    size_t nVertices = 0;
    auto fPrint = [&]([[maybe_unused]] lrs::lrs_dic<T> *P,
                      [[maybe_unused]] lrs::lrs_dat<T> *Q,
                      [[maybe_unused]] int const &col, T *out) -> void {
      for (int iCol = 0; iCol < nbCol; iCol++)
        os << " " << out[iCol];
      os << "\n";
      nVertices++;
    };
    lrs::Kernel_DualDescription(EXT, fPrint);
    os << "end\n";
    os << "*Total: nvertices=" << nVertices << "\n";
    return;
  }
  if (choice == "GAP") {
    os << "return [";
    bool IsFirst = true;
    auto fPrint = [&]([[maybe_unused]] lrs::lrs_dic<T> *P,
                      [[maybe_unused]] lrs::lrs_dat<T> *Q,
                      [[maybe_unused]] int const &col, T *out) -> void {
      if (!IsFirst)
        os << ",\n";
      os << "[";
      for (int iCol = shift; iCol < nbCol; iCol++) {
        if (iCol > shift)
          os << ",";
        os << out[iCol];
      }
      os << "]";
      IsFirst = false;
    };
    lrs::Kernel_DualDescription(EXT, fPrint);
    os << "];\n";
  }
  if (choice == "vertex_incidence") {
    std::vector<size_t> VertexIncd(nbRow, 0);
    T eScal;
    auto fUpdateIncd = [&]([[maybe_unused]] lrs::lrs_dic<T> *P,
                           [[maybe_unused]] lrs::lrs_dat<T> *Q,
                           [[maybe_unused]] int const &col, T *out) -> void {
      for (int iRow = 0; iRow < nbRow; iRow++) {
        eScal = 0;
        for (int iCol = 0; iCol < nbCol; iCol++)
          eScal += out[iCol] * EXT(iRow, iCol);
        if (eScal == 0)
          VertexIncd[iRow] += 1;
      }
    };
    lrs::Kernel_DualDescription(EXT, fUpdateIncd);
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
    auto fIncrement = [&]([[maybe_unused]] lrs::lrs_dic<T> *P,
                          [[maybe_unused]] lrs::lrs_dat<T> *Q,
                          [[maybe_unused]] int const &col,
                          [[maybe_unused]] T *out) -> void { nFacets++; };
    lrs::Kernel_DualDescription(EXT, fIncrement);
    os << "nFacets=" << nFacets << "\n";
    return;
  }
  if (choice == "qhull_incidence") {
    T eScal;
    auto fPrintIncd = [&]([[maybe_unused]] lrs::lrs_dic<T> *P,
                          [[maybe_unused]] lrs::lrs_dat<T> *Q,
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
    lrs::Kernel_DualDescription(EXT, fPrintIncd);
    return;
  }
  if (choice == "structure_vertex_facets") {
    std::vector<std::vector<size_t>> VertexIncd(nbRow);
    std::vector<Face> ListFace;
    size_t idx_facet = 0;
    auto f_insert = [&]([[maybe_unused]] lrs::lrs_dic<T> *P,
                        [[maybe_unused]] lrs::lrs_dat<T> *Q,
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
    lrs::Kernel_DualDescription(EXT, f_insert);
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
  HumanTime time1;
  try {
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
      std::cerr
          << "safe_integer  : integer arithmetic based on int64_t that fails\n";
      std::cerr << "    gracefully on overflow\n";
      std::cerr << "safe_rational : rational arithmetic based on int64_t that "
                   "fails\n";
      std::cerr << "    gracefully on overflow\n";
      std::cerr << "integer  : integer arithmetic on input\n";
      std::cerr << "rational : rational arithmetic on input\n";
      std::cerr << "Qsqrt2   : arithmetic over the field Q(sqrt(2))\n";
      std::cerr << "Qsqrt5   : arithmetic over the field Q(sqrt(5))\n";
      std::cerr
          << "RealAlgebraic=FileDesc  : For the real algebraic case of a ";
      std::cerr << "    field whose description is in FileDesc\n";
      return -1;
    }
    //
    std::string choice = argv[1];
    std::string arith = argv[2];
    std::string eFileI = argv[3];
    std::string eFileO = "stderr";
    if (argc == 5)
      eFileO = argv[4];
    auto call_lrs = [&](std::ostream &os) -> void {
      if (arith == "safe_integer") {
        using T = SafeInt64;
        return process<T>(eFileI, choice, os);
      }
      if (arith == "safe_rational") {
        using T = Rational<SafeInt64>;
        return process<T>(eFileI, choice, os);
      }
      if (arith == "integer") {
        using T = mpz_class;
        return process<T>(eFileI, choice, os);
      }
      if (arith == "rational") {
        using T = mpq_class;
        return process<T>(eFileI, choice, os);
      }
      if (arith == "double") {
        using T = double;
        return process<T>(eFileI, choice, os);
      }
      if (arith == "Qsqrt5") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 5>;
        return process<T>(eFileI, choice, os);
      }
      if (arith == "Qsqrt2") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 2>;
        return process<T>(eFileI, choice, os);
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
        return process<T>(eFileI, choice, os);
      }
      std::cerr << "Failed to find a matching field for arith=" << arith
                << "\n";
      std::cerr << "Available possibilities: rational, Qsqrt5, Qsqrt2, "
                   "RealAlgebraic\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(eFileO, call_lrs);
    std::cerr << "Normal termination of POLY_lrs\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_lrs\n";
    exit(e.eVal);
  }
  runtime(time1);
}
