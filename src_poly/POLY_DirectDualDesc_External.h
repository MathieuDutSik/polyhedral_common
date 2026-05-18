// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_DIRECTDUALDESC_EXTERNAL_H_
#define SRC_POLY_POLY_DIRECTDUALDESC_EXTERNAL_H_

// This file contains the code that drives the external dual-description
// programs (glrs, ppl_lcdd, lcdd_gmp, normaliz).  Each entry point writes the
// V-representation of EXT to a temporary file, invokes the chosen external
// binary, then parses the produced inequalities back into the requested form
// (incidence, inequality matrix, or face+inequality pair).

// clang-format off
#include "Basic_string.h"
#include "Basic_external_program.h"
#include "MAT_MatrixInt.h"
#include <algorithm>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef TIMINGS
#define TIMINGS_DUAL_DESC_EXTERNAL
#endif

#ifdef DEBUG
#define DEBUG_DUAL_DESC_EXTERNAL
#endif

#ifdef DISABLE_DEBUG_DUAL_DESC_EXTERNAL
#undef DEBUG_DUAL_DESC_EXTERNAL
#endif

template <typename T>
size_t GetShift([[maybe_unused]] MyMatrix<T> const &EXT,
                std::string const &eCommand) {
  size_t shift;
  if (eCommand == "normaliz") {
    shift = 0;
  } else {
    shift = 1;
  }
  return shift;
}

template <typename T, typename Finsert>
void DualDescExternalProgramGeneral(MyMatrix<T> const &EXT, Finsert f_insert,
                                    std::string const &eCommand,
                                    [[maybe_unused]] std::ostream &os) {
#ifdef TIMINGS_DUAL_DESC_EXTERNAL
  MicrosecondTime time;
#endif
  size_t n_row = EXT.rows();
  size_t n_col = EXT.cols();
  size_t DimEXT = n_col + 1;
  auto get_type_str = [&]() -> std::string {
    if (IsIntegralMatrix(EXT)) {
      return "integer";
    }
    return "rational";
  };
  std::string type_str = get_type_str();
  std::string rndStr = random_string(20);
  std::string prefix = "/tmp/";
  std::string suffix = "DD_" + std::to_string(n_row) + "_" +
                       std::to_string(n_col) + "_" + rndStr;
  std::string FileI, FileO, FileE;
  if (eCommand == "normaliz") {
    FileI = prefix + suffix + ".in";
    FileO = prefix + suffix + ".out";
    FileE = prefix + suffix + ".err";
  } else {
    FileI = prefix + suffix + ".ine";
    FileO = prefix + suffix + ".ext";
    FileE = prefix + suffix + ".err";
  }
  auto premature_ending = [&]() -> void {
    std::cerr << "DDD: FileI = " << FileI << "\n";
    std::cerr << "DDD: FileO = " << FileO << "\n";
    std::cerr << "DDD: FileE = " << FileE << "\n";
    if (IsExistingFile(FileE)) {
      for (auto const &line : ReadFullFile(FileE)) {
        std::cerr << line << "\n";
      }
    }
    throw TerminalException{1};
  };
  {
    std::ofstream osI(FileI);
    if (eCommand == "normaliz") {
      osI << "amb_space " << n_col << "\n";
      osI << "cone " << n_row << "\n";
      for (size_t i_row = 0; i_row < n_row; i_row++) {
        for (size_t i_col = 0; i_col < n_col; i_col++)
          osI << " " << EXT(i_row, i_col);
        osI << "\n";
      }
      osI << "SupportHyperplanes\n";
    } else {
      osI << "V-representation\n";
      osI << "begin\n";
      osI << n_row << " " << DimEXT << " " << type_str << "\n";
      for (size_t i_row = 0; i_row < n_row; i_row++) {
        osI << "0";
        for (size_t i_col = 0; i_col < n_col; i_col++)
          osI << " " << EXT(i_row, i_col);
        osI << "\n";
      }
      osI << "end\n";
    }
  }
#ifdef TIMINGS_DUAL_DESC_EXTERNAL
  os << "|DDD: FileWriting|=" << time << "\n";
#endif
#ifdef DEBUG_DUAL_DESC_EXTERNAL
  os << "DDD: FileO=" << FileO << " created\n";
#endif
  //
  // Now calling the external program
  //
#ifdef DEBUG_DUAL_DESC_EXTERNAL
  os << "DDD: program=" << eCommand << " FileI=" << FileI << " FileO=" << FileO
     << " FileE=" << FileE << "\n";
#endif
  int iret1;
  if (eCommand == "normaliz") {
    iret1 = RunExternalProgram(eCommand, {FileI}, std::nullopt, std::nullopt,
                               FileE);
  } else {
    iret1 = RunExternalProgram(eCommand, {FileI}, std::nullopt, FileO, FileE);
  }
#ifdef TIMINGS_DUAL_DESC_EXTERNAL
  os << "|DDD: glrs/ppl/cdd|=" << time << "\n";
#endif
#ifdef DEBUG_DUAL_DESC_EXTERNAL
  os << "DDD: External program terminated\n";
#endif
  if (iret1 != 0) {
    std::cerr << "DDD: iret1=" << iret1 << "\n";
    std::cerr << "DDD: The program has not terminated correctly\n";
    premature_ending();
  }
  size_t n_facet = 0;
  size_t n_insert = 0;
  auto check_consistency = [&]() -> void {
    if (n_insert != n_facet) {
      std::cerr << "DDD: Consistency error\n";
      std::cerr << "DDD: n_insert=" << n_insert << " n_facet=" << n_facet
                << "\n";
      premature_ending();
    }
  };
  //
  std::ifstream is(FileO);
  std::vector<T> LVal(DimEXT);
  T eScal;
  std::string line;
  size_t pos_wrt = 0;
  auto f_read = [&](const std::string &str) -> void {
    ParseScalar_inplace<T>(str, LVal[pos_wrt]);
    pos_wrt++;
  };
  auto process_line = [&]() -> void {
    STRING_Split_f(line, " ", f_read);
    pos_wrt = 0;
    f_insert(LVal);
    n_insert++;
  };
  size_t iLine = 0;
  // For this case we know at the beginning the number of inequalities
  if (eCommand == "normaliz") {
    bool has_n_facet = false;
    size_t iLineLimit = 0;
    while (std::getline(is, line)) {
#ifdef DEBUG_DUAL_DESC_EXTERNAL
      os << "DDD: iLine=" << iLine << " line=" << line << "\n";
#endif
      if (has_n_facet) {
        if (iLine < iLineLimit) {
          process_line();
        } else {
          break;
        }
      }
      if (!has_n_facet) {
        std::vector<std::string> LStr =
            STRING_Split(line, " support hyperplanes");
        if (LStr.size() == 2) {
          has_n_facet = true;
          n_facet = ParseScalar<size_t>(LStr[0]);
          iLineLimit = iLine + n_facet + 1;
        }
      }
      iLine++;
    }
    check_consistency();
  }
  if (eCommand == "ppl_lcdd" || eCommand == "lcdd_gmp") {
    size_t headersize, iLineLimit = 0;
    if (eCommand == "lcdd_gmp")
      headersize = 4;
    if (eCommand == "ppl_lcdd")
      headersize = 3;
    while (std::getline(is, line)) {
      if (iLine == headersize - 1) {
        // Determining the number of entries
        std::vector<std::string> LStr = STRING_Split(line, " ");
        n_facet = ParseScalar<size_t>(LStr[0]);
        iLineLimit = headersize + n_facet;
      }
      if (iLine >= headersize && (iLineLimit == 0 || iLine < iLineLimit))
        process_line();
      iLine++;
    }
    check_consistency();
  }
  if (eCommand == "glrs") {
    size_t headersize = 0;
    while (std::getline(is, line)) {
      if (line == "end")
        break;
      if (line == "begin")
        headersize = iLine + 2;
      if (headersize != 0 and iLine >= headersize)
        process_line();
      iLine++;
    }
  }
#ifdef TIMINGS_DUAL_DESC_EXTERNAL
  os << "|DDD: FileRead|=" << time << "\n";
#endif
#ifdef DEBUG_DUAL_DESC_EXTERNAL
  os << "DDD: FileI = " << FileI << "    FileO = " << FileO << "\n";
#endif
  if (n_insert == 0) {
    std::cerr << "DDD: We inserted zero entries\n";
    premature_ending();
  }
  RemoveFileIfExist(FileI);
  RemoveFileIfExist(FileO);
  RemoveFileIfExist(FileE);
}

template <typename T>
vectface DualDescExternalProgramIncidence(MyMatrix<T> const &EXT,
                                          std::string const &eCommand,
                                          std::ostream &os) {
  size_t n_row = EXT.rows();
  size_t n_col = EXT.cols();
  size_t DimEXT = n_col + 1;
  vectface vf(n_row);
  Face f(n_row);
  size_t shift = GetShift(EXT, eCommand);
  T eScal;
  auto f_insert = [&](std::vector<T> const &LVal) -> void {
    for (size_t i_row = 0; i_row < n_row; i_row++) {
      eScal = 0;
      for (size_t i = shift; i < DimEXT; i++)
        eScal += LVal[i] * EXT(i_row, i - shift);
      f[i_row] = static_cast<bool>(eScal == 0);
    }
    vf.push_back(f);
  };
  DualDescExternalProgramGeneral(EXT, f_insert, eCommand, os);
  return vf;
}

template <typename T, typename Tint>
std::vector<Tint> RescaleVec(std::vector<T> const &v) {
  static_assert(is_implementation_of_Q<T>::value,
                "Requires T to be an implementation of Q");
  static_assert(is_implementation_of_Q<Tint>::value ||
                    is_implementation_of_Z<Tint>::value,
                "Requires Tint to be an implementation of Z or Q");
  int cols = v.size();
  std::vector<Tint> dens(cols, 1);
  std::vector<Tint> vret = std::vector<Tint>(cols);
  for (int iCol = 0; iCol < cols; iCol++) {
    dens[iCol] = v[iCol].get_den();
  }
  Tint scale = LCMlist(dens);
  for (int iCol = 0; iCol < cols; iCol++) {
    vret[iCol] = (scale / v[iCol].get_den()) * v[iCol].get_num();
  }
  return vret;
}

template <>
vectface DualDescExternalProgramIncidence(MyMatrix<mpq_class> const &EXT,
                                          std::string const &eCommand,
                                          std::ostream &os) {
  using T = mpq_class;
  using Tint = mpz_class;
  using Tlong = int64_t;
  size_t n_row = EXT.rows();
  size_t n_col = EXT.cols();
  mpz_class n_col_mpz(n_col);
  size_t DimEXT = n_col + 1;
  size_t shift = GetShift(EXT, eCommand);
  Tint eScal;
  Tlong eScal_Tlong;
  vectface vf(n_row);
  Face f(n_row);
  MyMatrix<Tint> EXT_int = RescaleRows<T, Tint>(EXT);
  MyMatrix<Tlong> EXT_Tlong = MyMatrix<Tlong>(n_row, n_col);
  size_t max_bits = 0;
  auto f_bit = [&](mpz_class const &val) -> size_t {
    return mpz_sizeinbase(val.get_mpz_t(), 2);
  };
  for (size_t iRow = 0; iRow < n_row; iRow++) {
    for (size_t iCol = 0; iCol < n_col; iCol++) {
      max_bits = std::max(f_bit(EXT_int(iRow, iCol)), max_bits);
      EXT_Tlong(iRow, iCol) = EXT_int(iRow, iCol).get_si();
    }
  }
  max_bits += f_bit(n_col_mpz);
  std::vector<Tlong> LVal_Tlong(DimEXT, 0);

  auto f_insert = [&](std::vector<T> const &LVal) -> void {
    std::vector<Tint> LVal_int = RescaleVec<T, Tint>(LVal);
    size_t max_bits_LVal = 0;
    for (size_t i = shift; i < DimEXT; i++) {
      max_bits_LVal = std::max(max_bits_LVal, f_bit(LVal_int[i]));
      LVal_Tlong[i] = LVal_int[i].get_si();
    }

    if (max_bits + max_bits_LVal <= 60) {
      // safe to use long
      for (size_t i_row = 0; i_row < n_row; i_row++) {
        eScal_Tlong = 0;
        for (size_t i = shift; i < DimEXT; i++)
          eScal_Tlong += LVal_Tlong[i] * EXT_Tlong(i_row, i - shift);
        f[i_row] = static_cast<bool>(eScal_Tlong == 0);
      }
    } else {
      for (size_t i_row = 0; i_row < n_row; i_row++) {
        eScal = 0;
        for (size_t i = shift; i < DimEXT; i++)
          eScal += LVal_int[i] * EXT_int(i_row, i - shift);
        f[i_row] = static_cast<bool>(eScal == 0);
      }
    }
    vf.push_back(f);
  };
  DualDescExternalProgramGeneral(EXT, f_insert, eCommand, os);
  return vf;
}

template <typename T>
MyMatrix<T> DualDescExternalProgramIneq(MyMatrix<T> const &EXT,
                                        std::string const &eCommand,
                                        std::ostream &os) {
  size_t n_col = EXT.cols();
  size_t DimEXT = n_col + 1;
  size_t shift = GetShift(EXT, eCommand);
  int nbColRed = DimEXT - shift;
  MyVector<T> V(nbColRed);
  std::vector<MyVector<T>> ListVect;
  auto f_insert = [&](std::vector<T> const &LVal) -> void {
    for (int i = 0; i < nbColRed; i++) {
      V(i) = LVal[i + shift];
    }
    ListVect.push_back(V);
  };
  DualDescExternalProgramGeneral(EXT, f_insert, eCommand, os);
  return MatrixFromVectorFamily(ListVect);
}

template <typename T, typename Fprocess>
void DualDescExternalProgramFaceIneq(MyMatrix<T> const &EXT,
                                     std::string const &eCommand,
                                     Fprocess f_process, std::ostream &os) {
  size_t n_row = EXT.rows();
  size_t n_col = EXT.cols();
  size_t DimEXT = n_col + 1;
  size_t shift = GetShift(EXT, eCommand);
  int nbColRed = DimEXT - shift;
  std::pair<Face, MyVector<T>> pair{Face(n_row), MyVector<T>(nbColRed)};
  T eScal;
  auto f_insert = [&](std::vector<T> const &LVal) -> void {
    for (int i = 0; i < nbColRed; i++) {
      pair.second(i) = LVal[i + shift];
    }
    for (size_t i_row = 0; i_row < n_row; i_row++) {
      eScal = 0;
      for (size_t i = shift; i < DimEXT; i++)
        eScal += LVal[i] * EXT(i_row, i - shift);
      pair.first[i_row] = static_cast<bool>(eScal == 0);
    }
    f_process(pair);
  };
  DualDescExternalProgramGeneral(EXT, f_insert, eCommand, os);
}

// clang-format off
#endif  // SRC_POLY_POLY_DIRECTDUALDESC_EXTERNAL_H_
// clang-format on
