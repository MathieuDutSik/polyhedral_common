// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_DIRECTDUALDESC_H_
#define SRC_POLY_POLY_DIRECTDUALDESC_H_

// clang-format off
#include "Basic_string.h"
#include "POLY_c_cddlib.h"
#include "POLY_cddlib.h"
#include "POLY_lrslib.h"
#include "MAT_MatrixInt.h"
#include "POLY_DualDescription_PrimalDual.h"
#include "SmallPolytopes.h"
#include <string>
#include <vector>
// clang-format on

#ifdef TIMINGS
#define TIMINGS_DUAL_DESC
#endif

#ifdef DEBUG
#define DEBUG_DUAL_DESC
#endif

#ifdef KEY_VALUE
#define KEY_VALUE_DUAL_DESC
#endif

template <typename T> std::vector<size_t> Convert_T_To_Set(T const &val) {
  size_t pos = 0;
  std::vector<size_t> V;
  T valWork = val;
  T two = 2;
  while (true) {
    if (valWork == 0)
      break;
    T res = ResInt(valWork, two);
    if (res == 1) {
      V.push_back(pos);
    }
    valWork = (valWork - res) / two;
    pos++;
  }
  return V;
}

template <typename T> T Convert_Set_To_T(std::vector<size_t> const &V) {
  T ThePow = 1;
  size_t pos = 0;
  T retval = 0;
  for (auto &eV : V) {
    while (true) {
      if (pos == eV)
        break;
      ThePow *= 2;
      pos++;
    }
    retval += ThePow;
  }
  return retval;
}

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
                                    std::ostream &os) {
#ifdef TIMINGS_DUAL_DESC
  MicrosecondTime time;
#endif
  size_t n_row = EXT.rows();
  size_t n_col = EXT.cols();
  size_t DimEXT = n_col + 1;
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
      osI << n_row << " " << DimEXT << " integer\n";
      for (size_t i_row = 0; i_row < n_row; i_row++) {
        osI << "0";
        for (size_t i_col = 0; i_col < n_col; i_col++)
          osI << " " << EXT(i_row, i_col);
        osI << "\n";
      }
      osI << "end\n";
    }
  }
#ifdef TIMINGS_DUAL_DESC
  os << "|FileWriting|=" << time << "\n";
#endif
  //  os << "FileO=" << FileO << " created\n";
  //
  // Now calling the external program
  //
  std::string order;
  if (eCommand == "normaliz") {
    order = eCommand + " " + FileI;
  } else {
    order = eCommand + " " + FileI + " > " + FileO + " 2> " + FileE;
  }
  os << "order=" << order << "\n";
  int iret1 = system(order.c_str());
#ifdef TIMINGS_DUAL_DESC
  os << "|glrs/ppl/cdd|=" << time << "\n";
#endif
  os << "External program terminated\n";
  if (iret1 != 0) {
    os << "The program has not terminated correctly\n";
    os << "FileI = " << FileI << "\n";
    os << "FileO = " << FileO << "\n";
    os << "FileE = " << FileE << "\n";
    throw TerminalException{1};
  }
  size_t n_facet = 0;
  size_t n_insert = 0;
  auto check_consistency = [&]() -> void {
    if (n_insert != n_facet) {
      os << "Consistency error\n";
      os << "n_insert=" << n_insert << " n_facet=" << n_facet << "\n";
      os << "FileI = " << FileI << "\n";
      os << "FileO = " << FileO << "\n";
      os << "FileE = " << FileE << "\n";
      throw TerminalException{1};
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
      os << "iLine=" << iLine << " line=" << line << "\n";
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
#ifdef TIMINGS_DUAL_DESC
  os << "|FileRead|=" << time << "\n";
#endif
  os << "FileI = " << FileI << "    FileO = " << FileO << "\n";
  if (n_insert == 0) {
    std::cerr << "We inserted zero entries\n";
    os << "FileI = " << FileI << "\n";
    os << "FileO = " << FileO << "\n";
    os << "FileE = " << FileE << "\n";
    throw TerminalException{1};
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
#ifdef USE_ISINCD
    auto isincd = [&](size_t i_row) -> bool {
      eScal = 0;
      for (size_t i = shift; i < DimEXT; i++)
        eScal += LVal[i] * EXT(i_row, i - shift);
      return eScal == 0;
    };
    vf.InsertFaceRef(isincd);
#else
    for (size_t i_row = 0; i_row < n_row; i_row++) {
      eScal = 0;
      for (size_t i = shift; i < DimEXT; i++)
        eScal += LVal[i] * EXT(i_row, i - shift);
      f[i_row] = static_cast<bool>(eScal == 0);
    }
    vf.push_back(f);
#endif
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
  size_t n_row = EXT.rows();
  size_t n_col = EXT.cols();
  mpz_class n_col_mpz(n_col);
  size_t DimEXT = n_col + 1;
  size_t shift = GetShift(EXT, eCommand);
  Tint eScal;
  long eScal_long;
  vectface vf(n_row);
  Face f(n_row);
  MyMatrix<Tint> EXT_int = RescaleRows<T, Tint>(EXT);
  MyMatrix<long> EXT_long = MyMatrix<long>(n_row, n_col);
  size_t max_bits = 0;
  auto f_bit = [&](mpz_class const &val) -> size_t {
    return mpz_sizeinbase(val.get_mpz_t(), 2);
  };
  for (size_t iRow = 0; iRow < n_row; iRow++) {
    for (size_t iCol = 0; iCol < n_col; iCol++) {
      max_bits = std::max(f_bit(EXT_int(iRow, iCol)), max_bits);
      EXT_long(iRow, iCol) = EXT_int(iRow, iCol).get_si();
    }
  }
  max_bits += f_bit(n_col_mpz);
  std::vector<long> LVal_long(DimEXT, 0);

  auto f_insert = [&](std::vector<T> const &LVal) -> void {
    std::vector<Tint> LVal_int = RescaleVec<T, Tint>(LVal);
    size_t max_bits_LVal = 0;
    for (size_t i = shift; i < DimEXT; i++) {
      max_bits_LVal = std::max(max_bits_LVal, f_bit(LVal_int[i]));
      LVal_long[i] = LVal_int[i].get_si();
    }

    if (max_bits + max_bits_LVal <= 60) {
      // safe to use long
      for (size_t i_row = 0; i_row < n_row; i_row++) {
        eScal_long = 0;
        for (size_t i = shift; i < DimEXT; i++)
          eScal_long += LVal_long[i] * EXT_long(i_row, i - shift);
        f[i_row] = static_cast<bool>(eScal_long == 0);
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

template <typename T>
vectface DirectFacetComputationIncidence(MyMatrix<T> const &EXT,
                                         std::string const &ansProg,
                                         std::ostream &os) {
  std::string eProg;
  std::vector<std::string> ListProg;
  //
  eProg = "cdd_cbased";
  ListProg.push_back(eProg);
  if (ansProg == eProg) {
#ifdef USE_CDDLIB
    return cbased_cdd::DualDescription_incd(EXT);
#else
    std::cerr << "The code has been compiled without the CDDLIB library\n";
    throw TerminalException{1};
#endif
  }
  //
  if constexpr (is_ring_field<T>::value) {
    // CDD certainly requires the ring to be a field
    eProg = "cdd";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return cdd::DualDescription_incd(EXT);
    // If it is a field, then it makes sense to look at the internal ring
    eProg = "lrs_ring";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return lrs::DualDescription_incd_reduction(EXT);
  }
  //
  if constexpr (is_ring_field<T>::value &&
                has_reduction_subset_solver<T>::value) {
    // Small polytopes have special solutions
    eProg = "small_polytopes";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return SmallPolytope_Incidence(EXT);
  }
  // It applies to the field case or ring
  eProg = "lrs";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return lrs::DualDescription_incd(EXT);
  // It applies to the field case or ring
  eProg = "pd_lrs";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return POLY_DualDescription_PrimalDualIncidence(EXT, os);
  //
  //
  // The external programs are available only for rationl types
  //
  if constexpr (is_implementation_of_Q<T>::value) {
    eProg = "glrs";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return DualDescExternalProgramIncidence(EXT, "glrs", os);
    //
    eProg = "ppl_ext";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return DualDescExternalProgramIncidence(EXT, "ppl_lcdd", os);
    //
    eProg = "cdd_ext";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return DualDescExternalProgramIncidence(EXT, "lcdd_gmp", os);
    //
    eProg = "normaliz";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return DualDescExternalProgramIncidence(EXT, "normaliz", os);
  }
  //
  std::cerr << "ERROR in DirectFacetComputationIncidence\n";
  std::cerr << "No right program found with ansProg=" << ansProg
            << " or incorrect output\n";
  std::cerr << "ERROR: No right program found with ansProg=" << ansProg
            << " or incorrect output\n";
  std::cerr << "List of authorized programs :";
  bool IsFirst = true;
  for (auto &eP : ListProg) {
    if (!IsFirst)
      std::cerr << " ,";
    IsFirst = false;
    std::cerr << " " << eP;
  }
  std::cerr << "\n";
  throw TerminalException{1};
}

template <typename T>
MyMatrix<T> DirectFacetComputationInequalities(MyMatrix<T> const &EXT,
                                               std::string const &ansProg,
                                               std::ostream &os) {
  std::string eProg;
  std::vector<std::string> ListProg;
  //
  if constexpr (is_ring_field<T>::value) {
    // CDD certainly requires a field for working out.
    eProg = "cdd";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return cdd::DualDescription(EXT);
    // For lrs_ring, we certainly need to have a field for T
    eProg = "lrs_ring";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return lrs::DualDescription_reduction(EXT);
  }
  //
  if constexpr (is_ring_field<T>::value &&
                has_reduction_subset_solver<T>::value) {
    // Small polytopes have special solutions
    eProg = "small_polytopes";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return SmallPolytope_Ineq(EXT);
  }
  // lrs does not use divisions, so work even if not field.
  eProg = "lrs";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return lrs::DualDescription(EXT);
  // It applies to the field case or ring
  eProg = "pd_lrs";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return POLY_DualDescription_PrimalDualInequalities(EXT, os);
  //
  // The external programs are available only for rationl types
  //
  if constexpr (is_implementation_of_Q<T>::value) {
    eProg = "glrs";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return DualDescExternalProgramIneq(EXT, "glrs", os);
    //
    eProg = "ppl_ext";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return DualDescExternalProgramIneq(EXT, "ppl_lcdd", os);
    //
    eProg = "cdd_ext";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return DualDescExternalProgramIneq(EXT, "lcdd_gmp", os);
    //
    eProg = "normaliz";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return DualDescExternalProgramIneq(EXT, "normaliz", os);
  }
  //
  std::cerr << "ERROR in DirectFacetComputationInequalities\n";
  std::cerr << "No right program found with ansProg=" << ansProg
            << " or incorrect output\n";
  std::cerr << "List of authorized programs :";
  bool IsFirst = true;
  for (auto &eP : ListProg) {
    if (!IsFirst)
      std::cerr << " ,";
    IsFirst = false;
    std::cerr << " " << eP;
  }
  std::cerr << "\n";
  throw TerminalException{1};
}

template <typename T, typename Fprocess>
void DirectFacetComputationFaceIneq(MyMatrix<T> const &EXT,
                                    std::string const &ansProg,
                                    Fprocess f_process, std::ostream &os) {
  std::string eProg;
  std::vector<std::string> ListProg;
  //
  if constexpr (is_ring_field<T>::value) {
    // CDD requires for T to be a field
    eProg = "cdd";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return cdd::DualDescriptionFaceIneq(EXT, f_process);
    // For lrs_ring that computes in a subring, we need T to be a field.
    eProg = "lrs_ring";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return lrs::DualDescriptionFaceIneq_reduction(EXT, f_process);
  }
  // We need to make it work also for ase without reduction if that makes sense
  // which is not sure at all.
  if constexpr (is_ring_field<T>::value &&
                has_reduction_subset_solver<T>::value) {
    // Small polytopes can have special solutions
    eProg = "small_polytopes";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return SmallPolytope_FaceIneq(EXT, f_process);
  }
  // T can be a field or a ring here
  eProg = "lrs";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return lrs::DualDescriptionFaceIneq(EXT, f_process);
  // It applies to the field case or ring
  eProg = "pd_lrs";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return POLY_DualDescription_PrimalDualFaceIneq(EXT, f_process, os);
  //
  // The external programs are available only for rationl types
  //
  if constexpr (is_implementation_of_Q<T>::value) {
    eProg = "glrs";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return DualDescExternalProgramFaceIneq(EXT, "glrs", f_process, os);
    //
    eProg = "ppl_ext";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return DualDescExternalProgramFaceIneq(EXT, "ppl_lcdd", f_process, os);
    //
    eProg = "cdd_ext";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return DualDescExternalProgramFaceIneq(EXT, "lcdd_gmp", f_process, os);
    //
    eProg = "normaliz";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return DualDescExternalProgramFaceIneq(EXT, "normaliz", f_process, os);
  }
  //
  std::cerr << "ERROR in DirectFacetComputationFaceIneq\n";
  std::cerr << "No right program found with ansProg=" << ansProg
            << " or incorrect output\n";
  std::cerr << "List of authorized programs :";
  bool IsFirst = true;
  for (auto &eP : ListProg) {
    if (!IsFirst)
      std::cerr << " ,";
    IsFirst = false;
    std::cerr << " " << eP;
  }
  std::cerr << "\n";
  throw TerminalException{1};
}

template <typename T, typename Tgroup>
vectface DirectFacetOrbitComputation(MyMatrix<T> const &EXT, Tgroup const &GRP,
                                     std::string const &ansProg,
                                     std::ostream &os) {
#ifdef TIMINGS_DUAL_DESC
  MicrosecondTime time;
#endif
#ifdef KEY_VALUE_DUAL_DESC
  MicrosecondTime time_total;
#endif
  vectface ListIncd = DirectFacetComputationIncidence(EXT, ansProg, os);
#ifdef TIMINGS_DUAL_DESC
  os << "|DualDescription|=" << time << " |ListIncd|=" << ListIncd.size()
     << "\n";
#endif
  if (ListIncd.size() == 0) {
    std::cerr << "We found ListIncd to be empty. A clear error\n";
    throw TerminalException{1};
  }
  if (GRP.size() == 1) {
    return ListIncd;
  }
  vectface TheOutput = OrbitSplittingSet(ListIncd, GRP);
#ifdef KEY_VALUE_DUAL_DESC
  os << "KEY=(DirectFacetOrbitComputation_" << EXT.rows() << "_" << EXT.cols()
     << "_" << GRP.size() << "_" << ansProg << "_" << ListIncd.size() << "_"
     << TheOutput.size() << ") VALUE=(" << time_total << ")\n";
#endif
  return TheOutput;
}

template <typename T, typename Tgroup>
std::vector<std::pair<Face, MyVector<T>>>
DirectFacetIneqOrbitComputation(MyMatrix<T> const &EXT, Tgroup const &GRP,
                                std::string const &ansProg, std::ostream &os) {
#ifdef TIMINGS_DUAL_DESC
  MicrosecondTime time;
#endif
#ifdef KEY_VALUE_DUAL_DESC
  MicrosecondTime time_total;
#endif
  std::vector<std::pair<Face, MyVector<T>>> ListReturn;
  auto f_process = [&](std::pair<Face, MyVector<T>> const &pair_face) -> void {
    ListReturn.push_back(pair_face);
  };
  DirectFacetComputationFaceIneq(EXT, ansProg, f_process, os);
#ifdef TIMINGS_DUAL_DESC
  os << "|DualDescription|=" << time << "\n";
#endif
#ifdef DEBUG_DUAL_DESC
  os << "Found  |ListIncd|=" << ListReturn.size() << "\n";
#endif
  if (ListReturn.size() == 0) {
    std::cerr << "We found ListIncd to be empty. A clear error\n";
    throw TerminalException{1};
  }
  if (GRP.size() == 1) {
    return ListReturn;
  }
  std::vector<std::pair<Face, MyVector<T>>> TheOutput =
      OrbitSplittingMap(ListReturn, GRP);
#ifdef DEBUG_DUAL_DESC
  os << "|OrbitSplittingMap|=" << time << "\n";
#endif
#ifdef KEY_VALUE_DUAL_DESC
  os << "KEY=(DirectFacetIneqOrbitComputation_" << EXT.rows() << "_"
     << EXT.cols() << "_" << GRP.size() << "_" << ansProg << "_"
     << ListReturn.size() << "_" << TheOutput.size() << ") VALUE=("
     << time_total << ")\n";
#endif
  return TheOutput;
}

// clang-format off
#endif  // SRC_POLY_POLY_DIRECTDUALDESC_H_
// clang-format on
