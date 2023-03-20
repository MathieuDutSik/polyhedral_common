// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_DIRECTDUALDESC_H_
#define SRC_POLY_POLY_DIRECTDUALDESC_H_

#include "Basic_string.h"
#include "POLY_c_cddlib.h"
#include "POLY_cddlib.h"
#include "POLY_lrslib.h"
#include <string>
#include <vector>

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

template<typename T>
size_t GetShift([[maybe_unused]] MyMatrix<T> const &EXT, std::string const &eCommand) {
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
#ifdef TIMINGS
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
#ifdef TIMINGS
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
#ifdef TIMINGS
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
  auto check_consistency=[&]() -> void {
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
    size_t headersize = 7;
    while (std::getline(is, line)) {
      if (line == "end")
        break;
      if (iLine >= headersize)
        process_line();
      iLine++;
    }
  }
#ifdef TIMINGS
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
  auto f_insert=[&](std::vector<T> const& LVal) -> void {
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

template <typename T>
MyMatrix<T> DualDescExternalProgramIneq(MyMatrix<T> const &EXT,
                                        std::string const &eCommand,
                                        std::ostream &os) {
  //  size_t n_row = EXT.rows();
  size_t n_col = EXT.cols();
  size_t DimEXT = n_col + 1;
  size_t shift = GetShift(EXT, eCommand);
  int nbColRed = DimEXT - shift;
  MyVector<T> V(nbColRed);
  std::vector<MyVector<T>> ListVect;
  auto f_insert=[&](std::vector<T> const& LVal) -> void {
    for (int i=0; i<nbColRed; i++) {
      V(i) = LVal[i + shift];
    }
    ListVect.push_back(V);
  };
  DualDescExternalProgramGeneral(EXT, f_insert, eCommand, os);
  return MatrixFromVectorFamily(ListVect);
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
  eProg = "cdd";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return cdd::DualDescription_incd(EXT);
  //
  eProg = "lrs";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return lrs::DualDescription_incd(EXT);
  //
  eProg = "lrs_ring";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return lrs::DualDescription_incd_reduction(EXT);
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
  eProg = "cdd";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return cdd::DualDescription(EXT);
  //
  eProg = "lrs";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return lrs::DualDescription(EXT);
  //
  eProg = "lrs_ring";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return lrs::DualDescription_reduction(EXT);
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






template <typename T, typename Tgroup>
vectface DirectFacetOrbitComputation(MyMatrix<T> const &EXT, Tgroup const &GRP,
                                     std::string const &ansProg,
                                     std::ostream &os) {
#ifdef TIMINGS
  MicrosecondTime time;
#endif
  vectface ListIncd = DirectFacetComputationIncidence(EXT, ansProg, os);
#ifdef TIMINGS
  os << "|DualDescription|=" << time << " |ListIncd|=" << ListIncd.size()
     << "\n";
#endif
  if (ListIncd.size() == 0) {
    std::cerr << "We found ListIncd to be empty. A clear error\n";
    throw TerminalException{1};
  }
  vectface TheOutput = OrbitSplittingSet(ListIncd, GRP);
#ifdef TIMINGS
  os << "KEY=(OrbitSplitting_" << EXT.rows() << "_" << EXT.cols() << "_"
     << GRP.size() << "_" << ansProg << "_" << ListIncd.size() << "_"
     << TheOutput.size() << ") VALUE=" << time << "\n";
#endif
  return TheOutput;
}

template <typename T, typename Tgroup>
std::vector<std::pair<MyVector<T>,Face>> DirectFacetIneqOrbitComputation(MyMatrix<T> const &EXT, Tgroup const &GRP,
                                                                         std::string const &ansProg,
                                                                         std::ostream &os) {
#ifdef TIMINGS
  MicrosecondTime time;
#endif
  vectface ListIncd = DirectFacetComputationIncidence(EXT, ansProg, os);
#ifdef TIMINGS
  os << "|DualDescription|=" << time << " |ListIncd|=" << ListIncd.size()
     << "\n";
#endif
  if (ListIncd.size() == 0) {
    std::cerr << "We found ListIncd to be empty. A clear error\n";
    throw TerminalException{1};
  }
  vectface TheOutput = OrbitSplittingSet(ListIncd, GRP);
#ifdef TIMINGS
  os << "KEY=(OrbitSplitting_" << EXT.rows() << "_" << EXT.cols() << "_"
     << GRP.size() << "_" << ansProg << "_" << ListIncd.size() << "_"
     << TheOutput.size() << ") VALUE=" << time << "\n";
#endif
  return TheOutput;
}




// clang-format off
#endif  // SRC_POLY_POLY_DIRECTDUALDESC_H_
// clang-format on
