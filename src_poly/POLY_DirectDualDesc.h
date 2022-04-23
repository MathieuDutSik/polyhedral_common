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

template <typename T>
vectface DualDescExternalProgram(MyMatrix<T> const &EXT,
                                 std::string const &eCommand) {
#ifdef TIMINGS
  SingletonTime time1;
#endif
  size_t n_row = EXT.rows();
  size_t n_col = EXT.cols();
  size_t DimEXT = n_col + 1;
  std::string rndStr = random_string(20);
  std::string prefix = "/tmp/";
  std::string suffix = "DD_" + std::to_string(n_row) +
    "_" + std::to_string(n_col) + "_" + rndStr;
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
    std::ofstream os(FileI);
    if (eCommand == "normaliz") {
      os << "amb_space " << n_col << "\n";
      os << "cone " << n_row << "\n";
      for (size_t i_row = 0; i_row < n_row; i_row++) {
        for (size_t i_col = 0; i_col < n_col; i_col++)
          os << " " << EXT(i_row, i_col);
        os << "\n";
      }
      os << "SupportHyperplanes\n";
    } else {
      os << "V-representation\n";
      os << "begin\n";
      os << n_row << " " << DimEXT << " integer\n";
      for (size_t i_row = 0; i_row < n_row; i_row++) {
        os << "0";
        for (size_t i_col = 0; i_col < n_col; i_col++)
          os << " " << EXT(i_row, i_col);
        os << "\n";
      }
      os << "end\n";
    }
  }
#ifdef TIMINGS
  SingletonTime time2;
  std::cerr << "|FileWriting|=" << ms(time1, time2) << "\n";
#endif
  //  std::cerr << "FileO=" << FileO << " created\n";
  //
  // Now calling the external program
  //
  std::string order;
  if (eCommand == "normaliz") {
    order = eCommand + " " + FileI;
  } else {
    order = eCommand + " " + FileI + " > " + FileO + " 2> " + FileE;
  }
  std::cerr << "order=" << order << "\n";
  int iret1 = system(order.c_str());
#ifdef TIMINGS
  SingletonTime time3;
  std::cerr << "|glrs/ppl/cdd|=" << ms(time2, time3) << "\n";
#endif
  std::cerr << "External program terminated\n";
  if (iret1 != 0) {
    std::cerr << "The program has not terminated correctly\n";
    std::cerr << "FileO = " << FileO << "\n";
    throw TerminalException{1};
  }
  //  std::cerr << "iret1=" << iret1 << "\n";
  vectface ListFace(n_row);
  //
  std::ifstream is(FileO);
  size_t iLineLimit = 0;
  std::vector<T> LVal(DimEXT);
  T eScal;
  size_t shift;
  if (eCommand == "normaliz") {
    shift = 0;
  } else {
    shift = 1;
  }
#ifdef USE_ISINCD
  auto isincd = [&](size_t i_row) -> bool {
    eScal = 0;
    for (size_t i = shift; i < DimEXT; i++)
      eScal += LVal[i] * EXT(i_row, i - shift);
    return eScal == 0;
  };
#else
  Face f(n_row);
#endif
  size_t pos_wrt = 0;
  auto f_read = [&](const std::string &str) -> void {
    ParseScalar_inplace<T>(str, LVal[pos_wrt]);
    pos_wrt++;
  };
  std::string line;
  size_t iLine = 0;
  // For this case we know at the beginning the number of inequalities
  if (eCommand == "ppl_lcdd" || eCommand == "lcdd_gmp" || eCommand == "normaliz") {
    size_t headersize;
    if (eCommand == "lcdd_gmp")
      headersize = 4;
    if (eCommand == "ppl_lcdd")
      headersize = 3;
    if (eCommand == "normaliz")
      headersize = 2 + 5 + 6 + n_row + 6;
    std::cerr << "headersize=" << headersize << "\n";
    while (std::getline(is, line)) {
      std::cerr << "iLine=" << iLine << " line=" << line << "\n";
      if (iLine == headersize - 1) {
        std::cerr << "   Assigning iLineLimit\n";
        // Determining the number of entries
        std::vector<std::string> LStr = STRING_Split(line, " ");
        std::cerr << "LStr[0]=" << LStr[0] << "\n";
        iLineLimit = headersize + ParseScalar<size_t>(LStr[0]);
        std::cerr << "iLineLimit=" << iLineLimit << "\n";
      }
      if (iLine >= headersize && (iLineLimit == 0 || iLine < iLineLimit)) {
        STRING_Split_f(line, " ", f_read);
        pos_wrt = 0;
#ifdef USE_ISINCD
        ListFace.InsertFaceRef(isincd);
#else
        for (size_t i_row = 0; i_row < n_row; i_row++) {
          eScal = 0;
          for (size_t i = shift; i < DimEXT; i++)
            eScal += LVal[i] * EXT(i_row, i - shift);
          f[i_row] = static_cast<bool>(eScal == 0);
        }
        ListFace.push_back(f);
#endif
      }
      iLine++;
    }
  }
  if (eCommand == "glrs") {
    size_t headersize = 7;
    while (std::getline(is, line)) {
      if (line == "end")
        break;
      if (iLine >= headersize) {
        STRING_Split_f(line, " ", f_read);
        pos_wrt = 0;
#ifdef USE_ISINCD
        ListFace.InsertFaceRef(isincd);
#else
        for (size_t i_row = 0; i_row < n_row; i_row++) {
          eScal = 0;
          for (size_t i = shift; i < DimEXT; i++)
            eScal += LVal[i] * EXT(i_row, i - shift);
          f[i_row] = static_cast<bool>(eScal == 0);
        }
        ListFace.push_back(f);
#endif
      }
      iLine++;
    }
  }
#ifdef TIMINGS
  SingletonTime time4;
  std::cerr << "|FileRead|=" << ms(time3, time4) << "\n";
#endif
  std::cerr << "FileI = " << FileI << "    FileO = " << FileO << "\n";
  //  throw TerminalException{1};
  RemoveFileIfExist(FileI);
  RemoveFileIfExist(FileO);
  RemoveFileIfExist(FileE);
  return ListFace;
}

template <typename T>
vectface DirectFacetOrbitComputation_nogroup(MyMatrix<T> const &EXT,
                                             std::string const &ansProg) {
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
    return lrs::DualDescription_temp_incd(EXT);
  //
  eProg = "lrs_ring";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return lrs::DualDescription_temp_incd_reduction(EXT);
  //
  eProg = "glrs";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return DualDescExternalProgram(EXT, "glrs");
  //
  eProg = "ppl_ext";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return DualDescExternalProgram(EXT, "ppl_lcdd");
  //
  eProg = "cdd_ext";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return DualDescExternalProgram(EXT, "lcdd_gmp");
  //
  eProg = "normaliz";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return DualDescExternalProgram(EXT, "normaliz");
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
                                     std::string const &ansProg) {
#ifdef TIMINGS
  SingletonTime time1;
#endif
  vectface ListIncd = DirectFacetOrbitComputation_nogroup(EXT, ansProg);
#ifdef TIMINGS
  SingletonTime time2;
  std::cerr << "|DualDescription|=" << ms(time1, time2)
            << " |ListIncd|=" << ListIncd.size() << "\n";
#endif
  if (ListIncd.size() == 0) {
    std::cerr << "We found ListIncd to be empty. A clear error\n";
    throw TerminalException{1};
  }
  vectface TheOutput = OrbitSplittingSet(ListIncd, GRP);
#ifdef TIMINGS
  SingletonTime time3;
  std::cerr << "KEY=(OrbitSplitting_" << EXT.rows() << "_" << EXT.cols() << "_"
            << GRP.size()
            << "_" << ansProg << "_" << ListIncd.size() << "_"
            << TheOutput.size() << ") VALUE=" << ms(time2,time3) << "\n";
#endif
  return TheOutput;
}

#endif  // SRC_POLY_POLY_DIRECTDUALDESC_H_
