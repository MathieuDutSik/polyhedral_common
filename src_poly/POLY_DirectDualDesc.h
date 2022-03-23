#ifndef INCLUDE_POLY_DIRECTDUALDESC_H
#define INCLUDE_POLY_DIRECTDUALDESC_H

#include "Basic_string.h"
#include "POLY_c_cddlib.h"
#include "POLY_cddlib.h"
#include "POLY_lrslib.h"

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
  std::chrono::time_point<std::chrono::system_clock> time1 =
      std::chrono::system_clock::now();
#endif
  size_t n_row = EXT.rows();
  size_t n_col = EXT.cols();
  size_t DimEXT = n_col + 1;
  std::string rndStr = random_string(20);
  std::string prefix = "/tmp/";
  std::string FileO = prefix + "EXT_" + std::to_string(n_row) + "_" +
                      std::to_string(n_col) + "_" + rndStr + ".ext";
  std::string FileI = prefix + "INE_" + std::to_string(n_row) + "_" +
                      std::to_string(n_col) + "_" + rndStr + ".ext";
  std::string FileE = prefix + "INE_" + std::to_string(n_row) + "_" +
                      std::to_string(n_col) + "_" + rndStr + ".err";
  {
    std::ofstream os(FileO);
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
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 =
      std::chrono::system_clock::now();
  std::cerr << "|FileWriting|="
            << std::chrono::duration_cast<std::chrono::microseconds>(time2 -
                                                                     time1)
                   .count()
            << "\n";
#endif
  //  std::cerr << "FileO=" << FileO << " created\n";
  //
  // Now calling the external program
  //
  std::string order = eCommand + " " + FileO + " > " + FileI + " 2> " + FileE;
  std::cerr << "order=" << order << "\n";
  int iret1 = system(order.c_str());
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time3 =
      std::chrono::system_clock::now();
  std::cerr << "|glrs/ppl/cdd|="
            << std::chrono::duration_cast<std::chrono::microseconds>(time3 -
                                                                     time2)
                   .count()
            << "\n";
#endif
  std::cerr << "External program terminated\n";
  if (iret1 != 0) {
    std::cerr << "The program has not terminated correctly\n";
    std::cerr << "FileO=" << FileO << "\n";
    throw TerminalException{1};
  }
  //  std::cerr << "iret1=" << iret1 << "\n";
  vectface ListFace(n_row);
  //
  std::ifstream is(FileI);
  std::string line;
  size_t iLine = 0;
  size_t iLineLimit = 0;
  std::vector<T> LVal(DimEXT);
  T eScal;
#ifdef USE_ISINCD
  auto isincd = [&](size_t i_row) -> bool {
    eScal = 0;
    for (size_t i = 1; i < DimEXT; i++)
      eScal += LVal[i] * EXT(i_row, i - 1);
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
  if (eCommand == "ppl_lcdd" || eCommand == "lcdd_gmp") {
    size_t headersize;
    if (eCommand == "ppl_lcdd")
      headersize = 3;
    else
      headersize = 4;
    while (std::getline(is, line)) {
      //    std::cerr << "iLine=" << iLine << " line=" << line << "\n";
      if (iLine == headersize - 1) {
        std::vector<std::string> LStr = STRING_Split(line, " ");
        iLineLimit = headersize + ParseScalar<size_t>(LStr[0]);
        //      std::cerr << "iLineLimit=" << iLineLimit << "\n";
      }
      if (iLine >= headersize && (iLineLimit == 0 || iLine < iLineLimit)) {
        STRING_Split_f(line, " ", f_read);
        pos_wrt = 0;
#ifdef USE_ISINCD
        ListFace.InsertFaceRef(isincd);
#else
        for (size_t i_row = 0; i_row < n_row; i_row++) {
          eScal = 0;
          for (size_t i = 1; i < DimEXT; i++)
            eScal += LVal[i] * EXT(i_row, i - 1);
          f[i_row] = bool(eScal == 0);
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
          for (size_t i = 1; i < DimEXT; i++)
            eScal += LVal[i] * EXT(i_row, i - 1);
          f[i_row] = bool(eScal == 0);
        }
        ListFace.push_back(f);
#endif
      }
      iLine++;
    }
  }
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time4 =
      std::chrono::system_clock::now();
  std::cerr << "|FileRead|="
            << std::chrono::duration_cast<std::chrono::microseconds>(time4 -
                                                                     time3)
                   .count()
            << "\n";
#endif
  //  std::cerr << "FileI = " << FileI << "    FileO = " << FileO << "\n";
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
  std::chrono::time_point<std::chrono::system_clock> time1 =
      std::chrono::system_clock::now();
#endif
  vectface ListIncd = DirectFacetOrbitComputation_nogroup(EXT, ansProg);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 =
      std::chrono::system_clock::now();
  std::cerr << "|DualDescription|="
            << std::chrono::duration_cast<std::chrono::microseconds>(time2 -
                                                                     time1)
                   .count()
            << " |ListIncd|=" << ListIncd.size() << "\n";
#endif
  if (ListIncd.size() == 0) {
    std::cerr << "We found ListIncd to be empty. A clear error\n";
    throw TerminalException{1};
  }
  vectface TheOutput = OrbitSplittingSet(ListIncd, GRP);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time3 =
      std::chrono::system_clock::now();
  auto dur32 =
      std::chrono::duration_cast<std::chrono::microseconds>(time3 - time2)
          .count();
  std::cerr << "KEY=(OrbitSplitting," << EXT.rows() << "," << EXT.cols() << ","
            << GRP.size() << ","
            << "," << ansProg << "," << ListIncd.size() << ","
            << TheOutput.size() << ") VALUE=" << dur32 << "\n";
#endif
  return TheOutput;
}

#endif
