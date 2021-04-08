#ifndef TEMP_DIRECT_DUAL_DESC_H
#define TEMP_DIRECT_DUAL_DESC_H

#include "POLY_lrslib.h"
#include "POLY_cddlib.h"
#include "POLY_c_cddlib.h"
#include "Basic_string.h"



template<typename T>
std::vector<Face> CDD_PPL_ExternalProgram(MyMatrix<T> const& EXT, std::string const& eCommand)
{
  size_t n_row = EXT.rows();
  size_t n_col = EXT.cols();
  size_t DimEXT = n_col + 1;
  std::string rndStr = random_string(20);
  std::string prefix = "/tmp/";
  std::string FileO = prefix + "EXT_" + std::to_string(n_row) + "_" + std::to_string(n_col) + "_" + rndStr + ".ext";
  std::string FileI = prefix + "INE_" + std::to_string(n_row) + "_" + std::to_string(n_col) + "_" + rndStr + ".ext";
  std::string FileE = prefix + "INE_" + std::to_string(n_row) + "_" + std::to_string(n_col) + "_" + rndStr + ".err";
  std::ofstream os(FileO);
  os << "V-representation\n";
  os << "begin\n";
  os << n_row << " " << DimEXT << " integer\n";
  for (size_t i_row=0; i_row<n_row; i_row++) {
    os << "0";
    for (size_t i_col=0; i_col<n_col; i_col++)
      os << " " << EXT(i_row, i_col);
    os << "\n";
  }
  os << "end\n";
  os.close();
  //  std::cerr << "FileO=" << FileO << " created\n";
  //
  // Now calling the external program
  //
  std::string order = eCommand + " " + FileO + " > " + FileI + " 2> " + FileE;
  int iret1=system(order.c_str());
  if (iret1 != 0) {
    std::cerr << "The program has not terminated correctly\n";
    std::cerr << "FileO=" << FileO << "\n";
    throw TerminalException{1};
  }
  //  std::cerr << "iret1=" << iret1 << "\n";
  std::vector<Face> ListFace;
  //
  std::ifstream is(FileI);
  std::string line;
  size_t iLine = 0;
  size_t iLineLimit = 0;
  std::vector<T> LVal(DimEXT);
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
      std::vector<std::string> LStr = STRING_Split(line, " ");
      for (size_t i=0; i<DimEXT; i++)
        LVal[i] = ParseScalar<T>(LStr[i]);
      Face face(n_row);
      for (size_t i_row=0; i_row<n_row; i_row++) {
        T eScal=0;
        for (size_t i=1; i<DimEXT; i++)
          eScal += LVal[i] * EXT(i_row,i-1);
        if (eScal == 0)
          face[i_row] = 1;
      }
      ListFace.push_back(face);
    }
    iLine++;
  }
  //  std::cerr << "FileI = " << FileI << "    FileO = " << FileO << "\n";
  RemoveFileIfExist(FileI);
  RemoveFileIfExist(FileO);
  RemoveFileIfExist(FileE);
  return ListFace;
}






template<typename T, typename Tgroup>
std::vector<Face> DirectFacetOrbitComputation(MyMatrix<T> const& EXT, Tgroup const& GRP, std::string const& ansProg, std::ostream& os)
{
  int nbVert=EXT.rows();
  int nbCol=EXT.cols();
  bool WeAreDone=false;
  std::vector<Face> ListIncd;
#ifdef DEBUG_DIRECT_DUAL_DESC
  if (nbVert >= 56) {
    std::string FileSave="EXT_" + IntToString(nbVert);
    std::ofstream os_save(FileSave);
    WriteMatrix(os_save, EXT);
  }
#endif
  os << "DFOC prog=" << ansProg << " |EXT|=" << nbVert << " nbCol=" << nbCol << " |GRP|=" << GRP.size() << "\n";
  std::vector<std::string> ListProg;
  std::string eProg;
  //
  eProg = "cdd"; ListProg.push_back(eProg);
  if (ansProg == eProg) {
    ListIncd = cdd::DualDescription_incd(EXT);
    WeAreDone=true;
  }
  //
  eProg = "lrs"; ListProg.push_back(eProg);
  if (ansProg == eProg) {
    ListIncd = lrs::DualDescription_temp_incd(EXT);
    WeAreDone=true;
  }
  //
  eProg = "lrs_ring"; ListProg.push_back(eProg);
  if (ansProg == eProg) {
    ListIncd = lrs::DualDescription_temp_incd_reduction(EXT);
    WeAreDone=true;
  }
  //
  eProg = "ppl_ext"; ListProg.push_back(eProg);
  if (ansProg == eProg) {
    ListIncd = CDD_PPL_ExternalProgram(EXT, "ppl_lcdd");
    WeAreDone=true;
  }
  //
  eProg = "cdd_ext"; ListProg.push_back(eProg);
  if (ansProg == eProg) {
    ListIncd = CDD_PPL_ExternalProgram(EXT, "lcdd_gmp");
    WeAreDone=true;
  }
  //
  eProg = "cdd_cbased"; ListProg.push_back(eProg);
  if (ansProg == eProg) {
    ListIncd = cbased_cdd::DualDescription_incd(EXT);
    WeAreDone=true;
  }
  //
  if (!WeAreDone || ListIncd.size() == 0) {
    std::cerr << "ERROR: No right program found with ansProg=" << ansProg << " or incorrect output\n";
    std::cerr << "WeAreDone=" << WeAreDone << " |ListIncd|=" << ListIncd.size() << "\n";
    std::cerr << "List of authorized programs :";
    bool IsFirst=true;
    for (auto & eP : ListProg) {
      if (!IsFirst)
        std::cerr << " ,";
      IsFirst=false;
      std::cerr << " " << eP;
    }
    std::cerr << "\n";
    throw TerminalException{1};
  }
  std::vector<Face> TheOutput=OrbitSplittingSet(ListIncd, GRP);
  os << "DFOC |GRP|=" << GRP.size() << " |ListIncd|=" << ListIncd.size() << " |TheOutput|=" << TheOutput.size() << "\n";
  return TheOutput;
}


#endif
