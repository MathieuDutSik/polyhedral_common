#ifndef TEMP_DIRECT_DUAL_DESC_H
#define TEMP_DIRECT_DUAL_DESC_H

#include "POLY_lrslib.h"
#include "POLY_cddlib.h"
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
  //
  // Now calling the external program
  //
  std::string order = eCommand + " " + FileI + " " + FileO;
  std::vector<Face> ListFace;
  //
  std::ifstream is(FileI);
  std::string line;
  size_t iLine;
  size_t iLineLimit = 0;
  std::vector<T> LVal(DimEXT);
  while (std::getline(is, line)) {
    if (iLine == 1) {
      std::vector<std::string> LStr = STRING_Split(line, " ");
      iLineLimit = 2 + ParseScalar<size_t>(LStr[0]);
    }
    if (iLine >= 2 && (iLineLimit == 0 || iLine < iLineLimit)) {
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
  return ListFace;
}






template<typename T, typename Tgroup>
std::vector<Face> DirectFacetOrbitComputation(MyMatrix<T> const& EXT, Tgroup const& GRP, std::string const& ansProg, std::ostream& os)
{
  MyMatrix<T> EXTred=ColumnReduction(EXT);
  int nbVert=EXTred.rows();
  int nbCol=EXTred.cols();
  bool WeAreDone=false;
  std::vector<Face> ListIncd;
#ifdef DEBUG_DIRECT_DUAL_DESC
  if (nbVert >= 56) {
    std::string FileSave="EXT_" + IntToString(nbVert);
    std::ofstream os_save(FileSave);
    WriteMatrix(os_save, EXT);
  }
#endif
  os << "DFOC prog=" << ansProg << " |EXT|=" << nbVert << " nbCol=" << nbCol << "\n";
  if (ansProg == "cdd") {
    ListIncd=cdd::DualDescription_incd(EXTred);
    WeAreDone=true;
  }
  if (ansProg == "lrs") {
    ListIncd=lrs::DualDescription_temp_incd(EXTred);
    WeAreDone=true;
  }
  if (ansProg == "lrs_ring") {
    ListIncd=lrs::DualDescription_temp_incd_reduction(EXTred);
    WeAreDone=true;
  }
  if (ansProg == "ppl_ext") {
    ListIncd = CDD_PPL_ExternalProgram(EXTred, "ppl_lcdd");
    WeAreDone=true;
  }
  if (ansProg == "cdd_ext") {
    ListIncd = CDD_PPL_ExternalProgram(EXTred, "lcdd_gmp");
    WeAreDone=true;
  }
  if (!WeAreDone || ListIncd.size() == 0) {
    std::cerr << "No right program found with ansProg=" << ansProg << "\n";
    std::cerr << "Let us die\n";
    throw TerminalException{1};
  }
  std::vector<Face> TheOutput=OrbitSplittingSet(ListIncd, GRP);
  os << "DFOC |GRP|=" << GRP.size() << " |ListIncd|=" << ListIncd.size() << " |TheOutput|=" << TheOutput.size() << "\n";
  return TheOutput;
}


#endif
