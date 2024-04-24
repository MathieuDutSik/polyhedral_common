// Copyright (C) 2023 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_TSPACE_NAMELIST_H_
#define SRC_LATT_TSPACE_NAMELIST_H_

// clang-format off
#include "Namelist.h"
#include "MatrixGroupAverage.h"
// clang-format on

SingleBlock SINGLEBLOCK_Get_Tspace_Description() {
  std::map<std::string, std::string> ListStringValues1_doc;
  std::map<std::string, std::string> ListBoolValues1_doc;
  std::map<std::string, std::string> ListIntValues1_doc;
  ListStringValues1_doc["TypeTspace"] = "The possible type of T-space\n\
InvGroup: The space of matrices invariant under a group\n   \
RealQuad: The space corresponding to a real quadratic field\n   \
ImagQuad: The space corresponding to a real quadratic field\n   \
Raw     : the space given by matrices\n    \
File    : A file containing the Tspace";
  ListIntValues1_doc["RealImagDim"] = "Default: 0\n\
The dimension d of the space GL_d(R) for R a number ring";
  ListIntValues1_doc["RealImagSum"] = "Default: 0\n\
The sum of the two roots of the polynomial. Only relevant for RealQuad and ImagQuad";
  ListIntValues1_doc["RealImagProd"] = "Default: 0\n\
The product of the two roots of the polynomial. Only relevant for RealQuad and ImagQuad";
  ListStringValues1_doc["SuperMatMethod"] = "Default: NotNeeded\n\
NotNeeded: For InvGroup, RealQuad, ImagQuad, the supermat is computed internally and so this option needs to be selected.\n\
Compute: Compute from the basis of the T-space\n\
File: Use matrix read from the file";
  ListStringValues1_doc["FileInvGroup"] = "Default: unset\n\
The file containing the list of generators";
  ListStringValues1_doc["FileListMat"] = "Default: unset\n\
The file containing the list of matrices";
  ListStringValues1_doc["FileSuperMat"] = "Default: unset\n\
The file containing the SuperMat";
  ListStringValues1_doc["ListComm"] = "Default: trivial\n\
Trivial: This set it up to empty\n\
Use_realimag: Use the commuting from the real or imag case\n\
File: The ListComm will be obtained from the FileListComm";
  ListStringValues1_doc["FileListComm"] = "Default: unset\n\
The file containing the list of commutting matrices";
  ListStringValues1_doc["PtGroupMethod"] = "Default: Trivial\n\
Trivial: Use the group of order 2 generated by -Id_n\n\
Compute: Compute the pointwise group stabilizer by using Plesken Souvignier with the vector configuration obtained by the supermat\n\
InvGroupInit: Set the PtGroup from the group defining it\n\
File: Use the pointwise stabilizer obtained from the file";
  ListStringValues1_doc["FilePtGroupGenerator"] = "Default: unset\n\
unset: If unset then the file is not needed\n\
The filename for the list of generators of the pointwise stabilizer";
  ListStringValues1_doc["FileListSubspaces"] = "Default: unset\n\
unset: If unset then the file is not needed\n\
The filename for the list of subspaces that are preserved";
  ListStringValues1_doc["FileLinSpa"] = "Default: unset\n\
unset: If unset then the file is not needed\n\
The filename used for reading the whole T-space";
  SingleBlock BlockTSPACE;
  BlockTSPACE.setListStringValues(ListStringValues1_doc);
  BlockTSPACE.setListBoolValues(ListBoolValues1_doc);
  BlockTSPACE.setListIntValues(ListIntValues1_doc);
  return BlockTSPACE;
}

FullNamelist NAMELIST_GetOneTSPACE() {
  std::map<std::string, SingleBlock> ListBlock;
  ListBlock["TSPACE"] = SINGLEBLOCK_Get_Tspace_Description();
  return {ListBlock, "undefined"};
}

template<typename T>
LinSpaceMatrix<T> ReadLinSpaceFile(std::string const& eFile) {
  std::ifstream is(eFile);
  MyMatrix<T> SuperMat = ReadMatrix<T>(is);
  int n = SuperMat.rows();
  std::vector<MyMatrix<T>> ListMat = ReadListMatrix<T>(is);
  std::vector<std::vector<T>> ListLineMat;
  for (auto & eMat : ListMat) {
    std::vector<T> eV = GetLineVector(eMat);
    ListLineMat.push_back(eV);
  }
  std::vector<MyMatrix<T>> ListComm = ReadListMatrix<T>(is);
  std::vector<MyMatrix<T>> ListSubspaces = ReadListMatrix<T>(is);
  std::vector<MyMatrix<T>> PtStabGens = ReadListMatrix<T>(is);
  //
  MyMatrix<T> ListMatAsBigMat = GetListMatAsBigMat(ListMat);
  return {n, SuperMat, ListMat, ListLineMat, ListMatAsBigMat, ListComm, ListSubspaces, PtStabGens};
}

template<typename T>
void WriteLinSpace(std::ostream& os, LinSpaceMatrix<T> const& LinSpa) {
  WriteMatrix(os, LinSpa.SuperMat);
  WriteListMatrix(os, LinSpa.ListMat);
  WriteListMatrix(os, LinSpa.ListComm);
  WriteListMatrix(os, LinSpa.ListSubspaces);
  WriteListMatrix(os, LinSpa.PtStabGens);
}

template<typename T>
void WriteLinSpaceFile(std::string const& eFile, LinSpaceMatrix<T> const& LinSpa) {
  std::ofstream os(eFile);
  WriteLinSpace(os, LinSpa);
}

template<typename T, typename Tint>
LinSpaceMatrix<T> ReadTspace(SingleBlock const& Blk, std::ostream & os) {
  std::string TypeTspace = Blk.ListStringValues.at("TypeTspace");
  LinSpaceMatrix<T> LinSpaRet;
  auto set_paperwork=[&]() -> void {
    for (auto & eMat : LinSpaRet.ListMat) {
      std::vector<T> eV = GetLineVector(eMat);
      LinSpaRet.ListLineMat.push_back(eV);
    }
    int n_mat = LinSpaRet.ListMat.size();
    if (n_mat > 0) {
      LinSpaRet.ListMatAsBigMat = GetListMatAsBigMat(LinSpaRet.ListMat);
    }
  };
  auto set_listcomm=[&]() -> void {
    std::string ListComm = Blk.ListStringValues.at("ListComm");
    if (ListComm == "Trivial") {
      LinSpaRet.ListComm.clear();
      return;
    }
    if (ListComm == "Use_realimag") {
      if (TypeTspace != "RealQuad" && TypeTspace != "ImagQuad") {
        std::cerr << "We have TypeTspace=" << TypeTspace << "\n";
        std::cerr << "But only RealQuad and ImagQuad are allowed\n";
        throw TerminalException{1};
      }
      // Do nothing in that case
      return;
    }
    if (ListComm == "File") {
      std::string FileListComm = Blk.ListStringValues.at("FileListComm");
      LinSpaRet.ListComm = ReadListMatrixFile<T>(FileListComm);
      return;
    }
    std::cerr << "Failed to find an option for ListComm that suits\n";
    throw TerminalException{1};
  };
  auto set_subspaces=[&]() -> void {
    std::string FileListSubspaces = Blk.ListStringValues.at("FileListSubspaces");
    if (FileListSubspaces != "unset") {
      LinSpaRet.ListSubspaces = ReadListMatrixFile<T>(FileListSubspaces);
    }
  };
  auto set_pt_stab=[&]() -> void {
    std::string PtGroupMethod = Blk.ListStringValues.at("PtGroupMethod");
    if (PtGroupMethod == "Trivial") {
      MyVector<T> eGen = -IdentityMat<T>(LinSpaRet.n);
      LinSpaRet.PtStabGens = {eGen};
      return;
    }
    if (PtGroupMethod == "Compute") {
      std::vector<MyMatrix<Tint>> ListGens = ComputePointStabilizerTspace<T,Tint>(LinSpaRet.SuperMat, LinSpaRet.ListMat, os);
      for (auto & eGen : ListGens) {
        LinSpaRet.PtStabGens.push_back(UniversalMatrixConversion<T,Tint>(eGen));
      }
      return;
    }
    if (PtGroupMethod == "InvGroupInit") {
      std::string FileInvGroup = Blk.ListStringValues.at("FileInvGroup");
      LinSpaRet.PtStabGens = ReadListMatrixFile<T>(FileInvGroup);
      return;
    }
    if (PtGroupMethod == "File") {
      std::string FilePtGroupGenerator = Blk.ListStringValues.at("FilePtGroupGenerator");
      if (FilePtGroupGenerator == "unset") {
        std::cerr << "The FilePtGroupGenerator has not been set up, or set to unset\n";
        throw TerminalException{1};
      }
      LinSpaRet.PtStabGens = ReadListMatrixFile<T>(FilePtGroupGenerator);
      return;
    }
    std::cerr << "Failed to find an option for PtGroupMethod that suits\n";
    throw TerminalException{1};
  };
  auto set_supermat=[&]() -> void {
    std::string SuperMatMethod = Blk.ListStringValues.at("SuperMatMethod");
    if (TypeTspace != "RealQuad" && TypeTspace != "ImagQuad" && TypeTspace != "InvGroup") {
      if (SuperMatMethod == "NotNeeded") {
        std::cerr << "We have TypeTspace=" << TypeTspace << "\n";
        std::cerr << "For NotNeeded, the option needs to be RealQuad, ImagQuad or InvGroup\n";
        throw TerminalException{1};
      }
    } else {
      if (SuperMatMethod != "NotNeeded") {
        std::cerr << "For the options RealQuad, ImageQuad and InvGroup, the option has to be NotNeeded\n";
        throw TerminalException{1};
      }
    }
    if (TypeTspace == "InvGroup") {
      MyMatrix<T> eMat = IdentityMat<T>(LinSpaRet.n);
      LinSpaRet.SuperMat = OrbitBarycenterSymmetricMatrix(eMat, LinSpaRet.ListMat);
      return;
    }
    if (SuperMatMethod == "NotNeeded") {
      return;
    }
    if (SuperMatMethod == "Compute") {
      LinSpaRet.SuperMat = GetOnePositiveDefiniteMatrix<T,Tint>(LinSpaRet.ListMat, os);
      return;
    }
    if (SuperMatMethod == "File") {
      std::string FileSuperMat = Blk.ListStringValues.at("FileSuperMat");
      if (FileSuperMat == "unset") {
        std::cerr << "The FileSuperMat has not been set up, or set to unset\n";
        throw TerminalException{1};
      }
      LinSpaRet.SuperMat = ReadMatrixFile<T>(FileSuperMat);
      return;
    }
    std::cerr << "Failed to find an option for SuperMatMethod that suits\n";
    throw TerminalException{1};
  };
  if (TypeTspace == "RealQuad" || TypeTspace == "ImagQuad") {
    int n = Blk.ListIntValues.at("RealImagDim");
    int eSum = Blk.ListIntValues.at("RealImagSum");
    int eProd = Blk.ListIntValues.at("RealImagProd");
    if (TypeTspace == "RealQuad")
      LinSpaRet = ComputeRealQuadraticSpace<T>(n, eSum, eProd);
    if (TypeTspace == "ImagQuad")
      LinSpaRet = ComputeImagQuadraticSpace<T>(n, eSum, eProd);
    set_listcomm();
    set_subspaces();
    set_pt_stab();
    return LinSpaRet;
  }
  if (TypeTspace == "InvGroup") {
    std::string FileInvGroup = Blk.ListStringValues.at("FileInvGroup");
    std::vector<MyMatrix<T>> LGen = ReadListMatrixFile<T>(FileInvGroup);
    if (LGen.size() == 0) {
      std::cerr << "We have 0 matrices\n";
      throw TerminalException{1};
    }
    LinSpaRet.n = LGen[0].rows();
    LinSpaRet.ListMat = BasisInvariantForm(LinSpaRet.n, LGen);
    set_paperwork();
    set_supermat();
    set_listcomm();
    set_subspaces();
    set_pt_stab();
    return LinSpaRet;
  }
  if (TypeTspace == "Raw") {
    std::string FileListMat = Blk.ListStringValues.at("FileListMat");
    LinSpaRet.ListMat = ReadListMatrixFile<T>(FileListMat);
    if (LinSpaRet.ListMat.size() == 0) {
      std::cerr << "We have 0 matrices for ListMat\n";
      throw TerminalException{1};
    }
    LinSpaRet.n = LinSpaRet.ListMat[0].rows();
    set_paperwork();
    set_supermat();
    set_listcomm();
    set_subspaces();
    set_pt_stab();
    return LinSpaRet;
  }
  if (TypeTspace == "File") {
    std::string FileLinSpa = Blk.ListStringValues.at("FileLinSpa");
    return ReadLinSpaceFile<T>(FileLinSpa);
  }
  std::cerr << "Failed to find an option for TypeTspace that suits\n";
  throw TerminalException{1};
}



// clang-format off
#endif  // SRC_LATT_TSPACE_NAMELIST_H_
// clang-format on
