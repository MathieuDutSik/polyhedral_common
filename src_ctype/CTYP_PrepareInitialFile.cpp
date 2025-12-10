// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "CtypeMPI_types.h"
#include "Namelist.h"
#include "rational.h"
#include <netcdf>
// clang-format on

FullNamelist NAMELIST_InitialPreparation() {
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, std::string> ListStringValues1;
  ListIntValues1["n"] = 6;
  ListIntValues1["NprocInput"] = 100;
  ListIntValues1["NprocOutput"] = 100;
  ListBoolValues1["ApplyCanonicalization"] = true;
  ListStringValues1["OutputType"] = "int16";
  ListStringValues1["NatureInput"] = "text";
  ListStringValues1["InputFile"] = "InputFile";
  ListStringValues1["PrefixInput"] = "LOGinput_";
  ListStringValues1["PrefixOutput"] = "LOGoutput_";
  SingleBlock BlockDATA;
  BlockDATA.setListIntValues(ListIntValues1);
  BlockDATA.setListStringValues(ListStringValues1);
  BlockDATA.setListBoolValues(ListBoolValues1);
  ListBlock["DATA"] = BlockDATA;
  // Merging all data
  return FullNamelist(ListBlock);
}

template <typename Tinput, typename Toutput>
void AppendSingleCtype_T(MyMatrix<Tinput> const &M, int const &NbAdj,
                         size_t const &pos, netCDF::NcVar &var_Ctype,
                         netCDF::NcVar &var_NbAdj) {
  size_t n_vect = M.rows();
  size_t n = M.cols();
  std::vector<size_t> start1{pos};
  std::vector<size_t> count1{1};
  std::vector<size_t> start2{pos, 0, 0};
  std::vector<size_t> count2{1, n_vect, n};
  //
  int val = NbAdj;
  var_NbAdj.putVar(start1, count1, &val);
  //
  std::vector<Toutput> A(n_vect * n);
  int idx = 0;
  for (size_t i_vect = 0; i_vect < n_vect; i_vect++)
    for (size_t i = 0; i < n; i++) {
      A[idx] = M(i_vect, i);
      idx++;
    }
  var_Ctype.putVar(start2, count2, A.data());
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    using Tint = int;
    FullNamelist eFull = NAMELIST_InitialPreparation();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CTYP_MPI_Enumeration_c [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    //
    // Parsing the input
    //
    SingleBlock const &BlDATA = eFull.get_block("DATA");
    size_t n = BlDATA.get_int("n");
    size_t n_vect = std::pow(2, n) - 1;
    size_t NprocInput = BlDATA.get_int("NprocInput");
    size_t NprocOutput = BlDATA.get_int("NprocOutput");
    std::string PrefixInput = BlDATA.get_string("PrefixInput");
    std::string InputFile = BlDATA.get_string("InputFile");
    std::string PrefixOutput = BlDATA.get_string("PrefixOutput");
    std::string NatureInput = BlDATA.get_string("NatureInput");
    int posNature =
        PositionVect({std::string("text"), std::string("netcdf")}, NatureInput);
    if (posNature == -1) {
      std::cerr << "The NatureInput should be text or netcdf\n";
      std::cerr << "NatureInput=" << NatureInput << "\n";
      throw TerminalException{1};
    }
    std::string OutputType = BlDATA.get_string("OutputType");
    int posType = PositionVect({std::string("byte"), std::string("short"),
                                std::string("int"), std::string("int64")},
                               OutputType);
    if (posType == -1) {
      std::cerr << "The OutputType should be byte (8 bits), short (16 bits), "
                   "int (32 bits) or int64 (64-bits)\n";
      std::cerr << "OutputType=" << OutputType << "\n";
      throw TerminalException{1};
    }
    bool ApplyCanonicalization = BlDATA.get_bool("ApplyCanonicalization");
    //
    // Creating the netcdf output files.
    //
    std::vector<size_t> NbWrite(NprocOutput, 0);
    std::vector<netCDF::NcFile> ListNC(NprocOutput);
    ;
    std::vector<netCDF::NcVar> ListVar_Ctype;
    std::vector<netCDF::NcVar> ListVar_NbAdj;
    for (size_t iProcO = 0; iProcO < NprocOutput; iProcO++) {
      std::string FileO = PrefixOutput + std::to_string(iProcO) + ".nc";
      std::cerr << "Creating netcdf file FileO=" << FileO << "\n";
      ListNC[iProcO].open(FileO, netCDF::NcFile::replace, netCDF::NcFile::nc4);
      //
      netCDF::NcDim eDimNbCtype = ListNC[iProcO].addDim("number_ctype");
      netCDF::NcDim eDimN = ListNC[iProcO].addDim("n", n);
      netCDF::NcDim eDimNvect = ListNC[iProcO].addDim("n_vect", n_vect);
      //
      std::vector<std::string> LDim3{"number_ctype", "n_vect", "n"};
      std::vector<std::string> LDim1{"number_ctype"};
      //
      netCDF::NcVar varCtype =
          ListNC[iProcO].addVar("Ctype", OutputType, LDim3);
      varCtype.putAtt("long_name", "Ctype canonicalized coordinates");
      varCtype.putAtt("units", "nondimensional");
      //
      netCDF::NcVar varNbAdj =
          ListNC[iProcO].addVar("nb_adjacent", "int", LDim1);
      varNbAdj.putAtt("long_name", "number of adjacent Ctypes");
      varNbAdj.putAtt("units", "nondimensional");
      //
      ListVar_Ctype.push_back(varCtype);
      ListVar_NbAdj.push_back(varNbAdj);
    }
    //
    auto AppendSingleCtype = [&](MyMatrix<Tint> const &M, int const &NbAdj,
                                 size_t const &iProc) -> void {
      size_t pos = NbWrite[iProc];
      if (posType == 0)
        AppendSingleCtype_T<Tint, int8_t>(M, NbAdj, pos, ListVar_Ctype[iProc],
                                          ListVar_NbAdj[iProc]);
      if (posType == 1)
        AppendSingleCtype_T<Tint, int16_t>(M, NbAdj, pos, ListVar_Ctype[iProc],
                                           ListVar_NbAdj[iProc]);
      if (posType == 2)
        AppendSingleCtype_T<Tint, int32_t>(M, NbAdj, pos, ListVar_Ctype[iProc],
                                           ListVar_NbAdj[iProc]);
      if (posType == 3)
        AppendSingleCtype_T<Tint, int64_t>(M, NbAdj, pos, ListVar_Ctype[iProc],
                                           ListVar_NbAdj[iProc]);
      NbWrite[iProc] = pos + 1;
    };
    //
    auto InsertMatrix = [&](MyMatrix<Tint> const &M, int const &NbAdj) -> void {
      uint32_t seed = 0x1b873540;
      size_t e_hash = Matrix_Hash(M, seed);
      size_t iProc = e_hash % NprocOutput;
      AppendSingleCtype(M, NbAdj, iProc);
    };
    //
    auto InsertMatrixCan = [&](MyMatrix<Tint> const &M,
                               int const &NbAdj) -> void {
      if (M.rows() != static_cast<int>(n_vect) ||
          M.cols() != static_cast<int>(n)) {
        std::cerr << "We have |M|=" << M.rows() << " / " << M.cols() << "\n";
        std::cerr << "But n_vect=" << n_vect << " and n=" << n << "\n";
        throw TerminalException{1};
      }
      if (ApplyCanonicalization) {
        MyMatrix<Tint> Mcan =
            LinPolytopeAntipodalIntegral_CanonicForm<Tint>(M, std::cerr);
        InsertMatrix(Mcan, NbAdj);
      } else {
        InsertMatrix(M, NbAdj);
      }
    };
    //
    // Inserting data if it a text file
    //
    if (posNature == 0) {
      if (!IsExistingFile(InputFile)) {
        std::cerr << "The file InputFile=" << InputFile << " is missing\n";
        throw TerminalException{1};
      }
      std::ifstream is(InputFile);
      int nbType;
      is >> nbType;
      std::vector<size_t> ListHash(nbType);
      for (int iType = 0; iType < nbType; iType++) {
        std::cerr << "iType : " << iType << " / " << nbType << "\n";
        MyMatrix<Tint> eMat = ReadMatrix<Tint>(is);
        InsertMatrixCan(eMat, 0);
      }
    }
    //
    // Inserting data if it is a series of netcdf files.
    //
    if (posNature == 1) {
      for (size_t iProcI = 0; iProcI < NprocInput; iProcI++) {
        std::string FileI = PrefixInput + std::to_string(iProcI) + ".nc";
        if (!IsExistingFile(FileI)) {
          std::cerr << "The file FileI=" << FileI << " is missing\n";
          throw TerminalException{1};
        }
        netCDF::NcFile dataFile(FileI, netCDF::NcFile::read);
        if (dataFile.isNull()) {
          std::cerr << "Error while Netcdf opening of file=" << FileI << "\n";
          throw TerminalException{1};
        }
        netCDF::NcVar var_Ctype = dataFile.getVar("Ctype");
        if (var_Ctype.isNull()) {
          std::cerr << "Error while opening variable Ctype\n";
          throw TerminalException{1};
        }
        netCDF::NcVar var_NbAdj = dataFile.getVar("nb_adjacent");
        if (var_NbAdj.isNull()) {
          std::cerr << "Error while opening variable nb_adjacent\n";
          throw TerminalException{1};
        }
        int nbDim_Ctype = var_Ctype.getDimCount();
        int nbDim_NbAdj = var_NbAdj.getDimCount();
        if (nbDim_Ctype != 3 || nbDim_NbAdj != 1) {
          std::cerr
              << "Ctype should have dimension 3 and nb_adjacent dimension 1\n";
          throw TerminalException{1};
        }
        netCDF::NcDim eDim_nb_ctype = var_Ctype.getDim(0);
        int read_nb_ctype = eDim_nb_ctype.getSize();
        //
        netCDF::NcDim eDim_n_vect = var_Ctype.getDim(1);
        int read_n_vect = eDim_n_vect.getSize();
        if (read_n_vect != static_cast<int>(n_vect)) {
          std::cerr << "inconsistency in input of read_n_vect and n_vect\n";
          std::cerr << "read_n_vect=" << read_n_vect << " n_vect=" << n_vect
                    << "\n";
          throw TerminalException{1};
        }
        //
        netCDF::NcDim eDim_n = var_Ctype.getDim(2);
        int read_n = eDim_n.getSize();
        if (read_n != static_cast<int>(n)) {
          std::cerr << "inconsistency in input of read_n and n\n";
          std::cerr << "read_n=" << read_n << " n=" << n << "\n";
          throw TerminalException{1};
        }
        netCDF::NcType eType = var_Ctype.getType();
        //
        for (int iCtype = 0; iCtype < read_nb_ctype; iCtype++) {
          std::vector<size_t> start1{size_t(iCtype)};
          std::vector<size_t> count1{1};
          std::vector<size_t> start2{size_t(iCtype), 0, 0};
          std::vector<size_t> count2{1, n_vect, n};
          int NbAdj;
          var_NbAdj.getVar(start1, count1, &NbAdj);
          //
          MyMatrix<Tint> M(n_vect, n);
          if (eType == netCDF::NcType::nc_BYTE) {
            std::vector<int8_t> V(n_vect * n);
            var_Ctype.getVar(start2, count2, V.data());
            int idx = 0;
            for (size_t i_vect = 0; i_vect < n_vect; i_vect++)
              for (size_t i = 0; i < n; i++) {
                M(i_vect, i) = V[idx];
                idx++;
              }
          }
          if (eType == netCDF::NcType::nc_SHORT) {
            std::vector<short> V(n_vect * n);
            var_Ctype.getVar(start2, count2, V.data());
            int idx = 0;
            for (size_t i_vect = 0; i_vect < n_vect; i_vect++)
              for (size_t i = 0; i < n; i++) {
                M(i_vect, i) = V[idx];
                idx++;
              }
          }
          if (eType == netCDF::NcType::nc_INT) {
            std::vector<int> V(n_vect * n);
            var_Ctype.getVar(start2, count2, V.data());
            int idx = 0;
            for (size_t i_vect = 0; i_vect < n_vect; i_vect++)
              for (size_t i = 0; i < n; i++) {
                M(i_vect, i) = V[idx];
                idx++;
              }
          }
          if (eType == netCDF::NcType::nc_INT64) {
            std::vector<int64_t> V(n_vect * n);
            var_Ctype.getVar(start2, count2, V.data());
            int idx = 0;
            for (size_t i_vect = 0; i_vect < n_vect; i_vect++)
              for (size_t i = 0; i < n; i++) {
                M(i_vect, i) = V[idx];
                idx++;
              }
          }
          InsertMatrixCan(M, NbAdj);
        }
      }
    }
    std::cerr << "Normal termination of CTYP_PrepareInitialFile\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in CTYP_PrepareInitialFile\n";
    exit(e.eVal);
  }
  runtime(time);
}
