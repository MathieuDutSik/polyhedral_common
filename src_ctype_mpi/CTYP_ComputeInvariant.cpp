#include "CtypeMPI_types.h"
#include "NumberTheory.h"
#include "Namelist.h"
#include <boost/mpi.hpp>

#include <netcdf>

//#define DEBUG_INFO


template<typename T, typename Tint>
void NC_ReadMatrix_T(netCDF::NcVar & varCtype, MyMatrix<int> & M, size_t const& n_vect, size_t const& n, int const& pos)
{
  std::vector<size_t> start2{size_t(pos), 0, 0};
  std::vector<size_t> count2{1, n_vect, n};
  std::vector<T> V(n_vect * n);
  varCtype.getVar(start2, count2, V.data());
  int idx=0;
  for (size_t i_vect=0; i_vect<n_vect; i_vect++)
    for (size_t i=0; i<n; i++) {
      M(i_vect, i) = V[idx];
      idx++;
    }
}


int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  int irank_i;
  MPI_Comm_rank(MPI_COMM_WORLD,&irank_i);
  size_t irank=irank_i;
  //
  // Now reading the 
  //
  using Tint=int;
  std::string PrefixI=argv[1];
  std::string PrefixO=argv[2];
  std::string eFileI = PrefixI + std::to_string(irank) + ".nc";
  std::string eFileO = PrefixO + std::to_string(irank) + ".nc";
  if (!IsExistingFile(eFileI)) {
    std::cerr << "FileI=" << eFileI << " is missing\n";
    throw TerminalException{1};
  }
  if (IsExistingFile(eFileO)) {
    std::cerr << "FileO=" << eFileO << " should be missing\n";
    throw TerminalException{1};
  }
  //
  // The input file
  //
  std::cerr << "Reading eFileI=" << eFileI << "\n";
  netCDF::NcFile dataFileI(eFileI, netCDF::NcFile::read);
  netCDF::NcVar varCtypeI=dataFileI.getVar("Ctype");
  size_t n = varCtypeI.getDim(2).getSize();
  size_t n_vect = varCtypeI.getDim(1).getSize();
  size_t n_ctype = varCtypeI.getDim(0).getSize();
  netCDF::NcType eType=varCtypeI.getType();
  auto NC_ReadMatrix=[&](int const& pos) -> TypeCtypeExch<Tint> {
    MyMatrix<Tint> M(n_vect, n);
    if (eType == netCDF::NcType::nc_BYTE)
      NC_ReadMatrix_T<int8_t,Tint>(varCtypeI, M, n_vect, n, pos);
    if (eType == netCDF::NcType::nc_SHORT)
      NC_ReadMatrix_T<int16_t,Tint>(varCtypeI, M, n_vect, n, pos);
    if (eType == netCDF::NcType::nc_INT)
      NC_ReadMatrix_T<int32_t,Tint>(varCtypeI, M, n_vect, n, pos);
    if (eType == netCDF::NcType::nc_INT64)
      NC_ReadMatrix_T<int64_t,Tint>(varCtypeI, M, n_vect, n, pos);
    return {M};
  };
  netCDF::NcVar varNbAdjI=dataFileI.getVar("nb_adjacent");
  auto NC_GetNbAdjacent=[&](int const& pos) -> int {
    std::vector<size_t> start{size_t(pos)};
    std::vector<size_t> count{1};
    int nbAdjacent;
    varNbAdjI.getVar(start, count, &nbAdjacent);
    return nbAdjacent;
  };
  //
  // The output file
  //
  netCDF::NcFile dataFileO(eFileO, netCDF::NcFile::replace, netCDF::NcFile::nc4);
  netCDF::NcDim eDimNbCtype=dataFileO.addDim("number_ctype");
  netCDF::NcDim eDimN=dataFileO.addDim("n", n);
  netCDF::NcDim eDimNvect=dataFileO.addDim("n_vect", n_vect);
  //
  std::vector<std::string> LDim3{"number_ctype", "n_vect", "n"};
  std::vector<std::string> LDim1{"number_ctype"};
  //
  netCDF::NcVar varCtypeO = dataFileO.addVar("Ctype", "byte", LDim3);
  varCtypeO.putAtt("long_name", "Ctype canonicalized coordinates");
  //
  netCDF::NcVar varNbAdjO = dataFileO.addVar("nb_adjacent", "int", LDim1);
  varNbAdjO.putAtt("long_name", "number of adjacent Ctypes");
  //
  netCDF::NcVar varNbTripleO = dataFileO.addVar("nb_triple", "int", LDim1);
  varNbTripleO.putAtt("long_name", "number of triples");
  //
  netCDF::NcVar varNbIneqO = dataFileO.addVar("nb_ineq", "int", LDim1);
  varNbIneqO.putAtt("long_name", "number of inequalities from triples");
  //
  netCDF::NcVar varNbIneqAfterCritO = dataFileO.addVar("nb_ineq_after_crit", "int", LDim1);
  varNbIneqAfterCritO.putAtt("long_name", "number of inequalities after criterion reduction");
  //
  netCDF::NcVar varNbFreeO = dataFileO.addVar("nb_free", "int", LDim1);
  varNbFreeO.putAtt("long_name", "number of free vectors");
  //
  netCDF::NcVar varNbAutomO = dataFileO.addVar("nb_autom", "int", LDim1);
  varNbAutomO.putAtt("long_name", "number of automorphisms");
  //
  auto NC_WriteMatrix=[&](int const& pos, MyMatrix<Tint> const& eMat) -> void {
    std::vector<size_t> start2{size_t(pos), 0, 0};
    std::vector<size_t> count2{1, n_vect, n};
    std::vector<int8_t> V(n_vect * n);
    int idx=0;
    for (size_t i_vect=0; i_vect<n_vect; i_vect++)
      for (size_t i=0; i<n; i++) {
        V[idx] = eMat(i_vect, i);
        idx++;
      }
    varCtypeO.putVar(start2, count2, V.data());
  };
  auto NC_WriteEntry=[&](netCDF::NcVar & eVar, int const& pos, int value) -> void {
    std::vector<size_t> start{size_t(pos)};
    std::vector<size_t> count{1};
    eVar.putVar(start, count, &value);
  };
  //
  // Processing the data
  //
  std::string LogFileO="LOG_" + std::to_string(irank);
  std::ofstream log(LogFileO);
  log << "Beginning\n";
  for (size_t i_ctype=0; i_ctype<n_ctype; i_ctype++) {
#ifdef DEBUG_INFO
    std::cerr << "i_ctype=" << i_ctype << "\n";
#endif
    TypeCtypeExch<Tint> eType = NC_ReadMatrix(i_ctype);
#ifdef DEBUG_INFO
    std::cerr << "We have eType\n";
#endif
    int nb_adjacent = NC_GetNbAdjacent(i_ctype);
#ifdef DEBUG_INFO
    std::cerr << "We have nb_adjacent\n";
#endif
    //
    StructuralInfo info = CTYP_GetStructuralInfo(eType);
#ifdef DEBUG_INFO
    std::cerr << "We have info\n";
#endif
    //
    NC_WriteMatrix(i_ctype, eType.eMat);
    NC_WriteEntry(varNbAdjO, i_ctype, nb_adjacent);
    NC_WriteEntry(varNbTripleO, i_ctype, info.nb_triple);
    NC_WriteEntry(varNbIneqO, i_ctype, info.nb_ineq);
    NC_WriteEntry(varNbIneqAfterCritO, i_ctype, info.nb_ineq_after_crit);
    NC_WriteEntry(varNbFreeO, i_ctype, info.nb_free);
    NC_WriteEntry(varNbAutomO, i_ctype, info.nb_autom);
#ifdef DEBUG_INFO
    std::cerr << "Writes done\n";
    //
    size_t res = i_ctype % 1000;
    if (res == 0)
      std::cerr << "i_ctype=" << i_ctype << "/" << n_ctype << "\n";
#endif
    log << "i_ctype=" << i_ctype << "\n";
  }
  std::cerr << "Normal termination of the program\n";
  MPI_Finalize();
}
