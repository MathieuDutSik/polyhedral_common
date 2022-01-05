#include "NumberTheory.h"
#include "CtypeMPI_types.h"
#include "Namelist.h"
#include <unordered_map>

#include <boost/mpi.hpp>
#include "hash_functions.h"
#include <netcdf>
namespace mpi = boost::mpi;


FullNamelist NAMELIST_Timings()
{
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, std::string> ListStringValues1;
  ListIntValues1["n"]=6;
  ListStringValues1["TestFile"] = "WORK_6_0.nc";
  ListIntValues1["MaxRunTimeSecond"]=120;
  SingleBlock BlockDATA;
  BlockDATA.ListIntValues = ListIntValues1;
  BlockDATA.ListStringValues = ListStringValues1;
  BlockDATA.ListBoolValues = ListBoolValues1;
  ListBlock["DATA"]=BlockDATA;
  // Merging all data
  return {ListBlock, "undefined"};
}





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


template<typename T, typename Tint>
void NC_WriteMatrix_T(netCDF::NcVar & varCtype, MyMatrix<int> const& M, size_t const& n_vect, size_t const& n, int const& pos)
{
  std::vector<size_t> start2{size_t(pos), 0, 0};
  std::vector<size_t> count2{1, n_vect, n};
  std::vector<T> V(n_vect * n);
  int idx=0;
  for (size_t i_vect=0; i_vect<n_vect; i_vect++)
    for (size_t i=0; i<n; i++) {
      V[idx] = M(i_vect, i);
      idx++;
    }
  varCtype.putVar(start2, count2, V.data());
}






int main(int argc, char* argv[])
{
  //
  using Tint=int;
  //
  // The input file
  //
  FullNamelist eFull = NAMELIST_Timings();
  if (argc != 2) {
    std::cerr << "Number of argument is = " << argc << "\n";
    std::cerr << "This program is used as\n";
    std::cerr << "CTYP_MPI_Enumeration_c [file.nml]\n";
    std::cerr << "With file.nml a namelist file\n";
    NAMELIST_WriteNamelistFile(std::cerr, eFull);
    return -1;
  }
  std::string eFileName=argv[1];
  NAMELIST_ReadNamelistFile(eFileName, eFull);
  //
  // Parsing the input file
  //
  SingleBlock BlDATA = eFull.ListBlock["DATA"];
  //  int n=BlDATA.ListIntValues.at("n");
  int n = BlDATA.ListIntValues.at("n");
  std::string TestFile = BlDATA.ListStringValues.at("TestFile");
  int MaxRunTimeSecond = BlDATA.ListIntValues.at("MaxRunTimeSecond");
  //
  // The basic sizes
  //
  int n_vect = std::pow(2, n) - 1;
  //
  // The netcdf interface
  //
  netCDF::NcFile dataFile(TestFile, netCDF::NcFile::write);
  netCDF::NcVar varCtype=dataFile.getVar("Ctype");
  int n_read = varCtype.getDim(2).getSize();
  if (n_read != n) {
    std::cerr << "n_read=" << n_read << " n=" << n << "\n";
    return 0;
  }
  netCDF::NcType eType=varCtype.getType();
  netCDF::NcVar varNbAdj=dataFile.getVar("nb_adjacent");
  size_t nb_ent = varNbAdj.getDim(0).getSize();
  auto NC_ReadMatrix=[&](int const& pos) -> MyMatrix<Tint> {
    MyMatrix<Tint> M(n_vect, n);
    if (eType == netCDF::NcType::nc_BYTE)
      NC_ReadMatrix_T<int8_t,Tint>(varCtype, M, n_vect, n, pos);
    if (eType == netCDF::NcType::nc_SHORT)
      NC_ReadMatrix_T<int16_t,Tint>(varCtype, M, n_vect, n, pos);
    if (eType == netCDF::NcType::nc_INT)
      NC_ReadMatrix_T<int32_t,Tint>(varCtype, M, n_vect, n, pos);
    if (eType == netCDF::NcType::nc_INT64)
      NC_ReadMatrix_T<int64_t,Tint>(varCtype, M, n_vect, n, pos);
    return M;
  };
  //
  std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
  int nbAccess=0;
  while(true) {
    size_t e_rand = rand();
    size_t pos = e_rand % nb_ent;
    MyMatrix<Tint> eMat = NC_ReadMatrix(pos);
    //
    std::chrono::time_point<std::chrono::system_clock> aft = std::chrono::system_clock::now();
    int elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(aft - start).count();
    if (elapsed_seconds > MaxRunTimeSecond)
      break;
    nbAccess++;
  }
  double AvgTimeAccess = ((double)MaxRunTimeSecond) * ((double)1000000) / ((double)nbAccess);
  std::cerr << "nbAccess=" << nbAccess << " AvgTimeAccess=" << AvgTimeAccess << "\n";
}
