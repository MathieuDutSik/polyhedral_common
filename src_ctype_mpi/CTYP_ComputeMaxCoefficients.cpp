#include "NumberTheory.h"
#include "CtypeMPI_types.h"
#include "Namelist.h"
#include <netcdf>


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
  using Tint=int;
  //
  try {
    std::string Prefix=argv[1];
    //
    // Parsing the input
    //
    size_t iProc = 0;
    Tint MaxCoeff=0;
    while(true) {
      std::string eFile = Prefix + std::to_string(iProc) + ".nc";
      if (!IsExistingFile(eFile)) {
        std::cerr << "File=" << eFile << " is missing\n";
        throw TerminalException{1};
      }
      std::cerr << "Reading eFile=" << eFile << "\n";
      netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
      netCDF::NcVar varCtype=dataFile.getVar("Ctype");
      int n = varCtype.getDim(2).getSize();
      int n_vect = varCtype.getDim(1).getSize();
      size_t n_ctype = varCtype.getDim(0).getSize();
      netCDF::NcType eType=varCtype.getType();
      MyMatrix<Tint> M(n_vect, n);
      auto NC_ReadMatrix=[&](int const& pos) -> void {
        if (eType == netCDF::NcType::nc_BYTE)
          NC_ReadMatrix_T<int8_t,Tint>(varCtype, M, n_vect, n, pos);
        if (eType == netCDF::NcType::nc_SHORT)
          NC_ReadMatrix_T<int16_t,Tint>(varCtype, M, n_vect, n, pos);
        if (eType == netCDF::NcType::nc_INT)
          NC_ReadMatrix_T<int32_t,Tint>(varCtype, M, n_vect, n, pos);
        if (eType == netCDF::NcType::nc_INT64)
          NC_ReadMatrix_T<int64_t,Tint>(varCtype, M, n_vect, n, pos);
      };

      for (size_t i_ctype=0; i_ctype<n_ctype; i_ctype++) {
        NC_ReadMatrix(i_ctype);
        for (int i_vect=0; i_vect<n_vect; i_vect++)
          for (int i=0; i<n; i++) {
            Tint eVal = T_abs(M(i_vect,i));
            MaxCoeff = std::max(MaxCoeff, eVal);
          }
        size_t res = i_ctype % 1000;
        if (res == 0)
          std::cerr << "iProc=" << iProc << " i_ctype=" << i_ctype << "/" << n_ctype << " MaxCoeff=" << MaxCoeff << "\n";
      }
      iProc++;
    }
    std::cerr << "MaxCoeff=" << MaxCoeff << "\n";

  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
