// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRational.h"
#include "CtypeMPI_types.h"
#include "Namelist.h"
#include <boost/mpi.hpp>
#include <netcdf>
#include "Permutation.h"
#include "Group.h"
// clang-format on

template <typename T, typename Tint>
void NC_ReadMatrix_T(netCDF::NcVar &varCtype, MyMatrix<Tint> &M,
                     size_t const &n_vect, size_t const &n, int const &pos) {
  std::vector<size_t> start2{size_t(pos), 0, 0};
  std::vector<size_t> count2{1, n_vect, n};
  std::vector<T> V(n_vect * n);
  varCtype.getVar(start2, count2, V.data());
  int idx = 0;
  for (size_t i_vect = 0; i_vect < n_vect; i_vect++)
    for (size_t i = 0; i < n; i++) {
      M(i_vect, i) = V[idx];
      idx++;
    }
}

int main(int argc, char *argv[]) {
  try {
    if (argc != 2) {
      std::cerr << "CTYP_ComputeAverage [FileName]\n";
      throw TerminalException{1};
    }
    std::string FileName = argv[1];
    //
    // Now reading the
    //
    using Tint = int64_t;
    if (!IsExistingFile(FileName)) {
      std::cerr << "FileName=" << FileName << " is missing\n";
      throw TerminalException{1};
    }
    //
    // The input file
    //
    std::cerr << "Reading FileName=" << FileName << "\n";
    netCDF::NcFile dataFileI(FileName, netCDF::NcFile::read);
    netCDF::NcVar varCtypeI = dataFileI.getVar("Ctype");
    size_t n = varCtypeI.getDim(2).getSize();
    size_t n_vect = varCtypeI.getDim(1).getSize();
    size_t n_ctype = varCtypeI.getDim(0).getSize();
    netCDF::NcType ncType = varCtypeI.getType();
    auto NC_ReadMatrix = [&](int const &pos) -> TypeCtypeExch<Tint> {
      MyMatrix<Tint> M(n_vect, n);
      if (ncType == netCDF::NcType::nc_BYTE)
        NC_ReadMatrix_T<int8_t, Tint>(varCtypeI, M, n_vect, n, pos);
      if (ncType == netCDF::NcType::nc_SHORT)
        NC_ReadMatrix_T<int16_t, Tint>(varCtypeI, M, n_vect, n, pos);
      if (ncType == netCDF::NcType::nc_INT)
        NC_ReadMatrix_T<int32_t, Tint>(varCtypeI, M, n_vect, n, pos);
      if (ncType == netCDF::NcType::nc_INT64)
        NC_ReadMatrix_T<int64_t, Tint>(varCtypeI, M, n_vect, n, pos);
      return {M};
    };
    std::unordered_map<size_t, size_t> map1;
    size_t red_print = 1000;
    size_t pos = 0;
    size_t n_block = 0;
    for (size_t u = 0; u < n_ctype; u++) {
      TypeCtypeExch<Tint> eCtype = NC_ReadMatrix(u);
      size_t e_hash = std::hash<TypeCtypeExch<Tint>>()(eCtype);
      map1[e_hash] += 1;
      pos += 1;
      if (pos == red_print) {
        n_block += 1;
        std::cerr << "n_block=" << n_block << "\n";
        pos = 0;
      }
    }
    std::map<size_t, size_t> map2;
    for (auto &kv1 : map1) {
      map2[kv1.second] += 1;
    }
    for (auto &kv2 : map2) {
      std::cerr << "multiplicity=" << kv2.first << " attained " << kv2.second
                << " times\n";
    }
    std::cerr << "Normal termination in CTYP_ComputeHashStat\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in CTYP_ComputeHashStat\n";
    exit(e.eVal);
  }
}
