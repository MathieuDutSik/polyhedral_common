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

int main(int argc, char *argv[]) {
  try {
    if (argc != 4) {
      std::cerr << "CTYP_ComputeAverage [FileName] [VarName] [nCol]\n";
      throw TerminalException{1};
    }
    std::string FileName = argv[1];
    std::string VarName = argv[2];
    std::string nColStr = argv[3];
    int nCol = ParseScalar<int>(nColStr);
    //
    // Now reading the
    //
    if (!IsExistingFile(FileName)) {
      std::cerr << "FileName=" << FileName << " is missing\n";
      throw TerminalException{1};
    }
    //
    // The input file
    //
    std::cerr << "Reading FileName=" << FileName << "\n";
    netCDF::NcFile dataFileI(FileName, netCDF::NcFile::read);
    netCDF::NcVar var = dataFileI.getVar(VarName);
    size_t n_ctype = var.getDim(0).getSize();
    std::vector<int> V(n_ctype);
    var.getVar(V.data());
    //
    //
    int val_max = 0;
    int val_min = std::numeric_limits<int>::max();
    size_t sum = 0;
    std::map<int, size_t> map;
    for (size_t i = 0; i < n_ctype; i++) {
      int val = V[i];
      if (val < val_min) {
        val_min = val;
      }
      if (val > val_max) {
        val_max = val;
      }
      sum += static_cast<size_t>(val);
      map[val] += 1;
    }
    double sum_d = static_cast<double>(sum);
    double avg = sum_d / n_ctype;
    std::cerr << "val_min=" << val_min << " val_max=" << val_max
              << " avg=" << avg << "\n";

    std::vector<std::string> LStr;
    for (auto &kv : map) {
      std::string str =
          std::to_string(kv.first) + " & " + std::to_string(kv.second);
      LStr.push_back(str);
    }
    std::string strEmpty = "  &  ";
    std::string strSep = " & ";
    mpq_class nEnt = map.size();
    std::cerr << "nEnt=" << nEnt << "\n";
    mpq_class nCol_mpq = nCol;
    std::cerr << "nCol_mpq=" << nCol_mpq << "\n";
    mpq_class quot = nEnt / nCol_mpq;
    std::cerr << "quot=" << quot << "\n";
    mpq_class nRow_mpq = UniversalCeilScalarInteger<mpq_class, mpq_class>(quot);
    std::cerr << "nRow_mpq=" << nRow_mpq << "\n";
    int nRow = UniversalScalarConversion<int, mpq_class>(nRow_mpq);
    std::cerr << "nRow=" << nRow << "\n";
    std::vector<std::string> Lines(nRow);
    std::vector<size_t> LEnt(nRow, 0);
    int nTot = nRow * nCol;
    size_t pos = 0;
    for (int u = 0; u < nTot; u++) {
      if (LEnt[pos] > 0) {
        Lines[pos] += strSep;
      }
      if (pos < map.size()) {
        Lines[pos] += LStr[u];
      } else {
        Lines[pos] += strEmpty;
      }
      LEnt[pos] += 1;
      pos += 1;
      if (pos == static_cast<size_t>(nRow)) {
        pos = 0;
      }
    }
    for (int u = 0; u < nRow; u++) {
      std::cerr << Lines[u] << "\\\\\n";
    }
    std::cerr << "Normal termination of NC_ComputeAverage\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in NC_ComputeAverage\n";
    exit(e.eVal);
  }
}
