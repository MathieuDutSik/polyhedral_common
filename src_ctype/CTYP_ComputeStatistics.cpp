// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "CtypeMPI_types.h"
#include "Namelist.h"
#include <boost/mpi.hpp>
#include <netcdf>
// clang-format on

int main(int argc, char *argv[]) {
  //
  // Now reading the
  //
  std::string Prefix = argv[1];
  std::vector<std::string> LKey = {"nb_adjacent", "nb_triple",
                                   "nb_ineq",     "nb_ineq_after_crit",
                                   "nb_free",     "nb_autom"};
  //
  std::map<std::string, std::map<int, int>> RecMap;
  for (auto &eKey : LKey)
    RecMap[eKey] = {};
  //
  size_t irank = 0;
  while (true) {
    std::string eFile = Prefix + std::to_string(irank) + ".nc";
    if (!IsExistingFile(eFile)) {
      std::cerr << "eFile=" << eFile << " is missing\n";
      break;
    }
    netCDF::NcFile dataFileI(eFile, netCDF::NcFile::read);
    for (auto &eKey : LKey) {
      std::cerr << "eKey=" << eKey << " irank=" << irank << "\n";
      netCDF::NcVar var = dataFileI.getVar(eKey);
      size_t n_ctype = var.getDim(0).getSize();
      std::vector<int> V(n_ctype);
      std::vector<size_t> start{0};
      std::vector<size_t> count{n_ctype};
      var.getVar(start, count, V.data());
      std::map<int, int> &eMap = RecMap[eKey];
      for (size_t pos = 0; pos < n_ctype; pos++) {
        int value = V[pos];
        eMap[value]++;
      }
    }
    irank++;
  }
  //
  // Printing the data
  //
  for (auto &eKey : LKey) {
    std::cerr << "eKey=" << eKey << "\n";
    std::map<int, int> &eMap = RecMap[eKey];
    for (auto &kv : eMap)
      std::cerr << "  " << kv.first << " : " << kv.second << "\n";
  }
}
