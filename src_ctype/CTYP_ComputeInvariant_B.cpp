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

int main() {
  boost::mpi::environment env(boost::mpi::threading::serialized);
  if (env.thread_level() < boost::mpi::threading::serialized) {
    env.abort(-1);
  }
  boost::mpi::communicator world;
  int rank = world.rank();
  int size = world.size();
  std::cerr << "rank=" << rank << " size=" << size << "\n";
  size_t rank_s = static_cast<size_t>(rank);
  size_t size_s = static_cast<size_t>(size);
  //
  // Now reading the
  //
  using Tint = int64_t;
  std::string eFileI = "ctype_dim6.nc";
  std::string eFileO = "ctype_invariant6.nc";
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
  netCDF::NcVar varNbAdjI = dataFileI.getVar("nb_adjacent");
  auto NC_GetNbAdjacent = [&](int const &pos) -> int {
    std::vector<size_t> start{size_t(pos)};
    std::vector<size_t> count{1};
    int nbAdjacent;
    varNbAdjI.getVar(start, count, &nbAdjacent);
    return nbAdjacent;
  };
  std::vector<int> l_nb_adjacent, l_nb_triple, l_nb_ineq, l_nb_ineq_after_crit, l_nb_free, l_nb_autom;
  if (rank == 0) {
    l_nb_adjacent.resize(n_ctype);
    l_nb_triple.resize(n_ctype);
    l_nb_ineq.resize(n_ctype);
    l_nb_ineq_after_crit.resize(n_ctype);
    l_nb_free.resize(n_ctype);
    l_nb_autom.resize(n_ctype);
  }
  //
  // Processing the data
  //
  std::vector<size_t> l_idx;
  for (size_t i_ctype = 0; i_ctype < n_ctype; i_ctype++) {
    size_t res = i_ctype % size_s;
    if (res == rank_s) {
      l_idx.push_back(i_ctype);
    }
  }
  std::vector<int> loc_nb_adjacent, loc_nb_triple, loc_nb_ineq, loc_nb_ineq_after_crit, loc_nb_free, loc_nb_autom;
  std::string LogFileO = "LOG_" + std::to_string(rank);
  std::ofstream log(LogFileO);
  log << "Beginning\n";
  for (auto & i_ctype : l_idx) {
    log << "i_ctype=" << i_ctype << "\n";
    TypeCtypeExch<Tint> eType = NC_ReadMatrix(i_ctype);
    int nb_adjacent = NC_GetNbAdjacent(i_ctype);
    loc_nb_adjacent.push_back(nb_adjacent);
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    //    using Tint = mpz_class; // This is defined above
    using Tgroup = permutalib::Group<Telt, Tint>;
    StructuralInfo info =
        CTYP_GetStructuralInfo<Tint, Tgroup>(eType, std::cerr);
    //
    loc_nb_triple.push_back(info.nb_triple);
    loc_nb_ineq.push_back(info.nb_ineq);
    loc_nb_ineq_after_crit.push_back(info.nb_ineq_after_crit);
    loc_nb_free.push_back(info.nb_free);
    loc_nb_autom.push_back(info.nb_autom);
  }
  //
  // Collecting the data
  //
  if (rank == 0) {
    for (size_t u=0; u<l_idx.size(); u++) {
      size_t i_ctype = l_idx[u];
      l_nb_adjacent[i_ctype] = loc_nb_adjacent[u];
      l_nb_triple[i_ctype] = loc_nb_triple[u];
      l_nb_ineq[i_ctype] = loc_nb_ineq[u];
      l_nb_ineq_after_crit[i_ctype] = loc_nb_ineq_after_crit[u];
      l_nb_free[i_ctype] = loc_nb_free[u];
      l_nb_autom[i_ctype] = loc_nb_autom[u];
    }
    for (int irank=1; irank<size; irank++) {
      size_t irank_s = static_cast<size_t>(irank);
      world.recv(irank, 10, loc_nb_adjacent);
      world.recv(irank, 11, loc_nb_triple);
      world.recv(irank, 12, loc_nb_ineq);
      world.recv(irank, 13, loc_nb_ineq_after_crit);
      world.recv(irank, 14, loc_nb_free);
      world.recv(irank, 15, loc_nb_autom);
      size_t pos = 0;
      for (size_t i_ctype = 0; i_ctype < n_ctype; i_ctype++) {
        size_t res = i_ctype % size_s;
        if (res == irank_s) {
          l_nb_adjacent[i_ctype] = loc_nb_adjacent[pos];
          l_nb_triple[i_ctype] = loc_nb_triple[pos];
          l_nb_ineq[i_ctype] = loc_nb_ineq[pos];
          l_nb_ineq_after_crit[i_ctype] = loc_nb_ineq_after_crit[pos];
          l_nb_free[i_ctype] = loc_nb_free[pos];
          l_nb_autom[i_ctype] = loc_nb_autom[pos];
          pos += 1;
        }
      }
    }
  } else {
    world.send(0, 10, loc_nb_adjacent);
    world.send(0, 11, loc_nb_triple);
    world.send(0, 12, loc_nb_ineq);
    world.send(0, 13, loc_nb_ineq_after_crit);
    world.send(0, 14, loc_nb_free);
    world.send(0, 15, loc_nb_autom);
  }
  //
  // Writing the data
  //
  if (rank == 0) {
    netCDF::NcFile dataFileO(eFileO, netCDF::NcFile::replace,
                             netCDF::NcFile::nc4);
    netCDF::NcDim eDimNbCtype = dataFileO.addDim("number_ctype", n_ctype);
    netCDF::NcDim eDimN = dataFileO.addDim("n", n);
    netCDF::NcDim eDimNvect = dataFileO.addDim("n_vect", n_vect);
    //
    std::vector<std::string> LDim3{"number_ctype", "n_vect", "n"};
    std::vector<std::string> LDim1{"number_ctype"};
    //
    netCDF::NcVar varNbAdjO = dataFileO.addVar("nb_adjacent", "int", LDim1);
    varNbAdjO.putAtt("long_name", "number of adjacent Ctypes");
    varNbAdjO.putVar(l_nb_adjacent.data());
    //
    netCDF::NcVar varNbTripleO = dataFileO.addVar("nb_triple", "int", LDim1);
    varNbTripleO.putAtt("long_name", "number of triples");
    varNbTripleO.putVar(l_nb_triple.data());
    //
    netCDF::NcVar varNbIneqO = dataFileO.addVar("nb_ineq", "int", LDim1);
    varNbIneqO.putAtt("long_name", "number of inequalities from triples");
    varNbIneqO.putVar(l_nb_ineq.data());
    //
    netCDF::NcVar varNbIneqAfterCritO = dataFileO.addVar("nb_ineq_after_crit", "int", LDim1);
    varNbIneqAfterCritO.putAtt("long_name", "number of inequalities after criterion reduction");
    varNbIneqAfterCritO.putVar(l_nb_ineq_after_crit.data());
    //
    netCDF::NcVar varNbFreeO = dataFileO.addVar("nb_free", "int", LDim1);
    varNbFreeO.putAtt("long_name", "number of free vectors");
    varNbFreeO.putVar(l_nb_free.data());
    //
    netCDF::NcVar varNbAutomO = dataFileO.addVar("nb_autom", "int", LDim1);
    varNbAutomO.putAtt("long_name", "number of automorphisms");
    varNbAutomO.putVar(l_nb_autom.data());
  }
  //
  // Writing the data
  //
  if (rank == 0) {
    std::cerr << "Normal termination of the program\n";
  }
}
