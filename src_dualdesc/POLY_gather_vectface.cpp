// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include "POLY_RecursiveDualDesc_MPI.h"
// clang-format on



int main(int argc, char *argv[]) {
  boost::mpi::environment env(boost::mpi::threading::serialized);
  if (env.thread_level() < boost::mpi::threading::serialized) {
    env.abort(-1);
  }
  boost::mpi::communicator comm;
  int i_rank = comm.rank();
  int n_proc = comm.size();
  int i_proc_ret = 0;
  //
  int n = 10;
  int n_face = 5;
  vectface vf(n);
  for (int i_face=0; i_face<n_face; i_face++) {
    Face f = RandomFace(n);
    vf.push_back(f);
  }
  std::string FileLog =
    "log_" + std::to_string(n_proc) + "_" + std::to_string(i_rank);
  std::cerr << "We have moved. See the log in FileLog=" << FileLog << "\n";
  std::ofstream os(FileLog);

  os << "We have vf |vf|=" << vf.size() << " / " << vf.get_n() << "\n";
  vectface vf_tot = my_mpi_gather(comm, vf, i_proc_ret);
  os << "We have vf_tot |vf_tot|=" << vf_tot.size() << " / " << vf_tot.get_n() << " i_proc_ret=" << i_proc_ret << "\n";

}
