// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Group.h"
#include "Permutation.h"
#include "GenusEnumerationMPI.h"
#include <string>
#include <vector>
// clang-format on

/*
  GENUS_MPI_Enumerate: the enumeration of GENUS_Enumerate over MPI (see
  GenusEnumerationMPI.h), with the same namelist. The arithmetic is that of
  gmp. Every process reads the namelist and computes the same input; rank 0
  writes the result. An odd genus is enumerated through the genus of its
  even sublattices, each node then turning its even classes into odd ones.
 */

int main(int argc, char *argv[]) {
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator comm;
  maybe_install_gmp_pool();
  HumanTime time;
  int i_rank = comm.rank();
  // Only rank 0 reports; the other ranks write nothing.
  std::ostream null_stream(nullptr);
  std::ostream &os = (i_rank == 0) ? std::cerr : null_stream;
  try {
    using T = mpq_class;
    using Tint = mpz_class;
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt, mpz_class>;
    using Tmass = typename GenusEnumerationResult<T, Tgroup>::Tmass;
    FullNamelist eFull = NAMELIST_GetStandard_GENUS_ENUMERATE();
    if (argc != 2) {
      if (i_rank == 0) {
        std::cerr << "This program is used as\n";
        std::cerr << "mpirun -np [N] GENUS_MPI_Enumerate [file.nml]\n";
        std::cerr << "with file.nml the namelist of GENUS_Enumerate\n";
        eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      }
      return -1;
    }
    NAMELIST_ReadNamelistFile(argv[1], eFull);
    GenusEnumerationInput<T, Tgroup> input =
        ReadGenusEnumerationInput<T, Tint, Tgroup>(eFull, os);
    if (i_rank == 0) {
      std::cerr << "GENUS_MPI_Enumerate: " << comm.size() << " processes\n";
    }
    std::vector<MyMatrix<T>> ListSeedEnum = input.ListSeed;
    Tmass TotalMassEnum = input.TotalMass;
    if (!input.even) {
      // the genus of the even sublattices, of its own mass
      ListSeedEnum = {EvenSublatticeGram<T, Tint>(input.ListSeed[0])};
      TotalMassEnum = GenusMass<T, Tint>(ListSeedEnum[0]);
    }
    std::pair<GenusEnumerationResult<T, Tgroup>, bool> pair =
        GenusEnumerationMPI<T, Tint, Tgroup>(
            comm, ListSeedEnum, TotalMassEnum, input.prime, input.spec,
            input.max_runtime_second, os);
    GenusEnumerationResult<T, Tgroup> local = pair.first;
    bool complete = pair.second;
    if (!input.even) {
      // the odd classes over the even classes of this node
      local = GenusOddClassesOfEven<T, Tint, Tgroup>(
          input.ListSeed[0], input.TotalMass, local, input.spec.methods, os);
    }
    GenusEnumerationResult<T, Tgroup> result =
        GenusMPIGather<T, Tgroup>(comm, local, complete);
    if (i_rank == 0) {
      result.target_mass = input.TotalMass;
      result.complete = complete && result.accumulated_mass == input.TotalMass;
      int prime_out = (input.spec.strategy == NeighborStrategy::OutPrime)
                          ? input.prime
                          : 0;
      auto f = [&](std::ostream &os_out) -> void {
        WriteGenusEnumerationResult(os_out, result, input.OutFormat,
                                    prime_out);
      };
      FILE_PrintStderrStdoutFile(input.OutFile, f);
      std::cerr << "Normal termination of GENUS_MPI_Enumerate\n";
    }
  } catch (TerminalException const &e) {
    std::cerr << "Error in GENUS_MPI_Enumerate on rank " << i_rank << "\n";
    comm.abort(e.eVal);
  }
  if (i_rank == 0) {
    runtime(time);
  }
}
