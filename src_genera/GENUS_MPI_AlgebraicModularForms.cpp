// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Group.h"
#include "Permutation.h"
#include "AlgebraicModularFormsZ.h"
#include "MPI_functionality.h"
#include <boost/mpi.hpp>
#include <string>
#include <vector>
// clang-format on

/*
  GENUS_MPI_AlgebraicModularForms: the Hecke matrices of the algebraic
  modular forms of trivial weight of a genus of positive definite even
  lattices over Z, with or without a level structure at a prime N, from the
  list of its classes (AlgebraicModularFormsZ.h). The input is a namelist,
  the arithmetic that of gmp. With a single process it is the serial
  computation.

  The work is independent per class and per row, so it is split by class:
  the node i mod n_proc computes the automorphism group and the level
  orbits of the class i, then its rows at every prime. Between the two
  phases the orders (for the mass) and the level orbits (for the columns)
  are shared with every node. Every node builds the decision tree of the
  classes, which is deterministic; rank 0 assembles the matrices, checks
  them and writes the result.
 */

namespace {

bool AmfZIsPrime(int p) {
  if (p < 2) {
    return false;
  }
  for (int d = 2; d * d <= p; d++) {
    if (p % d == 0) {
      return false;
    }
  }
  return true;
}

} // namespace

int main(int argc, char *argv[]) {
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator comm;
  maybe_install_gmp_pool();
  HumanTime time;
  int i_rank = comm.rank();
  int n_proc = comm.size();
  // Only rank 0 reports; the other ranks write nothing.
  std::ostream null_stream(nullptr);
  std::ostream &os = (i_rank == 0) ? std::cerr : null_stream;
  try {
    using T = mpq_class;
    using Tint = mpz_class;
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using TintGroup = mpz_class;
    using Tgroup = permutalib::Group<Telt, TintGroup>;
    using Tmass = typename AmfZResult<TintGroup>::Tmass;
    FullNamelist eFull = NAMELIST_GetStandard_AMF_Z();
    if (argc != 2) {
      if (i_rank == 0) {
        std::cerr << "This program is used as\n";
        std::cerr << "mpirun -np [N] GENUS_MPI_AlgebraicModularForms "
                     "[file.nml]\n";
        std::cerr << "with file.nml a namelist file\n";
        eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      }
      return -1;
    }
    NAMELIST_ReadNamelistFile(argv[1], eFull);
    SingleBlock const &BlockDATA = eFull.get_block("DATA");
    SingleBlock const &BlockMETHOD = eFull.get_block("METHOD");
    SingleBlock const &BlockSYSTEM = eFull.get_block("SYSTEM");
    std::string ClassesFile = BlockDATA.get_string("ClassesFile");
    std::vector<int> ListPrime = BlockDATA.get_list_int("ListPrime");
    int level = BlockDATA.get_int("level");
    size_t max_level_lines =
        static_cast<size_t>(BlockDATA.get_int("max_level_lines"));
    std::string OutFormat = BlockSYSTEM.get_string("OutFormat");
    std::string OutFile = BlockSYSTEM.get_string("OutFile");
    GenusMethods methods;
    methods.automorphism = ParseLatticeMethodChoice(
        BlockMETHOD.get_string("automorphism_method"));
    methods.isomorphism = ParseLatticeMethodChoice(
        BlockMETHOD.get_string("isomorphism_method"));
    // The classes and the checks of the input.
    std::vector<MyMatrix<Tint>> ListGram = ReadListMatrixFile<Tint>(ClassesFile);
    size_t h = ListGram.size();
    if (h == 0) {
      std::cerr << "AMFZ: no class in " << ClassesFile << "\n";
      throw TerminalException{1};
    }
    int n = ListGram[0].rows();
    for (auto const &G : ListGram) {
      if (G.rows() != n || G.cols() != n || !IsSymmetricMatrix(G)) {
        std::cerr << "AMFZ: the classes are not symmetric matrices of the "
                     "same size\n";
        throw TerminalException{1};
      }
      for (int a = 0; a < n; a++) {
        if (ResidueOfInteger(G(a, a), 2) != 0) {
          std::cerr << "AMFZ: the classes must be even lattices\n";
          throw TerminalException{1};
        }
      }
    }
    MyMatrix<T> G0_T = UniversalMatrixConversion<T, Tint>(ListGram[0]);
    Tint det = UniversalScalarConversion<Tint, T>(DeterminantMat(G0_T));
    if (level != 1 && (!AmfZIsPrime(level) || ResidueOfInteger(det, level) == 0)) {
      std::cerr << "AMFZ: the level must be 1 or a prime not dividing the "
                   "determinant " << det << "\n";
      throw TerminalException{1};
    }
    for (int p : ListPrime) {
      if (!AmfZIsPrime(p) || ResidueOfInteger(det, p) == 0 || p == level) {
        std::cerr << "AMFZ: the Hecke prime " << p << " must be a prime not "
                  << "dividing the level and the determinant " << det
                  << "\n";
        throw TerminalException{1};
      }
    }
    if (i_rank == 0) {
      std::cerr << "GENUS_MPI_AlgebraicModularForms: " << n_proc
                << " processes, " << h << " classes of rank " << n
                << ", level " << level << "\n";
    }
    // Phase 1: the automorphism groups and the level orbits of the classes
    // of this node, then shared.
    std::vector<size_t> my_index;
    std::vector<TintGroup> my_order;
    std::vector<AmfZLevelData> my_level;
    std::map<size_t, std::vector<MyMatrix<Tint>>> my_gen;
    for (size_t i = static_cast<size_t>(i_rank); i < h;
         i += static_cast<size_t>(n_proc)) {
      MyMatrix<T> G_T = UniversalMatrixConversion<T, Tint>(ListGram[i]);
      LatticeAutInfo<T, Tint, Tgroup> info =
          GetLatticeAutInfo<T, Tint, Tgroup>(G_T, methods, os);
      my_index.push_back(i);
      my_order.push_back(info.order);
      my_level.push_back(AmfZLevelOrbits<Tint>(ListGram[i], info.ListGenMat,
                                               level, max_level_lines));
      my_gen[i] = std::move(info.ListGenMat);
    }
    std::vector<std::vector<size_t>> all_index;
    std::vector<std::vector<TintGroup>> all_order;
    std::vector<std::vector<AmfZLevelData>> all_level;
    boost::mpi::all_gather(comm, my_index, all_index);
    boost::mpi::all_gather(comm, my_order, all_order);
    boost::mpi::all_gather(comm, my_level, all_level);
    std::vector<TintGroup> ListAutOrder(h);
    std::vector<AmfZLevelData> ListLevel(h);
    for (size_t r = 0; r < all_index.size(); r++) {
      for (size_t k = 0; k < all_index[r].size(); k++) {
        ListAutOrder[all_index[r][k]] = all_order[r][k];
        ListLevel[all_index[r][k]] = std::move(all_level[r][k]);
      }
    }
    // The mass: the tree recognizes the classes of the list only, so the
    // list must be the whole genus.
    Tmass mass(0);
    for (auto const &order : ListAutOrder) {
      mass += Tmass(1) / UniversalScalarConversion<Tmass, TintGroup>(order);
    }
    Tmass mass_formula = GenusMass<T, Tint>(G0_T);
    if (mass != mass_formula) {
      std::cerr << "AMFZ: the classes have mass " << mass
                << " and the genus " << mass_formula
                << ": the list is not the whole genus\n";
      throw TerminalException{1};
    }
    if (i_rank == 0) {
      std::cerr << "GENUS_MPI_AlgebraicModularForms: mass " << mass
                << ", the list is the whole genus\n";
    }
    // Phase 2: the rows of the classes of this node, at every prime.
    InvariantFamily<Tint> family =
        get_invariant_family<Tint, Tgroup>(ListGram, os);
    // (i_prime, i, o, j, o2, value) and (i_prime, i, n_neighbor), flat
    std::vector<int64_t> my_entries;
    std::vector<int64_t> my_neighbor;
    for (size_t i : my_index) {
      for (size_t ip = 0; ip < ListPrime.size(); ip++) {
        AmfZRow row = AmfZComputeRow<T, Tint, Tgroup>(
            i, ListPrime[ip], family, my_gen[i], ListLevel, level, methods, os);
        for (auto const &[key, val] : row.entries) {
          my_entries.insert(my_entries.end(),
                            {static_cast<int64_t>(ip), static_cast<int64_t>(i),
                             static_cast<int64_t>(std::get<0>(key)),
                             static_cast<int64_t>(std::get<1>(key)),
                             static_cast<int64_t>(std::get<2>(key)), val});
        }
        my_neighbor.insert(my_neighbor.end(),
                           {static_cast<int64_t>(ip), static_cast<int64_t>(i),
                            static_cast<int64_t>(row.n_neighbor)});
      }
    }
    std::vector<std::vector<int64_t>> all_entries, all_neighbor;
    boost::mpi::gather(comm, my_entries, all_entries, 0);
    boost::mpi::gather(comm, my_neighbor, all_neighbor, 0);
    if (i_rank == 0) {
      AmfZResult<TintGroup> res;
      res.level = level;
      res.ListAutOrder = ListAutOrder;
      res.mass = mass;
      res.mass_formula = mass_formula;
      res.ListPrime = ListPrime;
      for (auto const &lv : ListLevel) {
        res.offset.push_back(res.n_state);
        res.ListNumberOrbit.push_back(lv.orbit_size.size());
        res.ListOrbitSize.push_back(lv.orbit_size);
        res.n_state += lv.orbit_size.size();
      }
      res.ListHecke.resize(ListPrime.size());
      // the number of neighbors depends on the genus only
      std::vector<int64_t> n_nb(ListPrime.size(), -1);
      for (auto const &vect : all_neighbor) {
        for (size_t k = 0; k < vect.size(); k += 3) {
          int64_t &ref = n_nb[vect[k]];
          if (ref != -1 && ref != vect[k + 2]) {
            std::cerr << "AMFZ: the classes do not have the same number of "
                         "neighbors at p=" << ListPrime[vect[k]] << "\n";
            throw TerminalException{1};
          }
          ref = vect[k + 2];
        }
      }
      for (auto const &vect : all_entries) {
        for (size_t k = 0; k < vect.size(); k += 6) {
          size_t row = res.offset[vect[k + 1]] + vect[k + 2];
          size_t col = res.offset[vect[k + 3]] + vect[k + 4];
          res.ListHecke[vect[k]][{row, col}] += vect[k + 5];
        }
      }
      // every state has all the neighbors of its class in its row
      for (size_t ip = 0; ip < ListPrime.size(); ip++) {
        res.ListNeighbor.push_back(static_cast<size_t>(n_nb[ip]));
        std::vector<int64_t> row_sum(res.n_state, 0);
        for (auto const &[key, val] : res.ListHecke[ip]) {
          row_sum[key.first] += val;
        }
        for (size_t s = 0; s < res.n_state; s++) {
          if (row_sum[s] != n_nb[ip]) {
            std::cerr << "AMFZ: the row of the state " << s << " at p="
                      << ListPrime[ip] << " sums to " << row_sum[s]
                      << " and not to the " << n_nb[ip] << " neighbors\n";
            throw TerminalException{1};
          }
        }
      }
      AmfZCheckResult(res);
      std::cerr << "GENUS_MPI_AlgebraicModularForms: " << res.n_state
                << " states\n";
      auto f = [&](std::ostream &os_out) -> void {
        WriteAmfZResult(os_out, res, OutFormat);
      };
      FILE_PrintStderrStdoutFile(OutFile, f);
      std::cerr << "Normal termination of GENUS_MPI_AlgebraicModularForms\n";
    }
  } catch (TerminalException const &e) {
    std::cerr << "Error in GENUS_MPI_AlgebraicModularForms on rank " << i_rank
              << "\n";
    comm.abort(e.eVal);
  }
  if (i_rank == 0) {
    runtime(time);
  }
}
