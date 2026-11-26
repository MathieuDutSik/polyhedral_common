// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Group.h"
#include "Permutation.h"
#ifdef ENABLE_FLINT_SUPPORT
# include "NumberTheoryFlint.h"
#endif
#include "AlgebraicModularFormsHerm.h"
#include "MPI_functionality.h"
#include <boost/mpi.hpp>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
// clang-format on

/*
  GENUS_MPI_AlgebraicModularForms: the Hecke matrices of the algebraic
  modular forms of trivial weight of a genus of positive definite even
  lattices over Z (AlgebraicModularFormsZ.h), or of positive definite
  Hermitian lattices over a CM field (Ring = Hermitian,
  AlgebraicModularFormsHerm.h), with or without a level structure at a
  prime, from the list of its classes. The input is a namelist, the
  arithmetic that of gmp. With a single process it is the serial
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

using T = mpq_class;
using Tint = mpz_class;
using Tidx = uint32_t;
using Telt = permutalib::SingleSidedPerm<Tidx>;
using TintGroup = mpz_class;
using Tgroup = permutalib::Group<Telt, TintGroup>;
using Tmass = typename AmfZResult<TintGroup>::Tmass;

// The classes of a genus as written by GENUS_HermitianEnumerate (CPP).
template <typename T> struct HermClasses {
  NumberField<T> K;
  int m;
  std::vector<MyMatrix<T>> h;
  std::vector<MyMatrix<T>> l_basis;
};

template <typename T>
HermClasses<T> ReadHermClassesFile(std::string const &file) {
  std::ifstream is(file);
  if (!is.good()) {
    std::cerr << "AMFH: the file " << file << " cannot be read\n";
    throw TerminalException{1};
  }
  HermClasses<T> res;
  res.K.d = 0;
  is >> res.K.deg;
  for (int k = 0; k <= res.K.deg; k++) {
    int64_t c;
    is >> c;
    res.K.Phi.push_back(c);
  }
  res.K.conj = ElemE<T>(res.K.deg, T(0));
  for (int k = 0; k < res.K.deg; k++) {
    std::string tok;
    is >> tok;
    res.K.conj[k] = ParseScalar<T>(tok);
  }
  is >> res.m;
  for (int k = 0; k < res.K.deg; k++) {
    res.h.push_back(ReadMatrix<T>(is));
  }
  size_t n_class;
  is >> n_class;
  for (size_t i = 0; i < n_class; i++) {
    res.l_basis.push_back(ReadMatrix<T>(is));
  }
  if (!is) {
    std::cerr << "AMFH: the file " << file << " ends too early\n";
    throw TerminalException{1};
  }
  return res;
}

// The data of phase 1, the same on every node after the exchange.
struct AmfPhase1 {
  std::vector<size_t> my_index;
  std::vector<TintGroup> ListAutOrder;
  std::vector<AmfZLevelData> ListLevel;
  std::map<size_t, std::vector<MyMatrix<Tint>>> my_gen;
};

// The orders and level orbits of the classes of this node, shared.
void AmfShare(boost::mpi::communicator &comm, size_t h,
              std::vector<TintGroup> const &my_order,
              std::vector<AmfZLevelData> const &my_level, AmfPhase1 &ph) {
  std::vector<std::vector<size_t>> all_index;
  std::vector<std::vector<TintGroup>> all_order;
  std::vector<std::vector<AmfZLevelData>> all_level;
  boost::mpi::all_gather(comm, ph.my_index, all_index);
  boost::mpi::all_gather(comm, my_order, all_order);
  boost::mpi::all_gather(comm, my_level, all_level);
  ph.ListAutOrder.resize(h);
  ph.ListLevel.resize(h);
  for (size_t r = 0; r < all_index.size(); r++) {
    for (size_t k = 0; k < all_index[r].size(); k++) {
      ph.ListAutOrder[all_index[r][k]] = all_order[r][k];
      ph.ListLevel[all_index[r][k]] = std::move(all_level[r][k]);
    }
  }
}

/*
  Phase 2, the end: the rows of every node gathered on rank 0, assembled,
  checked and written. res has the level, the orders, the masses and the
  operators (ListPrime, and ListPrimeName, ListAdjoint over a CM field).
 */
void AmfAssembleAndWrite(boost::mpi::communicator &comm,
                         AmfZResult<TintGroup> &res,
                         std::vector<AmfZLevelData> const &ListLevel,
                         std::vector<int64_t> const &my_entries,
                         std::vector<int64_t> const &my_neighbor,
                         std::string const &OutFormat,
                         std::string const &OutFile) {
  std::vector<std::vector<int64_t>> all_entries, all_neighbor;
  boost::mpi::gather(comm, my_entries, all_entries, 0);
  boost::mpi::gather(comm, my_neighbor, all_neighbor, 0);
  if (comm.rank() != 0) {
    return;
  }
  size_t n_op = res.ListPrime.size();
  for (auto const &lv : ListLevel) {
    res.offset.push_back(res.n_state);
    res.ListNumberOrbit.push_back(lv.orbit_size.size());
    res.ListOrbitSize.push_back(lv.orbit_size);
    res.n_state += lv.orbit_size.size();
  }
  res.ListHecke.resize(n_op);
  // the number of neighbors depends on the genus only
  std::vector<int64_t> n_nb(n_op, -1);
  for (auto const &vect : all_neighbor) {
    for (size_t k = 0; k < vect.size(); k += 3) {
      int64_t &ref = n_nb[vect[k]];
      if (ref != -1 && ref != vect[k + 2]) {
        std::cerr << "AMF: the classes do not have the same number of "
                     "neighbors for the operator " << vect[k] << "\n";
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
  for (size_t ip = 0; ip < n_op; ip++) {
    res.ListNeighbor.push_back(static_cast<size_t>(n_nb[ip]));
    std::vector<int64_t> row_sum(res.n_state, 0);
    for (auto const &[key, val] : res.ListHecke[ip]) {
      row_sum[key.first] += val;
    }
    for (size_t s = 0; s < res.n_state; s++) {
      if (row_sum[s] != n_nb[ip]) {
        std::cerr << "AMF: the row of the state " << s << " for the operator "
                  << ip << " sums to " << row_sum[s] << " and not to the "
                  << n_nb[ip] << " neighbors\n";
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
}

// One row in the flat arrays of AmfAssembleAndWrite.
void AmfAppendRow(size_t ip, size_t i, AmfZRow const &row,
                  std::vector<int64_t> &my_entries,
                  std::vector<int64_t> &my_neighbor) {
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

Tmass AmfMassOfOrders(std::vector<TintGroup> const &ListAutOrder) {
  Tmass mass(0);
  for (auto const &order : ListAutOrder) {
    mass += Tmass(1) / UniversalScalarConversion<Tmass, TintGroup>(order);
  }
  return mass;
}

void RunAmfZ(boost::mpi::communicator &comm, FullNamelist const &eFull,
             std::ostream &os) {
  int i_rank = comm.rank();
  int n_proc = comm.size();
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
                << "dividing the level and the determinant " << det << "\n";
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
  AmfPhase1 ph;
  std::vector<TintGroup> my_order;
  std::vector<AmfZLevelData> my_level;
  for (size_t i = static_cast<size_t>(i_rank); i < h;
       i += static_cast<size_t>(n_proc)) {
    MyMatrix<T> G_T = UniversalMatrixConversion<T, Tint>(ListGram[i]);
    LatticeAutInfo<T, Tint, Tgroup> info =
        GetLatticeAutInfo<T, Tint, Tgroup>(G_T, methods, os);
    ph.my_index.push_back(i);
    my_order.push_back(info.order);
    my_level.push_back(AmfZLevelOrbits<Tint>(ListGram[i], info.ListGenMat,
                                             level, max_level_lines));
    ph.my_gen[i] = std::move(info.ListGenMat);
  }
  AmfShare(comm, h, my_order, my_level, ph);
  // The mass: the tree recognizes the classes of the list only, so the
  // list must be the whole genus.
  Tmass mass = AmfMassOfOrders(ph.ListAutOrder);
  Tmass mass_formula = GenusMass<T, Tint>(G0_T);
  if (mass != mass_formula) {
    std::cerr << "AMFZ: the classes have mass " << mass << " and the genus "
              << mass_formula << ": the list is not the whole genus\n";
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
  for (size_t i : ph.my_index) {
    for (size_t ip = 0; ip < ListPrime.size(); ip++) {
      AmfZRow row = AmfZComputeRow<T, Tint, Tgroup>(
          i, ListPrime[ip], family, ph.my_gen[i], ph.ListLevel, level,
          methods, os);
      AmfAppendRow(ip, i, row, my_entries, my_neighbor);
    }
  }
  AmfZResult<TintGroup> res;
  res.level = level;
  res.ListAutOrder = ph.ListAutOrder;
  res.mass = mass;
  res.mass_formula = mass_formula;
  res.ListPrime = ListPrime;
  AmfAssembleAndWrite(comm, res, ph.ListLevel, my_entries, my_neighbor,
                      OutFormat, OutFile);
}

void RunAmfHerm(boost::mpi::communicator &comm, FullNamelist const &eFull,
                std::ostream &os) {
  int i_rank = comm.rank();
  int n_proc = comm.size();
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  SingleBlock const &BlockSYSTEM = eFull.get_block("SYSTEM");
  std::string ClassesFile = BlockDATA.get_string("ClassesFile");
  std::vector<int> ListPrime = BlockDATA.get_list_int("ListPrime");
  std::string LevelIdeal = BlockDATA.get_string("LevelIdeal");
  size_t max_level_lines =
      static_cast<size_t>(BlockDATA.get_int("max_level_lines"));
  std::string OutFormat = BlockSYSTEM.get_string("OutFormat");
  std::string OutFile = BlockSYSTEM.get_string("OutFile");
  HermClasses<T> cls = ReadHermClassesFile<T>(ClassesFile);
  size_t h = cls.l_basis.size();
  if (h == 0) {
    std::cerr << "AMFH: no class in " << ClassesFile << "\n";
    throw TerminalException{1};
  }
  // every node builds the same order (the computation is deterministic)
  CMOrder<T, Tint> ord = CMOrderOfField<T, Tint>(cls.K, true, os);
  RealifiedRing<T, Tint> ring = CMOrderRealifiedRing<T, Tint>(ord, cls.m, cls.h);
  HermLatticeData<T, Tint> L0 = HermLatticeOfBasis(ord, ring, cls.l_basis[0]);
  auto check_prime = [&](int p, std::string const &what) -> void {
    if (!AmfZIsPrime(p) || !HermIsUnimodularAt(L0, p)) {
      std::cerr << "AMFH: the " << what << " " << p << " must be a prime "
                << "where the genus is unimodular\n";
      throw TerminalException{1};
    }
  };
  // The level ideal N of E: the factors (ell,k) or (ell,k)^e, the k-th prime
  // above ell as HermIdealPrimesAbove lists them (their names).
  std::optional<HermLevel<Tint>> N;
  std::vector<int64_t> level_ell;
  if (LevelIdeal != "1") {
    std::vector<HermIdealPrime<Tint>> l_fact;
    std::vector<int> l_exp;
    for (std::string const &tok : STRING_Split(LevelIdeal, " ")) {
      int64_t ell = 0, k = 0, e = 1;
      char c1 = 0, c2 = 0, c3 = 0;
      std::istringstream is(tok);
      is >> c1 >> ell >> c2 >> k >> c3;
      if (!is || c1 != '(' || c2 != ',' || c3 != ')') {
        std::cerr << "AMFH: the factor " << tok << " of LevelIdeal is not of "
                     "the form (ell,k) or (ell,k)^e\n";
        throw TerminalException{1};
      }
      char c4 = 0;
      if (is >> c4) {
        if (c4 != '^' || !(is >> e) || e < 1) {
          std::cerr << "AMFH: the exponent of " << tok << " is not valid\n";
          throw TerminalException{1};
        }
      }
      check_prime(static_cast<int>(ell), "level prime");
      std::vector<HermIdealPrime<Tint>> l_prime = HermIdealPrimesAbove(ord, ell);
      if (i_rank == 0) {
        for (auto const &P : l_prime) {
          std::cerr << "GENUS_MPI_AlgebraicModularForms: above " << ell
                    << ", the prime " << P.name << " of norm " << P.q()
                    << "\n";
        }
      }
      if (k < 0 || static_cast<size_t>(k) >= l_prime.size()) {
        std::cerr << "AMFH: there is no prime " << tok << " above " << ell
                  << "\n";
        throw TerminalException{1};
      }
      l_fact.push_back(l_prime[k]);
      l_exp.push_back(static_cast<int>(e));
      level_ell.push_back(ell);
    }
    N = HermGetLevel<T, Tint>(ord, l_fact, l_exp);
  }
  // The Hecke operators: every prime of E above the primes of ListPrime.
  std::vector<HermNeighborPrime<T, Tint>> l_np;
  std::vector<int> ListPrimeOp;
  std::vector<std::string> ListPrimeName;
  std::vector<size_t> ListAdjoint;
  for (int p : ListPrime) {
    check_prime(p, "Hecke prime");
    for (int64_t ell : level_ell) {
      if (p == ell) {
        std::cerr << "AMFH: the Hecke prime " << p << " divides the level\n";
        throw TerminalException{1};
      }
    }
    std::vector<HermIdealPrime<Tint>> l_prime = HermIdealPrimesAbove(ord, p);
    size_t first = l_np.size();
    for (size_t k = 0; k < l_prime.size(); k++) {
      l_np.push_back(HermGetNeighborPrime(ord, l_prime, k));
      ListPrimeOp.push_back(p);
      ListPrimeName.push_back(l_prime[k].name);
      ListAdjoint.push_back(first + l_prime[k].conj_index);
    }
  }
  if (i_rank == 0) {
    std::cerr << "GENUS_MPI_AlgebraicModularForms: " << n_proc
              << " processes, " << h << " Hermitian classes of rank "
              << cls.m << " over a field of degree " << cls.K.deg << ", "
              << l_np.size() << " Hecke operators, level "
              << (N ? N->name : std::string("1")) << "\n";
  }
  // The line spaces of the level, for every class (needed for the columns).
  std::vector<HermLatticeData<T, Tint>> l_lat;
  std::vector<HermLevelSpace> level_ls;
  for (auto const &A : cls.l_basis) {
    l_lat.push_back(HermLatticeOfBasis(ord, ring, A));
    if (N) {
      level_ls.push_back(HermGetLevelSpace(l_lat.back(), *N));
      if (level_ls.back().n_line() > max_level_lines) {
        std::cerr << "AMFH: " << level_ls.back().n_line()
                  << " lines modulo the level, above max_level_lines\n";
        throw TerminalException{1};
      }
    }
  }
  // Phase 1: the unitary groups and the level orbits of the classes of
  // this node, then shared.
  AmfPhase1 ph;
  std::vector<TintGroup> my_order;
  std::vector<AmfZLevelData> my_level;
  for (size_t i = static_cast<size_t>(i_rank); i < h;
       i += static_cast<size_t>(n_proc)) {
    IntAutomorphismGroup<Tint, Tgroup> aut =
        ArithmeticAutomorphismGroupMultiple<T, Tint, Tgroup>(
            HermPairOfBasis(ring, cls.l_basis[i]), os);
    ph.my_index.push_back(i);
    my_order.push_back(aut.order);
    if (N) {
      my_level.push_back(
          HermAmfLevelOrbits(l_lat[i], level_ls[i], *N, aut.ListGen));
      if (my_level.back().orbit_size.empty()) {
        std::cerr << "AMFH: the class " << i << " has no level line modulo "
                  << N->name << "\n";
        throw TerminalException{1};
      }
    } else {
      my_level.push_back(AmfZLevelData{{0}, {1}});
    }
    ph.my_gen[i] = aut.ListGen;
  }
  AmfShare(comm, h, my_order, my_level, ph);
  // The mass of the list against that of the genus (HermitianGenusMass.h).
  Tmass mass = AmfMassOfOrders(ph.ListAutOrder);
  Tmass mass_formula = HermMassAnyField<T, Tint, Tgroup>(
      ord, ring, cls.l_basis[0], os);
  if (mass != mass_formula) {
    std::cerr << "AMFH: the classes have mass " << mass << " and the genus "
              << mass_formula << ": the list is not the whole genus\n";
    throw TerminalException{1};
  }
  if (i_rank == 0) {
    std::cerr << "GENUS_MPI_AlgebraicModularForms: mass " << mass
              << ", the list is the whole genus\n";
  }
  // Phase 2: the rows of the classes of this node, for every operator. The
  // classes are separated by the pairs (G, G Omega^T).
  std::vector<MyMatrix<Tint>> ListGram, ListSecond;
  for (auto const &A : cls.l_basis) {
    std::pair<MyMatrix<Tint>, MyMatrix<Tint>> GS = HermPairOfBasisInt(ring, A);
    ListGram.push_back(GS.first);
    ListSecond.push_back(GS.second);
  }
  InvariantFamily<Tint> family =
      get_invariant_family<Tint, Tgroup>(ListGram, os, 20000, ListSecond);
  std::vector<int64_t> my_entries;
  std::vector<int64_t> my_neighbor;
  for (size_t i : ph.my_index) {
    for (size_t ip = 0; ip < l_np.size(); ip++) {
      AmfZRow row = HermAmfComputeRow<T, Tint, Tgroup>(
          ord, ring, cls.l_basis, i, l_np[ip], family, ph.my_gen[i], N,
          level_ls, ph.ListLevel, os);
      AmfAppendRow(ip, i, row, my_entries, my_neighbor);
    }
  }
  AmfZResult<TintGroup> res;
  // the norm of the level ideal, and its factors
  res.level = N ? static_cast<int>(N->norm) : 1;
  res.level_name = N ? N->name : std::string("1");
  res.ListAutOrder = ph.ListAutOrder;
  res.mass = mass;
  res.mass_formula = mass_formula;
  res.ListPrime = ListPrimeOp;
  res.ListPrimeName = ListPrimeName;
  res.ListAdjoint = ListAdjoint;
  AmfAssembleAndWrite(comm, res, ph.ListLevel, my_entries, my_neighbor,
                      OutFormat, OutFile);
}

} // namespace

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
    std::string Ring = eFull.get_block("DATA").get_string("Ring");
    if (Ring == "Z") {
      RunAmfZ(comm, eFull, os);
    } else if (Ring == "Hermitian") {
      RunAmfHerm(comm, eFull, os);
    } else {
      std::cerr << "GENUS_MPI_AlgebraicModularForms: Ring is Z or "
                   "Hermitian, not " << Ring << "\n";
      throw TerminalException{1};
    }
    if (i_rank == 0) {
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
