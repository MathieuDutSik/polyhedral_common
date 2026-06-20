// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
#include "IsoDelaunayDomains.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

/*
  Random-walk search for an iso-Delaunay (L-type) domain in which every
  extreme ray, viewed as a vector in the t-space LinSpa.ListMat, defines a
  Gram matrix of full rank. Modelled on
  src_ctype/CTYP_LookForNoFreeVector.cpp: the objective being minimised
  is the number of non-full-rank rays of the current domain, computed via
  IsoDelaunayDomains.h::CountNonFullRankRays (the same rank logic used by
  WriteDetailedEntryGAP).
 */

template <typename T, typename Tint>
void process_A(FullNamelist const &eFull, int max_s, int n_try,
               int n_walk_steps, int max_iter) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;

  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  SingleBlock const &BlockSYSTEM = eFull.get_block("SYSTEM");
  SingleBlock const &BlockTSPACE = eFull.get_block("TSPACE");
  LinSpaceMatrix<T> LinSpa =
      ReadTspace<T, Tint, Tgroup>(BlockTSPACE, std::cerr);
  int dimEXT = LinSpa.n + 1;
  //
  std::string Prefix = BlockSYSTEM.get_string("Prefix");
  CreateDirectory(Prefix);
  //
  std::string FileDualDesc = BlockDATA.get_string("FileDualDescription");
  PolyHeuristicSerial<TintGroup> AllArr =
      Read_AllStandardHeuristicSerial_File<T, TintGroup>(FileDualDesc, dimEXT,
                                                         std::cerr);
  DataIsoDelaunayDomains<T, Tint, Tgroup> data =
      get_data_isodelaunay_domains<T, Tint, Tgroup>(eFull, AllArr, std::cerr);
  //
  for (int i_try = 0; i_try < n_try; i_try++) {
    std::cerr << "LATT_LookForFullRankRayDomain, i_try=" << i_try << "/"
              << n_try << "\n";
    LookForFullRankRayDomain<T, Tint, Tgroup>(data, Prefix, max_s,
                                              n_walk_steps, max_iter,
                                              std::cerr);
  }
}

void process_C(FullNamelist const &eFull, int max_s, int n_try,
               int n_walk_steps, int max_iter) {
  std::string arithmetic =
      GetNamelistStringEntry(eFull, "DATA", "arithmetic");
  if (arithmetic == "gmp") {
    using T = mpq_class;
    using Tint = mpz_class;
    return process_A<T, Tint>(eFull, max_s, n_try, n_walk_steps, max_iter);
  }
  if (arithmetic == "gmp_boost") {
    using T = boost::multiprecision::mpq_rational;
    using Tint = boost::multiprecision::mpz_int;
    return process_A<T, Tint>(eFull, max_s, n_try, n_walk_steps, max_iter);
  }
  if (arithmetic == "multi_boost") {
    using T = boost::multiprecision::cpp_rational;
    using Tint = boost::multiprecision::cpp_int;
    return process_A<T, Tint>(eFull, max_s, n_try, n_walk_steps, max_iter);
  }
  if (arithmetic == "safe") {
    using T = Rational<SafeInt64>;
    using Tint = SafeInt64;
    return process_A<T, Tint>(eFull, max_s, n_try, n_walk_steps, max_iter);
  }
  std::cerr << "LATT_LookForFullRankRayDomain: Failed to find a matching "
               "entry for arithmetic="
            << arithmetic << "\n";
  std::cerr << "Available types: gmp, gmp_boost, multi_boost, safe\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    FullNamelist eFull =
        NAMELIST_GetStandard_COMPUTE_LATTICE_IsoDelaunayDomains();
    if (argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_LookForFullRankRayDomain [file.nml] [n_try] [max_s] "
                   "[n_walk_steps] [max_iter]\n";
      std::cerr << "\n";
      std::cerr << "with:\n";
      std::cerr << "file.nml      : namelist describing the t-space (see "
                   "LATT_SerialLattice_IsoDelaunayDomain). SYSTEM.Prefix is "
                   "used as the output directory.\n";
      std::cerr << "n_try         : the number of independent restarts\n";
      std::cerr << "max_s         : Gram matrices with at most this many "
                   "non-full-rank rays are written to disk. Use 0 to keep "
                   "only successful hits.\n";
      std::cerr << "n_walk_steps  : number of random adjacency jumps used "
                   "to escape a local minimum (e.g. 50)\n";
      std::cerr << "max_iter      : per-restart cap on outer loop iterations. "
                   "Guarantees finite runtime even when the target "
                   "curr_count=0 is unreachable (n<=5 in the Classic t-space).\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    unsigned seed = get_random_seed();
    std::cerr << "seed=" << seed << "\n";
    srand(seed);
    std::string eFileName = argv[1];
    int n_try = ParseScalar<int>(argv[2]);
    int max_s = ParseScalar<int>(argv[3]);
    int n_walk_steps = ParseScalar<int>(argv[4]);
    int max_iter = ParseScalar<int>(argv[5]);
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    process_C(eFull, max_s, n_try, n_walk_steps, max_iter);
    //
    std::cerr << "Normal termination of LATT_LookForFullRankRayDomain\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_LookForFullRankRayDomain\n";
    exit(e.eVal);
  }
  runtime(time);
}
