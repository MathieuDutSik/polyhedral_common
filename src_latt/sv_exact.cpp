/* sv.c  simple driver for shvec                              */
/* Version July 11, 2005                                      */
/* Copyright: Frank Vallentin 2005, frank.vallentin@gmail.com */

// clang-format off
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
# include "NumberTheoryBoostGmpInt.h"
#else
# include "NumberTheory.h"
#endif
#include "ShortestUniversal.h"
#include "Shvec_exact.h"
// clang-format on

[[noreturn]] void die_sv(std::string const &last_words) {
  std::cout << "sv.c: " << last_words << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
    using T = boost::multiprecision::mpq_rational;
    using Tint = boost::multiprecision::mpz_int;
#else
    using T = mpq_class;
    using Tint = mpz_class;
#endif

    bool coset = false;
    bool NeedBound = false;
    int i, j, mode;
    T bound = 0;
    mode = TempShvec_globals::TEMP_SHVEC_MODE_UNDEF;
    int c;
    while ((c = getopt(argc, argv, "hb:s:t:mMTcgevl")) != -1)
      switch (c) {
      case 'h':
        printf("Usage: sv [options] <file\n");
        printf("-h  show this help\n");
        printf("-bN compute vectors v with (v, v) <= N\n");
        printf("-tN compute the first N coefficients of the theta-series\n");
        printf("-l determine the vectors with (v-c, v-c) <= N with N defined "
               "later\n");
        printf("-v determine the vectors with (v-c, v-c) = N with N defined "
               "later\n");
        printf("-m  determine the minimum\n");
        printf("-M  compute minimal vectors\n");
        printf(
            "-T  compute one vector of fixed length (quaetsion of Han Tran)\n");
        printf("-c  find shortest vectors in a coset\n");
        printf("-e  do some additional checks\n");
        return 0;
      case 'b':
        mode = TempShvec_globals::TEMP_SHVEC_MODE_BOUND;
        NeedBound = 1;
        break;
      case 'M':
        mode = TempShvec_globals::TEMP_SHVEC_MODE_SHORTEST_VECTORS;
        break;
      case 'T':
        mode = TempShvec_globals::TEMP_SHVEC_MODE_HAN_TRAN;
        coset = false;
        NeedBound = true;
        break;
      case 'm':
        mode = TempShvec_globals::TEMP_SHVEC_MODE_MINIMUM;
        break;
      case 'v':
        mode = TempShvec_globals::TEMP_SHVEC_MODE_VINBERG_ALGO;
        coset = true;
        NeedBound = true;
        break;
      case 'l':
        mode = TempShvec_globals::TEMP_SHVEC_MODE_LORENTZIAN;
        coset = true;
        NeedBound = true;
        break;
      case 'c':
        coset = true;
        break;
      case 'e':
        break;
      default:
        die_sv("invalid option\nTry 'sv -h' for more information.\n");
      }
    std::cerr << "main mode=" << mode << "\n";
    //
    // First reading data
    //
    int dim;
    std::cin >> dim;
    MyMatrix<T> gram_matrix(dim, dim);
    T eT;
    for (i = 0; i < dim; i++)
      for (j = 0; j <= i; j++) {
        std::cin >> eT;
        gram_matrix(i, j) = eT;
        gram_matrix(j, i) = eT;
      }
    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
        std::cerr << " " << gram_matrix(i, j);
      }
      std::cerr << "\n";
    }
    MyVector<T> cosetVect = ZeroVector<T>(dim);
    if (coset) {
      for (i = 0; i < dim; i++) {
        std::cin >> eT;
        cosetVect(i) = eT;
      }
    }
    std::cerr << "coset=" << coset << " eV =";
    for (i = 0; i < dim; i++)
      std::cerr << " " << cosetVect(i);
    std::cerr << "\n";
    if (NeedBound) {
      std::cin >> eT;
      bound = eT;
    }
    std::cerr << "NeedBound=" << NeedBound << " bound=" << bound << "\n";

    //
    // Defining info and computing with it
    //
    T_shvec_request<T> request =
        initShvecReq(gram_matrix, cosetVect, bound, mode);
    std::cerr << "Before computeShvec mode=" << mode << "\n";
    T_shvec_info<T, Tint> info =
        T_computeShvec<T, Tint>(request, mode, std::cerr);
    int nbVect = info.short_vectors.size();
    std::cerr << "After computeShvec |V|=" << nbVect << "\n";
    //
    // Checking central symmetry
    //
#ifdef CHECK_SHVEC
    bool IsCentrallySymmetric = true;
    for (i = 0; i < dim; i++) {
      T twoT = 2 * eT;
      if (!IsInteger(twoT)) {
        IsCentrallySymmetric = false;
      }
    }
    if (IsCentrallySymmetric) {
      MyVector<T> TwoCoset_T = -2 * cosetVect;
      MyVector<Tint> TwoCoset = UniversalVectorConversion<Tint, T>(TwoCoset_T);
      std::unordered_set<MyVector<Tint>> set;
      for (auto &V : info.short_vectors) {
        set.insert(V);
      }
      for (auto &V : info.short_vectors) {
        MyVector<Tint> Vdiff = TwoCoset - V;
        if (set.count(Vdiff) == 0) {
          std::cerr << "The problem is centrally symmetric\n";
          std::cerr << "cosetVect=" << StringVector(cosetVect) << "\n";
          std::cerr << "But the vector V=" << StringVector(V)
                    << " is in the solution\n";
          std::cerr << "while Vdiff=" << StringVector(Vdiff)
                    << " is not contained\n";
          throw TerminalException{1};
        }
      }
    }
#endif
    //
    // Data output
    //
#ifdef CHECK_SHVEC
    std::map<T, size_t> map_norm;
#endif
    std::cout << nbVect << "\n";
    for (i = 0; i < nbVect; i++) {
      MyVector<Tint> const &V = info.short_vectors[i];
      std::cout << "[ unset ]: ";
      for (int iCol = 0; iCol < dim; iCol++)
        std::cout << " " << V(iCol);
      std::cout << "\n";
#ifdef CHECK_SHVEC
      MyVector<T> V_T = UniversalVectorConversion<T, Tint>(V);
      if (coset && NeedBound) {
        T eNorm(0);
        for (int i = 0; i < dim; i++) {
          for (int j = 0; j < dim; j++) {
            eNorm += (V_T(i) + cosetVect(i)) * (V_T(j) + cosetVect(j)) *
                     gram_matrix(i, j);
          }
        }
        if (eNorm > bound) {
          std::cerr << "eNorm=" << eNorm << " bound=" << bound << "\n";
          throw TerminalException{1};
        }
        map_norm[eNorm] += 1;
      }
#endif
    }
#ifdef CHECK_SHVEC
    std::cerr << "map_norm =";
    for (auto &kv : map_norm) {
      std::cerr << " (" << kv.first << "," << kv.second << ")";
    }
    std::cerr << "\n";
#endif
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in sv_exact\n";
    exit(e.eVal);
  }
  runtime(time);
}
