// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "NumberTheoryRealField.h"
#ifdef QUADRATIC_FIELDS
#include "NumberTheoryQuadField.h"
#endif
#include "Group.h"
#include "Permutation.h"
#include "LatticeStabEquiCan.h"
#include "InvariantVectorFamily.h"
#include "SignatureSymmetric.h"
#include "LatticePleskenSouvignier.h"
#include "LatticeAutomorphismVinberg.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void ComputeAutomorphism(std::string const &FileListMat,
                         std::string const &OutFormat,
                         std::string const &method, std::ostream &os) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  std::vector<MyMatrix<T>> ListMat = ReadListMatrixFile<T>(FileListMat);
  if (ListMat.empty()) {
    std::cerr << "LATT_Automorphism: The input matrix list in " << FileListMat
              << " is empty\n";
    throw TerminalException{1};
  }
  if (!IsSymmetricMatrix(ListMat[0]) ||
      !IsPositiveDefinite(ListMat[0], std::cerr)) {
    std::cerr << "LATT_Automorphism: The first input Gram matrix in "
              << FileListMat << " is not symmetric positive definite\n";
    throw TerminalException{1};
  }
  if (method != "graph" && method != "plesken_souvignier" &&
      method != "vinberg" && method != "weyl_chamber") {
    std::cerr << "LATT_Automorphism: unknown method=" << method << "\n";
    std::cerr << "The methods are graph, plesken_souvignier, vinberg and "
                 "weyl_chamber\n";
    throw TerminalException{1};
  }
  std::vector<MyMatrix<Tint>> ListGen;
  // The order, when the method computes it on the way.
  std::optional<mpz_class> opt_order;
  if (method == "graph") {
    ListGen = ArithmeticAutomorphismGroupMultiple<T, Tint, Tgroup>(ListMat,
                                                                   std::cerr);
    if (OutFormat == "GAP_order") {
      // The action on the canonic family is faithful, the family being of
      // full rank.
      CanonicVectorFamily<Tint> fam =
          GetCanonicVectorFamily<T, Tint>(ListMat[0], std::cerr);
      opt_order = OrderFromGens<Tint, Tgroup>(fam.get_full(), ListGen);
    }
  } else {
    // Plesken-Souvignier and Vinberg work with the lattice coordinates over
    // the rationals; the algebraic fields only have the graph method.
    if constexpr (is_implementation_of_Q<T>::value) {
      if (method == "plesken_souvignier") {
        PleskenSouvignierAutomResult<Tint> result =
            PleskenSouvignierLatticeAutomorphism<T, Tint>(ListMat, std::cerr);
        ListGen = result.ListGen;
        opt_order =
            PleskenSouvignierGroupOrder<mpz_class>(result.ListOrbitSize);
      } else {
        // The reflections in the roots of norm 2 are integral only for an
        // integral form, and the decomposition is the one of a single form.
        if (ListMat.size() != 1 || !IsIntegralMatrix(ListMat[0])) {
          std::cerr << "LATT_Automorphism: the " << method
                    << " method needs a single integral Gram matrix\n";
          throw TerminalException{1};
        }
        VinbergAutom<Tint> vin =
            (method == "vinberg")
                ? ComputeAutomorphismVinberg<T, Tint, Tgroup>(ListMat[0],
                                                              std::cerr)
                : ComputeAutomorphismWeylChamber<T, Tint, Tgroup>(ListMat[0],
                                                                  std::cerr);
        ListGen = vin.ListGen;
        opt_order = vin.order;
      }
    } else {
      std::cerr << "LATT_Automorphism: the method " << method
                << " is only available for the rational arithmetics\n";
      throw TerminalException{1};
    }
  }
  if (OutFormat == "GAP_order") {
    os << "return " << *opt_order << ";\n";
    return;
  }
  if (OutFormat == "GAP") {
    os << "return ";
    WriteListMatrixGAP(os, ListGen);
    os << ";\n";
    return;
  }
  if (OutFormat == "Oscar") {
    os << ListGen.size() << "\n";
    for (auto &eMat : ListGen) {
      WriteMatrix(os, eMat);
    }
    return;
  }
  std::cerr << "Failed to find a matching type for OutFormat=" << OutFormat
            << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 3 && argc != 5 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_Automorphism [arith] [ListMat] [OutFormat] "
                << "[OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_Automorphism [arith] [ListMat]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_Automorphism [arith] [ListMat] [OutFormat] "
                << "[OutFile] [method]\n";
      std::cerr << "\n";
      std::cerr << "arith values:\n";
      std::cerr << "  gmp         : mpq_class / mpz_class (default choice)\n";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << "  gmp_boost   : the boost bindings to the gmp types\n";
#endif
#ifdef ENABLE_BOOST_TYPES
      std::cerr << "  multi_boost : the boost multiprecision types\n";
#endif
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "  flint       : fmpq_class / fmpz_class\n";
#endif
#ifdef QUADRATIC_FIELDS
      std::cerr << "  Qsqrt2      : the form over the field Q(sqrt(2))\n";
#endif
#ifdef QUADRATIC_FIELDS
      std::cerr << "  Qsqrt5      : the form over the field Q(sqrt(5))\n";
#endif
      std::cerr << "  RealAlgebraic=FileDesc : the form over the real\n";
      std::cerr << "              algebraic field described in FileDesc\n";
      std::cerr << "              (the lattice is Z^n in every case, so the\n";
      std::cerr << "              automorphisms are integral matrices)\n";
      std::cerr << "OutFormat values:\n";
      std::cerr << "  GAP       : ListGen returned as a GAP-readable list of\n";
      std::cerr << "              integral matrix generators (default)\n";
      std::cerr << "  Oscar     : ListGen returned in the Oscar format\n";
      std::cerr << "  GAP_order : the order |Aut(GramMat)| only, computed "
                << "via\n";
      std::cerr << "              the permutation action on a full-rank\n";
      std::cerr << "              invariant vector family (skips the matrix\n";
      std::cerr << "              lift, much cheaper)\n";
      std::cerr << "method values:\n";
      std::cerr << "  graph              : the automorphisms of a vector "
                   "family as a\n";
      std::cerr << "                       graph (default)\n";
      std::cerr << "  plesken_souvignier : the backtrack of Plesken and "
                   "Souvignier\n";
      std::cerr << "                       (rational arithmetics)\n";
      std::cerr << "  vinberg            : Aut(L) = W(R) x| Aut(L, rho) "
                   "for a single\n";
      std::cerr << "                       integral Gram matrix (rational "
                   "arithmetics)\n";
      std::cerr << "  weyl_chamber       : as vinberg, Aut(L, rho) from the "
                   "Weyl-chamber\n";
      std::cerr << "                       family, which spans the lattice\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileListMat = argv[2];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    std::string method = "graph";
    if (argc >= 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    if (argc == 6) {
      method = argv[5];
    }
    //
    auto prt = [&](std::ostream &os) -> void {
      // The lattice is Z^n whatever field the form takes its values in, so
      // Tint stays mpz_class: only the Gram matrices leave the rationals.
      ArithmeticDispatch dispatch;
      auto f_arith = [&]<typename T>() -> void {
        ComputeAutomorphism<T>(FileListMat, OutFormat, method, os);
      };
      add_rational_arithmetics_by_library(dispatch, f_arith);
      add_algebraic_field_arithmetics(dispatch, f_arith);
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(OutFile, prt);
    std::cerr << "Normal termination of LATT_Automorphism\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_Automorphism\n";
    exit(e.eVal);
  }
  runtime(time);
}
