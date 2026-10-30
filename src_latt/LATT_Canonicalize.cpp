// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Group.h"
#include "Permutation.h"
#include "LatticeStabEquiCan.h"
#include "LatticeOrbitSplitCanonical.h"
#include "SignatureSymmetric.h"
#include "LatticeRootDecomposition.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void ComputeCanonical(std::string const &FileI, std::string const &OutFormat,
                      std::string const &method, std::ostream &os) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  MyMatrix<T> eMat = ReadMatrixFile<T>(FileI);
  if (!IsSymmetricMatrix(eMat) || !IsPositiveDefinite(eMat, std::cerr)) {
    std::cerr << "LATT_Canonicalize: The input Gram matrix in " << FileI
              << " is not symmetric positive definite\n";
    throw TerminalException{1};
  }
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  auto f_canonical = [&]() -> MyMatrix<Tint> {
    if (method == "graph") {
      return ComputeCanonicalForm<T, Tint, Tgroup>(eMat, std::cerr);
    }
    if (method == "plesken_souvignier") {
      // Orbit splitting over the Plesken-Souvignier automorphism group:
      // nothing quadratic in the family, so it reaches lattices like Leech
      // whose vector families the graph method cannot handle. The method
      // can fail (surviving orbit tie, straddled sign components); the
      // failure is a property of the isometry class, so the forms it does
      // canonicalize remain comparable with each other.
      std::optional<MyMatrix<Tint>> opt =
          OrbitSplitCanonicalizeGram<T, Tint, Tgroup>(eMat, std::cerr);
      if (!opt) {
        std::cerr << "LATT_Canonicalize: the orbit splitting could not "
                     "select a canonical prefix, no canonical form for this "
                     "input\n";
        throw TerminalException{1};
      }
      return *opt;
    }
    if (method == "vinberg") {
      // The roots are the vectors of norm 2, which is the decomposition of
      // an integral form.
      if (!IsIntegralMatrix(eMat)) {
        std::cerr << "LATT_Canonicalize: the vinberg method needs an "
                     "integral Gram matrix\n";
        throw TerminalException{1};
      }
      return ComputeCanonicalFormRootDecomposed<T, Tint, Tgroup>(eMat,
                                                                 std::cerr);
    }
    if (method == "weyl") {
      // As vinberg, but with the simple roots of one positive system in
      // place of the whole root system, so the Weyl group leaves the
      // configuration.
      if (!IsIntegralMatrix(eMat)) {
        std::cerr << "LATT_Canonicalize: the weyl method needs an integral "
                     "Gram matrix\n";
        throw TerminalException{1};
      }
      return ComputeCanonicalFormWeylReduced<T, Tint, Tgroup>(eMat,
                                                              std::cerr);
    }
    if (method == "weyl_chamber") {
      // The characteristic vectors in the fundamental Weyl chamber and the
      // simple roots: a family spanning L, so no glue step.
      if (!IsIntegralMatrix(eMat)) {
        std::cerr << "LATT_Canonicalize: the weyl_chamber method needs an "
                     "integral Gram matrix\n";
        throw TerminalException{1};
      }
      return ComputeCanonicalFormWeylChamber<T, Tint, Tgroup>(eMat,
                                                              std::cerr);
    }
    std::cerr << "LATT_Canonicalize: unknown method=" << method << "\n";
    std::cerr << "The methods are graph, plesken_souvignier, vinberg, weyl "
                 "and weyl_chamber\n";
    throw TerminalException{1};
  };
  MyMatrix<Tint> B = f_canonical();
  MyMatrix<T> B_T = UniversalMatrixConversion<T,Tint>(B);
  MyMatrix<T> eMat_red = B_T * eMat * B_T.transpose();
  if (OutFormat == "CPP") {
    WriteMatrix(os, eMat_red);
    return;
  }
  if (OutFormat == "PYTHON") {
    os << "{\"Basis\":" << StringMatrixPYTHON(B)
       << ", \"eG\":" << StringMatrixPYTHON(eMat_red) << "}\n";
    return;
  }
  if (OutFormat == "GAP") {
    os << "return ";
    WriteMatrixGAP(os, eMat_red);
    os << ";\n";
    return;
  }
  if (OutFormat == "GAP_full") {
    os << "return rec(Basis:=";
    WriteMatrixGAP(os, B);
    os << ", eG:=";
    WriteMatrixGAP(os, eMat_red);
    os << ");\n";
    return;
  }
  std::cerr << "Failed to find a matching entry\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 3 && argc != 5 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_Canonicalize [arith] [GramI]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_Canonicalize [arith] [GramI] [OutFormat] [OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_Canonicalize [arith] [GramI] [OutFormat] [OutFile] "
                   "[method]\n";
      std::cerr << "\n";
      std::cerr << "GramI (input) : The gram matrix on input\n";
      std::cerr << "OutFile: The filename of the data in output\n";
      std::cerr << "method values:\n";
      std::cerr << "  graph   : from the characteristic vector family "
                   "(default)\n";
      std::cerr << "  plesken_souvignier : orbit splitting over the "
                   "automorphism group\n";
      std::cerr << "              computed by the Plesken-Souvignier "
                   "backtrack; viable on\n";
      std::cerr << "              large families (Leech), can decline on "
                   "some inputs\n";
      std::cerr << "  vinberg : along the root decomposition, for an "
                   "integral Gram matrix\n";
      std::cerr << "  weyl    : as vinberg, with only the simple roots in "
                   "the family\n";
      std::cerr << "  weyl_chamber : the characteristic vectors in the Weyl "
                   "chamber and the\n";
      std::cerr << "              simple roots, a family spanning the "
                   "lattice\n";
      std::cerr << "The canonical forms of two methods are not comparable "
                   "with each other.\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileI = argv[2];
    std::string OutFormat = "CPP";
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
    auto f = [&](std::ostream &os) -> void {
      ArithmeticDispatch dispatch;
      add_rational_arithmetics_by_library(
          dispatch, [&]<typename T>() -> void {
            ComputeCanonical<T>(FileI, OutFormat, method, os);
          });
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    //
    std::cerr << "Normal termination of LATT_Canonicalize\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_Canonicalize\n";
    exit(e.eVal);
  }
  runtime(time);
}
