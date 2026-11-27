// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "CtypeValidSimplices.h"
#include "Permutation.h"
#include "Group.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T, typename Tint>
void process(std::string const &FileM, std::string const &OutFormat,
             std::ostream &os_out) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tint_grp = mpz_class;
  using Tgroup = permutalib::Group<Telt, Tint_grp>;
  using Tidx_ctype = uint16_t;
  MyMatrix<T> GramMat = ReadMatrixFile<T>(FileM);
  if (!IsSymmetricMatrix(GramMat) || !IsPositiveDefinite(GramMat, std::cerr)) {
    std::cerr << "The matrix in FileM=" << FileM
              << " should be symmetric positive definite\n";
    throw TerminalException{1};
  }
  HumanTime time;
  ctype_valid::ResultValidSimplices<T, Tint> res =
      ctype_valid::ComputeValidSimplices<T, Tint, Tgroup, Tidx_ctype>(
          GramMat, std::cerr);
  std::map<int, size_t> map_rank;
  for (auto &R : res.cone.ListRay) {
    map_rank[RankMat(R)]++;
  }
  std::cerr << "|ListIneq|=" << res.cone.ListIneq.rows()
            << " |ListRay|=" << res.cone.ListRay.size() << " ranks:";
  for (auto &kv : map_rank) {
    std::cerr << " " << kv.first << ":" << kv.second;
  }
  std::cerr << "\n";
  size_t n_simplex = 0, n_valid = 0;
  for (auto &cell : res.ListCell) {
    n_simplex += cell.ListSimplex.size();
    for (bool v : cell.ListValid) {
      if (v) {
        n_valid++;
      }
    }
  }
  std::cerr << "|ListCell|=" << res.ListCell.size()
            << " candidate simplices=" << n_simplex << " valid=" << n_valid
            << " (counted per orbit of Delaunay polytopes)\n";
  std::cerr << "|ComputeValidSimplices|=" << time << "\n";
  if (OutFormat == "GAP") {
    ctype_valid::WriteValidSimplicesGAP(os_out, res);
    return;
  }
  if (OutFormat == "None") {
    return;
  }
  std::cerr << "No type available for OutFormat=" << OutFormat << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_CtypeValidSimplices [arith] [FileM] [OutFormat] "
                   "[OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_CtypeValidSimplices [arith] [FileM]\n";
      std::cerr << "FileM is a Gram matrix in the interior of a generic "
                   "iso-edge domain C.\n";
      std::cerr << "The program lists the simplices spanned by its Delaunay "
                   "polytopes that are valid on the whole of C.\n";
      std::cerr << "OutFormat is GAP or None\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileM = argv[2];
    std::string OutFormat = "None";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    auto f = [&](std::ostream &os) -> void {
      ArithmeticDispatch dispatch;
      add_rational_integer_arithmetics_by_library(
          dispatch, [&]<typename T, typename Tint>() -> void {
            process<T, Tint>(FileM, OutFormat, os);
          });
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    std::cerr << "Normal termination of LATT_CtypeValidSimplices\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_CtypeValidSimplices\n";
    exit(e.eVal);
  }
  runtime(time1);
}
