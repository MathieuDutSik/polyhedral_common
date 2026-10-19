// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#ifdef QUADRATIC_FIELDS
#include "NumberTheoryQuadField.h"
#endif
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "GRP_GroupFile.h"
#include "Group.h"
#include "POLY_Kskeletton.h"
#include "Permutation.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T, typename Tgroup>
void MainFunctionFaceLattice_A(FullNamelist const &eFull, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
#ifdef DEBUG_POLY_KSKELETTON
  os << "SKEL: Reading PROC\n";
#endif
  SingleBlock const &BlockPROC = eFull.get_block("PROC");
  std::string const &FACfile = BlockPROC.get_string("FACfile");
  MyMatrix<T> FAC = ReadMatrixFile<T>(FACfile);
#ifdef DEBUG_POLY_KSKELETTON
  os << "SKEL: |FAC|=" << FAC.rows() << " / " << FAC.cols() << "\n";
#endif
  if (size_t(FAC.rows()) > size_t(std::numeric_limits<Tidx>::max())) {
    std::cerr << "SKEL: We have |FAC|=" << FAC.rows() << "\n";
    std::cerr << "SKEL: But <Tidx>::max()="
              << size_t(std::numeric_limits<Tidx>::max()) << "\n";
    throw TerminalException{1};
  }
  if (RankMat(FAC) != FAC.cols()) {
    std::cerr << "The matrix FAC should be of full rank\n";
    throw TerminalException{1};
  }
  //
  std::string method_spann = BlockPROC.get_string("method_spann");
  std::string method_final = BlockPROC.get_string("method_final");
  bool ComputeTotalNumberFaces = BlockPROC.get_bool("ComputeTotalNumberFaces");
#ifdef DEBUG_POLY_KSKELETTON
  os << "SKEL: method_final=" << method_final << "\n";
  os << "SKEL: method_spann=" << method_spann << "\n";
  os << "SKEL: ComputeTotalNumberFaces=" << ComputeTotalNumberFaces << "\n";
#endif
  //
  MyMatrix<T> EXT;
  if (method_spann == "ExtremeRays" ||
      method_spann == "ExtremeRaysNonSimplicial") {
    std::string EXTfile = BlockPROC.get_string("EXTfile");
    EXT = ReadMatrixFile<T>(EXTfile);
    if (FAC.cols() != EXT.cols()) {
      std::cerr << "The dimension of EXT and FAC should be the same\n";
      throw TerminalException{1};
    }
  }
  //
  std::string GRPfile = BlockPROC.get_string("GRPfile");
  Tgroup GRP = ReadGroupFile<Tgroup>(GRPfile);
#ifdef DEBUG_POLY_KSKELETTON
  os << "SKEL: |GRP|=" << GRP.size() << "\n";
#endif
  //
  int LevSearch = BlockPROC.get_int("LevSearch");
#ifdef DEBUG_POLY_KSKELETTON
  os << "SKEL: LevSearch=" << LevSearch << "\n";
#endif
  if (LevSearch == -1) {
    int nbCol = FAC.cols();
    LevSearch = nbCol - 2;
  }
  //
  std::string OutFile = BlockPROC.get_string("OutFile");
  std::string OutFormat = BlockPROC.get_string("OutFormat");
#ifdef DEBUG_POLY_KSKELETTON
  os << "SKEL: OutFile=" << OutFile << " OutFormat=" << OutFormat << "\n";
#endif
  //
  std::vector<vectface> TheOutput =
      EnumerationFaces(GRP, FAC, EXT, LevSearch, method_spann, method_final,
                       ComputeTotalNumberFaces, os);
  //
  auto f_print = [&](std::ostream &os_out) -> void {
    OutputFaces_stream(TheOutput, os_out, OutFormat);
  };
  FILE_PrintStderrStdoutFile(OutFile, f_print);
  //
  SingleBlock const &BlockGROUP = eFull.get_block("GROUP");
  bool ComputeAutGroup = BlockGROUP.get_bool("ComputeAutGroup");
  if (ComputeAutGroup) {
    //    using Tgr = GraphBitset;
    using Tgr = GraphListAdj;
    Tgroup GRPfull =
        ComputeGroupFromOrbitFaces<Tgroup, Tgr>(TheOutput, GRP, os);
#ifdef DEBUG_POLY_KSKELETTON
    os << "SKEL: |GRPfull|=" << GRPfull.size() << "\n";
#endif
    std::string FileGroup = BlockGROUP.get_string("FileGroup");
    std::string OutFormat = BlockGROUP.get_string("OutFormat");
#ifdef DEBUG_POLY_KSKELETTON
    os << "SKEL: FileGroup=" << FileGroup << "\n";
    os << "SKEL: OutFormat=" << OutFormat << "\n";
#endif
    WriteGroupFormat(FileGroup, OutFormat, GRPfull);
  }
}

template <typename Tgroup>
void MainFunctionFaceLattice(FullNamelist const &eFull) {
  SingleBlock const &BlockPROC = eFull.get_block("PROC");
  std::string const &arith = BlockPROC.get_string("Arithmetic");
  ArithmeticDispatch dispatch;
  auto f_arith = [&]<typename T>() -> void {
    MainFunctionFaceLattice_A<T, Tgroup>(eFull, std::cerr);
  };
  add_rational_arithmetics(dispatch, f_arith);
  add_algebraic_field_arithmetics(dispatch, f_arith);
  dispatch.run(arith);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    //
    FullNamelist eFull = NAMELIST_GetStandard_FaceLattice();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_FaceLatticeGen [file.nml]\n";
      std::cerr << "with file.nml a namelist\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string FileNML = argv[1];
    NAMELIST_ReadNamelistFile(FileNML, eFull);
    //
    MainFunctionFaceLattice<Tgroup>(eFull);
    std::cerr << "Normal termination of POLY_FaceLatticeGen\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_FaceLatticeGen\n";
    exit(e.eVal);
  }
  runtime(time1);
}
