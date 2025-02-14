// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include "Group.h"
#include "POLY_RecursiveDualDesc.h"
#include "POLY_RecursiveDualDesc_MPI.h"
#include "Permutation.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    using T = mpq_class;
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    if (argc != 7) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_EvaluateBalinski [FileFAC] [FileGRP] [VFformat] "
                   "[vf_undone] [max_iter] [max_depth]\n";
      std::cerr << "\n";
      return -1;
    }
    std::string FileFAC = argv[1];
    std::string FileGRP = argv[2];
    std::string VFformat = argv[3];
    std::string FileUndone = argv[4];
    std::string str_max_iter = argv[5];
    std::string str_max_depth = argv[6];
    //
    MyMatrix<T> FAC = ReadMatrixFile<T>(FileFAC);
    Tgroup GRP = ReadGroupFile<Tgroup>(FileGRP);
    int n = FAC.rows();
    vectface vf_undone_orbit = ReadFacets<T, Tgroup>(VFformat, FileUndone, n);
    size_t max_iter = ParseScalar<size_t>(str_max_iter);
    size_t max_depth = ParseScalar<size_t>(str_max_depth);
    //
    vectface vf_undone = GenerateOrbits(vf_undone_orbit, GRP);
    MyMatrix<T> EXT_undone = GetVertexSet_from_vectface(FAC, vf_undone);
    size_t n_iter = 0;
    auto f_recur = [&](const std::pair<size_t, Face> &pfr) -> bool {
      n_iter++;
      std::cerr << "BAL:  f_recur n_iter=" << n_iter << "\n";
      if (n_iter == max_iter)
        return false;
      if (pfr.first > max_depth)
        return false;
      return true;
    };
    bool test = EvaluationConnectednessCriterion_Kernel(
        FAC, GRP, EXT_undone, vf_undone, f_recur, std::cerr);
    std::cerr << "Obtained result=" << test << "\n";
    std::cerr << "Normal termination of POLY_EvaluateBalinski\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_EvaluateBalinski\n";
    exit(e.eVal);
  }
  runtime(time1);
}
