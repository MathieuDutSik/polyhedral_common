// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include "NumberTheoryRealField.h"
#include "NumberTheoryQuadField.h"
#include "POLY_RecursiveDualDesc.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

int main(int argc, char *argv[]) {
  try {
    if (argc != 5) {
      std::cerr << "POLY_ReadPartialEnum [FileGRP] [DatabaseInput] [OutFormat] [OutFile]\n";
      return -1;
    }
    std::string FileGRP = argv[1];
    std::string DatabaseI = argv[2];
    std::string OutFormat = argv[3];
    std::string OutFile = argv[4];
    //
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using T = mpq_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    std::ifstream GRPfs(FileGRP);
    Tgroup GRP = ReadGroup<Tgroup>(GRPfs);
    std::map<Tidx, int> LFact = GRP.factor_size();
    Tidx n_act = GRP.n_act();
    std::pair<size_t, size_t> ep = get_delta(LFact, n_act);
    size_t delta = ep.second;
    using Torbsize = uint32_t;
    DataFaceOrbitSize<Torbsize, Tgroup> dfo(GRP);
    //
    std::string eFileFN = DatabaseI + ".nb";
    std::string eFileFB = DatabaseI + ".fb";
    std::string eFileFF = DatabaseI + ".ff";
    std::string OFileEXT = DatabaseI + ".ext";
    std::string OFileMethod = DatabaseI + ".method";
    bool overwrite = false;
    FileNumber fn(eFileFN, overwrite);
    size_t n_orbit = fn.getval();
    FileBool fb(eFileFB, n_orbit);
    FileFace ff(eFileFF, delta, n_orbit);
    std::map<size_t, size_t> map_incidence;
    std::map<Tint, size_t> map_stabsize;
    vectface TheOutput(n_act);
    for (size_t i_orbit=0; i_orbit<n_orbit; i_orbit++) {
      Face f = ff.getface(i_orbit);
      std::pair<Face, Tint> pair = dfo.ConvertFace(f);
      Tint stabSize = GRP.size() / pair.second;
      Face f_red = pair.first;
      TheOutput.push_back(f_red);
      size_t cnt = f_red.count();
      map_incidence[cnt] += 1;
      map_stabsize[stabSize] += 1;
    }
    for (auto & kv : map_incidence) {
      std::cerr << "incidence " << kv.first << " attained " << kv.second << " times\n";
    }
    for (auto & kv : map_stabsize) {
      std::cerr << "stbsize " << kv.first << " attained " << kv.second << " times\n";
    }
    MyMatrix<T> EXT;
    OutputFacets(EXT, GRP, TheOutput, OutFile, OutFormat, std::cerr);
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something went wrong\n";
    exit(e.eVal);
  }
}
