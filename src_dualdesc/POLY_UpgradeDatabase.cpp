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
  HumanTime time;
  try {
    using T = mpq_class;
    //
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    if (argc != 2) {
      std::cerr << "POLY_UpgradeDatabase [ThePrefix]\n";
      return -1;
    }
    std::string ThePrefix = argv[1];
    //
    std::string eFileEXT = ThePrefix + ".ext";
    std::string eFileGRP = ThePrefix + ".grp";
    std::string eFileOrbitSize = ThePrefix + ".orbitsize";
    std::string eFileNB = ThePrefix + ".nb";
    std::string eFileFF = ThePrefix + ".ff";
    IsExistingFileDie(eFileEXT);
    IsExistingFileDie(eFileGRP);
    IsExistingFileDie(eFileNB);
    IsExistingFileDie(eFileFF);
    //
    MyMatrix<T> EXT = ReadMatrixFile<T>(eFileEXT);
    size_t n_row = EXT.rows();
    //
    std::ifstream is_grp(eFileGRP);
    Tgroup GRP;
    is_grp >> GRP;
    //
    FileNumber fn(eFileNB, false);
    size_t n_orbit = fn.getval();
    auto read_vf = [&]() -> vectface {
      FileFace ff(eFileFF, n_row, n_orbit);
      vectface vf(n_row);
      for (size_t i_orbit = 0; i_orbit < n_orbit; i_orbit++) {
        Face f = ff.getface(i_orbit);
        vf.push_back(f);
      }
      return vf;
    };
    vectface vf = read_vf();
    //
    std::map<Tidx, int> LFact = GRP.factor_size();
    std::pair<size_t, size_t> ep = get_delta(LFact, n_row);
    size_t n_bit_orbsize = ep.first;
    size_t delta = ep.second;
    std::vector<Tint> ListPossOrbsize = GetAllPossibilities<Tidx, Tint>(LFact);
    std::unordered_map<Tint, size_t> OrbSize_map;
    for (size_t i_poss = 0; i_poss < ListPossOrbsize.size(); i_poss++) {
      OrbSize_map[ListPossOrbsize[i_poss]] = i_poss;
    }
    //
    vectface vfo(delta);
    for (size_t i_orbit = 0; i_orbit < n_orbit; i_orbit++) {
      Face f = vf[i_orbit];
      Tint orbSize = GRP.OrbitSize_OnSets(f);
      size_t idx_orb = OrbSize_map[orbSize];
      Face f_full(delta);
      for (size_t i = 0; i < n_row; i++)
        f_full[i] = f[i];
      size_t work_idx = idx_orb;
      size_t i_acc = n_row;
      for (size_t i = 0; i < n_bit_orbsize; i++) {
        bool val = work_idx % 2;
        f_full[i_acc] = val;
        i_acc++;
        work_idx = work_idx / 2;
      }
      vfo.push_back(f_full);
    }
    //
    // Now writing the new data
    //
    std::ofstream os_orbitsize(eFileOrbitSize);
    os_orbitsize << ListPossOrbsize;
    //
    RemoveFileIfExist(eFileFF);
    FileFace ff(eFileFF, delta);
    ff.direct_write(vfo.serial_get_std_vector_uint8_t());
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something went wrong\n";
    exit(e.eVal);
  }
  runtime(time);
}
