// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
//#include "NumberTheory.h"
#include "Group.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include "POLY_RecursiveDualDesc.h"
#include "Permutation.h"

int main(int argc, char *argv[]) {
  SingletonTime time1;
  try {
    if (argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_DatabaseRestructuration FileGRP] [DatabaseInput] [NprocInput] [DatabaseOutput] [NprocOutput]\n";
      std::cerr << "Be careful when using it, it depends on so many aspects of the code\n";
      return -1;
    }
    std::string FileGRP = argv[1];
    std::string DatabaseI = argv[2];
    int NprocI = ParseScalar<int>(argv[3]);
    std::string DatabaseO = argv[4];
    int NprocO = ParseScalar<int>(argv[5]);
    //
    // Reading the group
    //
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    std::ifstream GRPfs(FileGRP);
    Tgroup GRP = ReadGroup<Tgroup>(GRPfs);
    std::map<Tidx,int> LFact = GRP.factor_size();
    Tidx n_act = GRP.n_act();
    std::vector<size_t,size_t> ep = get_delta(LFact, n_act);
    size_t delta = ep.second;
    //
    // Now the streams
    //
    std::vector<size_t> List_shift(NprocO, 0);
    std::vector<FileNumber*> List_FN(NprocO, nullptr);
    std::vector<FileBool*> List_FB(NprocO, nullptr);
    std::vector<FileFace*> List_FF(NprocO, nullptr);
    for (int iProc=0; iProc<NprocO; iProc++) {
      std::string MainPrefix = DatabaseO + "_todo";
      std::string eFileFN = MainPrefix + ".nb";
      std::string eFileFB = MainPrefix + ".fb";
      std::string eFileFF = MainPrefix + ".ff";
      List_FN[iProc] = new FileNumber(eFileNB, true);
      List_FB[iProc] = new FileBool(eFileFB);
      List_FF[iProc] = new FileFace(eFileFF, delta);
    }
    //
    // Now reading and preprocessing it
    //
    for (int iProc=0; iProc<NprocI; iProc++) {
      std::string MainPrefix = DatabaseI + "_todo";
      std::string eFileFN = MainPrefix + ".nb";
      std::string eFileFB = MainPrefix + ".fb";
      std::string eFileFF = MainPrefix + ".ff";
      FileNumber fn(eFileFN, false);
      size_t n_orbit = fn.getval();
      FileBool fb(eFileFB, n_orbit);
      FileFace ff(eFileFF, delta, n_orbit);
      for (size_t pos=0; pos<n_orbit; pos++) {
        Face f = ff.getface(pos);
        bool val = fb.getbit(pos);
        Face f_res(n_act);
      }
    }

    //
    // Now writing the n_orbit to the fileNB
    //
    for (int iProc=0; iProc<NprocO; iProc++) {
      List_FN[iProc]->setval(List_shift[iProc]);
    }
    //
    // deleting the object which is critical for writing down all data to file.
    //
    for (int iProc=0; iProc<NprocO; iProc++) {
      delete List_FN[iProc];
      delete List_FB[iProc];
      delete List_FF[iProc];
    }
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const &e) {
    std::cerr << "Something went wrong in the computation, please debug\n";
    exit(e.eVal);
  }
  runtime(time1);
}
