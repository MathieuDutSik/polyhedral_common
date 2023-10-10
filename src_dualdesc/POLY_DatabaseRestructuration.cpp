// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "Group.h"
#include "TypeConversion.h"
#include "TypeConversionFinal.h"
#include "NumberTheoryGmp.h"
#include "POLY_RecursiveDualDesc.h"
#include "POLY_RecursiveDualDesc_MPI.h"
#include "Permutation.h"
// clang-format on

int main(int argc, char *argv[]) {
  SingletonTime time1;
  try {
    if (argc != 8) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_DatabaseRestructuration [FileGRP] [DatabaseInput] "
                   "[NprocInput] [DatabaseOutput] [NprocOutput] [IgnoreLastI] [IgnoreLastO]\n";
      std::cerr << "Be careful when using it, it depends on so many aspects of "
                   "the code\n";
      return -1;
    }
    std::string FileGRP = argv[1];
    std::string DatabaseI = argv[2];
    int NprocI = ParseScalar<int>(argv[3]);
    std::string DatabaseO = argv[4];
    int NprocO = ParseScalar<int>(argv[5]);
    bool IgnoreLastI = (bool)ParseScalar<int>(argv[6]);
    bool IgnoreLastO = (bool)ParseScalar<int>(argv[7]);
    assert(!IgnoreLastI || NprocI > 1);
    assert(!IgnoreLastO || NprocO > 1);
    //
    // Reading the group
    //
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    std::ifstream GRPfs(FileGRP);
    Tgroup GRP = ReadGroup<Tgroup>(GRPfs);
    std::map<Tidx, int> LFact = GRP.factor_size();
    Tidx n_act = GRP.n_act();
    std::pair<size_t, size_t> ep = get_delta(LFact, n_act);
    size_t delta = ep.second;
    std::string dbnameI = (NprocI>1) ? "database" : "D_"+std::to_string(n_act);
    std::string dbnameO = (NprocO>1) ? "database" : "D_"+std::to_string(n_act);
    std::cerr << NprocI << " " << dbnameI << std::endl;
    std::cerr << NprocO << " " << dbnameO << std::endl;
    //
    // Reading the EXT and method file
    //
    std::string eFileEXT; 
    std::string eFileMethod;
    {
      std::string eDir = DatabaseI;
      if(NprocI > 1)
        update_path_using_nproc_iproc(eDir, NprocI, 0);
      eFileEXT = eDir + dbnameI + ".ext";
      eFileMethod = eDir + dbnameI + ".method";
    }
    MyMatrix<mpq_class> EXT = ReadMatrixFile<mpq_class>(eFileEXT);
    //
    // Now the streams
    //
    std::vector<size_t> List_shift(NprocO - IgnoreLastO, 0);
    std::vector<FileNumber *> List_FN(NprocO - IgnoreLastO, nullptr);
    std::vector<FileBool *> List_FB(NprocO - IgnoreLastO, nullptr);
    std::vector<FileFace *> List_FF(NprocO - IgnoreLastO, nullptr);
    for (int iProc = 0; iProc < NprocO - IgnoreLastO; iProc++) {
      std::string eDir = DatabaseO;
      if( NprocO > 1)
        update_path_using_nproc_iproc(eDir, NprocO, iProc);
      CreateDirectory(eDir);
      std::string eFileFN = eDir + dbnameO + ".nb";
      std::string eFileFB = eDir + dbnameO + ".fb";
      std::string eFileFF = eDir + dbnameO + ".ff";
      std::string OFileEXT = eDir + dbnameO + ".ext";
      std::string OFileMethod = eDir + dbnameO + ".method";
      List_FN[iProc] = new FileNumber(eFileFN, true);
      List_FB[iProc] = new FileBool(eFileFB);
      List_FF[iProc] = new FileFace(eFileFF, delta);
      WriteMatrixFile(OFileEXT, EXT);
      CopyOperation(eFileMethod, OFileMethod);
    }
    //
    // Now reading process by process and remapping via hash.
    //
    for (int iProc = 0; iProc < NprocI - IgnoreLastI; iProc++) {
      std::string eDir = DatabaseI;
      if(NprocI > 1)
        update_path_using_nproc_iproc(eDir, NprocI, iProc);
      std::string eFileFN = eDir + dbnameI + ".nb";
      std::string eFileFB = eDir + dbnameI + ".fb";
      std::string eFileFF = eDir + dbnameI + ".ff";
      FileNumber fn(eFileFN, false);
      size_t n_orbit = fn.getval();
      FileBool fb(eFileFB, n_orbit);
      FileFace ff(eFileFF, delta, n_orbit);
      for (size_t pos = 0; pos < n_orbit; pos++) {
        Face f = ff.getface(pos);
        bool val = fb.getbit(pos);
        size_t e_hash = mpi_get_hash(f, n_act);
        size_t iProcO = e_hash % size_t(NprocO - IgnoreLastO);
        size_t &shift = List_shift[iProcO];
        List_FB[iProcO]->setbit(shift, val);
        List_FF[iProcO]->setface(shift, f);
        shift++;
      }
    }
    //
    // Now writing the n_orbit to the fileNB
    //
    for (int iProc = 0; iProc < NprocO - IgnoreLastO; iProc++) {
      List_FN[iProc]->setval(List_shift[iProc]);
    }
    //
    // deleting the object which is critical for writing down all data to file.
    //
    for (int iProc = 0; iProc < NprocO - IgnoreLastO; iProc++) {
      delete List_FN[iProc];
      delete List_FB[iProc];
      delete List_FF[iProc];
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something went wrong in the computation, please debug\n";
    exit(e.eVal);
  }
  runtime(time1);
}
