// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Group.h"
#include "NumberTheoryGmp.h"
#include "POLY_RecursiveDualDesc.h"
#include "POLY_RecursiveDualDesc_MPI.h"
#include "Permutation.h"

int main(int argc, char *argv[]) {
  SingletonTime time1;
  try {
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    if (argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_FindNeededBalinski [FileEXT] [FileGRP] [DatabaseInput] [NprocInput] [MinIncd] [DirectoryOutput]\n";
      std::cerr << "\n";
      std::cerr << "        FileEXT : it is the file of the polytope as input\n";
      std::cerr << "        FileGRP : it is the file of the group used for the \n";
      std::cerr << "  DatabaseInput : it is the database where the data is located\n";
      std::cerr << "     NprocInput : it is the number of processors or it can be just the string serial\n";
      std::cerr << "        MinIncd : The minimum incidence of the data to be considered\n";
      std::cerr << "DirectoryOutput : The directory where the examples are done\n";
      return -1;
    }
    std::string FileEXT = argv[1];
    MyMatrix<T> EXT = ReadMatrixFile<T>(FileEXT);
    //
    std::string FileGRP = argv[2];
    Tgroup GRP = ReadGroupFile<Tgroup>(FileGRP);
    std::vector<Telt> LGen = GRP.GeneratorsOfGroup();
    std::map<Tidx,int> LFact = GRP.factor_size();
    Tidx n_act = GRP.n_act();
    std::pair<size_t,size_t> ep = get_delta(LFact, n_act);
    size_t delta = ep.second;
    //
    std::string DatabaseI = argv[3];
    //
    int NprocI = -1;
    std::string NprocI_str = argv[4];
    if (NprocI_str == "serial") {
      NprocI = ParseScalar<int>(NprocI_str);
    }
    //
    int MinIncd = ParseScalar<int>(argv[5]);
    std::string DatabaseI = argv[6];
    //
    // Reading the group
    //
    std::ifstream GRPfs(FileGRP);
    Tgroup GRP = ReadGroup<Tgroup>(GRPfs);
    std::map<Tidx,int> LFact = GRP.factor_size();
    Tidx n_act = GRP.n_act();
    std::pair<size_t,size_t> ep = get_delta(LFact, n_act);
    size_t delta = ep.second;
    //
    // The function for outputting
    //
    size_t iPolytope = 0;
    auto write_polytope=[&](Face const& f) -> void {
      MyMatrix<T> EXTred = SelectRow(EXT, f);
      using Tidx_value = uint16_t;
      WeightMatrix<true, T, Tidx_value> WMat = GetWeightMatrix<T, Tidx_value>(EXT);
      Tgroup GRPred = GetStabilizerWeightMatrix<T, Tgr, Tgroup, Tidx_value>(WMat);
      //
      std::string FileGRP_out = DirectoryOutput + "GRP_" + std::to_stinrg(iPolytope);
      std::ofstream osGRP(FileGRP_out);
      WriteGroup(osGRP, GRPred);
      //
      std::string FileEXT_out = DirectoryOutput + "EXT_" + std::to_stinrg(iPolytope);
      std::ofstream osEXT(FileEXT_out);
      WriteMatrix(osEXT, EXTred);
      //
      iPolytope++;
    };
    auto f_process_database=[&](std::string const& eDir) -> void {
      std::string eFileNB = ePrefix + ".nb";
      std::string eFileFB = ePrefix + ".fb";
      std::string eFileFF = ePrefix + ".ff";
      FileNumber* fn = new FileNumber(eFileNB, false);
      size_t n_orbit = fn->getval();
      FileBool* fb = new FileBool(eFileFB, n_orbit);
      FileFace* ff = new FileFace(eFileFF, delta, n_orbit);

      for (size_t i_orbit=0; i_orbit<n_orbit; i_orbit++) {
        bool status = fb->getbit(i_orbit);
        if (!status) {
          Face f = ff->getface(i_orbit);
          int eIncd = f.count();
          if (eIncd >= MinIncd) {
            vectface vf_orbit = OrbitFace(f, LGen);
            MyMatrix<T> FAC_undone = GetVertexSet_from_vectface(bb.EXT, vf_orbit);
            bool test = EvaluationConnectednessCriterion_PreKernel(EXT, GRP, FAC_undone, vf_orbit, std::cerr);
            if (!test) {
              std::cerr << "i_orbit=" << i_orbit << " needs to be done eIncd=" << eIncd << "\n";
              if (DirectoryOutput != "unset") {
                write_polytope(f);
              }
            }
          }
        }
      }
      delete fn;
      delete fb;
      delete ff;
    };
    //
    // Now reading process by process and remapping via hash.
    //
    if (NprocI == -1) {
      std::string ePrefix = DatabaseI;
      f_process_database(ePrefix);
    } else {
      for (int iProc=0; iProc<NprocI; iProc++) {
        std::string ePrefix = DatabaseI;
        update_path_using_nproc_iproc(ePrefix, NprocI, iProc);
        ePrefix += "database";
        f_process_database(ePrefix);
      }
    }
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const &e) {
    std::cerr << "Something went wrong in the computation, please debug\n";
    exit(e.eVal);
  }
  runtime(time1);
}
