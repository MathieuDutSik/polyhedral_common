#include "NumberTheory.h"
#include "coxeter_dynkin.h"


int main(int argc, char* argv[])
{
  try {
    if (argc != 5) {
      std::cerr << "COXDYN_FindExtensionsCoxMat [FileCoxMat] [opt_sph_eucl] [opt_lorentzian] [OutFile]\n";
      std::cerr << "Used for debugging the enumeration of extensions of Coxeter-Dynkin diagrams\n";
      throw TerminalException{1};
    }
    using T = mpq_class;
    std::string FileCoxMat = argv[1];
    //
    std::string opt_se = argv[2];
    bool only_spherical = true;
    bool setval = false;
    if (opt_se == "spherical") {
      only_spherical = true;
      setval = true;
    }
    if (opt_se == "euclidean") {
      only_spherical = false;
      setval = true;
    }
    if (!setval) {
      std::cerr << "We have opt_se=" << opt_se << "\n";
      std::cerr << "Only allowed values are \"spherical\" and \"euclidean\"\n";
      throw TerminalException{1};
    }
    //
    std::string opt_lor = argv[3];
    bool only_lorentzian = true;
    setval = false;
    if (opt_lor == "lorentzian") {
      only_lorentzian = true;
      setval = true;
    }
    if (opt_lor == "general") {
      only_lorentzian = false;
      setval = true;
    }
    if (!setval) {
      std::cerr << "We have opt_lor=" << opt_lor << "\n";
      std::cerr << "Only allowed values are \"lorentzian\" and \"general\"\n";
      throw TerminalException{1};
    }
    //
    DiagramSelector DS;
    DS.OnlySimplyLaced = false;
    DS.OnlyLorentzianAdmissible = only_lorentzian;
    DS.OnlySpherical = only_spherical;
    //
    MyMatrix<T> M = ReadMatrixFile<T>(FileCoxMat);
    std::vector<std::string> LOut;
    std::string symb = coxdyn_matrix_to_string(M);
    std::cerr << "symb=" << symb << "\n";
    
    std::vector<MyVector<T>> LVect = FindDiagramExtensions(M, DS);
    MyMatrix<T> Mtot = MatrixFromVectorFamily(LVect);
    std::string OutFile = argv[4];
    std::ofstream os(OutFile);
    os << "return " << StringMatrixGAP(Mtot) << ";\n";
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Something went wrong\n";
    exit(e.eVal);
  }
}
