#include "NumberTheory.h"
#include "coxeter_dynkin.h"


int main(int argc, char* argv[])
{
  try {
    if (argc != 3) {
      std::cerr << "COXDYN_FindExtensions [strin] [option]\n";
      throw TerminalException{1};
    }
    using T = mpq_class;
    std::string str = argv[1];
    std::string option = argv[2];
    bool only_spherical = true;
    bool setval = false;
    if (option == "spherical") {
      only_spherical = true;
      setval = true;
    }
    if (option == "euclidean") {
      only_spherical = false;
      setval = true;
    }
    if (!setval) {
      std::cerr << "We have option=" << option << "\n";
      std::cerr << "Only allowed values are \"spherical\" and \"euclidean\"\n";
      throw TerminalException{1};
    }
    DiagramSelector DS;
    DS.OnlySimplyLaced = false;
    DS.OnlyLorentzianAdmissible = true;
    DS.OnlySpherical = only_spherical;
    //
    MyMatrix<T> M = string_to_coxdyn_matrix<T>(str);
    std::vector<std::string> LOut;
    for (auto & eV : FindDiagramExtensions(M, DS)) {
      MyMatrix<T> Mext = ExtendMatrix(M, eV);
      std::string stro = coxdyn_matrix_to_string(Mext);
      LOut.push_back(stro);
    }
    size_t pos_found=0;
    std::cerr << "input = " << str << " only_spherical=" << only_spherical << "\n";
    for (auto & stro : LOut) {
      std::cerr << "pos=" << pos_found << " diagram=" << stro << "\n";
      pos_found++;
    }
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
