#include "NumberTheory.h"
#include "coxeter_dynkin.h"


int main(int argc, char* argv[])
{
  try {
    if (argc != 4) {
      std::cerr << "COXDYN_FindExtensions [strin] [opt_sph_eucl] [opt_lorentzian]\n";
      throw TerminalException{1};
    }
    using T = mpq_class;
    std::string str = argv[1];
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
    std::string opt_lor = argv[2];
    bool only_lorentzian = true;
    setval = false;
    if (opt_lor == "spherical") {
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
