// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "coxeter_dynkin.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 4) {
      std::cerr
          << "COXDYN_FindExtensions [strin] [opt_sph_eucl] [opt_lorentzian]\n";
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
    MyMatrix<T> M = string_to_coxdyn_matrix<T>(str);
    std::vector<std::string> LOut;
    std::vector<MyVector<T>> LVect = FindDiagramExtensions(M, DS);
    for (auto &eV : LVect) {
      MyMatrix<T> Mext = ExtendMatrix(M, eV);
      std::string stro = coxdyn_matrix_to_string(Mext);
      LOut.push_back(stro);
    }
    size_t len = LVect.size();
    std::cerr << "input = " << str << " only_spherical=" << only_spherical
              << "\n";
    for (size_t i = 0; i < len; i++)
      std::cerr << "i=" << i << " diag=" << LOut[i]
                << " v=" << StringVector(LVect[i]) << "\n";
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something went wrong\n";
    exit(e.eVal);
  }
  runtime(time);
}
