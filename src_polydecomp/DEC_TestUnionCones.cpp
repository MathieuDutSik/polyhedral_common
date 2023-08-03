// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "Decompositions.h"
// clang-format on

int main(int argc, char *argv[]) {
  try {
    if (argc != 5 && argc != 4) {
      std::cerr << "DEC_TestUnionCones [FileI] [FileO]\n";
      std::cerr << "or\n";
      std::cerr << "DEC_TestUnionCones [FileI]\n";
      std::cerr << "\n";
      std::cerr << "  ------ FileI -------\n";
      std::cerr << "\n";
      std::cerr << "FileI is the input file containing \n";
      std::cerr << "\n";
      std::cerr << "  ------ FileO -------\n";
      std::cerr << "\n";
      std::cerr << "FileO is the output file. If absent then the std::cerr is used\n";
      throw TerminalException{1};
    }
    using T = mpq_class;
    //
    std::string FileI = argv[1];
    std::ifstream is(FileI);
    //
    // The polyhedral cones.
    std::vector<ConeSimpDesc<T>> l_cones;
    //
    size_t n_domain;
    is >> n_domain;
    for (size_t i = 0; i < n_domain; i++) {
      std::cerr << "i=" << i << " / " << n_domain << "\n";
      MyMatrix<T> EXT = ReadMatrix<T>(is);
      std::cerr << "We have read EXT, |EXT|=" << EXT.rows() << "/" << EXT.cols()
                << "\n";
      MyMatrix<T> FAC = ReadMatrix<T>(is);
      std::cerr << "We have read FAC, |FAC|=" << FAC.rows() << "/" << FAC.cols()
                << "\n";
      ConeSimpDesc<T> eCone{EXT, FAC};
      l_cones.push_back(eCone);
    }
    //
    std::optional<ConeSimpDesc<T>> opt = TestPolyhedralPartition(l_cones);
    if (!opt) {
      std::cerr << "Failed to find the full polyhedral cone\n";
      throw TerminalException{1};
    }
    ConeSimpDesc<T> const& cone = *opt;
    //
    auto do_print = [&](std::ostream &os) -> void {
      WriteMatrix(os, cone.EXT);
      WriteMatrix(os, cone.FAC);
    };
    //
    if (argc == 4) {
      do_print(std::cerr);
    } else {
      std::string FileO = argv[4];
      std::ofstream os(FileO);
      do_print(os);
    }
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
