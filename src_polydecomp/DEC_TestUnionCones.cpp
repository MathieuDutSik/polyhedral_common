// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "Decompositions.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 2 && argc != 3) {
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
    std::string FileO = "stderr";
    if (argc == 3)
      FileO = argv[2];
    //
    // The polyhedral cones.
    //
    std::ifstream is(FileI);
    size_t n_domain;
    is >> n_domain;
    std::vector<ConeSimpDesc<T>> l_cones;
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
    if (FileO == "stderr") {
      do_print(std::cerr);
    } else {
      if (FileO == "stderr") {
        do_print(std::cout);
      } else {
        std::ofstream os(FileO);
        do_print(os);
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something went wrong\n";
    exit(e.eVal);
  }
  runtime(time1);
}
