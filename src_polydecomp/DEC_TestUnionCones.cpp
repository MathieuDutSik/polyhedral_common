// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "Decompositions.h"
// clang-format on

template<typename T>
std::vector<ConeSimpDesc<T>> ReadFamilyCones(std::string const& FileI) {
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
  return l_cones;
}



int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    using T = mpq_class;
    //
    FullNamelist eFull = NAMELIST_TestUnionCones();
    if (argc != 2) {
      std::cerr << "DEC_TestUnionCones [FileNML]\n";
      NAMELIST_WriteNamelistFile(std::cerr, eFull, true);
      return -1;
    }
    std::string FileNML = argv[1];
    NAMELIST_ReadNamelistFile(FileNML, eFull);
    SingleBlock BlockPROC = eFull.ListBlock.at("PROC");
    //
    std::string FileI = BlockPROC.ListStringValues.at("FileI");
    std::string FileO = BlockPROC.ListStringValues.at("FileO");
    bool TestPairwiseIntersection = BlockPROC.ListBoolValues.at("TestPairwiseIntersection");
    //
    // The polyhedral cones.
    //
    std::vector<ConeSimpDesc<T>> l_cones = ReadFamilyCones<T>(FileI);
    //
    // Processing data
    //
    std::optional<ConeSimpDesc<T>> opt = TestPolyhedralPartition(TestPairwiseIntersection, l_cones);
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
