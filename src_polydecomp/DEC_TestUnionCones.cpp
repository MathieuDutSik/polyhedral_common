// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Decompositions.h"
// clang-format on

template <typename T>
std::vector<ConeSimpDesc<T>> ReadFamilyCones(std::string const &FileI) {
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
    bool TestPairwiseIntersection =
        BlockPROC.ListBoolValues.at("TestPairwiseIntersection");
    bool BreakConnectedComponents =
        BlockPROC.ListBoolValues.at("BreakConnectedComponents");
    //
    // Reading the polyhedral cones.
    //
    using Tcone = std::vector<ConeSimpDesc<T>>;
    Tcone l_cones = ReadFamilyCones<T>(FileI);
    //
    // Splitting by connected components
    //
    std::vector<Tcone> ll_cones;
    if (BreakConnectedComponents) {
      std::vector<std::vector<size_t>> eListList =
          ConnectedComponentsPolyhedral(l_cones);
      for (auto &eList : eListList) {
        Tcone sing_cone;
        for (auto &idx : eList) {
          sing_cone.push_back(l_cones[idx]);
        }
        ll_cones.push_back(sing_cone);
      }
    } else {
      ll_cones.push_back(l_cones);
    }
    std::cerr << "Working with |ll_cones|=" << ll_cones.size() << "\n";
    //
    // Processing data
    //
    std::vector<ConeSimpDesc<T>> l_big_cone;
    for (auto &u_cone : ll_cones) {
      std::cerr << "|u_cone|=" << u_cone.size() << "\n";
      std::optional<ConeSimpDesc<T>> opt =
        TestPolyhedralPartition(TestPairwiseIntersection, u_cone, std::cerr);
      if (!opt) {
        std::cerr << "Failed to find the full polyhedral cone\n";
        throw TerminalException{1};
      }
      ConeSimpDesc<T> cone = *opt;
      std::cerr << "|EXT|=" << cone.EXT.rows() << " |FAC|=" << cone.FAC.rows()
                << "\n";
      l_big_cone.push_back(cone);
    }
    //
    auto do_print = [&](std::ostream &os) -> void {
      os << l_big_cone.size() << "\n";
      for (auto &cone : l_big_cone) {
        WriteMatrix(os, cone.EXT);
        WriteMatrix(os, cone.FAC);
      }
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
