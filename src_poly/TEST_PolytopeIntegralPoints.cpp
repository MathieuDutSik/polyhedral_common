#include "NumberTheory.h"
#include "POLY_PolytopeInt.h"

template <typename T>
MyMatrix<T> ReordListPoint(const std::vector<MyVector<T>> &ListPoint) {
  std::set<MyVector<T>> e_set;
  for (auto &ePt : ListPoint)
    e_set.insert(ePt);
  std::vector<MyVector<T>> ListPt;
  for (auto &ePt : e_set)
    ListPt.push_back(ePt);
  return MatrixFromVectorFamily(ListPt);
}

int main(int argc, char *argv[]) {
  try {
    if (argc != 2 && argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "TEST_PolytopeIntegralPoints [FAC]\n";
      std::cerr << "or\n";
      std::cerr << "TEST_PolytopeIntegralPoints [FAC] [opt]\n";
      std::cerr << "\n";
      std::cerr << "FAC: list of facets\n";
      return -1;
    }
    //
    //  std::cerr << "Reading input\n";
    //
    std::ifstream isFAC(argv[1]);
    std::string opt = "all";
    if (argc == 3)
      opt = argv[2];
    if (opt != "all" && opt != "lp" && opt != "iter") {
      std::cerr << "The allowed options are all, lp or iter\n";
      throw TerminalException{1};
    }
    using T = mpq_class;
    using Tint = int;
    MyMatrix<T> FAC = ReadMatrix<T>(isFAC);
    //
    MyMatrix<Tint> ListIntPoint1, ListIntPoint2;
    //
    if (opt == "all" || opt == "iter") {
      MyMatrix<T> EXT = cdd::DualDescription(FAC);
      std::chrono::time_point<std::chrono::system_clock> time1 =
          std::chrono::system_clock::now();
      MyMatrix<Tint> ListIntPoint1 =
          ReordListPoint(GetListIntegralPoint<T, Tint>(FAC, EXT));
      std::chrono::time_point<std::chrono::system_clock> time2 =
          std::chrono::system_clock::now();
      std::cerr << "|GetListIntegralPoint|    = "
                << std::chrono::duration_cast<std::chrono::microseconds>(time2 -
                                                                         time1)
                       .count()
                << "\n";
    }
    //
    if (opt == "all" || opt == "lp") {
      std::chrono::time_point<std::chrono::system_clock> time1 =
          std::chrono::system_clock::now();
      MyMatrix<Tint> ListIntPoint2 =
          ReordListPoint(GetListIntegralPoint_LP<T, Tint>(FAC));
      std::chrono::time_point<std::chrono::system_clock> time2 =
          std::chrono::system_clock::now();
      std::cerr << "|GetListIntegralPoint_LP| = "
                << std::chrono::duration_cast<std::chrono::microseconds>(time2 -
                                                                         time1)
                       .count()
                << "\n";
    }
    //
    if (opt == "all") {
      std::cerr << "|ListIntPoint1|=" << ListIntPoint1.rows() << "\n";
      std::cerr << "|ListIntPoint2|=" << ListIntPoint2.rows() << "\n";
      if (ListIntPoint1 != ListIntPoint2) {
        std::cerr << "We found a different set of points. Please debug\n";
        throw TerminalException{1};
      }
    }
    //
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
