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
    bool allow_euclidean = true;
    bool setval = false;
    if (option == "spherical") {
      allow_euclidean = false;
      setval = true;
    }
    if (option == "euclidean") {
      allow_euclidean = true;
      setval = true;
    }
    if (!setval) {
      std::cerr << "We have option=" << option << "\n";
      std::cerr << "Only allowed values are spherical or euclidean\n";
      throw TerminalException{1};
    }

    MyMatrix<T> M = string_to_coxdyn_matrix<T>(str);
    std::vector<std::string> LOut;
    for (auto & eV : FindDiagramExtensions(M, allow_euclidean)) {
      MyMatrix<T> Mext = ExtendMatrix(M, eV);
      std::string stro = coxdyn_matrix_to_string(Mext);
      LOut.push_back(stro);
    }
    size_t pos_found=0;
    std::cerr << "input = " << str << " allow_euclidean=" << allow_euclidean << "\n";
    for (auto & stro : LOut) {
      std::cerr << "pos=" << pos_found << " diagram=" << stro << "\n";
      pos_found++;
    }
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
