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
    if (option == "spherical")
      allow_euclidean = false;
    if (option == "euclidean")
      allow_euclidean = true;

    MyMatrix<T> M = string_to_coxdyn_matrix<T>(str);
    std::vector<std::string> LOut;
    for (auto & eV : FindDiagramExtensions(M, allow_euclidean)) {
      MyMatrix<T> Mext = ExtendMatrix(M, eV);
      std::string stro = coxdyn_matrix_to_string(Mext);
      LOut.push_back(stro);
    }
    size_t pos_found=0;
    for (auto & stro : LOut) {
      std::cerr << "pos=" << pos_found << " diagram=" << stro << "\n";
      pos_found++;
    }
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
