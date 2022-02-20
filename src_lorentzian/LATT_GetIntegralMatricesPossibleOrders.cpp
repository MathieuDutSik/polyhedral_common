#include "NumberTheory.h"
#include "lorentzian_linalg.h"

int main(int argc, char* argv[])
{
  try {
    if (argc != 2) {
      std::cerr << "LATT_GetIntegralMatricesPossibleOrders [dim]\n";
      throw TerminalException{1};
    }
    using T = mpq_class;
    std::string dim_str = argv[1];
    T dim = ParseScalar<T>(dim_str);
    //
    std::vector<T> l_order = GetIntegralMatricesPossibleOrders(dim);
    std::cerr << "dim=" << dim << " l_order =";
    for (auto & order : l_order)
      std::cerr << " " << order;
    std::cerr << "\n";
    //
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Something went wrong\n";
    exit(e.eVal);
  }
}
