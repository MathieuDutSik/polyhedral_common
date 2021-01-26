#include "NumberTheory.h"
#include "SimulDiophantApprox.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_DiophantApprox [FileIn] [epsilon]\n";
      std::cerr << "\n";
      std::cerr << "FileIn  : The vector file on on input\n";
      std::cerr << "epsilon : The epsilon value on input\n";
      return -1;
    }
    //    using T = mpq_class;
    using T = double;
    using Tint = mpz_class;
    //
    std::ifstream is(argv[1]);
    MyVector<T> V=ReadVector<T>(is);
    //
    std::string epsilon_str = argv[2];
    std::stringstream s(epsilon_str);
    T epsilon;
    s >> epsilon;
    //
    DiophantResult<Tint> eRes = SimultaneousDiophantineApproximation<T,Tint>(V, epsilon);
    //
    std::cerr << "Proposed solution = " << GapStringDiophantineApprox(eRes) << "\n";
    T penalty = ComputeDiophantinePenalty(V, eRes);
    std::cerr << "penalty=" << penalty << "\n";
    //
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
