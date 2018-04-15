#include "POLY_LinearProgramming.h"
#include "mpreal_related.h"
#include "NumberTheory.h"


int main(int argc, char *argv[])
{
  try {
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_process_csdp [nbDigit] [CsdpRationalIn] [CsdpRealOut]\n";
      std::cerr << "\n";
      std::cerr << "CsdpRationalIn : The file containing the rational input file\n";
      std::cerr << "CsdpRealOut    : The file containing the real output file\n";
      return -1;
    }
    //
    std::cerr << "Reading input\n";
    int nbDigit;
    sscanf(argv[1], "%d", &nbDigit);
    std::string CsdpRationalIn = argv[2];
    std::string CsdpRealOut = argv[3];
    //
    mpfr::mpreal::set_default_prec(nbDigit);
    //
    std::vector<std::string> ListLine = ReadFullFile(CsdpRationalIn);
    //
    std::ofstream os(CsdpRealOut);
    os << std::setprecision(nbDigit);
    os << ListLine[0] << "\n";
    os << ListLine[1] << "\n";
    os << ListLine[2] << "\n";
    os << ListLine[3] << "\n";
    int nbLine=ListLine.size();
    for (int iLine=4; iLine<nbLine; iLine++) {
      std::string eLine = ListLine[iLine];
      std::vector<std::string> LStr = STRING_Split(eLine, " ");
      os << LStr[0] << " " << LStr[1] << " " << LStr[2] << " " << LStr[3];
      mpq_class eVal_q;
      std::istringstream(LStr[4]) >> eVal_q;
      mpfr::mpreal eVal_mpreal = UniversalTypeConversion<mpfr::mpreal,mpq_class>(eVal_q);
      os << " " << eVal_mpreal << "\n";
    }
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
