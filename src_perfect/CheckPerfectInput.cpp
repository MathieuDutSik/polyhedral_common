#include "Permlib_specific.h"
#include "Temp_PerfectForm.h"
#include "Permutation.h"
#include "Group.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 2 && argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CheckPerfectInput sizCheck [perfectinput]\n";
      std::cerr << "CheckPerfectInput sizCheck < stdin\n";
      throw TerminalException{1};
    }
    using T=mpq_class;
    int n;
    sscanf(argv[1], "%d", &n);
    //
    auto check_stream=[](std::istream& is) -> void {
      MyMatrix<T> G = ReadMatrix<T>(is);
      int n = G.rows();
      MyMatrix<T> SHV = ReadMatrix<T>(is);
      int n_SHV = SHV.rows();
      int pDim = (n * (n+1)) / 2;
      MyMatrix<T> PerfectCone(n_SHV, pDim);
      for (int i_SHV=0; i_SHV<n_SHV; i_SHV++) {
        size_t pos = 0;
        for (int i=0; i<n; i++)
          for (int j=i; j<n; j++) {
            PerfectCone(i_SHV, pos) = SHV(i_SHV,i) * SHV(i_SHV,j);
            pos++;
          }
      }

      size_t n_facet;
      is >> n_facet;
      size_t n_facet_work = n_facet;
      if (n != -1)
        n_facet_work = n;
      //
      for (size_t i_facet=0; i_facet<n_facet_work; i_facet++) {
        std::cerr << "i_facet=" << i_facet << " / " << n_facet << " / " << n_facet_work << "\n";
        T val;
        is >> val;
        std::vector<size_t> V = Convert_T_To_Set(val);
        T val2 = Convert_Set_To_T<T>(V);
        if (val != val2) {
          std::cerr << "The conversion val to val2 is bugged\n";
          throw TerminalException{1};
        }
        Face f(n_SHV);
        for (auto & eV : V)
          f[eV] = 1;
        TestFacetness(PerfectCone, f);
      }
    };


    if (argc == 2) {
      check_stream(std::cin);
    } else {
      std::string eFileName=argv[2];
      std::ifstream is(eFileName);
      check_stream(is);
    }
    std::cerr << "Completion of the program\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Something went wrong\n";
    exit(e.eVal);
  }
}
