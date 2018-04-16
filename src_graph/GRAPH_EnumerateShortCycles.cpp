#include "Basic_file.h"
#include "GRAPH_GraphicalFunctions.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRAPH_EnumerateShortCycles [DataGraph]\n";
      std::cerr << "\n";
      std::cerr << "DataGraph : The file containing the graph\n";
      return -1;
    }
    //
    std::ifstream GRAfs(argv[1]);
    GraphBitset eGR=GRAPH_Read<GraphBitset>(GRAfs);
    //
    std::vector<std::vector<int>> ListCycles = GRAPH_FindAllCycles(eGR);
    int nbCycle=ListCycles.size();
    std::cerr << "nbCycle=" << nbCycle << "\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
