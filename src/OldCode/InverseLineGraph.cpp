#include "Basic_file.h"
#include "GraphicalFunctions.h"

int main(int argc, char *argv[])
{
  if (argc != 2) {
    std::cerr << "Number of argument is = " << argc << "\n";
    std::cerr << "This program is used as\n";
    std::cerr << "InverseLineGraph [DataGraph]\n";
    std::cerr << "\n";
    std::cerr << "DataGraph : The file containing the graph\n";
    return -1;
  }
  //
  fprintf(stderr, "Reading input\n");
  //
  std::ifstream GRAfs;
  std::string eFile=argv[1];
  if (IsExistingFile(eFile) == false) {
    std::cerr << "Missing file\n";
    std::cerr << "eFile=" << eFile << "\n";
    exit(1);
  }
  GRAfs.open(eFile);
  GraphBitset eGR=GRAPH_Read<GraphBitset>(GRAfs);
  GRAfs.close();
  //
  std::vector<std::vector<int> > PartListLabel=InverseLineGraphConnected(eGR);
  for (auto & eVect : PartListLabel) {
    for (auto & eVal : eVect)
      std::cerr << " " << eVal;
    std::cerr << "\n";
  }
  std::cerr << "Completion of the program\n";
}
