#include "Basic_file.h"
#include "GRAPH_GraphicalFunctions.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "EnumerateDistanceSembedding [DataGraph] [S] [MaxIter] [FileOut]\n";
      std::cerr << "\n";
      std::cerr << "DataGraph : The file containing the graph\n";
      std::cerr << "S         : The value of k used for the enumeration\n";
      std::cerr << "MaxIter   : Number of iteration in the backtracking\n";
      std::cerr << "FileOut   : The file where data is saved\n";
      std::cerr << "It returns the found distance s L1-embeddings\n";
      return -1;
    }
    std::cerr << "Reading input\n";
    //
    std::ifstream GRAfs(argv[1]);
    GraphBitset eGR=GRAPH_Read<GraphBitset>(GRAfs);
    //
    int sDist;
    sscanf(argv[2], "%d", &sDist);
    //
    long MaxIter;
    sscanf(argv[3], "%ld", &MaxIter);
    //
    long iter;
    std::vector<MyMatrix<int> > ListEmbedding=GRAPH_S_Embedding(eGR, sDist, MaxIter, iter);
    //
    std::ofstream OUTfs(argv[4]);
    OUTfs << "return rec(iter:=" << iter << ", ListEmbedding:=[";
    int nbEmbed=ListEmbedding.size();
    std::cerr << "nbEmbed=" << nbEmbed << "\n";
    for (int iEmbed=0; iEmbed<nbEmbed; iEmbed++) {
      if (iEmbed>0)
	OUTfs << ",\n";
      WriteMatrixGAP(OUTfs, ListEmbedding[iEmbed]);
    }
    OUTfs << "]);\n";
    //
    std::cerr << "Completion of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
