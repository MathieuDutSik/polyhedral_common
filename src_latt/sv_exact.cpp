/* sv.c  simple driver for shvec                              */
/* Version July 11, 2005                                      */
/* Copyright: Frank Vallentin 2005, frank.vallentin@gmail.com */

#include "ShortestUniversal.h"
//#include "Shvec_double.h"
#include "Shvec_exact.h"


[[ noreturn]] void die_sv(std::string const& last_words)
{
  std::cout << "sv.c: " << last_words << "\n";
  throw TerminalException{1};
}


int main(int argc, char *argv[])
{
  using T=mpq_class;
  using Tint=mpz_class;
  try {
    int i, j, mode, coset;
    mpq_class bound = 0;
    mode = TempShvec_globals::TEMP_SHVEC_MODE_UNDEF;
    coset = 0;
    int NeedBound=0;
    int c;
    while ((c = getopt (argc, argv, "hb:s:t:mMcgev")) != -1)
      switch (c) 
	{
	case 'h':
	  printf("Usage: sv [options] <file\n");
	  printf("-h  show this help\n");
	  printf("-bN compute vectors v with (v, v) <= N\n");
	  printf("-tN compute the first N coefficients of the theta-series\n");
	  printf("-v determine the vectors with (v-c, v-c) <= N with N defined later\n");
	  printf("-m  determine the minimum\n");
	  printf("-M  compute minimal vectors\n");
	  printf("-c  find shortest vectors in a coset\n");
	  printf("-e  do some additional checks\n");
	  return 0;
	case 'b':
	  mode = TempShvec_globals::TEMP_SHVEC_MODE_BOUND;
	  NeedBound=1;
	  break;
	case 'M':
	  mode = TempShvec_globals::TEMP_SHVEC_MODE_SHORTEST_VECTORS;
	  break;
	case 'm':
	  mode = TempShvec_globals::TEMP_SHVEC_MODE_MINIMUM;
	  break;
	case 'v':
	  mode = TempShvec_globals::TEMP_SHVEC_MODE_VINBERG;
	  coset=1;
	  NeedBound=1;
	  break;
	case 'c':
	  coset = 1;
	  break;
	case 'e':
	  break;
	default:
	  die_sv("invalid option\nTry 'sv -h' for more information.\n");
	}
    std::cerr << "main mode=" << mode << "\n";
    //
    // First reading data
    //
    int dim;
    std::cin >> dim;
    MyMatrix<T> gram_matrix(dim,dim);
    mpq_class eT;
    for (i = 0; i < dim; i++)
      for (j = 0; j <= i; j++) {
	std::cin >> eT;
	gram_matrix(i,j)=eT;
	gram_matrix(j,i)=eT;
      }
    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
	std::cerr << " " << gram_matrix(i,j);
      }
      std::cerr << "\n";
    }
    
    MyVector<T> cosetVect(dim);
    if (coset) {
      for (i = 0; i < dim; i++)
	{
	  std::cin >> eT;
	  cosetVect(i) = eT;
	}
    }
    std::cerr << "coset=" << coset << " eV =";
    for (i = 0; i < dim; i++)
      std::cerr << " " << cosetVect(i);
    std::cerr << "\n";
    if (NeedBound) {
      std::cin >> eT;
      bound=eT;
    }
    std::cerr << "NeedBound=" << NeedBound << " bound=" << bound << "\n";

    //
    // Defining info and computing with it
    //
    T_shvec_info<T,Tint> info;
    initShvecReq<T,Tint>(dim, gram_matrix, info);
    info.request.bound = bound;
    info.request.mode = mode;
    info.request.coset = cosetVect;
    info.minimum = -44;
    std::cerr << "Before computeShvec mode=" << info.request.mode << "\n";
    int result=T_computeShvec(info);
    int nbVect=info.short_vectors.size();
    std::cerr << "After computeShvec |V|=" << nbVect << "\n";
    //
    // Data output
    //
    std::cout << nbVect << "\n";
    for (i = 0; i < nbVect; i++) {
      std::cout << "[ unset ]: ";
      for (int iCol=0; iCol<dim; iCol++)
	std::cout << " " << info.short_vectors[i](iCol);
      std::cout << "\n";
    }
    return result;
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
