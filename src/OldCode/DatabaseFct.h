#if !defined DATABASE_INCLUDE_SYMPOL
#include "GroupFct.h"
#include "gmpxx.h"
#include <vector>
#include <set>
using namespace std;

typedef struct {
  int nbOrbit;
  int nbOrbitUndone;
  vector<vector<int> > ListOrbit;
  vector<int> ListStatus;
  vector<int> ListIncidence;
  mpz_class nbElementUndone;
  vector<mpz_class> ListOrbitSize;
  TheGroupFormat GRP;
} DatabaseOrbit;
#endif
#define DATABASE_INCLUDE_SYMPOL
