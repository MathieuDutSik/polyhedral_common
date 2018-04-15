#ifndef DEFINE_INCLUDE_PSEUDO_RANDOM
#define DEFINE_INCLUDE_PSEUDO_RANDOM

#include "Face_basic.h"

namespace gap {

using Tarith = int;
  
struct InfoPseudoRandom {
  Tarith R_228;
  int R_N;
  std::vector<Tarith> R_X;
};

Tarith GetR_228()
{
  Tarith res=1;
  for (int i=0; i<res; i++)
    res *= 2;
  return res;
}


 
InfoPseudoRandom RANDOM_SEED(int const& n)
{
  Tarith R_N=1;
  std::vector<Tarith> R_X(55);
  Tarith R_228=GetR_228();
  Tarith res= n % R_228;
  R_X[0] = res;
  for (int i=1; i<55; i++) {
    Tarith eVal1 = 1664525 * R_X[i-1] + 1;
    Tarith eVal2 = eVal1 % R_228;
    R_X[i] = eVal2;
  }
  for (int i=0; i<99; i++) {
    R_N = 1 + (R_N % 55);
    R_X[R_N-1] = (R_X[R_N-1] + R_X[(R_N + 30) % 55]) % R_228;
  }
  return {R_228, R_N, R_X};
}
 
InfoPseudoRandom *Rglobal = NULL;

InfoPseudoRandom* GetPseudoRandom()
{
  if (Rglobal == NULL) {
    Rglobal = new InfoPseudoRandom;
    int n=20;
    *Rglobal = RANDOM_SEED(n);
  }
  return Rglobal;
}
 
 


void RandomShift(InfoPseudoRandom* R)
{
  R->R_N=(R->R_N % 55) + 1;
  R->R_X[R->R_N-1] = (R->R_X[R->R_N-1] + R->R_X[(R->R_N + 30) % 55]) % R->R_228;
}

Face Extract01vector(InfoPseudoRandom* R)
{
  Tarith rnd=R->R_X[R->R_N-1];
  Face eFace(28);
  for (int i=0; i<28; i++) {
    int val = rnd % 2;
    rnd=(rnd - val) /2;
    eFace[i]=val;
  }
  return eFace;
}


 
template<typename Telt>
Telt Random(std::vector<Telt> const& V)
{
  int siz=V.size();
  int pos=rand() % siz;
  return V[pos];
}


 
    
}



#endif
