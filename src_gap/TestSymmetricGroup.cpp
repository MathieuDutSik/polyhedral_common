#include "list.h"
#include "Permutation.h"
#include "StabChainMain.h"
#include "NumberTheory.h"

int main(int argc, char *argv[])
{
  try {
    using Telt = perm::DoubleSidedPerm;
    using Tint = mpz_class;
    int n=10;
    std::vector<int> ePermV1(n);
    for (int i=0; i<n; i++) {
      int iNext=i+1;
      if (iNext == n)
	iNext=0;
      ePermV1[i] = iNext;
    }
    Telt ePerm1(ePermV1);
    //
    std::vector<int> ePermV2(n);
    for (int i=0; i<n; i++)
      ePermV2[i] = i;
    ePermV2[1]=0;
    ePermV2[0]=1;
    Telt ePerm2(ePermV2);
    //
    std::vector<Telt> LGen{ePerm1, ePerm2};
    //
    gap::StabChain<Telt> S = gap::MinimalStabChain<Telt,Tint>(LGen);
    std::cerr << "S=" << S << "\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
  std::cerr << "Completion of the program\n";
}
