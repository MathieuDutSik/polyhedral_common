#ifndef DEFINE_STAB_CHAIN_MAIN_H
#define DEFINE_STAB_CHAIN_MAIN_H


#include "StabChain.h"
#include "stbcrand.h"


namespace gap {

  
// The main function
// Right now we do not implement the PCGS algorithm
//
template<typename Telt, typename Tint>
StabChain<Telt> StabChainOp(std::vector<Telt> const& Lgen, StabChainOptions<Tint> const& options)
{
  int degree = LargestMovedPoint( Lgen );
  if (degree > 100) {
    std::cerr << "SEARCH : Before call to StabChainRandomPermGroup\n";
    Telt TheId(degree);
    return StabChainRandomPermGroup(Lgen, TheId, options);
  }
  std::cerr << "SEARCH : Doing the ordinary Schreier Sims\n";
  int n=Lgen[0].size();
  StabChain<Telt> Stot = EmptyStabChain<Telt>(n);
  if (!IsTrivial(Lgen)) {
    int eLev=0;
    Stot.UseCycle=true;
    StabChainStrong(Stot, eLev, Lgen, options );
  }
  if (!options.reduced && options.base.size() > 0) {
    int TheLev=0;
    ExtendStabChain(Stot, TheLev, options.base);
  }
  /*
    The business with StabChainOptions look eminently dangerous and a reliable replacement
    has to be found.
    It is a record of the options chosen for the stabilizer chain that is outside of the
    variable itself!
  if (options.random > 0) {
        if IsBound( StabChainOptions( Parent( G ) ).random )  then
            options.random := Minimum( StabChainOptions( Parent( G ) ).random,
                                      options.random );
        fi;
        StabChainOptions( G ).random := options.random;
	fi;*/  
  return Stot;
}



template<typename Telt, typename Tint>
StabChain<Telt> MinimalStabChain(std::vector<Telt> const& LGen)
{
  int largMov=LargestMovedPoint(LGen);
  StabChainOptions<Tint> options = GetStandardOptions<Tint>();
  options.base = ClosedInterval(0, largMov);
  return StabChainOp(LGen, options);
}





}


#endif
