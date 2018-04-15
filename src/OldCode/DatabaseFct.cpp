/* Mathieu Dutour Sikiric and Thomas Rehn */
#include "DatabaseFct.h"

void FuncInsert(vector<int> eList, DatabaseOrbit & OneData)
{
  set<int> eSet, fSet;
  int nb, nbOrbit, iOrbit;
  int test, eOrbitSize;
  TheGroupFormat TheStab;
  nb=eList.size();
  nbOrbit=OneData.nbOrbit;
  for (iOrbit=0; iOrbit<nbOrbit; iOrbit++)
    {
      test=TestEquivalence(OneData.GRP, eList, OneData.ListOrbit[iOrbit]);
      if (test == 1)
	return;
    }
  OneData.nbOrbit++;
  OneData.nbOrbitUndone++;
  GetStabilizer(OneData.GRP, eList, TheStab);
  eOrbitSize=(OneData.GRP.group->order())/(TheStab.group->order());
  OneData.ListStatus.push_back(1);
  OneData.ListOrbit.push_back(eList);
  OneData.ListOrbitSize.push_back(eOrbitSize);
  OneData.ListIncidence.push_back(eList.size());
  OneData.nbElementUndone=OneData.nbElementUndone + eOrbitSize;
}

void SetOrbitAsDone(int i, DatabaseOrbit & OneData)
{
  OneData.ListStatus[i]=0;
  OneData.nbElementUndone=OneData.nbElementUndone - OneData.ListOrbitSize[i];
  OneData.nbOrbitUndone--;
}

void GetUndoneOrbitMinIncidence(int &iOrbitSelect, DatabaseOrbit & OneData)
{
  int nbOrbit, iOrbit, MinIncidence;
  int IsFirst;
  nbOrbit=OneData.nbOrbit;
  IsFirst=1;
  MinIncidence=-30;
  for (iOrbit=0; iOrbit<nbOrbit; iOrbit++)
    if (OneData.ListStatus[iOrbit] == 1)
      {
	if (IsFirst == 1)
	  {
	    MinIncidence=OneData.ListIncidence[iOrbit];
	    iOrbitSelect=iOrbit;
	  }
	else
	  {
	    if (OneData.ListIncidence[iOrbit] < MinIncidence)
	      {
		MinIncidence=OneData.ListIncidence[iOrbit];
		iOrbitSelect=iOrbit;
	      }
	  }
	IsFirst=0;
      }
}

