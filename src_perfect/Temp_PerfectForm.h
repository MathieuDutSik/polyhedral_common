#ifndef TEMP_PERFECT_FORM_INCLUDE
#define TEMP_PERFECT_FORM_INCLUDE

#include "POLY_ThreadDualDescription.h"
#include "Temp_Positivity.h"
#include "Temp_Tspace_General.h"
#include "MatrixGroup.h"
#include "Parallel_Classes.h"



template<typename T, typename Tint>
struct NakedPerfect {
  MyMatrix<T> eGram;
  MyMatrix<Tint> SHV;
  MyMatrix<Tint> SHVred;
  MyMatrix<T> PerfDomEXT;
  std::vector<std::vector<int>> ListBlock;
  std::vector<int> ListPos;
};



template<typename T, typename Tint>
MyMatrix<T> GetNakedPerfectConeClassical(MyMatrix<Tint> const& M)
{
  int nbRow=M.rows();
  int n=M.cols();
  int dimSymm=n*(n+1)/2;
  MyMatrix<T> RetMat(nbRow, dimSymm);
  for (int iRow=0; iRow<nbRow; iRow++) {
    MyVector<Tint> V = GetMatrixRow(M, iRow);
    MyMatrix<T> M(n,n);
    for (int u=0; u<n; u++)
      for (int v=0; v<n; v++)
	M(u,v) = V(u) * V(v);
    MyVector<T> Vm=SymmetricMatrixToVector(M);
    AssignMatrixRow(RetMat, iRow, Vm);
  }
  return RetMat;
}



template<typename T, typename Tint>
NakedPerfect<T,Tint> GetNakedPerfectCone(LinSpaceMatrix<T> const&LinSpa, MyMatrix<T> const& eGram, Tshortest<T,Tint> const& RecSHV)
{
  int nbSHV=RecSHV.SHV.rows();
  std::vector<int> ListPos(nbSHV);
  int nbMat=LinSpa.ListMat.size();
  int n=eGram.rows();
  MyMatrix<T> RyshkovLoc(nbSHV, nbMat);
  for (int iSHV=0; iSHV<nbSHV; iSHV++) {
    MyVector<Tint> eVect=GetMatrixRow(RecSHV.SHV, iSHV);
    for (int iMat=0; iMat<nbMat; iMat++) {
      T eSum=EvaluationQuadForm<T,Tint>(LinSpa.ListMat[iMat], eVect);
      RyshkovLoc(iSHV, iMat)=eSum;
    }
  }
  //  std::cerr << "RyshkovLoc=\n";
  //  WriteMatrix(std::cerr, RyshkovLoc);
  int nbBlock=0;
  std::vector<std::vector<int>> ListBlock;
  for (int iSHV=0; iSHV<nbSHV; iSHV++) {
    MyVector<T> eVec1=GetMatrixRow(RyshkovLoc, iSHV);
    bool IsMatch=false;
    for (int iBlock=0; iBlock<nbBlock; iBlock++)
      if (!IsMatch) {
	int jSHV=ListBlock[iBlock][0];
	MyVector<T> eVec2=GetMatrixRow(RyshkovLoc, jSHV);
	bool test=IsVectorPositiveMultiple(eVec1, eVec2);
	if (test) {
	  IsMatch=true;
	  ListBlock[iBlock].push_back(iSHV);
	  ListPos[iSHV]=iBlock;
	}
      }
    if (!IsMatch) {
      ListBlock.push_back({iSHV});
      ListPos[iSHV]=nbBlock;
      nbBlock++;
    }
  }
  std::cerr << "nbBlock=" << nbBlock << "\n";
  MyMatrix<T> PerfDomEXT(nbBlock, nbMat);
  MyMatrix<Tint> SHVred(nbBlock, n);
  for (int iBlock=0; iBlock<nbBlock; iBlock++) {
    int iSHV=ListBlock[iBlock][0];
    MyVector<Tint> eVect=GetMatrixRow(RecSHV.SHV, iSHV);
    AssignMatrixRow(SHVred, iBlock, eVect);
    for (int iMat=0; iMat<nbMat; iMat++)
      PerfDomEXT(iBlock, iMat)=RyshkovLoc(iSHV, iMat);
  }
  return {eGram, RecSHV.SHV, SHVred, PerfDomEXT, ListBlock, ListPos};
}


struct PerfEquivInfo {
  int iOrbit;
  MyMatrix<int> eMatEquiv;
  Face eInc;
};


template<typename T, typename Tint, typename Tgroup>
struct SinglePerfect {
  MyMatrix<T> eGram;
  MyMatrix<Tint> SHV;
  MyMatrix<T> PerfDomEXT;
  Tgroup PerfDomGRP;
  std::vector<std::vector<int>> ListBlock;
  std::vector<PerfEquivInfo> ListEquivInfo;
  int eStatus;
  int eCons;
};



template<typename T, typename Tint, typename Tgroup>
Tgroup MapLatticeGroupToConeGroup(NakedPerfect<T,Tint> const& eNaked, Tgroup const& GRPshv)
{
  using Telt=typename Tgroup::Telt;
  int nbBlock=eNaked.ListBlock.size();
  std::vector<Telt> ListGen;
  std::vector<Telt> LGen = GRPshv.GeneratorsOfGroup();
  for (auto & eGen : LGen) {
    std::vector<int> v(nbBlock);
    for (int iBlock=0; iBlock<nbBlock; iBlock++) {
      int iSHV=eNaked.ListBlock[iBlock][0];
      int jSHV=OnPoints(iSHV, eGen);
      int jBlock=eNaked.ListPos[jSHV];
      v[iBlock]=jBlock;
    }
    ListGen.push_back(Telt(v));
  }
  return Tgroup(ListGen, nbBlock);
}




template<typename T,typename Tint, typename Tgroup>
SinglePerfect<T,Tint,Tgroup> GetPerfectCone(LinSpaceMatrix<T> const&LinSpa, MyMatrix<T> const& eGram, Tshortest<T,int> const& RecSHV)
{
  NakedPerfect<T,Tint> eNaked=GetNakedPerfectCone(LinSpa, eGram, RecSHV);
  Tgroup TheGRPshv=PERF_Automorphism(LinSpa, eGram, RecSHV.SHV);
  Tgroup PerfDomGRP=MapLatticeGroupToConeGroup(eNaked, TheGRPshv);
  return {eGram, RecSHV.SHV, eNaked.PerfDomEXT, PerfDomGRP, eNaked.ListBlock, {}, 0, 0};
}


template<typename T, typename Tint, typename Tgroup>
EquivTest<MyMatrix<Tint>> PERF_TestEquivalence(LinSpaceMatrix<T> const&LinSpa,
				       MyMatrix<T> const&ePerf1, MyMatrix<T> const&ePerf2,
				       MyMatrix<Tint> const& SHV1, MyMatrix<Tint> const& SHV2)
{
  using Telt=typename Tgroup::Telt;
  MyMatrix<T> T_SHV1=ConvertMatrixUniversal<T,Tint>(SHV1);
  MyMatrix<T> T_SHV2=ConvertMatrixUniversal<T,Tint>(SHV2);
  WeightMatrix<std::vector<T>> WMat1=GetWeightMatrix_ListComm(T_SHV1, ePerf1, LinSpa.ListComm);
  WeightMatrix<std::vector<T>> WMat2=GetWeightMatrix_ListComm(T_SHV2, ePerf2, LinSpa.ListComm);
  EquivTest<Telt> eResEquiv=GetEquivalenceAsymmetricMatrix<std::vector<T>,T,Telt>(WMat1, WMat2);
  if (!eResEquiv.TheReply) {
    return {false, {}};
  }
  MyMatrix<T> M3=RepresentVertexPermutation(T_SHV1, T_SHV2, eResEquiv.TheEquiv);
  if (IsIntegralMatrix(M3)) {
    MyMatrix<Tint> eMatEquiv=ConvertMatrixUniversal<Tint,T>(M3);
    return {true, eMatEquiv};
  }
  std::cerr << "Need to write some code here\n";
  throw TerminalException{1};
}






struct QueryEquivInfo {
  bool result;
  int nbOrbit;
  PerfEquivInfo eEquiv;
};




template<typename T, typename Tint, typename Tgroup>
struct ListPerfectForm {
public:
  QueryEquivInfo IsPresentNoLock(int const& iOrbitStart,
				 LinSpaceMatrix<T> const&LinSpa,
				 MyMatrix<T> const&eGram,
				 MyMatrix<int> const& SHV)
  {
    int nbOrbit=ListPerfect.size();
    Face eTrivFace;
    for (int iOrbit=iOrbitStart; iOrbit<nbOrbit; iOrbit++) {
      MyMatrix<T> oldGram=GetGram(iOrbit);
      MyMatrix<int> oldSHV=GetSHV(iOrbit);
      EquivTest<MyMatrix<int>> eLinSpaRes=PERF_TestEquivalence<T,Tint,Tgroup>(LinSpa, eGram, oldGram, SHV, oldSHV);
      if (eLinSpaRes.TheReply) {
	PerfEquivInfo eEquiv{iOrbit, eLinSpaRes.TheEquiv, eTrivFace};
	return {true, nbOrbit, eEquiv};
      }
    }
    return {false, nbOrbit, {-1, {}, eTrivFace}};
  }
  PerfEquivInfo InsertForm(LinSpaceMatrix<T> const&LinSpa,
			   MyMatrix<T> const&eGram)
  {
    std::vector<PerfEquivInfo> eListEquivInfo;
    Face eTrivFace;
    int n=LinSpa.n;
    Tshortest<T,int> RecSHV=T_ShortestVector<T,int>(eGram);
    QueryEquivInfo eQuery=IsPresentNoLock(0, LinSpa, eGram, RecSHV.SHV);
    if (eQuery.result)
      return eQuery.eEquiv;
    std::lock_guard<std::mutex> lk(mul);
    QueryEquivInfo fQuery=IsPresentNoLock(eQuery.nbOrbit, LinSpa, eGram, RecSHV.SHV);
    if (fQuery.result)
      return fQuery.eEquiv;
    SinglePerfect<T,Tint,Tgroup> eSing=GetPerfectCone(LinSpa, eGram, RecSHV);
    ListPerfect.push_back(eSing);
    int len=ListPerfect.size();
    MyMatrix<int> eMatEquiv=IdentityMat<int>(n);
    return {len, eMatEquiv, eTrivFace};
  }
  std::vector<int> GetListOrbitWorkEff(std::vector<int> const & ListOrbWork)
  {
    std::vector<int> ListOrbitWorkEff;
    std::lock_guard<std::mutex> lk(mul);
    for (auto eOrb : ListOrbWork) {
      if (ListPerfect[eOrb].eCons == 0) {
	ListPerfect[eOrb].eCons=1;
	ListOrbitWorkEff.push_back(eOrb);
      }
    }
    return ListOrbitWorkEff;
  }
  void SetStatus(int const& eOrb, int const& eVal)
  {
    ListPerfect[eOrb].eStatus=eVal;
  }
  MyMatrix<int> GetSHV(int eOrb) const
  {
    return ListPerfect[eOrb].SHV;
  }
  MyMatrix<T> GetPerfDomEXT(int const& eOrb) const
  {
    return ListPerfect[eOrb].PerfDomEXT;
  }
  Tgroup GetPerfDomGRP(int const& eOrb) const
  {
    return ListPerfect[eOrb].PerfDomGRP;
  }
  MyMatrix<T> GetGram(int eOrb) const
  {
    return ListPerfect[eOrb].eGram;
  }
  void InsertEquivInfo(int eOrb, PerfEquivInfo const& eEquiv)
  {
    ListPerfect[eOrb].ListEquivInfo.push_back(eEquiv);
  }
  int GetNbPerf() const
  {
    int len=ListPerfect.size();
    return len;
  }
  std::vector<std::vector<int>> VoronoiAlgo_THR_GetPartition(
           int NbThr, int & nbConsTodo,
	   int & IsFinished)
  {
    int iOrbit, iThr;
    int nbOrbit=ListPerfect.size();
    std::vector<std::vector<int>> ThePartition;
    nbConsTodo=0;
    if (NbThr == 0) {
      std::cerr << "Error in NbThr\n";
      std::cerr << "NbThr=" << NbThr << "\n";
      throw TerminalException{1};
    }
    std::vector<int> eNewList;
    for (iThr=0; iThr<NbThr; iThr++)
      ThePartition.push_back(eNewList);
    iThr=0;
    IsFinished=1;
    for (iOrbit=0; iOrbit<nbOrbit; iOrbit++) {
      if (ListPerfect[iOrbit].eCons == 0) {
	nbConsTodo++;
	ThePartition[iThr].push_back(iOrbit);
	iThr++;
	if (iThr == NbThr)
	  iThr=0;
      }
      if (ListPerfect[iOrbit].eStatus == 0)
	IsFinished=0;
    }
    return ThePartition;
  }
private:
  std::mutex mul;
  std::vector<SinglePerfect<T,Tint,Tgroup>> ListPerfect;
};



template<typename T, typename Tint>
struct RecShort {
  std::function<bool(MyMatrix<T> const&)> IsAdmissible;
  std::function<Tshortest<T,Tint>(MyMatrix<T> const&)> ShortestFunction;
};

template<typename Tint>
bool TestInclusionSHV(MyMatrix<Tint> const& TheSHVbig, MyMatrix<Tint> const& TheSHVsma)
{
  int nbRowSma=TheSHVsma.rows();
  int nbRowBig=TheSHVbig.rows();
  int n=TheSHVsma.cols();
  for (int iRowSma=0; iRowSma<nbRowSma; iRowSma++) {
    bool WeMatch=false;
    for (int iRowBig=0; iRowBig<nbRowBig; iRowBig++)
      if (!WeMatch) {
	int SumErr=0;
	for (int i=0; i<n; i++) {
	  Tint eVal=TheSHVbig(iRowBig, i);
	  Tint fVal=TheSHVsma(iRowSma, i);
	  if (eVal != fVal)
	    SumErr++;
	}
	if (SumErr == 0)
	  WeMatch=true;
      }
    if (!WeMatch)
      return false;
  }
  return true;
}


#undef DEBUG_FLIP
//#define DEBUG_FLIP
template<typename T, typename Tint>
std::pair<MyMatrix<T>,Tshortest<T,Tint>> Kernel_Flipping_Perfect(RecShort<T, Tint> const& eRecShort, MyMatrix<T> const& eMatIn, MyMatrix<T> const& eMatDir)
{
  std::vector<MyMatrix<T>> ListMat;
  std::vector<Tshortest<T,Tint>> ListShort;
  // Memoization procedure
  auto RetriveShortestDesc=[&](MyMatrix<T> const& eMat) -> Tshortest<T,Tint> {
    int len=ListMat.size();
    for (int i=0; i<len; i++)
      if (ListMat[i] == eMat)
        return ListShort[i];
    Tshortest<T,Tint> RecSHV=eRecShort.ShortestFunction(eMat);
    ListMat.push_back(eMat);
    ListShort.push_back(RecSHV);
    return RecSHV;
  };
  Tshortest<T,Tint> const RecSHVperf=RetriveShortestDesc(eMatIn);
#ifdef DEBUG_FLIP
  std::cerr << "Kernel_Flipping_Perfect : SHVinformation=\n";
  int nbSHV=RecSHVperf.SHV.rows();
  int n=RecSHVperf.SHV.cols();
  int nbZero_sumMat=0;
  std::vector<int> ListScal(nbSHV);
  for (int iSHV=0; iSHV<nbSHV; iSHV++) {
    MyVector<Tint> V = GetMatrixRow(RecSHVperf.SHV, iSHV);
    T sumMatIn  = EvaluationQuadForm(eMatIn, V);
    T sumMatDir = EvaluationQuadForm(eMatDir, V);
    int eVal=1;
    if (sumMatDir == 0) {
      nbZero_sumMat++;
      eVal=0;
    }
    ListScal[iSHV]=eVal;
    std::cerr << "iSHV=" << iSHV << " V=";
    for (int i=0; i<n; i++)
      std::cerr << V(i) << " ";
    std::cerr << "sumMatIn=" << sumMatIn << " sumMatDir=" << sumMatDir << "\n";
  }
  MyMatrix<Tint> SHVface(nbZero_sumMat, n);
  int idx=0;
  for (int iSHV=0; iSHV<nbSHV; iSHV++) {
    if (ListScal[iSHV] == 0) {
      for (int i=0; i<n; i++)
	SHVface(idx, i) = RecSHVperf.SHV(iSHV,i);
      idx++;
    }
  }
  std::cerr << "SHVface=\n";
  WriteMatrix(std::cerr, SHVface);
  std::cerr << "nbZero_sumMat=" << nbZero_sumMat << "\n";
  MyMatrix<T> ConeClassicalInput = GetNakedPerfectConeClassical<T>(RecSHVperf.SHV);
  std::cerr << "RankMat(ConeClassicalInput)=" << RankMat(ConeClassicalInput) << "\n";
  MyMatrix<T> ConeClassicalFace = GetNakedPerfectConeClassical<T>(SHVface);
  std::cerr << "RankMat(ConeClassicalFace)=" << RankMat(ConeClassicalFace) << "\n";
#endif
  T TheUpperBound=1;
  T TheLowerBound=0;
#ifdef DEBUG_FLIP
  int iterLoop=0;
#endif
  while(true) {
    MyMatrix<T> Qupp = eMatIn + TheUpperBound*eMatDir;
    bool test=eRecShort.IsAdmissible(Qupp);
#ifdef DEBUG_FLIP
    iterLoop++;
    std::cerr << "iterLoop=" << iterLoop << " TheUpperBound=" << TheUpperBound << " TheLowerBound=" << TheLowerBound << " test=" << test << "\n";
#endif
    if (!test)
      TheUpperBound=(TheUpperBound + TheLowerBound)/2;
    else {
      Tshortest<T,Tint> RecSHV=RetriveShortestDesc(Qupp);
#ifdef DEBUG_FLIP
      std::cerr << "ITER: RecSHV.eMin=" << RecSHV.eMin << "\n";
      std::cerr << "ITER: RecSHV.SHV=\n";
      WriteMatrix(std::cerr, RecSHV.SHV);
#endif
      if (RecSHV.eMin == RecSHVperf.eMin) {
	T nLow=TheUpperBound;
	T nUpp=2*TheUpperBound;
	TheLowerBound=nLow;
	TheUpperBound=nUpp;
      }
      else
	break;
    }
  }
#ifdef DEBUG_FLIP
  std::cerr << "FIRST LOOP FINISHED TheUpperBound=" << TheUpperBound << " TheLowerBound=" << TheLowerBound << "\n";
#endif
  while(true) {
#ifdef DEBUG_FLIP
    std::cerr << "Now TheUpperBound=" << TheUpperBound << " TheLowerBound=" << TheLowerBound << "\n";
#endif
    MyMatrix<T> Qlow = eMatIn + TheLowerBound*eMatDir;
    MyMatrix<T> Qupp = eMatIn + TheUpperBound*eMatDir;
    Tshortest<T,Tint> RecSHVlow=RetriveShortestDesc(Qlow);
    Tshortest<T,Tint> RecSHVupp=RetriveShortestDesc(Qupp);
#ifdef DEBUG_FLIP
    std::cerr << "RecSHVupp.eMin=" << RecSHVupp.eMin << " RecSHVlow.eMin=" << RecSHVlow.eMin << "\n";
    MyMatrix<T> ConeClassicalLow = GetNakedPerfectConeClassical<T>(RecSHVlow.SHV);
    std::cerr << "RankMat(ConeClassicalLow)=" << RankMat(ConeClassicalLow) << "\n";
    MyMatrix<T> ConeClassicalUpp = GetNakedPerfectConeClassical<T>(RecSHVupp.SHV);
    std::cerr << "RankMat(ConeClassicalUpp)=" << RankMat(ConeClassicalUpp) << "\n";
#endif
    bool test1=RecSHVupp.eMin == RecSHVperf.eMin;
    bool test2=TestInclusionSHV(RecSHVperf.SHV, RecSHVlow.SHV);
    if (test1) {
#ifdef DEBUG_FLIP
      std::cerr << "Return Qupp\n";
#endif
      return {std::move(Qupp), std::move(RecSHVupp)};
    }
    if (!test2) {
#ifdef DEBUG_FLIP
      std::cerr << "Qperf=\n";
      WriteMatrix(std::cerr, eMatIn);
      std::cerr << "Qlow=\n";
      WriteMatrix(std::cerr, Qlow);
      //
      std::cerr << "RecSHVperf.SHV=\n";
      WriteMatrix(std::cerr, RecSHVperf.SHV);
      std::cerr << "RecSHVlow.SHV=\n";
      WriteMatrix(std::cerr, RecSHVlow.SHV);
      std::cerr << "Return Qlow\n";
#endif
      return {std::move(Qlow),std::move(RecSHVlow)};
    }
    T TheGamma=(TheLowerBound + TheUpperBound)/2;
    MyMatrix<T> Qgamma = eMatIn + TheGamma*eMatDir;
#ifdef DEBUG_FLIP
    std::cerr << "Qgamma=\n";
    WriteMatrix(std::cerr, Qgamma);
#endif
    Tshortest<T,Tint> RecSHVgamma=RetriveShortestDesc(Qgamma);
#ifdef DEBUG_FLIP
    std::cerr << "|RecSHVgamma.SHV|=" << RecSHVgamma.SHV.rows() << "\n";
    WriteMatrix(std::cerr, RecSHVgamma.SHV);
    MyMatrix<T> ConeClassicalGamma = GetNakedPerfectConeClassical<T>(RecSHVgamma.SHV);
    std::cerr << "RankMat(ConeClassicalGamma)=" << RankMat(ConeClassicalGamma) << "\n";
#endif
    if (RecSHVgamma.eMin >= RecSHVperf.eMin) {
      TheLowerBound=TheGamma;
    }
    else {
#ifdef DEBUG_FLIP
      std::cerr << "Assigning TheUpperBound to TheGamma=" << TheGamma << "\n";
#endif
      TheUpperBound=TheGamma;
      int nbRowGamma=RecSHVgamma.SHV.rows();
      for (int iRowGamma=0; iRowGamma<nbRowGamma; iRowGamma++) {
	MyVector<Tint> eVectShort=GetMatrixRow(RecSHVgamma.SHV, iRowGamma);
#ifdef DEBUG_FLIP
	std::cerr << "iRowGamma=" << iRowGamma << " / " << nbRowGamma << " eVectShort=";
	WriteVector(std::cerr, eVectShort);
#endif
	T rVal=EvaluationQuadForm<T,Tint>(eMatDir, eVectShort);
	T qVal=EvaluationQuadForm<T,Tint>(eMatIn, eVectShort);
	if (rVal < 0) {
	  T TheVal=(RecSHVperf.eMin - qVal)/rVal;
	  if (TheVal < TheUpperBound) {
#ifdef DEBUG_FLIP
	    std::cerr << "iRowGamma=" << iRowGamma << " Assigning TheUpperBound to TheVal=" << TheVal << "\n";
	    std::cerr << "rVal=" << rVal << " qVal=" << qVal << " RecSHVperf.eMin=" << RecSHVperf.eMin << "\n";
#endif
	    TheUpperBound=TheVal;
	  }
	}
      }
    }
  }
}


template<typename T, typename Tint>
std::pair<MyMatrix<T>,Tshortest<T,Tint>> Flipping_Perfect(MyMatrix<T> const& eMatIn, MyMatrix<T> const& eMatDir)
{
  std::function<bool(MyMatrix<T> const&)> IsAdmissible=[](MyMatrix<T> const& eMat) -> bool {
    return IsPositiveDefinite<T>(eMat);
  };
  std::function<Tshortest<T,Tint>(MyMatrix<T> const&)> ShortestFunction=[](MyMatrix<T> const& eMat) -> Tshortest<T,Tint> {
    return T_ShortestVector<T,Tint>(eMat);
  };
  RecShort<T,Tint> eRecShort{IsAdmissible, ShortestFunction};
  return Kernel_Flipping_Perfect(eRecShort, eMatIn, eMatDir);
}



template<typename T, typename Tint, typename Tgroup>
void VoronoiAlgo_THR_BlockTreatment(MainProcessor &MProc, int TheId,
				    DataBank<PolyhedralEntry<T,Tgroup>> &TheBank,
				    LinSpaceMatrix<T> const& LinSpa,
				    ListPerfectForm<T,Tint,Tgroup> &ListPerf,
				    PolyHeuristic<mpz_class> const& AllArr,
				    std::vector<int> const& ListOrbWork)
{
  std::vector<int> ListOrbitWorkEff=ListPerf.GetListOrbitWorkEff(ListOrbWork);
  for (auto& eOrb : ListOrbitWorkEff) {
    MyMatrix<T> eGram=ListPerf.GetGram(eOrb);
    MyMatrix<T> PerfDomEXT=ListPerf.GetPerfDomEXT(eOrb);
    Tgroup PerfDomGRP=ListPerf.GetPerfDomGRP(eOrb);
    vectface TheOutput=DUALDESC_THR_AdjacencyDecomposition(MProc,
             TheId, TheBank,
	     PerfDomEXT, PerfDomGRP,
	     AllArr);
    for (auto& eOrbB : TheOutput) {
      MyVector<T> eVectOrb=FindFacetInequality(PerfDomEXT, eOrbB);
      MyMatrix<T> DirMat=LINSPA_GetMatrixInTspace(LinSpa, eVectOrb);
      MyMatrix<T> NewPerf=Flipping_Perfect<T,int>(eGram, DirMat).first;
      PerfEquivInfo eEquiv=ListPerf.InsertForm(LinSpa, NewPerf);
      eEquiv.eInc=eOrbB;
      ListPerf.InsertEquivInfo(eOrb, eEquiv);
    }
    ListPerf.SetStatus(eOrb, 1);
  }
}









template<typename T>
MyMatrix<T> GetOnePerfectForm(LinSpaceMatrix<T> const& LinSpa)
{
  int nbMat=LinSpa.ListMat.size();
  MyMatrix<T> ThePerfMat=LinSpa.SuperMat;
  while(true) {
    Tshortest<T,int> RecSHV=T_ShortestVector<T,int>(ThePerfMat);
    int nbShort=RecSHV.SHV.rows();
    MyMatrix<T> ScalMat(nbShort, nbMat);
    for (int iShort=0; iShort<nbShort; iShort++) {
      MyVector<int> eVectShort=RecSHV.SHV.row(iShort);
      for (int iMat=0; iMat<nbMat; iMat++) {
	T eNorm=EvaluationQuadForm<T,int>(LinSpa.ListMat[iMat], eVectShort);
	ScalMat(iShort, iMat)=eNorm;
      }
    }
    SelectionRowCol<T> eSelect=TMat_SelectRowCol(ScalMat);
    int TheRank=eSelect.TheRank;
    if (TheRank == nbMat)
      break;
    MyVector<T> eVect=eSelect.NSP.row(0);
    MyMatrix<T> DirMat=LINSPA_GetMatrixInTspace(LinSpa, eVect);
    MyMatrix<T> eMatRet=Flipping_Perfect<T,int>(ThePerfMat, DirMat).first;
    ThePerfMat=eMatRet;
  }
  return ThePerfMat;
}








template<typename T, typename Tint, typename Tgroup>
Tgroup PERF_Automorphism(LinSpaceMatrix<T> const& LinSpa,
				 MyMatrix<T> const& ePerf,
				 MyMatrix<Tint> const& SHV)
{
  MyMatrix<T> T_SHV=ConvertMatrixUniversal<T,Tint>(SHV);
  WeightMatrix<std::vector<T>> WMat=GetWeightMatrix_ListComm(T_SHV, ePerf, LinSpa.ListComm);
  return GetStabilizerAsymmetricMatrix<std::vector<T>, Tgroup>(WMat);
}








template<typename T, typename Tint, typename Tgroup>
void VoronoiAlgo_THR_EnumeratePerfectForm(MainProcessor &MProc, int TheId,
					  DataBank<PolyhedralEntry<T,Tgroup>> &TheBank,
					  LinSpaceMatrix<T> const& LinSpa,
					  PolyHeuristic<mpz_class> const& AllArr,
					  ListPerfectForm<T,Tint,Tgroup> &ListPerf)
{
  int IsFinished, NewId;
  int nbConsTodo, iThr, NbThr;
  MyMatrix<T> ePerf=GetOnePerfectForm<T>(LinSpa);
  PerfEquivInfo eEquiv=ListPerf.InsertForm(LinSpa, ePerf);
  while(true) {
    NbThr=MProc.MPU_NumberFree();
    std::vector<std::vector<int>> ThePartition=ListPerf.VoronoiAlgo_THR_GetPartition(NbThr+1, nbConsTodo, IsFinished);
    if (IsFinished == 1)
      break;
    for (iThr=1; iThr<NbThr; iThr++)
      if (ThePartition[iThr].size() > 0) {
	NewId=MProc.MPU_GetId();
	if (NewId != -1) {
	  std::thread my_thread(VoronoiAlgo_THR_BlockTreatment<T>, std::ref(MProc), NewId, std::ref(TheBank), std::ref(LinSpa), std::ref(ListPerf), std::ref(AllArr), std::ref(ThePartition[iThr]));
	  my_thread.detach();
	  MProc.MPU_Terminate(NewId);
	}
      }
    if (ThePartition[0].size() > 0)
      VoronoiAlgo_THR_BlockTreatment<T>(MProc, TheId, TheBank, LinSpa, ListPerf, AllArr, ThePartition[0]);
  }
}

template<typename T, typename Tint, typename Tgroup>
void VoronoiAlgo_PrintListMatrix(std::ostream &os,
				 LinSpaceMatrix<T> const& LinSpa,
				 ListPerfectForm<T,Tint,Tgroup> const& ListPerf)
{
  int nbGram=ListPerf.GetNbPerf();
  os << nbGram << "\n";
  for (int iGram=0; iGram<nbGram; iGram++) {
    os << iGram << " " << nbGram << "\n";
    MyMatrix<T> eGram=ListPerf.GetGram(iGram);
    MyMatrix<int> eSHV=ListPerf.GetSHV(iGram);
    WriteMatrix(os, eGram);
    WriteMatrix(os, eSHV);
  }
}

template<typename T,typename Tint>
struct SimplePerfect {
  MyMatrix<T> Gram;
};

template<typename T,typename Tint>
std::istream& operator>>(std::istream& is, SimplePerfect<T,Tint>& obj)
{
  MyMatrix<T> eG=ReadMatrix<T>(is);
  obj={eG};
  return is;
}


template<typename T,typename Tint>
  std::ostream& operator<<(std::ostream& os, SimplePerfect<T,Tint> const& obj)
{
  WriteMatrix(os, obj.Gram);
  return os;
}


template<typename T>
struct SimplePerfectInv {
  std::vector<T> ListSingleInv;
};

template<typename T>
bool operator==(SimplePerfectInv<T> const& x, SimplePerfectInv<T> const& y)
{
  return x.ListSingleInv == y.ListSingleInv;
}

template<typename T>
std::istream& operator>>(std::istream& is, SimplePerfectInv<T>& obj)
{
  std::vector<T> ListSingleInv;
  is >> ListSingleInv;
  obj={ListSingleInv};
  return is;
}
template<typename T>
std::ostream& operator<<(std::ostream& os, SimplePerfectInv<T> const& obj)
{
  os << obj.ListSingleInv;
  return os;
}

template<typename T>
bool operator<(SimplePerfectInv<T> const& x, SimplePerfectInv<T> const& y)
{
  return x.ListSingleInv < y.ListSingleInv;
}

template <typename T,typename Tint>
  struct invariant_info<SimplePerfect<T,Tint>> {
  typedef SimplePerfectInv<T> invariant_type;
};

template <typename T,typename Tint>
  struct equiv_info<SimplePerfect<T,Tint>> {
  typedef MyMatrix<Tint> equiv_type;
};

template<typename T,typename Tint>
SimplePerfectInv<T> ComputeInvariantSimplePerfect(MyMatrix<T> const& eGram)
{
  int n=eGram.rows();
  Tshortest<T,Tint> RecSHV=T_ShortestVector<T,Tint>(eGram);
  MyMatrix<T> eG = eGram / RecSHV.eMin;
  int nbSHV=RecSHV.SHV.size();
  T eDet=DeterminantMat(eG);
  Tint PreIndex=Int_IndexLattice(RecSHV.SHV);
  Tint eIndex=T_abs(PreIndex);
  WeightMatrix<T> WMat(nbSHV);
  for (int i=0; i<nbSHV-1; i++)
    for (int j=i+1; j<nbSHV; j++) {
      MyVector<Tint> V1(n);
      MyVector<Tint> V2(n);
      for (int iCol=0; iCol<n; iCol++) {
	V1(i)=RecSHV.SHV(i,iCol);
	V2(i)=RecSHV.SHV(j,iCol);
      }
      T eScal=ScalarProductQuadForm<T,Tint>(eGram, V1, V2);
      WMat.Update(i,j,eScal);
    }
  T ePolyInv_T=GetInvariantWeightMatrix(WMat);
  std::vector<T> ListSingleInv{nbSHV, eDet, eIndex, ePolyInv_T};
  return {ListSingleInv};
}

template<typename T>
struct DataLinSpa {
  LinSpaceMatrix<T> LinSpa;
  bool SavingPerfect;
  bool FullDataInMemory;
  std::string PrefixPerfect;
  std::string PrefixPolyhedral;
  bool ReturnAll;
  mpz_class UpperLimitMethod4;
  bool NeedCheckStabilization;
};


template<typename T, typename Tint, typename Tgroup>
EquivTest<MyMatrix<Tint>> SimplePerfect_TestEquivalence(
		      DataLinSpa<T> const&eData,
		      MyMatrix<T> const& Gram1,
		      MyMatrix<T> const& Gram2)
{
  using Telt=typename Tgroup::Telt;
  Tshortest<T,Tint> RecSHV1=T_ShortestVector<T,Tint>(Gram1);
  Tshortest<T,Tint> RecSHV2=T_ShortestVector<T,Tint>(Gram2);
  MyMatrix<T> T_SHV1=ConvertMatrixUniversal<T,Tint>(RecSHV1.SHV);
  MyMatrix<T> T_SHV2=ConvertMatrixUniversal<T,Tint>(RecSHV2.SHV);
  WeightMatrix<std::vector<T>> WMat1=GetWeightMatrix_ListComm(T_SHV1, Gram1, eData.LinSpa.ListComm);
  WeightMatrix<std::vector<T>> WMat2=GetWeightMatrix_ListComm(T_SHV2, Gram2, eData.LinSpa.ListComm);
  EquivTest<Telt> eResEquiv=GetEquivalenceAsymmetricMatrix<std::vector<T>, Telt>(WMat1, WMat2);
  if (!eResEquiv.TheReply) {
    return {false, {}};
  }
  MyMatrix<T> M3=RepresentVertexPermutation(T_SHV1, T_SHV2, eResEquiv.TheEquiv);
  if (IsIntegralMatrix(M3)) {
    MyMatrix<Tint> eMatEquiv=ConvertMatrixUniversal<Tint,T>(M3);
    return {true, eMatEquiv};
  }
  std::vector<MyVector<T>> ListMatVect;
  for (auto & eMat : eData.LinSpa.ListMat) {
    MyVector<T> eVect=SymmetricMatrixToVector(eMat);
    ListMatVect.push_back(eVect);
  }
  MyMatrix<T> ListMatVectB=MatrixFromVectorFamily(ListMatVect);
  std::function<bool(MyMatrix<T>)> IsMatrixCorrect=[&](MyMatrix<T> const& M) -> bool {
    if (!IsIntegralMatrix(M))
      return false;
    if (eData.NeedCheckStabilization) {
      for (auto & eMat : eData.LinSpa.ListMat) {
	MyMatrix<T> eProd=M*eMat*TransposedMat(M);
	MyVector<T> eVect=SymmetricMatrixToVector(eProd);
	SolMatResult<T> eSol=SolutionMat(ListMatVectB, eVect);
	if (!eSol.result)
	  return false;
      }
    }
    return true;
  };
  auto ConvertEquiv=[](EquivTest<MyMatrix<T>> const& eEq) -> EquivTest<MyMatrix<Tint>> {
    if (!eEq.TheReply)
      return {false, {}};
    MyMatrix<Tint> eMat_I=ConvertMatrixUniversal<Tint,T>(eEq.TheEquiv);
    return {true, eMat_I};
  };
  Tgroup GRP1=GetStabilizerAsymmetricMatrix<std::vector<T>, Tgroup>(WMat1);
  if (GRP1.size() < eData.UpperLimitMethod4) {
    return ConvertEquiv(LinPolytopeIntegral_Isomorphism_Method4(T_SHV1, T_SHV2, GRP1, eResEquiv.TheEquiv, IsMatrixCorrect));
  }
  else {
    EquivTest<MyMatrix<T>> fResEquiv=LinPolytopeIntegral_Isomorphism_Method8(T_SHV1, T_SHV2, GRP1, eResEquiv.TheEquiv);
    if (!fResEquiv.TheReply)
      return {false, {}};
    if (IsMatrixCorrect(fResEquiv.TheEquiv))
      return ConvertEquiv(fResEquiv);
    return ConvertEquiv(LinPolytopeIntegral_Isomorphism_Method4(T_SHV1, T_SHV2, GRP1, eResEquiv.TheEquiv, IsMatrixCorrect));
  }
}



template<typename T, typename Tint, typename Tgroup>
Tgroup SimplePerfect_Stabilizer(DataLinSpa<T> const& eData, MyMatrix<T> const& Gram, Tshortest<T,Tint> const& RecSHV)
{
  using Telt=typename Tgroup::Telt;
  //
  // Functionality for checking quality of equivalences
  //
  std::vector<MyVector<T>> ListMatVect;
  for (auto & eMat : eData.LinSpa.ListMat) {
    MyVector<T> eVect=SymmetricMatrixToVector(eMat);
    ListMatVect.push_back(eVect);
  }
  MyMatrix<T> ListMatVectB=MatrixFromVectorFamily(ListMatVect);
  MyMatrix<T> T_SHV=ConvertMatrixUniversal<T,Tint>(RecSHV.SHV);
  std::function<bool(MyMatrix<T>)> IsMatrixCorrect=[&](MyMatrix<T> const& M) -> bool {
    if (!IsIntegralMatrix(M))
      return false;
    if (eData.NeedCheckStabilization) {
      for (auto & eMat : eData.LinSpa.ListMat) {
	MyMatrix<T> eProd=M*eMat*TransposedMat(M);
	MyVector<T> eVect=SymmetricMatrixToVector(eProd);
	SolMatResult<T> eSol=SolutionMat(ListMatVectB, eVect);
	if (!eSol.result)
	  return false;
      }
    }
    return true;
  };
  auto IsCorrectGroup=[&](Tgroup const& g) -> bool {
    std::vector<Telt> LGen = g.GeneratorsOfGroup();
    for (auto & eGen : LGen) {
      MyMatrix<T> M=RepresentVertexPermutation(T_SHV, T_SHV, eGen);
      if (!IsMatrixCorrect(M))
	return false;
    }
    return true;
  };
  //
  // Now the computation itself
  //
  WeightMatrix<std::vector<T>> WMat=GetWeightMatrix_ListComm(T_SHV, Gram, eData.LinSpa.ListComm);
  Tgroup GRPshv1=GetStabilizerAsymmetricMatrix<std::vector<T>, Tgroup>(WMat);
  if (IsCorrectGroup(GRPshv1))
    return GRPshv1;
  if (GRPshv1.size() < eData.UpperLimitMethod4) {
    return LinPolytopeIntegral_Stabilizer_Method4(T_SHV, GRPshv1, IsMatrixCorrect);
  }
  else {
    Tgroup GRPshv2=LinPolytopeIntegral_Stabilizer_Method8(T_SHV, GRPshv1);
    if (IsCorrectGroup(GRPshv2))
      return GRPshv2;
    return LinPolytopeIntegral_Stabilizer_Method4(T_SHV, GRPshv2, IsMatrixCorrect);
  }
}






template<typename T, typename Tint, typename Tgroup>
  std::vector<SimplePerfect<T,Tint>> EnumerationPerfectMatrices(
	   MainProcessor &MProc, int const& TheId,
	   DataBank<PolyhedralEntry<T,Tgroup>> &TheBank,
	   DataLinSpa<T> const& eData,
	   PolyHeuristic<mpz_class> const& AllArr)
{
  std::function<bool(SimplePerfectInv<T> const&,SimplePerfectInv<T> const&)> CompFCT=[](SimplePerfectInv<T> const& x, SimplePerfectInv<T> const& y) -> bool {
    return x < y;
  };
  std::function<void(TrivialBalinski &,SimplePerfect<T,Tint> const&,SimplePerfectInv<T> const&,std::ostream&)> UpgradeBalinskiStat=[](TrivialBalinski const& eStat, SimplePerfect<T,Tint> const& eEnt, SimplePerfectInv<T> const& eInv, std::ostream&os) -> void {
  };
  std::function<EquivTest<MyMatrix<Tint>>(SimplePerfect<T,Tint> const&, SimplePerfect<T,Tint> const&)> fEquiv=[&](SimplePerfect<T,Tint> const& x, SimplePerfect<T,Tint> const& y) -> EquivTest<MyMatrix<Tint>> {
    return SimplePerfect_TestEquivalence<T,Tint,Tgroup>(eData, x.Gram, y.Gram);
  };
  NewEnumerationWork<SimplePerfect<T,Tint>> ListOrbit(AllArr.Saving, AllArr.eMemory, eData.PrefixPerfect, CompFCT, UpgradeBalinskiStat, fEquiv, MProc.GetO(TheId));
  auto FuncInsert=[&](SimplePerfect<T,Tint> const& x, std::ostream&os) -> int {
    SimplePerfectInv<T> eInv=ComputeInvariantSimplePerfect<T,Tint>(x.Gram);
    return ListOrbit.InsertEntry({x, eInv}, os);
  };
  int nbPresentOrbit=ListOrbit.GetNbEntry();
  if (nbPresentOrbit == 0) {
    MyMatrix<T> eGram=GetOnePerfectForm(eData.LinSpa);
    int RetVal=FuncInsert({eGram}, MProc.GetO(TheId));
    if (RetVal != -1) {
      std::cerr << "RetVal=" << RetVal << "\n";
      std::cerr << "The first orbit should definitely be new\n";
      std::cerr << "Otherwise, we have a clear bug\n";
      throw TerminalException{1};
    }
  }
  MProc.GetO(TheId) << "We have inserted eInc\n";
  int nbSpannThread=0;
  std::condition_variable cv;
  std::mutex mtx_cv;
  auto WaitStuck=[&](int const& MyId) -> void {
    std::unique_lock<std::mutex> lk(mtx_cv);
    cv.wait(lk, [&]{return ListOrbit.IsStuck() == false;});
  };
  auto WaitComplete=[&](int const& MyId) -> void {
    std::unique_lock<std::mutex> lk(mtx_cv);
    cv.wait(lk, [&]{return nbSpannThread == 0;});
  };
  auto CompSimplePerfectAdjacency=[&](int const& MyId, std::ostream& os) -> void {
    int nbWork=0;
    nbSpannThread++;
    while(true) {
      bool testStuck=ListOrbit.IsStuck();
      if (testStuck) {
	WaitStuck(MyId);
      }
      int eEntry=ListOrbit.GetNonTreatedOrbit(os);
      bool IsComplete=ListOrbit.GetCompleteStatus();
      if (eEntry == -1)
	break;
      if (IsComplete)
	break;
      SimplePerfect<T,Tint> ePERF=ListOrbit.GetRepresentative(eEntry);
      SimplePerfectInv<T> eInv=ListOrbit.GetInvariant(eEntry);
      Tshortest<T,Tint> RecSHV=T_ShortestVector<T,Tint>(ePERF.Gram);
      NakedPerfect<T,Tint> eNaked=GetNakedPerfectCone(eData.LinSpa, ePERF.Gram, RecSHV);
      Tgroup GRPshv=SimplePerfect_Stabilizer<T,Tint,Tgroup>(eData, ePERF.Gram, RecSHV);
      Tgroup PerfDomGRP=MapLatticeGroupToConeGroup(eNaked, GRPshv);
      CondTempDirectory eDir(AllArr.Saving, eData.PrefixPolyhedral + "ADM" + IntToString(eEntry) + "/");
      vectface TheOutput=DUALDESC_THR_AdjacencyDecomposition(MProc, MyId, TheBank, eNaked.PerfDomEXT, PerfDomGRP, AllArr, eDir.str(), 0);
      for (auto& eOrbB : TheOutput) {
	MyVector<T> eVectOrb=FindFacetInequality(eNaked.PerfDomEXT, eOrbB);
	MyMatrix<T> DirMat=LINSPA_GetMatrixInTspace(eData.LinSpa, eVectOrb);
	MyMatrix<T> NewPerf=Flipping_Perfect<T,int>(ePERF.Gram, DirMat).first;
	int eVal=FuncInsert({NewPerf}, os);
	if (eVal == -1)
	  cv.notify_one();
      }
      ListOrbit.SetEntryAsDone(eEntry, os);
      nbWork++;
    }
    if (MyId != TheId)
      MProc.MPU_Terminate(MyId);
    nbSpannThread--;
    if (nbWork > 0) {
      cv.notify_all();
    }
  };
  int NbThr=MProc.MPU_NumberFree();
  std::vector<std::thread> ListThreads(NbThr);
  MProc.GetO(TheId) << "We have NbThr=" << NbThr << "\n";
  while(true) {
    bool IsCompleteSpann=ListOrbit.GetCompleteStatus();
    MProc.GetO(TheId) << "IsCompleteSpann=" << IsCompleteSpann << "\n";
    if (IsCompleteSpann)
      WaitComplete(TheId);
    else {
      for (int iThr=0; iThr<NbThr; iThr++) {
	int NewId=MProc.MPU_GetId();
	MProc.GetO(TheId) << "NewId=" << NewId << "\n";
	if (NewId != -1) {
	  MProc.GetO(NewId) << "Before call to my_thread\n";
	  ListThreads[iThr]=std::thread(CompSimplePerfectAdjacency, NewId, std::ref(MProc.GetO(NewId)));
	  MProc.GetO(NewId) << "After call to my_thread\n";
	  ListThreads[iThr].detach();
	}
      }
      CompSimplePerfectAdjacency(TheId, MProc.GetO(TheId));
    }
    if (nbSpannThread == 0)
      break;
  }
  bool IsComplete=ListOrbit.GetCompleteStatus();
  MProc.GetO(TheId) << "IsComplete=" << IsComplete << "\n";
  if (!IsComplete) {
    std::cerr << "Major error in the code. We should be complete now\n";
    throw TerminalException{1};
  }
  std::vector<SimplePerfect<T,Tint>> ReturnData;
  if (eData.ReturnAll) {
    ReturnData=ListOrbit.GetListRepr();
    ListOrbit.FuncClear();
  }
  return ReturnData;
}




FullNamelist NAMELIST_GetStandard_COMPUTE_PERFECT()
{
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListStringValues1["LinSpaFile"]="unset.gram";
  ListStringValues1["OUTfile"]="unset.out";
  ListBoolValues1["SavingPerfect"]=false;
  ListBoolValues1["NeedCheckStabilization"]=false;
  ListBoolValues1["FullDataInMemory"]=true;
  ListBoolValues1["ReturnAll"]=false;
  ListStringValues1["PrefixPerfect"]="/irrelevant/";
  SingleBlock BlockDATA;
  BlockDATA.ListIntValues=ListIntValues1;
  BlockDATA.ListBoolValues=ListBoolValues1;
  BlockDATA.ListDoubleValues=ListDoubleValues1;
  BlockDATA.ListStringValues=ListStringValues1;
  BlockDATA.ListListStringValues=ListListStringValues1;
  ListBlock["DATA"]=BlockDATA;
  // METHOD
  std::map<std::string, int> ListIntValues2;
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, double> ListDoubleValues2;
  std::map<std::string, std::string> ListStringValues2;
  std::map<std::string, std::vector<std::string>> ListListStringValues2;
  ListIntValues2["NPROC"]=1;
  ListIntValues2["UpperLimitMethod4"]=10000;
  ListBoolValues2["Saving"]=false;
  ListStringValues2["PrefixPolyhedral"]="/irrelevant/";
  ListBoolValues2["FullDataInMemory"]=true;
  ListStringValues2["SplittingHeuristicFile"]="unset.heu";
  ListStringValues2["AdditionalSymmetryHeuristicFile"]="unset.heu";
  ListStringValues2["DualDescriptionHeuristicFile"]="unset.heu";
  ListStringValues2["StabEquivFacetHeuristicFile"]="unset.heu";
  ListStringValues2["MethodInitialFacetSetFile"]="unset.heu";
  SingleBlock BlockMETHOD;
  BlockMETHOD.ListIntValues=ListIntValues2;
  BlockMETHOD.ListBoolValues=ListBoolValues2;
  BlockMETHOD.ListDoubleValues=ListDoubleValues2;
  BlockMETHOD.ListStringValues=ListStringValues2;
  BlockMETHOD.ListListStringValues=ListListStringValues2;
  ListBlock["METHOD"]=BlockMETHOD;
  // BANK
  std::map<std::string, int> ListIntValues3;
  std::map<std::string, bool> ListBoolValues3;
  std::map<std::string, double> ListDoubleValues3;
  std::map<std::string, std::string> ListStringValues3;
  std::map<std::string, std::vector<std::string>> ListListStringValues3;
  ListStringValues3["BankSaveHeuristicFile"]="unset.heu";
  ListStringValues3["Prefix"]="./unset/";
  ListBoolValues3["Saving"]=false;
  ListBoolValues3["FullDataInMemory"]=true;
  SingleBlock BlockBANK;
  BlockBANK.ListIntValues=ListIntValues3;
  BlockBANK.ListBoolValues=ListBoolValues3;
  BlockBANK.ListDoubleValues=ListDoubleValues3;
  BlockBANK.ListStringValues=ListStringValues3;
  BlockBANK.ListListStringValues=ListListStringValues3;
  ListBlock["BANK"]=BlockBANK;
  // Merging all data
  return {ListBlock, "undefined"};
}







template<typename T, typename Tint,typename Tgroup>
void TreatPerfectLatticesEntry(FullNamelist const& eFull)
{
  SingleBlock BlockBANK=eFull.ListBlock.at("BANK");
  SingleBlock BlockDATA=eFull.ListBlock.at("DATA");
  SingleBlock BlockMETHOD=eFull.ListBlock.at("METHOD");
  //
  bool BANK_IsSaving=BlockBANK.ListBoolValues.at("Saving");
  bool BANK_Memory=BlockBANK.ListBoolValues.at("FullDataInMemory");
  std::string BANK_Prefix=BlockBANK.ListStringValues.at("Prefix");
  CreateDirectory(BANK_Prefix);
  FctsDataBank<PolyhedralEntry<T,Tgroup>> recFct=GetRec_FctsDataBank<T,Tgroup>();
  DataBank<PolyhedralEntry<T,Tgroup>> TheBank(BANK_IsSaving, BANK_Memory, BANK_Prefix, recFct);
  //
  std::cerr << "Reading DATA\n";
  std::string LINSPAfile=BlockDATA.ListStringValues.at("LinSpaFile");
  IsExistingFileDie(LINSPAfile);
  std::cerr << "LINSPAfile=" << LINSPAfile << "\n";
  std::string OUTfile=BlockDATA.ListStringValues.at("OUTfile");
  std::cerr << "OUTfile=" << OUTfile << "\n";
  //
  std::ifstream LINSPAfs(LINSPAfile);
  LinSpaceMatrix<T> LinSpa;
  LINSPAfs >> LinSpa;
  bool NeedCheckStabilization=BlockDATA.ListBoolValues.at("NeedCheckStabilization");
  bool SavingPerfect=BlockDATA.ListBoolValues.at("SavingPerfect");
  bool FullDataInMemory=BlockDATA.ListBoolValues.at("FullDataInMemory");
  bool ReturnAll=BlockDATA.ListBoolValues.at("ReturnAll");
  std::string PrefixPerfect=BlockDATA.ListStringValues.at("PrefixPerfect");
  CreateDirectory(PrefixPerfect);
  std::string PrefixPolyhedral=BlockMETHOD.ListStringValues.at("PrefixPolyhedral");
  CreateDirectory(PrefixPolyhedral);
  int UpperLimitMethod4=BlockMETHOD.ListIntValues.at("UpperLimitMethod4");
  std::string CVPmethod=BlockMETHOD.ListStringValues.at("CVPmethod");
  mpz_class UpperLimitMethod4_z=UpperLimitMethod4;
  DataLinSpa<T> eData;
  eData.LinSpa=LinSpa;
  eData.SavingPerfect=SavingPerfect;
  eData.FullDataInMemory=FullDataInMemory;
  eData.PrefixPerfect=PrefixPerfect;
  eData.PrefixPolyhedral=PrefixPolyhedral;
  eData.ReturnAll=ReturnAll;
  eData.UpperLimitMethod4=UpperLimitMethod4_z;
  eData.NeedCheckStabilization=NeedCheckStabilization;
  //
  std::cerr << "Creating MPROC\n";
  int NbThr=BlockMETHOD.ListIntValues.at("NPROC");
  MainProcessor MProc(NbThr);
  int TheId=MProc.MPU_GetId();
  //
  PolyHeuristic<mpz_class> AllArr=AllStandardHeuristic<mpz_class>();
  //
  std::string HeuSplitFile=BlockMETHOD.ListStringValues.at("SplittingHeuristicFile");
  if (HeuSplitFile != "unset.heu") {
    IsExistingFileDie(HeuSplitFile);
    std::ifstream SPLITfs(HeuSplitFile);
    AllArr.Splitting=ReadHeuristic<mpz_class>(SPLITfs);
  }
  //
  std::string AddiSymmFile=BlockMETHOD.ListStringValues.at("AdditionalSymmetryHeuristicFile");
  if (AddiSymmFile != "unset.heu") {
    IsExistingFileDie(AddiSymmFile);
    std::ifstream SYMMfs(AddiSymmFile);
    AllArr.AdditionalSymmetry=ReadHeuristic<mpz_class>(SYMMfs);
  }
  //
  std::string DualDescFile=BlockMETHOD.ListStringValues.at("DualDescriptionHeuristicFile");
  if (DualDescFile != "unset.heu") {
    IsExistingFileDie(DualDescFile);
    std::ifstream DDfs(DualDescFile);
    AllArr.DualDescriptionProgram=ReadHeuristic<mpz_class>(DDfs);
  }
  //
  std::string StabEquivFile=BlockMETHOD.ListStringValues.at("StabEquivFacetHeuristicFile");
  if (StabEquivFile != "unset.heu") {
    IsExistingFileDie(StabEquivFile);
    std::ifstream STEQfs(StabEquivFile);
    AllArr.StabEquivFacet=ReadHeuristic<mpz_class>(STEQfs);
  }
  //
  std::string MethodFacetFile=BlockMETHOD.ListStringValues.at("MethodInitialFacetSetFile");
  if (MethodFacetFile != "unset.heu") {
    IsExistingFileDie(MethodFacetFile);
    std::ifstream MIFSfs(MethodFacetFile);
    AllArr.InitialFacetSet=ReadHeuristic<mpz_class>(MIFSfs);
  }
  //
  std::string BankSaveFile=BlockBANK.ListStringValues.at("BankSaveHeuristicFile");
  if (BankSaveFile != "unset.heu") {
    IsExistingFileDie(BankSaveFile);
    std::ifstream BANKfs(BankSaveFile);
    AllArr.BankSave=ReadHeuristic<mpz_class>(BANKfs);
  }
  //
  bool DD_Saving=BlockMETHOD.ListBoolValues.at("Saving");
  bool DD_Memory=BlockMETHOD.ListBoolValues.at("FullDataInMemory");
  AllArr.Saving=DD_Saving;
  AllArr.eMemory=DD_Memory;
  //
  std::vector<SimplePerfect<T,Tint>> LPerf=EnumerationPerfectMatrices<T,Tint>(MProc, TheId,
									      TheBank, eData, AllArr);
  std::cerr << "We now have ListDel\n";
  //
  if (ReturnAll) {
    std::ofstream OUTfs(OUTfile);
    int nbPerf=LPerf.size();
    OUTfs << "nbPerf=" << nbPerf << "\n";
    for (int iPerf=0; iPerf<nbPerf; iPerf++) {
      OUTfs << "iPerf=" << iPerf << "/" << nbPerf << "\n";
      WriteMatrix(OUTfs, LPerf[iPerf].Gram);
    }
  }
}




#endif
