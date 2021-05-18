#ifndef INCLUDE_PERM_FCT_H
#define INCLUDE_PERM_FCT_H



template<typename T, typename Tidx>
std::vector<Tidx> SortingPerm(std::vector<T> const & ListV)
{
  struct PairData {
    Tidx i;
    T x;
  };
  Tidx len=ListV.size();
  std::vector<PairData> ListPair(len);
  for (Tidx i=0; i<len; i++) {
    PairData ePair{i, ListV[i]};
    ListPair[i]=ePair;
  }
  sort(ListPair.begin(), ListPair.end(),
       [](PairData const & x1, PairData const& x2) -> bool {
	 if (x1.x < x2.x)
	   return true;
	 if (x2.x < x1.x)
	   return false;
	 return x1.i< x2.i;
       });
  std::vector<Tidx> v(len);
  for (Tidx i=0; i<len; i++)
    v[i]=ListPair[i].i;
  return v;
}




template<typename T, typename Telt>
MyMatrix<T> RepresentVertexPermutation(MyMatrix<T> const& EXT1, MyMatrix<T> const& EXT2, Telt const& ePerm)
{
  SelectionRowCol<T> eSelect=TMat_SelectRowCol(EXT1);
  std::vector<int> ListRowSelect=eSelect.ListRowSelect;
  MyMatrix<T> M1=SelectRow(EXT1, ListRowSelect);
  MyMatrix<T> M1inv=Inverse(M1);
  size_t nbRow=ListRowSelect.size();
  std::vector<int> ListRowSelectImg(nbRow);
  for (size_t iRow=0; iRow<nbRow; iRow++)
    ListRowSelectImg[iRow]=ePerm.at(iRow);
  MyMatrix<T> M2=SelectRow(EXT2, ListRowSelectImg);
  return M1inv*M2;
}



template<typename T, typename Telt>
Telt GetPermutationOnVectors(MyMatrix<T> const& EXT1, MyMatrix<T> const& EXT2)
{
  using Tidx = typename Telt::Tidx;
  size_t nbVect=EXT1.rows();
  std::vector<MyVector<T>> EXTrow1(nbVect), EXTrow2(nbVect);
  for (size_t iVect=0; iVect<nbVect; iVect++) {
    EXTrow1[iVect]=GetMatrixRow(EXT1, iVect);
    EXTrow2[iVect]=GetMatrixRow(EXT2, iVect);
  }
  Telt ePerm1=Telt(SortingPerm<MyVector<T>,Tidx>(EXTrow1));
  Telt ePerm2=Telt(SortingPerm<MyVector<T>,Tidx>(EXTrow2));
  Telt ePermRet=(~ePerm1) * ePerm2;
#ifdef DEBUG
  for (size_t iVect=0; iVect<nbVect; iVect++) {
    size_t jVect=ePermRet.at(iVect);
    if (EXTrow2[jVect] != EXTrow1[iVect]) {
      std::cerr << "iVect=" << iVect << " jVect=" << jVect << "\n";
      std::cerr << "Id:";
      for (size_t k=0; k<nbVect; k++)
	std::cerr << " " << k;
      std::cerr << "\n";
      std::cerr << " p:";
      for (size_t k=0; k<nbVect; k++)
	std::cerr << " " << ePermRet.at(k);
      std::cerr << "\n";

      std::cerr << "perm1:";
      for (size_t k=0; k<nbVect; k++)
	std::cerr << " " << ePerm1.at(k);
      std::cerr << "\n";
      std::cerr << "perm2:";
      for (size_t k=0; k<nbVect; k++)
	std::cerr << " " << ePerm2.at(k);
      std::cerr << "\n";


      std::cerr << "EXTrow1[iVect]=";
      WriteVector(std::cerr, EXTrow1[iVect]);
      std::cerr << "EXTrow2[jVect]=";
      WriteVector(std::cerr, EXTrow2[jVect]);
      std::cerr << "EXTrow1=\n";
      WriteMatrix(std::cerr, EXT1);
      std::cerr << "EXTrow2=\n";
      WriteMatrix(std::cerr, EXT2);
      std::cerr << "Error in GetPermutationOnVectors\n";
      throw TerminalException{1};
    }
  }
#endif
  return ePermRet;
}









#endif
