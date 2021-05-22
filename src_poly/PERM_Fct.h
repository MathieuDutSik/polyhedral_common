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
  std::vector<int> const& ListRowSelect=eSelect.ListRowSelect;
  MyMatrix<T> M1=SelectRow(EXT1, ListRowSelect);
  MyMatrix<T> M1inv=Inverse(M1);
  size_t nbRow=ListRowSelect.size();
  std::vector<int> ListRowSelectImg(nbRow);
  for (size_t iRow=0; iRow<nbRow; iRow++)
    ListRowSelectImg[iRow]=ePerm.at(iRow);
  MyMatrix<T> M2=SelectRow(EXT2, ListRowSelectImg);
  return M1inv*M2;
}


template<typename T, typename Tfield, typename Telt>
EquivTest<MyMatrix<Tfield>> RepresentVertexPermutationTest(MyMatrix<T> const& EXT1, MyMatrix<T> const& EXT2, Telt const& ePerm)
{
  static_assert(is_ring_field<Tfield>::value, "Requires Tfield to be a field in DivideVector");
  size_t nbRow=EXT1.rows();
  size_t nbCol=EXT1.cols();
  SelectionRowCol<T> eSelect=TMat_SelectRowCol(EXT1); // needs another version that do not use the type conversion
  std::vector<int> const& ListRowSelect=eSelect.ListRowSelect;
  MyMatrix<T> M1=SelectRow(EXT1, ListRowSelect);
  MyMatrix<Tfield> M1_field=ConvertMatrixUniversal<Tfield,T>(M1);
  MyMatrix<Tfield> M1inv_field=Inverse(M1_field);
  size_t nbRow_select=ListRowSelect.size();
  if (nbRow_select != nbCol) {
    std::cerr << "nbRow_select should be equal to nbCol\n";
    throw TerminalException{1};
  }
  std::vector<int> ListRowSelectImg(nbRow);
  for (size_t iRow=0; iRow<nbRow_select; iRow++)
    ListRowSelectImg[iRow]=ePerm.at(iRow);
  MyMatrix<T> M2=SelectRow(EXT2, ListRowSelectImg);
  MyMatrix<T> M2_field=ConvertMatrixUniversal<Tfield,T>(M2);
  MyMatrix<Tfield> EqMat = M1inv_field * M2_field;
  // Now testing that we have EXT1 EqMat = EXT2
  for (size_t iRow=0; iRow<nbRow; iRow++) {
    size_t iRowImg = ePerm.at(iRow);
    for (size_t iCol=0; iCol<nbCol; iCol++) {
      Tfield eSum = -EXT2(iRowImg, iCol);
      for (size_t jRow=0; jRow<nbCol; jRow++) {
        eSum += EqMat(jRow,iCol) * EXT1(iRow, jRow);
      }
      if (eSum != 0) {
        return {false, {}};
      }
    }
  }
  return {true, EqMat};
}


template<typename T, typename Tfield, typename Telt>
EquivTest<Telt> RepresentVertexPermutationTest(MyMatrix<T> const& EXT1, MyMatrix<T> const& EXT2, MyMatrix<Tfield> const& P)
{
  using Tidx = typename Telt::Tidx;
  size_t n_rows = EXT1.rows();
  size_t n_cols = EXT1.cols();
  MyMatrix<T> VectorContain(1,n_cols);
  ContainerMatrix<T> Cont(EXT2);
  Cont.SetPtr(&VectorContain);
  //
  // We are testing if EXT1 P = perm(EXT2) 
  std::vector<Tidx> V(n_rows);
  for (size_t i_row=0; i_row<n_rows; i_row++) {
    for (size_t i_col=0; i_col<n_cols; i_col++) {
      Tfield eSum1 = 0;
      for (size_t j_row=0; j_row<n_cols; j_row++)
        eSum1 += EXT1(i_row, j_row) * P(j_row, i_col);
      T eSum2 = UniversalTypeConversion<T,Tfield>(eSum1);
      Tfield eSum3 = UniversalTypeConversion<Tfield,T>(eSum2);
      if (eSum1 != eSum3) {
        return {false,{}}; // We fail because the image is not integral.
      }
      VectorContain(0, i_col) = eSum2;
    }
    std::pair<bool, size_t> epair = Cont.GetIdx();
    if (!epair.first) {
      return {false,{}}; // We fail because the image does not belong to EXT2
    }
    V[i_row] = epair.second;
  }
  return {true, std::move(Telt(V))};
}







template<typename T, typename Telt>
EquivTest<MyMatrix<T>> RepresentVertexPermutation(MyMatrix<T> const& EXT1, MyMatrix<T> const& EXT2, Telt const& ePerm)
{
  SelectionRowCol<T> eSelect=TMat_SelectRowCol(EXT1);
  std::vector<int> const& ListRowSelect=eSelect.ListRowSelect;
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






template<typename Tgroup>
Tgroup GetGroupListGen(std::vector<std::vector<unsigned int>> const& ListGen, size_t const& nbVert)
{
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  std::vector<Telt> generatorList;
  std::vector<Tidx> gList(nbVert);
  for (auto & eGen : ListGen) {
    for (size_t iVert=0; iVert<nbVert; iVert++)
      gList[iVert]=eGen[iVert];
    generatorList.push_back(Telt(gList));
  }
  return Tgroup(generatorList, nbVert);
}


template<typename Tgr>
bool CheckListGenerators(std::vector<std::vector<unsigned int>> const& ListGen, Tgr const& eGR)
{
  size_t nbVert = eGR.GetNbVert();
  for (auto & eGen : ListGen) {
    for (size_t iVert=0; iVert<nbVert; iVert++) {
      int eColor = eGR.GetColor(iVert);
      int iVert_img = eGen[iVert];
      int fColor = eGR.GetColor(iVert_img);
      if (eColor != fColor) {
        return false;
      }
      //
      for (auto & jVert : eGR.Adjacency(iVert)) {
        int jVert_img = eGen[jVert];
        bool test = eGR.IsAdjacent(iVert_img, jVert_img);
        if (!test)
          return false;
      }
    }
  }
  return true;
}


template<typename Tgr, typename Tgroup>
void PrintStabilizerGroupSizes(std::ostream& os, Tgr const& eGR)
{
  bliss::Graph g=GetBlissGraphFromGraph(eGR);
  size_t nbVert=eGR.GetNbVert();
  std::vector<std::vector<unsigned int>> ListGen1 = BLISS_GetListGenerators(eGR);
  std::vector<std::vector<unsigned int>> ListGen2 = TRACES_GetListGenerators(eGR);
  auto siz1 = GetGroupListGen<Tgroup>(ListGen1, nbVert).size();
  auto siz2 = GetGroupListGen<Tgroup>(ListGen2, nbVert).size();
  bool test1 = CheckListGenerators(ListGen1, eGR);
  bool test2 = CheckListGenerators(ListGen2, eGR);
  os << "|GRP bliss|=" << siz1 << " |GRP traces|=" << siz2 << " test1=" << test1 << " test2=" << test2 << "\n";
}





#endif
