#include "MAT_Matrix.h"
#include "Namelist.h"
#include "NumberTheory.h"

FullNamelist NAMELIST_GetStandardAnalysis()
{
  std::map<std::string, SingleBlock> ListBlock;
  // PROC
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string> > ListListStringValues1;
  std::map<std::string, std::vector<int> > ListListIntValues1;
  std::map<std::string, std::vector<double> > ListListDoubleValues1;
  ListStringValues1["InputFile"]="unset";
  ListDoubleValues1["thrEigenvalueZero"]=0.0001;
  ListIntValues1["MaximumValueIntSearch"]=500;
  ListBoolValues1["SearchKernel"]=false;
  ListBoolValues1["SearchOrthKernel"]=false;
  ListBoolValues1["ShowOrthKernel"]=false;
  ListBoolValues1["ShowKernel"]=false;
  ListBoolValues1["CanonicalizeByMinValue"]=false;
  ListBoolValues1["CheckPrimalDualCancellation"]=false;
  ListStringValues1["PrimalDualPairFile"]="unset";
  SingleBlock BlockPROC;
  BlockPROC.ListIntValues=ListIntValues1;
  BlockPROC.ListBoolValues=ListBoolValues1;
  BlockPROC.ListDoubleValues=ListDoubleValues1;
  BlockPROC.ListStringValues=ListStringValues1;
  BlockPROC.ListListStringValues=ListListStringValues1;
  BlockPROC.ListListIntValues=ListListIntValues1;
  BlockPROC.ListListDoubleValues=ListListDoubleValues1;
  ListBlock["PROC"]=BlockPROC;
  // Final part
  return {ListBlock, "undefined"};
}





int main(int argc, char *argv[])
{
  std::cerr << std::setprecision(20);
  try {
    FullNamelist eFull=NAMELIST_GetStandardAnalysis();
    if (argc != 2) {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "AnalysisSingMat [Namelist]\n");
      NAMELIST_WriteNamelistFile(std::cerr, eFull);
      return -1;
    }
    std::string eFileName=argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    SingleBlock BlPROC = eFull.ListBlock.at("PROC");
    std::string InputFile = BlPROC.ListStringValues.at("InputFile");
    double thrEigenvalueZero = BlPROC.ListDoubleValues.at("thrEigenvalueZero");
    double MaximumValueIntSearch = BlPROC.ListIntValues.at("MaximumValueIntSearch");
    bool SearchKernel = BlPROC.ListBoolValues.at("SearchKernel");
    bool SearchOrthKernel = BlPROC.ListBoolValues.at("SearchOrthKernel");
    bool ShowOrthKernel = BlPROC.ListBoolValues.at("ShowOrthKernel");
    bool ShowKernel = BlPROC.ListBoolValues.at("ShowKernel");
    bool CanonicalizeByMinValue = BlPROC.ListBoolValues.at("CanonicalizeByMinValue");
    bool CheckPrimalDualCancellation = BlPROC.ListBoolValues.at("CheckPrimalDualCancellation");
    std::string PrimalDualPairFile = BlPROC.ListStringValues.at("PrimalDualPairFile");
    // reading the matrix
    if (!IsExistingFile(InputFile)) {
      std::cerr << "The InputFile=" << InputFile << " is missing\n";
      throw TerminalException{1};
    }
    std::ifstream INmat(InputFile);
    MyMatrix<double> eMat=ReadMatrix<double>(INmat);
    int nbRow=eMat.rows();
    if (CheckPrimalDualCancellation) {
      /*
      if (!IsExistingFile(PrimalDualPairFile)) {
	std::cerr << "The PrimalDualPairFile=" << PrimalDualPairFile << " is missing\n";
	throw TerminalException{1};
	}*/
      std::ifstream PAIRmat(PrimalDualPairFile);
      std::cerr << "PrimalDualPairFile=" << PrimalDualPairFile << "\n";
      MyMatrix<double> PairMat=ReadMatrix<double>(PAIRmat);
      MyMatrix<double> eProd = eMat * PairMat;
      double eScal = eProd.trace();
      std::cerr << "PairDual signature = " << eScal << "\n";
    }

    
    std::cerr << "Read eMat rows/cols = " << eMat.rows() << " / " << eMat.cols() << "\n";
    //
    double sumCoeff=0;
    for (int iRow=0; iRow<nbRow; iRow++) {
      for (int iCol=0; iCol<nbRow; iCol++) {
	sumCoeff += fabs(eMat(iRow, iCol));
      }
    }
    std::cerr << "sumCoeff = " << sumCoeff << "\n";
    // creating the RecSparse operators.
    Eigen::SelfAdjointEigenSolver<MyMatrix<double>> eig(eMat);
    MyVector<double> ListEig=eig.eigenvalues();
    MyMatrix<double> ListVect=eig.eigenvectors();
    std::cerr << "ListEigenvalues =";
    for (int iRow=0; iRow<nbRow; iRow++)
      std::cerr << " " << ListEig(iRow);
    std::cerr << "\n";
    int DimKernel=0;
    std::vector<MyVector<double>> ListKernel;
    std::vector<MyVector<double>> ListOrthKernel;
    WriteVector(std::cerr, ListEig);
    auto TheCan=[&](MyVector<double> const& V) -> MyVector<double> {
      if (!CanonicalizeByMinValue)
	return V;
      double eVal = V.cwiseAbs().minCoeff();
      std::cerr << "eVal=" << eVal << "\n";
      return V / eVal;
    };
    for (int iRow=0; iRow<nbRow; iRow++) {
      double eEig=ListEig(iRow);
      MyVector<double> eCol=GetMatrixColumn(ListVect, iRow);
      MyVector<double> TheDiff = eEig * eCol - eMat * eCol;
      double err=TheDiff.cwiseAbs().maxCoeff();
      //      std::cerr << "iRow=" << iRow << " eEig=" << eEig << " err=" << err << "\n";
      if (fabs(eEig) < thrEigenvalueZero) {
	DimKernel++;
	ListKernel.push_back(eCol);
	if (ShowKernel)
	  WriteVector(std::cerr, TheCan(eCol));
      }
      else {
	ListOrthKernel.push_back(eCol);
	if (ShowOrthKernel) {
	  WriteVector(std::cerr, TheCan(eCol));
	}
      }
    }
    using Tint = mpz_class;
    auto Reduction=[&](std::vector<MyVector<double>> const& inpListVect) -> std::pair<double,MyMatrix<Tint>> {
      MyMatrix<double> inpMat = MatrixFromVectorFamily(inpListVect);
      MyMatrix<double> RedMat = CanonicalizeBasisVectorSpace(inpMat);
      std::vector<MyVector<Tint>> LVect;
      int nbRow=RedMat.rows();
      double SumErr=0;
      for (int iRow=0; iRow<nbRow; iRow++) {
	MyVector<double> eRow = GetMatrixRow(RedMat, iRow);
	std::pair<double,MyVector<Tint>> PairRed = FindBestIntegerApproximation<Tint,double>(eRow, MaximumValueIntSearch);
	LVect.push_back(PairRed.second);
	SumErr += PairRed.first;
      }
      return {SumErr,MatrixFromVectorFamily(LVect)};
    };
    if (SearchKernel) {
      std::pair<double,MyMatrix<Tint>> PairRed = Reduction(ListKernel);
      std::cerr << "Approximate kernel found with error=" << PairRed.first << "\n";
      WriteMatrix(std::cerr, PairRed.second);
    }
    if (SearchOrthKernel) {
      std::pair<double,MyMatrix<Tint>> PairRed = Reduction(ListOrthKernel);
      std::cerr << "Approximate orthogonal of kernel found with error=" << PairRed.first << "\n";
      WriteMatrix(std::cerr, PairRed.second);
    }
    std::cerr << "DimKernel=" << DimKernel << " nbRow=" << nbRow << "\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
