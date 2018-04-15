template<typename T>
struct MultiArray {
public:
  MultiArray()
  {
    dim=0;
    ListElt=nullptr;
  }
  MultiArray(std::vector<size_t> const& inpVectDim)
  {
    dim=eVectDim.size();
    int eProd=1;
    for (int iDim=0; iDim<dim; iDim++) {
      int eDim=inpVectDim[iDim];
      eProd=eProd*eDim;
      if (eDim < 0) {
	std::cerr << "iDim=" << iDim << " val=" << eVectDim[iDim] << "\n";
	std::cerr << "(zero dimension allowed for things like nullspace)\n";
	exit(1);
      }
    }
    eVectDim=inpVectDim;
    ListElt=new T[eProd];
  }
  MultiArray(MultiArray<T> const& eArr)
  {
    eVectDim=eArr.GetDims();
    dim=eVectDim.size();
    int eProd=1;
    for (int iDim=0; iDim<dim; iDim++) {
      int eDim=eVectDim[iDim];
      eProd=eProd*eDim;
    }
    ListElt=new T[eProd];
    for (int ielt=0; ielt<eProd; ielt++) {
      T eVal=eArr.getelt(ielt);
      ListElt[ielt]=eVal;
    }
  }
  T getelt(int idx) const
  {
    return ListElt[idx];
  }
  void setelt(int idx, T eVal)
  {
    ListElt[idx]=eVal;
  }
  MultiArray<T> operator=(MultiArray<T> const& eArr)
  {
    eVectDim=eArr.GetDims();
    int eProd=1;
    for (int iDim=0; iDim<dim; iDim++) {
      int eDim=eVectDim[iDim];
      eProd=eProd*eDim;
    }
    ListElt=new T[eProd];
    for (int ielt=0; ielt<eProd; ielt++) {
      T eVal=eArr.getelt(ielt);
      ListElt[ielt]=eVal;
    }
    return *this;
  }
  ~MultiArray()
  {
    delete [] ListElt;
  }
  // easier parts of the job
  std::vector<size_t> GetDims(void) const
  {
    return eVectDim;
  }
  void a(std::vector<size_t> const& ePos, T const& eVal)
  {
    int idx=0;
    int eProd=1;
    for (int iDim=0; iDim<dim; iDim++) {
      int jDim=dim-1-iDim;
      if (ePos[jDim] < 0 || ePos[jDim] >= eVectDim[jDim]) {
	std::cerr << "ePos is incorrect\n";
	std::cerr << "jDim=" << jDim << " eDim=" << eVectDim[jDim] << "\n";
	std::cerr << "ePos=" << ePos[jDim] << "\n";
	exit(1);
      }
      idx=idx*eVectDim[jDim] + ePos[jDim];
      eProd=eProd*eVectDim[jDim];
    }
    if (idx >= eProd) {
      std::cerr << "Wrong position for idx\n";
      std::cerr << "idx=" << idx << "\n";
    }
    ListElt[idx]=eVal;
  }
  T g(std::vector<size_t> const& ePos) const
  {
    int idx=0;
    int eProd=1;
    for (int iDim=0; iDim<dim; iDim++) {
      int jDim=dim-1-iDim;
      if (ePos[jDim] < 0 || ePos[jDim] >= eVectDim[jDim]) {
	std::cerr << "ePos is incorrect\n";
	std::cerr << "jDim=" << jDim << " eDim=" << eVectDim[jDim] << "\n";
	std::cerr << "ePos=" << ePos[jDim] << "\n";
	exit(1);
      }
      idx=idx*eVectDim[jDim] + ePos[jDim];
      eProd=eProd*eVectDim[jDim];
    }
    if (idx >= eProd) {
      std::cerr << "Wrong position for idx\n";
      std::cerr << "idx=" << idx << "\n";
    }
    return ListElt[idx];
  }
  MultiArray<T> DimensionExtraction(size_t const& iDim, size_t const& eDim) const
  {
    std::vector<size_t> NewVectDim;
    for (size_t i=0; i<dim; i++)
      if (i != iDim)
	NewVectDim.push_back(eVectDim[i]);
    MultiArray<T> NewArr(NewVectDim);
    size_t eProdLow=1;
    for (size_t i=0; i<iDim; i++)
      eProdLow=eProdLow*eVectDim[i];
    size_t eProdUpp=1;
    for (size_t i=iDim+1; i<dim; i++)
      eProdUpp=eProdUpp*eVectDim[i];
    size_t eShift=eProdLow*eVectDim[iDim];
    for (size_t iUpp=0; iUpp<eProdUpp; iUpp++)
      for (size_t iLow=0; iLow<eProdLow; iLow++) {
	size_t idx1=iLow + eProdLow*eDim + eShift*iUpp;
	size_t idx2=iLow + eProdLow*iUpp;
	T eVal=ListElt[idx1];
	NewArr.setelt(idx2, eVal);
      }
    return NewArr;
  }
private:
  size_t dim;
  std::vector<size_t> eVectDim;
  T *ListElt;
};

