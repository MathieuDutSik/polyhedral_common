#ifndef DEFINE_MATRIX_SERIALIZATION_H
#define DEFINE_MATRIX_SERIALIZATION_H


namespace boost { namespace serialization {

    //
    // MyMatrix data type
    //

    
    template<class Archive, typename T>
      inline void serialize(Archive & ar,
			    MyMatrix<T> & matrix,
			    const unsigned int version)
      {
	int rows = matrix.rows();
	int cols = matrix.cols();
	ar & make_nvp("rows", rows);
	ar & make_nvp("cols", cols);
	matrix.resize(rows, cols); // no-op if size does not change!
	
	// always save/load row-major
	for(int r = 0; r < rows; ++r)
	  for(int c = 0; c < cols; ++c)
	    ar & make_nvp("val", matrix(r,c));
      }

    //
    // MySparseMatrix data type
    //

    template<class Archive, typename T>
      inline void load(Archive & ar,
		       MySparseMatrix<T> & val,
		       const unsigned int version)
      {
	std::cerr << "load(MySparseMatrix<T>), step 1\n";
	int nbRow, nbCol, nnz;
	ar & make_nvp("rows", nbRow);
	ar & make_nvp("cols", nbCol);
	ar & make_nvp("nnz", nnz);
	using T2=Eigen::Triplet<T>;
	std::vector<T2> tripletList(nnz);
	for (int iNNZ=0; iNNZ<nnz; iNNZ++) {
	  int iRow, iCol;
	  T eVal;
	  ar & make_nvp("iRow", iRow);
	  ar & make_nvp("iCol", iCol);
	  ar & make_nvp("eVal", eVal);
	  tripletList[iNNZ]=T2(iRow,iCol,eVal);
	}
	MySparseMatrix<T> SpMat(nbRow, nbCol);
	SpMat.setFromTriplets(tripletList.begin(), tripletList.end());
	val=SpMat;
	std::cerr << "load(MySparseMatrix<T>), step 2\n";
      }
    
    template<class Archive, typename T>
      inline void save(Archive & ar,
		       MyMatrix<T> const& val,
		       const unsigned int version)
      {
	std::cerr << "save(MySparseMatrix<T>), step 1\n";
	int nbRow=val.rows();
	int nbCol=val.cols();
	int nnz=val.nonZeros();
	ar & make_nvp("rows", nbRow);
	ar & make_nvp("cols", nbCol);
	ar & make_nvp("nnz", nnz);
	for (int k=0; k<val.outerSize(); ++k)
	  for (typename MySparseMatrix<T>::InnerIterator it(val,k); it; ++it) {
	    T eVal=it.value();
	    int iRow=it.row();   // row index
	    int iCol=it.col();   // col index (here it is equal to k)
	    ar & make_nvp("iRow", iRow);
	    ar & make_nvp("iCol", iCol);
	    ar & make_nvp("eVal", eVal);
	  }
	std::cerr << "save(MySparseMatrix<T>), step 2\n";
      }
    
    template<class Archive, typename T>
      inline void serialize(Archive & ar,
			    MySparseMatrix<T> & val,
			    const unsigned int version)
      {
	std::cerr << "split_free(MySparseMatrix<T>), step 1\n";
	split_free(ar, val, version);
	std::cerr << "split_free(MySparseMatrix<T>), step 2\n";
      }
}}    




#endif
