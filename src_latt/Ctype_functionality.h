#ifndef INCLUDE_CTYPE_FUNCTIONALITY
#define INCLUDE_CTYPE_FUNCTIONALITY



struct triple {
  int8_t i;
  int8_t j;
  int8_t k;
};



template<typename T>
MyMatrix<T> CTYP_TheFlipping(MyMatrix<T> const& TheCtype, std::vector<triple> const& TheInfo)
{
  size_t n_rows = TheCtype.rows();
  size_t n_cols = TheCtype.cols();
  Face ListIchange(n_rows);
  for (auto & e_triple : TheInfo)
    ListIchange[e_triple.i] = 1;
  MyMatrix<T> RetMat(n_rows, n_cols);
  size_t idx=0;
  for (auto & e_triple : TheInfo) {
    int8_t j = triple.j;
    int8_t k = triple.k;
    //
    for (size_t i_col=0; i_col<n_cols; i_col++)
      RetMat(idx, i_col) = -TheCtype(j, i_col) + TheCtype(k, i_col);
    idx++;
    //
    for (size_t i_col=0; i_col<n_cols; i_col++)
      RetMat(idx, i_col) =  TheCtype(j, i_col) - TheCtype(k, i_col);
    idx++;
  }
  for (size_t i_row=0; i_row<n_rows; i_row++) {
    if (ListIchange(i_row) == 0) {
      for (size_t i_col=0; i_col<n_cols; i_col++)
        RetMat(idx, i_col) =  TheCtype(i_row, i_col);
      idx++;
    }
  }
  return RetMat;
}






#endif


