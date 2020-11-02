#ifndef INCLUDE_CTYPE_FUNCTIONALITY
#define INCLUDE_CTYPE_FUNCTIONALITY



template<typename T>
std::vector<MyMatrix<T>> CTYP_GetBasis(int n)
{
  std::vector<MyMatrix<T>> ListSymmMat;
  for (int i=0; i<n; i++) {
    MyMatrix<T> eMat = ZeroMatrix<T>(n,n);
    eMat(i,i) = 1;
    ListSymmMat.push_back(eMat);
  }
  for (int i=0; i<n; i++)
    for (int j=i+1; j<n; j++) {
      MyMatrix<T> eMat = ZeroMatrix<T>(n,n);
      eMat(i,j) = 1;
      eMat(j,i) = 1;
      ListSymmMat.push_back(eMat);
    }
  return ListSymmMat;
}



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

template<typename T>
std::vector<triple> CTYP_GetListTriple(MyMatrix<T> const& TheCtype)
{
  int8_t n_edge = TheCtype.rows();
  int8_t n_cols = TheCtype.cols();
  std::vector<triple> ListTriples;
  auto get_position=[&](MyVector<T> const& eV, int8_t start_idx) -> int8_t {
    auto get_nature=[&](int8_t pos) -> bool {
      for (int8_t i_col=0; i_col<n_cols; i_col++)
        if (TheCtype(pos, i_col) != eV(i_col))
          return false;
      return true;
    };
    for (int8_t k=start_idx+1; k<n_edge; k++) {
      if (get_nature(k))
        return k
    }
    return -1;
  };
  for (int8_t i=0; i<n_edge; i++)
    for (int8_t j=i+1; j<n_edge; j++) {
      MyVector<T> eDiff(n_cols);
      for (int8_t i_col=0; i_col<n_cols; i_col++)
        eDiff(i_col) = - TheCtype(i, i_col) - TheCtype(j,i_col);
      int8_t pos = get_position(eDiff, j);
      if (pos != -1)
        ListTriples.push_back({i,j,pos});
    }
  return ListTriples;
}

template<typename T>
struct CtypeCombData {
  MyMatrix<T> ListInequalities;
  std::vector<std::vector<triple>> ListInformations;
};


template<typename T>
std::vector<MyMatrix<T>> CTYP_GetInequalitiesOfCtype(MyMatrix<T> const& TheCtype)
{
  std::vector<triple> ListTriples = CTYP_GetListTriples(TheCtype);
  int8_t n = TheCtype.cols();
  int8_t tot_dim = n*(n+1) / 2;
  auto ComputeInequality=[&](MyVector<T> const& V1, MyVector<T> const& V2) -> MyVector<T> {
    MyVector<T> TheVector(tot_dim);
    for (int8_t i=0; i<n; i++) {
      TheVector(idx) = V1(i) * V1(i) - V2(i) * V2(i);
      idx++;
    }
    for (int8_t i=0; i<n; i++)
      for (int8_t j=i+1; j<n; j++) {
        // Factor 2 removed for simplification and faster code.
        TheVector(idx) = V1(i) * V1(j) - V2(i) * V2(j);
        idx++;
      }
    return TheVector;
  };
  std::unordered_map<MyVector<T>, std::vector<triple>> Tot_map;
  auto FuncInsertInequality=[&](int8_t i, int8_t j, int8_t k) -> void {
    MyVector<T> V1(n), V2(n);
    for (int8_t i_col=0; i_col<n; i_col++) {
      V1(i_col) = 2 * TheCtype(k, i_col) + TheCtype(i, i_col);
      V2(i_col) = TheCtype(i, i_col);
    }
    MyVector<T> TheVector = ComputeInequality(V1, V2);
    triple TheInfo = {i,j,k};
    std::vector<triple>& list_trip = Tot_map.find(TheVector);
    list_trip.push_back(TheInfo);
  };
  for (auto & e_triple : ListTriples) {
    FuncInsertInequality(e_triple.i, e_triple.j, e_triple.k);
    FuncInsertInequality(e_triple.j, e_triple.k, e_triple.i);
    FuncInsertInequality(e_triple.k, e_triple.i, e_triple.j);
  }
  size_t n_ineq = Tot_map.size();
  MyMatrix<T> ListInequalities(n_ineq, tot_dim);
  std::vector<std::vector<triple>> ListInformations;
  size_t i_ineq=0;
  for (auto & kv : Tot_map) {
    for (int8_t i_col=0; i_col<tot_dim; i_col++)
      ListInequalities(i_ineq, i_col) = kv.first(i_col);
    i_ineq++;
    ListInformations.push_back(std::move(kv.second));
  }
  // Reducing by redundancy
  std::vector<int> ListIrred = cdd::RedundancyReductionClarkson(ListInequalities);
  std::vector<MyMatrix<T>> ListCtype;
  for (auto & e_int : ListIrred) {
    ListCtype.push_back(CTYP_TheFlipping(TheCtype, ListInformations[e_int]));
  }
  return ListCtype;
}







#endif


