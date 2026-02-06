// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_COMPLEX_AND_CONTRACTION_H_
#define SRC_PERFECT_COMPLEX_AND_CONTRACTION_H_

template<typename Tint>
struct BoundEntry {
  int iOrb;
  int sign;
  MyMatrix<Tint> M;
};

template<typename Tint>
struct ListBoundEntry {
  std::vector<BoundEntry<Tint>> l_bound;
};

template<typename Tint>
struct FullBoundary {
  std::vector<ListBoundEntry<Tint>> ll_bound;
};


template<typename Tint>
struct FaceEntry {
  Tint value;
  int iOrb;
  MyMatrix<Tint> M;
};






// Compute the boundary d(c) for c a part of the full complex.
template<typename T, typename Tint, typename Tgroup>
std::vector<FaceEntry<Tint>> compute_boundary(std::vector<FaceEntry<Tint>> const& c_idim, int const& idim, PerfectComplexTopDimInfo<T,Tint,Tgroup> const& pctdi, ResultStepEnumeration<T, Tint, Tgroup> const& rse, std::ostream& os) {
  std::unordered_map<MyMatrix<Tint>, size_t> map_ext_set;
  std::vector<FaceEntry<Tint>> chain_ret;
  auto tot_set=[&](MyMatrix<Tint> const& EXTin) -> MyMatrix<Tint> {
    std::set<MyVector<Tint>> set;
    for (int i_row=0; i_row<EXTin.rows(); i_row++) {
      MyVector<Tint> V = GetMatrixRow(EXTin, i_row);
      set.insert(V);
    }
    MyMatrix<Tint> M(EXTin.rows(), EXTin.cols());
    int pos = 0;
    for (auto & V: set) {
      AssignMatrixRow(M, pos, V);
      pos += 1;
    }
    return M;
  };
  auto f_insert=[&](Tint const& value, int const& iOrb, MyMatrix<Tint> const& M) -> void {
    MyMatrix<Tint> const& EXT1 = rse.level[idim + 1][iOrb].EXT;
    OrientationInfo const& or_info = rse.level[idim + 1][iOrb].or_info;
    std::vector<MyMatrix<T>> const& ListMat = pctdi.ListMat;
    MyMatrix<Tint> EXT2 = EXT1 * M;
    MyMatrix<Tint> EXT3 = tot_set(EXT2);
    int& idx = map_ext_set[EXT3];
    if (idx == 0) {
      idx = map_ext_set.size();
      FaceEntry<Tint> fe{value, iOrb, M};
      chain_ret.emplace_back(std::move(fe));
    } else {
      int sign = 1;
      if (rse.level[idim + 1][iOrb].is_orientable) {
        MyMatrix<Tint> t = Inverse(M) * chain_ret[idx - 1].M;
        sign = get_face_orientation(EXT1, ListMat, or_info, t);
      }
      chain_ret[idx - 1].value += sign * value;
    }
  };
  if (!rse.boundary) {
    std::cerr << "COMPCONT: The boundary should have been constructed\n";
    throw TerminalException{1};
  }
  FullBoundary<Tint> const& bnd = *rse.boundary;
  for (auto & fe: c_idim) {
    int iOrb = fe.iOrb;
    for (auto & ebnd: bnd.ll_bound[idim].l_bound[iOrb]) {
      Tint value = fe.value * ebnd.sign;
      MyMatrix<Tint> M = ebnd.M * fe.M;
      f_insert(value, ebnd.iOrb, M);
    }
  }
  std::vector<FaceEntry<Tint>> chain_ret_red;
  for (auto & fe: chain_ret) {
    if (fe.multiplicity != 0) {
      chain_ret_red.emplace_back(std::move(fe));
    }
  }
  return chain_ret_red;
}

// clang-format off
#endif  // SRC_PERFECT_COMPLEX_AND_CONTRACTION_H_
// clang-format on
