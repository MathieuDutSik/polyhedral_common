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
  int iOrb;
  Tint value;
  MyMatrix<Tint> M;
};

template<typename T, typename Tint, typename Tgroup>
struct ChainBuilder {
private:
  std::unordered_map<MyMatrix<Tint>, size_t> map_ext_set;
  int idim;
  PerfectComplexTopDimInfo<T,Tint,Tgroup> const& pctdi;
  ResultStepEnumeration<T, Tint, Tgroup> const& rse;
  std::ostream& os;
public:
  std::vector<FaceEntry<Tint>> chain;
  ChainBuilder(int _idim, PerfectComplexTopDimInfo<T,Tint,Tgroup> const& _pctdi, ResultStepEnumeration<T, Tint, Tgroup> const& _rse, std::ostream& _os) : pctdi(_pctdi), rse(_rse), os(_os) {
  }
  void f_insert(Tint const& value, int const& iOrb, MyMatrix<Tint> const& M) {
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
    MyMatrix<Tint> const& EXT1 = rse.level[idim][iOrb].EXT;
    OrientationInfo const& or_info = rse.level[idim][iOrb].or_info;
    std::vector<MyMatrix<T>> const& ListMat = pctdi.ListMat;
    MyMatrix<Tint> EXT2 = EXT1 * M;
    MyMatrix<Tint> EXT3 = tot_set(EXT2);
    int& idx = map_ext_set[EXT3];
    if (idx == 0) {
      idx = map_ext_set.size();
      FaceEntry<Tint> fe{iOrb, value, M};
      chain.emplace_back(std::move(fe));
    } else {
      int sign = 1;
      if (rse.level[idim][iOrb].is_orientable) {
        MyMatrix<Tint> t = Inverse(M) * chain_ret[idx - 1].M;
        sign = get_face_orientation(EXT1, ListMat, or_info, t);
      }
      chain[idx - 1].value += sign * value;
    }
  }
  std::vector<FaceEntry<Tint>> get_faces() const {
    std::vector<FaceEntry<Tint>> chain_ret;
    for (auto & fe: chain_ret) {
      if (fe.multiplicity != 0) {
        chain_ret.push_back(fe);
      }
    }
    return chain_ret;
  }
};


template<typename T, typename Tint, typename Tgroup>
bool is_equal_chain(std::vector<FaceEntry<Tint>> const& chain1, std::vector<FaceEntry<Tint>> const& chain2, int const& idim, PerfectComplexTopDimInfo<T,Tint,Tgroup> const& pctdi, ResultStepEnumeration<T, Tint, Tgroup> const& rse, std::ostream& os) {
  ChainBuilder<T,Tint,Tgroup> chain_builder(idim, pctdi, rse, os);
  for (auto & fe: chain1) {
    chain_builder(fe.value, fe.iOrb, fe.M);
  }
  for (auto & fe: chain1) {
    chain_builder(-fe.value, fe.iOrb, fe.M);
  }
  std::vector<FaceEntry<Tint>> l_faces = chain_builder.get_faces();
  if (l_faces.size() == 0) {
    return true;
  } else {
    return false;
  }
}




// Compute the boundary d(c) for c a part of the full complex.
template<typename T, typename Tint, typename Tgroup>
std::vector<FaceEntry<Tint>> compute_boundary(std::vector<FaceEntry<Tint>> const& c_idim, int const& idim, PerfectComplexTopDimInfo<T,Tint,Tgroup> const& pctdi, ResultStepEnumeration<T, Tint, Tgroup> const& rse, std::ostream& os) {
  ChainBuilder<T,Tint,Tgroup> chain_builder(idim+1, pctdi, rse, os);
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
      chain_builder.f_insert(value, ebnd.iOrb, M);
    }
  }
  return chain_builder.get_faces();
}

// We should have d_{i+1}(d_i(c)) = 0 consistently.
template<typename T, typename Tint, typename Tgroup>
bool is_product_zero(int const& idim, PerfectComplexTopDimInfo<T,Tint,Tgroup> const& pctdi, ResultStepEnumeration<T, Tint, Tgroup> const& rse, std::ostream& os) {
  int dim = rse.level[idim].size();
  int n = pctdi.LinSpa.n;
  for (int iOrb=0; iOrb<dim; iOrb++) {
    FaceEntry<Tint> fe{iOrb, Tint(1), IdentityMat<Tint>(n)};
    std::vector<FaceEntry<Tint>> chain1{fe};
    std::vector<FaceEntry<Tint>> chain2 = compute_boundary(chain1, idim, pctdi, rse, os);
    std::vector<FaceEntry<Tint>> chain3 = compute_boundary(chain2, idim+1, pctdi, rse, os);
    if (chain3.size() > 0) {
      return false;
    }
  }
}






template<typename T, typename Tint, typename Tgroup>
std::vector<FaceEntry<Tint>> contracting_homotopy(std::vector<FaceEntry<Tint>> const& c_idim, int const& idim, PerfectComplexTopDimInfo<T,Tint,Tgroup> const& pctdi, ResultStepEnumeration<T, Tint, Tgroup> const& rse, std::ostream& os) {
  

}








// clang-format off
#endif  // SRC_PERFECT_COMPLEX_AND_CONTRACTION_H_
// clang-format on
