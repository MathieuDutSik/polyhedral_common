// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLYDECOMP_TRIPLES_H_
#define SRC_POLYDECOMP_TRIPLES_H_


#ifdef DEBUG
#define DEBUG_TRIPLE
#endif

/*
  When working with polyhedral complexes, one needs often to
  deal with a situation where we do not have simple equivalence
  for the lower dimensional cells.

  One example is perfect cones for which we have a nice equivalence
  theory, but for which the equivalence of vectors is not simple.
  ---
  The strategy is then to keep track of all the ways those data
  structures are contained in top dimensional cells. We keep track
  of all of those.

  The advantages:
  * We can compute the stabilizer of cells by having a generating
    set.
  * We have a relatively efficient algorithm which does not depend
    on so many underlying structures.
  * This was already done in
    + Lorentzian.g for GetFullComplex / SaturationAndStabilizer
    + ProjectiveSystem.g for SaturationAndStabilizer

  Decisions:
  * In that code we use the symmetries in order to compute
    the possible mappings. This is different from Lorentzian.g
    and ProjectiveSystem.g where we generate the full set of
    facets, without trying to exploit the symmetries.
    It is unsure whether this is the right strategy.
  * We want that code to be used whenever useful. Of course we
    want to avoid code duplication.
  * But at the same time the types must contain the relevant
    information but also contain other relevant information.
    The solution to that is to have minimal types which are
    declared before hand and accessed via subscript.
  * 

  The minimal type are going to be

  In order to model polyhedral complex, one way to keep track of
  it is to 
 */

template <typename Tint> struct sing_adj {
  size_t jCone;
  Face f_ext;
  MyMatrix<Tint> eMat;
};

// The minimal types used for 
template<typename Tint, typename Tgroup>
struct TopConeMin {
  using Telt = typename Tgroup::Telt;
  // EXT is used for checks, mapping representation from Tgroup to MyMatrix and finding matching elements.
  MyMatrix<Tint> EXT;
  Tgroup GRP_ext;
  std::vector<sing_adj<Tint>> l_sing_adj;
  MyMatrix<Tint> find_matrix(Telt const& x, [[maybe_unused]] std::ostream& os) const {
    return FindTransformation(EXT, EXT, x);
  }
};


// A face of the cellular complex is determined by all the ways in which
// f_ext is the subset of the corresponding face
template <typename Tint>
struct triple {
  size_t iCone;
  Face f_ext;
  MyMatrix<Tint> eMat;
};

template <typename Ttopcone>
std::optional<MyMatrix<typename Ttopcone::Tint>>
test_equiv_triple(std::vector<Ttopcone> const &l_cones,
                  triple<typename Ttopcone::Tint> const &ef1,
                  triple<typename Ttopcone::Tint> const &ef2,
                  std::ostream& os) {
  using Tgroup = typename Ttopcone::Tgroup;
  using Telt = typename Tgroup::Telt;
  using Tint = typename Ttopcone::Tint;
  if (ef1.iCone != ef2.iCone)
    return {};
  size_t iC = ef1.iCone;
  const Ttopcone &eC = l_cones[iC];
  std::optional<Telt> test =
      eC.GRP_ext.RepresentativeAction_OnSets(ef1.f_ext, ef2.f_ext);
  if (!test)
    return {};
  MyMatrix<Tint> eMat = eC.find_matrix(*test, os);
  return Inverse(ef1.eMat) * eMat * ef2.eMat;
}


template <typename Ttopcone>
triple<typename Ttopcone::Tint>
canonicalize_triple(std::vector<Ttopcone> const &l_cones,
                    triple<typename Ttopcone::Tint> const &t1,
                    std::ostream& os) {
  using Tgroup = typename Ttopcone::Tgroup;
  using Telt = typename Tgroup::Telt;
  using Tint = typename Ttopcone::Tint;
  size_t iC = t1.iCone;
  const Ttopcone &eC = l_cones[iC];
  Face f_can = eC.GRP_ext.OptCanonicalImage(t1.f_ext);
  // It is somewhat suboptimal to have to compute the equivalence.
  // It seems possible to modify NewCanonicImage in nsi.h
  // of permutalib in order to get the canonicalizing element
  // directly. But there is some work required for that.
  std::optional<Telt> test =
      eC.GRP_ext.RepresentativeAction_OnSets(t1.f_ext, f_can);
#ifdef SANITY_CHECK_TRIPLE
  if (!test) {
    std::cerr << "TRIP: Failed to find the equivalence\n";
    throw TerminalException{1};
  }
#endif
  Telt const& x = *test;
  MyMatrix<Tint> x_mat = eC.find_matrix(x, os);
  MyMatrix<Tint> eMat = t1.eMat * x_mat;
  triple<Tint> t2{iC, f_can, eMat};
#ifdef SANITY_CHECK_TRIPLE
  auto get_ext=[&](triple<Tint> const& t) -> std::set<MyVector<Tint>> {
    std::set<MyVector<Tint>> set;
    for (auto &ePt : FaceToVector<int>(t.f_ext)) {
      MyVector<Tint> V = GetMatrixRow(l_cones[t.iCone].EXT, ePt);
      MyVector<Tint> Vimg = t.eMat.transpose() * V;
      set.insert(Vimg);
    }
    return set;
  };
  if (get_ext(t1) != get_ext(t2)) {
    std::cerr << "TRIP: The canonicalization failed\n";
    throw TerminalException{1};
  }
#endif
  return t2;
}

template <typename Ttopcone>
std::optional<MyMatrix<typename Ttopcone::Tint>>
test_triple_in_listtriple(std::vector<Ttopcone> const &l_cones,
                          std::vector<triple<typename Ttopcone::Tint>> const &lt1,
                          triple<typename Ttopcone::Tint> const &ef2,
                          std::ostream& os) {
  using Tint = typename Ttopcone::Tint;
  for (auto &ef1: lt1) {
    std::optional<MyMatrix<Tint>> opt = test_equiv_triple(l_cones, ef1, ef2, os);
    if (opt) {
      MyMatrix<Tint> const& M = *opt;
      return M;
    }
  }
  return {};
}





/*
  Generate the list of entries in the face and the list of stabilizer generators
 */
template <typename Ttopcone>
std::pair<std::vector<triple<typename Ttopcone::Tint>>, std::vector<MyMatrix<typename Ttopcone::Tint>>>
get_spanning_list_triple(
    std::vector<Ttopcone> const &l_cones,
    const triple<typename Ttopcone::Tint> &ef_input, std::ostream& os) {
  //  os << "Beginning of get_spanning_list_triple\n";
  using Tgroup = typename Ttopcone::Tgroup;
  using Telt = typename Tgroup::Telt;
  using Tint = typename Ttopcone::Tint;
  std::vector<MyMatrix<Tint>> ListMatrGen;
  std::vector<MyVector<Tint>> EXTinner;
  // That value of dim should be overwritten later
  for (auto &ePt : FaceToVector<int>(ef_input.f_ext)) {
    MyVector<Tint> V = GetMatrixRow(l_cones[ef_input.iCone].EXT, ePt);
    // In GAP Vimg = V A and in transpose we get Vimg^T = A^T V^T
    MyVector<Tint> Vimg = ef_input.eMat.transpose() * V;
    EXTinner.push_back(std::move(Vimg));
  }
#ifdef DEBUG_TRIPLE
  os << "TRIP: |EXTinner|=" << EXTinner.size() << "\n";
#endif
  auto f_insert_generator = [&](const MyMatrix<Tint> &eMatrGen) -> void {
    ListMatrGen.push_back(eMatrGen);
#ifdef DEBUG_TRIPLE
    for (auto &eV : EXTinner) {
      MyVector<Tint> Vimg = eMatrGen.transpose() * eV;
      if (std::find(EXTinner.begin(), EXTinner.end(), Vimg) == EXTinner.end()) {
        std::cerr
            << "Error: The generator does not preserve the face globally\n";
        throw TerminalException{1};
      }
    }
#endif
  };
  std::vector<triple<Tint>> l_triple;
  auto f_insert = [&](const triple<Tint> &ef_A) -> void {
    for (const auto &ef_B : l_triple) {
      std::optional<MyMatrix<Tint>> equiv_opt =
        test_equiv_triple(l_cones, ef_A, ef_B, os);
      if (equiv_opt) {
        f_insert_generator(*equiv_opt);
        return;
      }
    }
    l_triple.emplace_back(std::move(ef_A));
    const Ttopcone &uC = l_cones[ef_A.iCone];
    Tgroup stab = uC.GRP_ext.Stabilizer_OnSets(ef_A.f_ext);
    MyMatrix<Tint> eInv = Inverse(ef_A.eMat);
    for (auto &eGen : stab.SmallGeneratingSet()) {
      MyMatrix<Tint> eMatGen = uC.find_matrix(eGen, os);
      MyMatrix<Tint> TransGen = eInv * eMatGen * ef_A.eMat;
      f_insert_generator(TransGen);
    }
  };
  f_insert(ef_input);
  size_t curr_pos = 0;
  while (true) {
    size_t len = l_triple.size();
    if (curr_pos == len) {
      break;
    }
#ifdef DEBUG_TRIPLE
    os << "TRIP: curr_pos=" << curr_pos << " len=" << len << "\n";
#endif
    for (size_t i = curr_pos; i < len; i++) {
      // We cannot use const& for triple because the underlying array
      // of std::vector can be freed and reallocated.
      triple<Tint> ef = l_triple[i];
      const Ttopcone &eC = l_cones[ef.iCone];
#ifdef DEBUG_TRIPLE
      os << "TRIP: i=" << i << " iCone=" << ef.iCone << "\n";
      size_t n_facet = 0;
#endif
      for (auto &e_sing_adj : l_cones[ef.iCone].l_sing_adj) {
        size_t jCone = e_sing_adj.jCone;
#ifdef DEBUG_TRIPLE
        os << "TRIP: |        ef.f_ext|=" << SignatureFace(        ef.f_ext) << "\n";
        os << "TRIP: |e_sing_adj.f_ext|=" << SignatureFace(e_sing_adj.f_ext) << "\n";
        int n_act = eC.GRP_ext.n_act();
        os << "TRIP: n_act=" << n_act << " jCone=" << jCone << "\n";
#endif
        std::vector<std::pair<Face, Telt>> l_pair =
            FindContainingOrbit(eC.GRP_ext, e_sing_adj.f_ext, ef.f_ext);
#ifdef DEBUG_TRIPLE
        os << "TRIP: |l_pair|=" << l_pair.size() << "\n";
#endif
#ifdef DEBUG_TRIPLE
        vectface vfo =
            OrbitFace(e_sing_adj.f_ext, eC.GRP_ext.GeneratorsOfGroup());
        n_facet += vfo.size();
        Tgroup stab = eC.GRP_ext.Stabilizer_OnSets(ef.f_ext);
        os << "TRIP: |vfo|=" << vfo.size() << " |stab|=" << stab.size() << "\n";
        size_t n_match = 0;
        vectface vfcont(vfo.get_n());
        for (auto &eFace : vfo) {
          if (is_subset(ef.f_ext, eFace)) {
            os << "TRIP: V=" << StringFace(eFace) << "\n";
            vfcont.push_back(eFace);
            n_match++;
          }
        }
        os << "TRIP: |vfcont|=" << vfcont.size() << "\n";
        vectface vfs = OrbitSplittingSet(vfcont, stab);
        os << "TRIP: |l_pair|=" << l_pair.size()
           << " |e_sing_adj.f_ext|=" << SignatureFace(e_sing_adj.f_ext)
           << " |ef.f_ext|=" << SignatureFace(ef.f_ext)
           << " n_match=" << n_match << " |vfs|=" << vfs.size() << "\n";
        if (vfs.size() != l_pair.size()) {
          std::cerr << "We have a size error\n";
          throw TerminalException{1};
        }
#endif
        for (auto &e_pair : l_pair) {
#ifdef DEBUG_TRIPLE
          os << "TRIP: face=" << StringFace(e_pair.first) << " elt=" << e_pair.second << "\n";
#endif
          MyMatrix<Tint> eMat1 = eC.find_matrix(e_pair.second, os);
          const Ttopcone &fC = l_cones[jCone];
          MyMatrix<Tint> eMatAdj = e_sing_adj.eMat * eMat1 * ef.eMat;
          MyMatrix<Tint> EXTimg = fC.EXT * eMatAdj;
#ifdef DEBUG_TRIPLE
          os << "TRIP: EXTinner=\n";
          WriteMatrix(os, MatrixFromVectorFamily(EXTinner));
          os << "TRIP: EXTimg=\n";
          WriteMatrix(os, EXTimg);
          std::set<MyVector<Tint>> setEXTinner;
          for (auto & V: EXTinner) {
            setEXTinner.insert(V);
          }
          std::vector<MyVector<Tint>> EXTinner_int_EXTimg;
          for (int iV=0; iV<EXTimg.rows(); iV++) {
            MyVector<Tint> V = GetMatrixRow(EXTimg, iV);
            if (setEXTinner.count(V) == 1) {
              EXTinner_int_EXTimg.push_back(V);
            }
          }
          os << "TRIP: EXTinner_int_EXTimg=\n";
          WriteMatrix(os, MatrixFromVectorFamilyDim(EXTimg.cols(), EXTinner));
#endif
          ContainerMatrix<Tint> Cont(EXTimg);
          Face faceNew(EXTimg.rows());
          for (auto &e_line : EXTinner) {
            std::optional<size_t> opt = Cont.GetIdx_v(e_line);
            if (!opt) {
              std::cerr << "TRIP: The vector is not in the image. Clear bug\n";
              throw TerminalException{1};
            }
            size_t idx = *opt;
            faceNew[idx] = 1;
          }
          triple<Tint> efNew{jCone, faceNew, eMatAdj};
          f_insert(efNew);
        }
      }
#ifdef DEBUG_TRIPLE
      os << "TRIP: iCone=" << ef.iCone
         << " |ef.f_ext|=" << SignatureFace(ef.f_ext)
         << " n_facet=" << n_facet << "\n";
#endif
    }
#ifdef DEBUG_TRIPLE
    os << "TRIP: Now |l_triple|=" << l_triple.size() << "\n";
#endif
    curr_pos = len;
  }
#ifdef DEBUG_TRIPLE
  os << "TRIP: |l_triple|=" << l_triple.size()
     << " |ListMatrGen|=" << ListMatrGen.size() << "\n";
  os << "l_triple.iCon =";
  for (auto &e_ent : l_triple)
    os << " " << e_ent.iCone;
  os << "\n";
#endif
  return {std::move(l_triple), std::move(ListMatrGen)};
}

// clang-format off
#endif  // SRC_POLYDECOMP_TRIPLES_H_
// clang-format on
