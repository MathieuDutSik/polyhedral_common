// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_MATRIXGROUPPLESKEN_H_
#define SRC_GROUP_MATRIXGROUPPLESKEN_H_

// This was an attempt at devising a scheme based on Plesken
// Souvignier kind of algorithms.

/*
There are several challenges for this implementation.
---We need an invariant vector family which we may obtain
   from a positive definite forms.
---We have already the one from the Lorentzian family or
   from the Q_inv. But it is unlikely to be of full rank.
---If the configuration is of full rank we win. If not, it is
   coming from an isotropic vector.
---From the fact that the group of interest is finite, we know
   that there exist a positive definite quadratic form that need
   to be found.
---From the Qinv, we can get a quadratic form preserving the
   n-1 dimensional space. That partial quadratic form, then has to be
   extended.
   In essence we ask that any isometry preserving the n-1 dimensional
   space also would preserve that one.
   We need to obtain that form from computation not involving the
   group.
---If we limit ourselves to the case of just preserving v^{perp} then
   we do not have finiteness. Let us take the form -x0^2 + x1^2 + x2^2
   The space v^{perp} is of dimension 2 and the restricted quadratic
   form is a_u^2, that is positive semidefinite but not definite. As a
   conclusion the space of isometries is
   (a_u, a_v) -> (a_u + C a_v, a_v).
   It is an isometry for all C because v is isotropic. All those isometries
   are extendible to full isometries of the 3-dim space by the known lemma.
---Under this infinite group there is a single orbit of isotropic vectors
   different from v.
---Now we have the additional finiteness requirement. That additional
   requirement is encapsulated into the preservation of the Qinv form.
   The vector v of v^{perp} is preserved. The orthogonal H of v in v^{perp}
   for Qinv is preserved as well.
---Now, we need to find another isotropic vector. The previous argument
   suggests that there is a single vector to be found.
   We can work with H^{perp} for the Lorentzian scalar product. It is
   a two dimensional space of signature (1,1). It has two isotropic
   vectors in it.
   One vector would be v, the other the one we are looking for.
---From this we can build a full dimensional vector system and nice
   quadratic form. And so with the Shortest vector problem, we can
   get vectors.
---From ListMatr, we can compute the space of invariant quadratic forms
   and this can get us good efficient signatures.
 */

/*
There are several challenges for this implementation.
---We need an invariant vector family which we may obtain
   from a positive definite forms.
---We have already the one from the Lorentzian family or
   from the Q_inv. But it is unlikely to be of full rank.
---If the configuration is of full rank we win. If not, it is
   coming from an isotropic vector.
---From the fact that the group of interest is finite, we know
   that there exist a positive definite quadratic form that need
   to be found.
---From the Qinv, we can get a quadratic form preserving the
   n-1 dimensional space. That partial quadratic form, then has to be
   extended.
   In essence we ask that any isometry preserving the n-1 dimensional
   space also would preserve that one.
   We need to obtain that form from computation not involving the
   group.
---If we limit ourselves to the case of just preserving v^{perp} then
   we do not have finiteness. Let us take the form -x0^2 + x1^2 + x2^2
   The space v^{perp} is of dimension 2 and the restricted quadratic
   form is a_u^2, that is positive semidefinite but not definite. As a
   conclusion the space of isometries is
   (a_u, a_v) -> (a_u + C a_v, a_v).
   It is an isometry for all C because v is isotropic. All those isometries
   are extendible to full isometries of the 3-dim space by the known lemma.
---Under this infinite group there is a single orbit of isotropic vectors
   different from v.
---Now we have the additional finiteness requirement. That additional
   requirement is encapsulated into the preservation of the Qinv form.
   The vector v of v^{perp} is preserved. The orthogonal H of v in v^{perp}
   for Qinv is preserved as well.
---Now, we need to find another isotropic vector. The previous argument
   suggests that there is a single vector to be found.
   We can work with H^{perp} for the Lorentzian scalar product. It is
   a two dimensional space of signature (1,1). It has two isotropic
   vectors in it.
   One vector would be v, the other the one we are looking for.
---From this we can build a full dimensional vector system and nice
   quadratic form. And so with the Shortest vector problem, we can
   get vectors.
---From ListMatr, we can compute the space of invariant quadratic forms
   and this can get us good efficient signatures.
 */

template <typename T, typename Telt>
std::vector<MyMatrix<T>>
GetListQuadraticForms(FiniteIsotropicMatrixGroupHelper<T, Telt> const &helper) {
  int n = helper.n;
  int n_vect = helper.EXTfaithful.rows();
  MyMatrix<T> BasisSp = RowReduction(helper.EXTfaithful);
  std::optional<MyMatrix<T>> opt1 =
      ListSolutionMat(BasisSp, helper.EXTfaithful);
#ifdef SANITY_CHECK_MATRIX_GROUP
  if (!opt1) {
    std::cerr << "opt1 : Failed to solution the linear systems\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<T> const &EXTfaithful_red = *opt1;
  std::optional<MyVector<T>> opt2 = SolutionMat(BasisSp, helper.Visotrop);
#ifdef SANITY_CHECK_MATRIX_GROUP
  if (!opt2) {
    std::cerr << "opt2 : Failed to solution the linear systems\n";
    throw TerminalException{1};
  }
#endif
  MyVector<T> const &Visotrop_red = *opt2;
  MyMatrix<T> Qinv = GetQmatrix(EXTfaithful_red);
  MyMatrix<T> eProd1 = MatrixFromVectorFamily(Qinv * Visotrop_red);
  MyMatrix<T> PosDefSpace = NullspaceMat(eProd1);
  MyMatrix<T> eProd2 = helper.G * PosDefSpace.transpose();
  MyMatrix<T> Sign11space = NullspaceMat(eProd2);
  MyMatrix<T> G11 = Sign11space * helper.G * Sign11space.transpose();
  auto get_second_isotropic = [&]() -> MyVector<T> {
    std::vector<MyVector<T>> LVect = GetBasisIsotropicVectors(G11);
    for (auto &V : LVect) {
      if (!IsVectorMultiple(V, helper.Visotrop)) {
        MyVector<T> Vret = Sign11space.transpose() * V;
        return Vret;
      }
    }
    std::cerr << "Find the other vector\n";
    throw TerminalException{1};
  };
  MyVector<T> VisotropSecond = get_second_isotropic();

  MyMatrix<T> FullBasis(n, n);
  for (int i_vect = 0; i_vect < n_vect; i_vect++) {
    for (int i = 0; i < n; i++)
      FullBasis(i_vect, i) = BasisSp(i_vect, i);
  }
  for (int i = 0; i < n; i++)
    FullBasis(n - 1, i) = VisotropSecond(i);
  //
  MyMatrix<T> prePosDefMat1 = ZeroMatrix<T>(n, n);
  for (int i = 0; i < n - 1; i++)
    for (int j = 0; j < n - 1; j++)
      prePosDefMat1(i, j) = Qinv(i, j);
  MyMatrix<T> prePosDefMat2 = ZeroMatrix<T>(n, n);
  prePosDefMat2(n - 1, n - 1) = 1;
  //
  MyMatrix<T> InvFullBasis = Inverse(FullBasis);
  MyMatrix<T> PosDefMat1 =
      InvFullBasis * prePosDefMat1 * InvFullBasis.transpose();
  MyMatrix<T> PosDefMat2 =
      InvFullBasis * prePosDefMat2 * InvFullBasis.transpose();

  return {PosDefMat1, PosDefMat2};
}

template <typename T, typename Telt>
std::vector<MyMatrix<T>>
GetListQuadraticForms(FiniteMatrixGroupHelper<T, Telt> const &helper) {
  MyMatrix<T> Qinv = GetQmatrix(helper.EXTfaithful);
  return {Qinv};
}

template <typename T, typename Tmod, typename Tgroup, typename Thelper>
std::optional<std::vector<MyMatrix<T>>>
PleskenSouvignier_Subspace_Stabilizer(std::vector<MyMatrix<T>> const &ListMatr,
                                      Thelper const &helper) {
  using Tint = typename underlying_ring<T>::ring_type;
  using Tidx = typename Tgroup::Telt::Tidx;
  int n = helper.n;
  std::vector<MyMatrix<T>> BasisSymmMat = BasisInvariantForm(n, ListMatr);
  std::vector<MyMatrix<T>> ListPosDef = GetListQuadraticForms(helper);
  std::vector<MyVector<T>> ListVect;
  std::vector<T> Vdiag;
  T idx = 0;
  for (auto &ePosDef : ListPosDef) {
    if (ePosDef.rows() > 1) {
      MyMatrix<T> eNSP = NullspaceIntMat(ePosDef);
      MyMatrix<Tint> eNSP_i = UniversalMatrixConversion<Tint, T>(eNSP);
      CanSolIntMat<Tint> eCan = ComputeCanonicalFormFastReduction(eNSP);
      MyMatrix<T> eGnsp = eNSP * ePosDef * eNSP.transpose();
      std::vector<MyMatrix<Tint>> ListMatrRed;
      for (auto &eMatr : ListMatr) {
        MyMatrix<Tint> eMatr_i = UniversalMatrixConversion<Tint, T>(eMatr);
        MyMatrix<Tint> eProd = eNSP_i * eMatr_i;
        std::optional<MyMatrix<Tint>> opt = CanSolutionIntMatMat(eCan, eProd);
#ifdef SANITY_CHECK_MATRIX_GROUP
        if (!opt) {
          std::cerr << "The NSP space is not preserved\n";
          throw TerminalException{1};
        }
#endif
        ListMatrRed.push_back(*opt);
      }
      MyMatrix<Tint> SHV_e =
          ExtractInvariantBreakingVectorFamily(eGnsp, ListMatrRed);
      MyMatrix<Tint> SHVbreak = SHV_e * eNSP_i;
      std::vector<T> Vdiag(SHV_e.rows(), idx);
      // Getting a full dimensional family from the other vectors.
      T jdx = 0;
      for (auto &fPosDef : ListPosDef) {
        if (idx != jdx) {
          MyMatrix<T> fNSP = NullspaceIntMat(fPosDef);
          MyMatrix<Tint> fNSP_i = UniversalMatrixConversion<Tint, T>(fNSP);
          MyMatrix<T> fGnsp = fNSP * fPosDef * fNSP.transpose();
          MyMatrix<Tint> SHV_f =
              ExtractInvariantVectorFamilyFullRank<T, Tint>(fGnsp);
          MyMatrix<Tint> fProd = SHV_f * fNSP_i;
          SHVbreak = Concatenation(SHVbreak, fProd);
          std::vector<T> V(SHV_f.rows(), jdx);
          Vdiag.insert(Vdiag.end(), V.begin(), V.end());
        }
        jdx++;
      }
      //
#ifdef SANITY_CHECK_MATRIX_GROUP
      if (RankMat(SHVbreak) != n) {
        std::cerr << "The matrix SHVbreak is not of full rank\n";
        throw TerminalException{1};
      }
#endif
      //
      const bool use_scheme = true;
      using Tfield = typename overlying_field<T>::field_type;
      std::vector<std::vector<Tidx>> ListListIdx =
          GetListGenAutomorphism_ListMat_Vdiag<T, Tfield, Tidx, use_scheme>(
              SHVbreak, BasisSymmMat, Vdiag);
      std::vector<MyMatrix<T>> NewListMatr;
      for (auto &eListIdx : ListListIdx) {
        std::optional<MyMatrix<T>> opt =
            FindMatrixTransformationTest(SHVbreak, SHVbreak, eListIdx);
#ifdef SANITY_CHECK_MATRIX_GROUP
        if (!opt) {
          std::cerr << "Failed to represent the permutation\n";
          throw TerminalException{1};
        }
#endif
        NewListMatr.push_back(*opt);
      }
    }
    idx++;
  }
  std::cerr << "Failed to find a breaking invariant family\n";
  throw TerminalException{1};
}

// clang-format off
#endif  // SRC_GROUP_MATRIXGROUP_H_
// clang-format on
