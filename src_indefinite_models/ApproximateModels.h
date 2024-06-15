// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_MODELS_APPROXIMATEMODELS_H_
#define SRC_INDEFINITE_MODELS_APPROXIMATEMODELS_H_


template<typename T>
struct ApproximateModel {
  std::function<std::vector<MyVector<T>>(void)> GetApproximateGroup;
  std::function<void(std::vector<MyMatrix<T>>)> SetListClassesOrbitwise;
  std::function<std::vector<MyVector<T>>(T)> GetCoveringOrbitRepresentatives;
  std::function<MyVector<T>(T)> GetOneOrbitRepresentative;
};

template<typename T>
bool INDEF_FORM_IsEven(MyMatrix<T> const& Qmat) {
  if (!IsIntegralMatrix(Qmat)) {
    return false;
  }
  T two(2);
  for (int u=0; u<Qmat.rows(); u++) {
    T res = ResInt(Qmat(u,u), two);
    if (res != 0) {
      return false;
    }
  }
  return true;
}


template<typename T>
struct InternalEichler {
  MyMatrix<T> Qmat;
  std::vector<MyVector<T>> ListClasses;
  std::vector<MyVector<T>> GRP;
};

/*
  Based on paragraph 10 of Eichler book
  Quadratische Formen und orthogonale gruppen
  Also used the Scattone thesis as basis.

  Note: There are other works based on Eichler work that are analogous:
  ---Wall CTC, On The Orthogonal Groups of Unimodular Quadratic Forms
  ---James D.G. Indefinite Quadratic Forms of Determinant pm2p
  They also require the same hypothesis of having two hyperplanes.
*/
template<typename T, typename Tgroup>
ApproximateModel<T> INDEF_FORM_EichlerCriterion_TwoHyperplanesEven(MyMatrix<T> const& Qmat) {
  int n = Qmat.rows();
  std::vector<int> LIdx{1,0,3,2};
  for (int iRow=0; iRow<4; iRow++) {
    int iColCrit = LIdx[iRow];
    for (int iCol=0; iCol<4; iCol++) {
      if (iCol == iColCrit) {
        if (Qmat(iRow,iCol) != 1) {
          std::cerr << "We do not have the structure of 2U, A\n";
          throw TerminalException{1};
        }
      } else {
        if (Qmat(iRow,iCol) != 0) {
          std::cerr << "We do not have the structure of 2U, B\n";
          throw TerminalException{1};
        }
      }
    }
  }
  for (int iRow=0; iRow<4; iRow++) {
    for (int iCol=4; iCol<n; iCol++) {
      if (Qmat(iRow,iCol) != 0) {
        std::cerr << "We do not have the structure of 2U, C\n";
        throw TerminalException{1};
      }
    }
  }
  if (!INDEF_FORM_IsEven(Qmat)) {
    std::cerr << "MODIND: The lattice is not even\n";
    throw TerminalException{1};
  }
  // ComputeClasses
  std::vector<MyVector<T>> ListClasses = GetTranslationClasses(Gmat);
  std::vector<MyMatrix<T>> GRPstart;
  InternalEichler<T> ie{Qmat, ListClasses, GRPstart};
  std::shared_ptr<InternalEichler<T>> shr_ptr=make_shared<InternalEichler<T>>(ie);
  std::function<void(std::vector<MyMatrix<T>>)> SetListClassesOrbitwise = [=](std::vector<MyMatrix<T>> const& GRPmatr) -> void {
    shr_ptr->GRPstart = GRPmatr;
    MyMatrix<T> const& Qmat = shr_ptr->Qmat;
    MyVector<T> ZerV = ZeroVector<T>(4);
    std::vector<MyVector<T>> ListClassesExt;
    for (auto & eV : ListClasses) {
      MyVector<T> NewV = Concatenate(ZerV, eV);
      ListClassesExt.emplace_back(std::move(NewV));
    }
    auto GetPosition=[&](MyVector<T> const& eVect) -> size_t {
      for (size_t i=0; i<ListClassesExt.size(); i++) {
        MyVector<T> eDiff = ListClassesExt[i] - eVect;
        std::optional<MyVector<T>> opt = SolutionIntMat(Qmat, eDiff);
        if (opt) {
          return i;
        }
      }
      std::cerr << "Failed to find a matching coset in GetPosition\n";
      throw TerminalException{1};
    };
    using Telt = typename Tgroup::Telt;
    using Tidx = typename Telt::Tidx;
    std::vector<Telt> ListPermGens;
    for (auto & eMatrGen : GRPmatr) {
      std::vector<Tidx> eList;
      for (auto & eClassExt : ListClassesExt) {
        MyVector<T> x_eM = eMatrGen.transpose() * eClassExt;
        eList.push_back(GetPosition(x_eM));
      }
      Telt ePerm(eList);
      ListPermGens.push_back(ePerm);
    }
    Tgroup GRPperm(ListPermGens);
    std::vector<size_t> ListRepr = DecomposeOrbitPoint_FullRepr<Tgroup>(GRPperm);
    std::vector<MyVector<T>> ListClasses;
    for (auto & iRepr : ListRepr) {
      MyVector<T> eV = shr_ptr->Listclasses[iRepr];
      ListClasses.push_back(eV);
    }
    shr_ptr->ListClasses = ListClasses;
  };
  std::function<std::vector<MyVector<T>>(void)> GetApproximateGroup = [=]() -> std::vector<MyVector<T>> {
    MyMatrix<T> Qmat = shr_ptr->Qmat;
    int n = shr_ptr->Qmat.rows();
    // The quadratic form 2 x1 x2 + 2 x3 x4
    // We formulate it as (1/2) det U(x1,x2,x3,x4)
    // with U(x1,x2,x3,x4) = [ [x1 , -x4] , [x3 , x2] ]
    // If we multiply by a matrix in SL(2,Z) on the left and right
    auto GetLeftMultiplication = [&](MyMatrix<T> const& P) -> MyMatrix<T> {
      // We compute the product [[a, b],[c,d]] U(x1,x2,x3,x4)
      // This gives us
      // [[a x1 + b x3 , -a x4 + b x2] , [ c x1 + d x3 , -c x4 + d x2]]
      // So Y = XA,
      // Y1 =  a x1 + b x3
      // Y2 = -c x4 + d x2
      // Y3 =  c x1 + d x3
      // Y4 =  a x4 - b x2
      T a = P(0,0);
      T b = P(0,1);
      T c = P(1,0);
      T d = P(1,1);
      MyMatrix<T> RetMat = ZeroMatrix<T>(4,4);
      RetMat(0,0) = a;
      RetMat(0,2) = c;
      RetMat(1,1) = d;
      RetMat(1,3) = -b;
      RetMat(2,0) = b;
      RetMat(2,2) = d;
      RetMat(3,1) = -c;
      RetMat(3,3) = a;
      return RetMat;
    };
    auto GetRightMultiplication = [&](MyMatrix<T> const& P) -> MyMatrix<T> {
      // We compute the product [ [x1 , -x4] , [x3 , x2] ]      [[a, b],[c,d]]
      // This gives us
      // [[a x1 - c x4 , b x1 - d x4] , [a x3 + c x2 , b x3 + d x2]]
      // So Y = XA,
      // Y1 = a x1 - c x4
      // Y2 = b x3 + d x2
      // Y3 = a x3 + c x2
      // Y4 = -b x1 + d x4
      T a = P(0,0);
      T b = P(0,1);
      T c = P(1,0);
      T d = P(1,1);
      MyMatrix<T> RetMat = ZeroMatrix<T>(4,4);
      RetMat(0,0) = a;
      RetMat(0,3) = -b;
      RetMat(1,1) = d;
      RetMat(1,2) = c;
      RetMat(2,1) = b;
      RetMat(2,2) = a;
      RetMat(3,0) = -c;
      RetMat(3,3) = d;
      return RetMat;
    };
    auto ExpandGenerator=[&](MyMatrix<T> const& eGen) -> MyMatrix<T> {
      MyMatrix<T> BigMat = IdentityMat<T>(n);
      for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
          BigMat(i,j) = eGen(i,j);
        }
      }
      return BigMat;
    };
    std::vector<MyMatrix<T>> ListGenerators;
    auto FuncInsert=[&](MyMatrix<T> const& TheMat) -> void {
      if (!IsIntegralMatrix(TheMat)) {
        std::cerr << "The matrix TheMat is not integral\n";
        throw TerminalException{1};
      }
    };

        ListGenerators:=[];
        FuncInsert:=function(TheMat)
            if IsIntegralMat(TheMat)=false then
                Error("
            fi;
            if TheMat * Qmat * TransposedMat(TheMat) <> Qmat then
                Error("The matrix is not preserving the quadratic form");
            fi;
            Add(ListGenerators, TheMat);
        end;
        for eGen in GeneratorsOfGroup(GRPstart)
        do
            FuncInsert(eGen);
        od;
        # Generators of SL2(Z)
        eGenS:=[[0,-1],[1,0]];
        eGenT:=[[1,1],[0,1]];
        for eGen in [eGenS, eGenT]
        do
            FuncInsert(ExpandGenerator(GetLeftMultiplication(eGen)));
            FuncInsert(ExpandGenerator(GetRightMultiplication(eGen)));
        od;
        # Now looking at the isotropic vectors, generating the Eichler transvections
        for i in [1..4]
        do
            eVect:=ListWithIdenticalEntries(n,0);
            eVect[i]:=1;
            eProd:=eVect * Qmat;
            BasisOrth:=NullspaceIntMat(TransposedMat([eProd]));
            for eOrth in BasisOrth
            do
                FuncInsert(INDEF_FORM_Eichler_Transvection(Qmat, eVect, eOrth));
            od;
        od;
        return Group(ListGenerators);
    end;
    # Notion of divisor
    # For an element y in a lattice L, define div(y) the integer such that
    # (L,y) = div(y) Z.
    # For each v in L* we can find smallest d such that w = dv in L.
    # And we then have div(w) = d
    EnumerateVectorOverDiscriminant:=function(X)
        # For each vector v in M* / M
        # We want (x, v) divided by d.
        #
        # We should have e_I.x = d alpha_I
        # We should have (x,x) = X
        # For the first coordinates, we can write them as (0 , 0 , d , du)
        # By changing the value of u, we can shift the value of u which shifts the value
        # by 2 d^2, + or -.
        # The lattice is even so the 2 is not a problem.
        #
        # So, if we have some X, we want to find u\in Z and v in Z^{n-4} such that X = 2d^2 u + A[v]
        # with A the Gram matrix of M.
        # v must actually be of the form v = v0 + d A^{-1} h.
        # We are only interested on the value modulo 2d^2 of A on those vectors. So, it is a finite problem.
        # But it is a hard one as the search space might be very large.
        #
        if X mod 2 = 1 then
            return [];
        fi;
        Xdiv2 := X / 2;
        # The two hyperbolic planes are unimodular so for the discriminant, we only need
        # to consider the Gmat part.
        Ginv:=Inverse(Gmat);
        Print("Det(Gmat)=", DeterminantMat(Gmat), "\n");
        ListSolution:=[];
        for eClass1 in ListClasses
        do
            # The vector eClass1 is defined modulo an element of Gmat Z^n
            # Thus eClass2 is defined modulo an element of Z^n. This is an elements of L*
            # Thus eClass3 is defined modulo an element of d Z^n
            # So, we can write v = eClass3 + d w
            # which gets us A[v] = A[eClass3] + 2d w^T Gmat eClass3 + d^2 A[w]
            # So, the problem is actually linear.
            eClass2:=eClass1 * Ginv;
            DivD_frac:=RemoveFractionPlusCoef(eClass2).TheMult;
            DivD:=NumeratorRat(DivD_frac);
            eClass3:=eClass2 * DivD;
            #
            X_res1:=Xdiv2 mod DivD;
            X_res2:=Xdiv2 mod (DivD * DivD);
            Aclass3_norm:=(eClass3 * Gmat * eClass3) / 2;
            Aclass3_res1:=Aclass3_norm mod DivD;
            Aclass3_res2:=Aclass3_norm mod (DivD * DivD);
            if Aclass3_res1=X_res1 then # if not we cannot find a solution
                # and so
                # (X - A[eClass3])/2 = 2d w^T Gmat eClass3 + ....
                diff:=(X_res2 - Aclass3_res2) / DivD;
                eProd:=eClass3 * Gmat;
                eRec:=GcdVector(eProd);
                if eRec.TheGcd=0 then
                    quot:=0;
                else
                    quot:=diff / eRec.TheGcd;
                fi;
                if IsInt(quot) then
                    w:=quot * eRec.ListCoef;
                    eClass4:=eClass3 + DivD * w;
                    Aclass4_norm:=(eClass4 * Gmat * eClass4) / 2;
                    Aclass4_res2:=Aclass4_norm mod (DivD * DivD);
                    if Aclass4_res2<>X_res2 then
                        Error("A bug to resolve");
                    fi;
                    u:=(Xdiv2 - Aclass4_norm) / (DivD * DivD);
                    eSolution:=Concatenation([0,0,DivD, DivD*u], eClass4);
                    eNorm:=eSolution*Qmat*eSolution;
                    if eNorm<>X then
                        Error("A bug to resolve");
                    fi;
                    Add(ListSolution, eSolution);
                fi;
            fi;
        od;
        return ListSolution;
    end;
    GetCoveringOrbitRepresentatives:=function(X)
        local ListDiv, ListSolution, eDiv, eSol;
        if X = 0 then
            return EnumerateVectorOverDiscriminant(0);
        fi;
        ListDiv:=GetSquareDivisors(X);
        ListSolution:=[];
        for eDiv in ListDiv
        do
            for eSol in EnumerateVectorOverDiscriminant(X / (eDiv*eDiv))
            do
                Add(ListSolution, eDiv * eSol);
            od;
        od;
        return ListSolution;
    end;
    GetOneOrbitRepresentative:=function(X)
        local Xdiv;
        if X mod 2 = 1 then
            return fail;
        fi;
        Xdiv:=X / 2;
        return Concatenation([1,Xdiv], ListWithIdenticalEntries(n-2, 0));
    end;
    return rec(ApproximateGroup:=GetApproximateGroup(),
               SetListClassesOrbitwise:=SetListClassesOrbitwise,
               GetCoveringOrbitRepresentatives:=GetCoveringOrbitRepresentatives,
               GetOneOrbitRepresentative:=GetOneOrbitRepresentative);


}



// clang-format off
#endif  // SRC_INDEFINITE_MODELS_APPROXIMATEMODELS_H_
// clang-format on
