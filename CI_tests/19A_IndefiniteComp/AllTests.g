Read("../common.g");
Read("../access_points.g");

# The programs are built with ENABLE_FLINT_SUPPORT=1 and run on flint,
# the faster arithmetic for these computations.
IndefOptions:=rec(arith:="flint");;
Print("Beginning of AllTests\n");

TestStab:=function(eList)
    local eProg, eGram, TmpDir, FileIn, FileOut, FileErr, TheCommand, U, eP;
    Print("Running TestStab with INDEF_FORM_AutomorphismGroup\n");
    eGram:=GetGramMatrixFromList(eList);
    GRP:=INDEF_FORM_Stabilizer(eGram, IndefOptions);
    if GRP="program failure" then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    for eP in GeneratorsOfGroup(GRP)
    do
        if eP * eGram * TransposedMat(eP) <> eGram then
            Print("TestEquiv: Failed at U * eGram * TransposedMat(U) <> eGram\n");
            return false;
        fi;
    od;
    return true;
end;

TestEquiv:=function(eList)
    local TmpDir, eProg, eGram1, eP, eGram2, FileM1, FileM2, FileOut, FileErr, TheCommand, U;
    Print("Running TestEquiv with INDEF_FORM_TestEquivalence\n");
    eGram1:=GetGramMatrixFromList(eList);
    eP:=RandomIntegralUnimodularMatrix(Length(eGram1));
    eGram2:=eP * eGram1 * TransposedMat(eP);
    U:=INDEF_FORM_TestEquivalence(eGram1, eGram2, IndefOptions);
    if U="program failure" then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    if U = fail then
        Print("eGram1=\n");
        PrintArray(eGram1);
        Print("eGram2=\n");
        PrintArray(eGram2);
        Print("TestEquiv: Failed at U = fail\n");
        return rec(fct:="TestEquiv", eGram1:=eGram1, eGram2:=eGram2, eP:=eP);
    fi;
    if U * eGram1 * TransposedMat(U) <> eGram2 then
        Print("TestEquiv: Failed at U * eGram1 * TransposedMat(U) <> eGram2\n");
        return false;
    fi;
    return true;
end;

GetNrOrbitRepresentative:=function(eGram, Xnorm)
    local U;
    Print("Running GetNrOrbitRepresentative with INDEF_FORM_GetOrbitRepresentative\n");
    U:=INDEF_FORM_GetOrbitRepresentative(eGram, Xnorm, IndefOptions);
    if U="program failure" then
        Print("The output file is not existing. That qualifies as a fail\n");
        return rec(is_correct:=false);
    fi;
    return rec(is_correct:=true, nr_orbit:=Length(U));
end;

TestOrbitRepresentative:=function(eList)
    local eGram1, eP, eGram2, Xnorm, res1, res2, nr_orbit1, nr_orbit2;
    eGram1:=GetGramMatrixFromList(eList);
    eP:=RandomIntegralUnimodularMatrix(Length(eGram1));
    eGram2:=eP * eGram1 * TransposedMat(eP);
    for Xnorm in [0, 2]
    do
        res1:=GetNrOrbitRepresentative(eGram1, Xnorm);
        if res1.is_correct=false then
            Print("TestOrbitRepresentative: Failed at res1.is_correct\n");
            return false;
        fi;
        nr_orbit1:=res1.nr_orbit;
        #
        res2:=GetNrOrbitRepresentative(eGram2, Xnorm);
        if res2.is_correct=false then
            Print("TestOrbitRepresentative: Failed at res2.is_correct\n");
            return false;
        fi;
        nr_orbit2:=res2.nr_orbit;
        #
        if nr_orbit1<>nr_orbit2 then
            Print("TestOrbitRepresentative: Failed at nr_orbi1<>nr_orbit2\n");
            return false;
        fi;
    od;
    return true;
end;

GetNrOrbitIsotropic:=function(eGram, kDim, choice)
    local TmpDir, eProg, FileM, FileOut, FileErr, TheCommand, U;
    Print("Running GetNrOrbitIsotropic with INDEF_FORM_GetOrbit_IsotropicKplane\n");
    U:=INDEF_FORM_GetOrbitIsotropic(eGram, kDim, choice, IndefOptions);
    if U="program failure" then
        Print("The output file is not existing. That qualifies as a fail\n");
        return rec(is_correct:=false);
    fi;
    return rec(is_correct:=true, nr_orbit:=Length(U));
end;

TestOrbitIsotropic:=function(eList, kDim)
    local eGram1, eP, eGram2, Xnorm, choice, res1, res2, nr_orbit1, nr_orbit2;
    eGram1:=GetGramMatrixFromList(eList);
    eP:=RandomIntegralUnimodularMatrix(Length(eGram1));
    eGram2:=eP * eGram1 * TransposedMat(eP);
    Xnorm:=2;
    for choice in ["plane", "flag"]
    do
        res1:=GetNrOrbitIsotropic(eGram1, kDim, choice);
        if res1.is_correct=false then
            Print("TestOrbitIsotropic: Failed at res1.is_correct\n");
            return false;
        fi;
        nr_orbit1:=res1.nr_orbit;
        #
        res2:=GetNrOrbitIsotropic(eGram2, kDim, choice);
        if res2.is_correct=false then
            Print("TestOrbitIsotropic: Failed at res2.is_correct\n");
            return false;
        fi;
        nr_orbit2:=res2.nr_orbit;
        #
        if nr_orbit1<>nr_orbit2 then
            Print("TestOrbitIsotropic: Failed at nr_orbit1<>nr_orbit2\n");
            return false;
        fi;
    od;
    return true;
end;

# The number of positive eigenvalues of a nondegenerate symmetric matrix:
# its characteristic polynomial has only real roots, so Descartes' rule of
# signs counts them exactly.
NumberPositiveEigenvalues:=function(eGram)
    local coeffs, signs;
    coeffs:=CoefficientsOfUnivariatePolynomial(CharacteristicPolynomial(eGram));
    signs:=Filtered(List(coeffs, SignInt), x -> x <> 0);
    return Number([1..Length(signs)-1], i -> signs[i] <> signs[i+1]);
end;

# Q and -Q are equivalent when the signature (p,q) is balanced, p = q, for
# the lattices tested here (U + 2U is isometric to its negative), and never
# when p <> q.
TestEquivNegative:=function(eList)
    local eGram, U, p, n;
    Print("Running TestEquivNegative with INDEF_FORM_TestEquivalence\n");
    eGram:=GetGramMatrixFromList(eList);
    n:=Length(eGram);
    p:=NumberPositiveEigenvalues(eGram);
    U:=INDEF_FORM_TestEquivalence(eGram, -eGram, IndefOptions);
    if U="program failure" then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    if 2 * p <> n then
        if U <> fail then
            Print("TestEquivNegative: Q and -Q found equivalent, signature (", p, ",", n - p, ")\n");
            return false;
        fi;
        return true;
    fi;
    if U = fail then
        Print("TestEquivNegative: Q and -Q not found equivalent, signature (", p, ",", n - p, ")\n");
        return false;
    fi;
    if U * eGram * TransposedMat(U) <> -eGram then
        Print("TestEquivNegative: Failed at U * eGram * TransposedMat(U) <> -eGram\n");
        return false;
    fi;
    return true;
end;

# The orbits of vectors of norm Xnorm for a definite form and for its
# negative (norm -Xnorm), compared with the known number.
TestDefiniteOrbits:=function(eRec)
    local eGram, eSign, res;
    eGram:=GetGramMatrixFromList(eRec.eList);
    for eSign in [1, -1]
    do
        res:=GetNrOrbitRepresentative(eSign * eGram, eSign * eRec.Xnorm);
        if res.is_correct=false then
            Print("TestDefiniteOrbits: Failed at res.is_correct\n");
            return false;
        fi;
        if res.nr_orbit <> eRec.nr_orbit then
            Print("TestDefiniteOrbits: eList=", eRec.eList, " sign=", eSign, " nr_orbit=", res.nr_orbit, " but expected ", eRec.nr_orbit, "\n");
            return false;
        fi;
    od;
    return true;
end;

FullTest:=function()
    local eRec, ListRec, eList, k, result;
    ListRec:=[];
    Add(ListRec, rec(eList:=["U", "2U"], k:=2));
    Add(ListRec, rec(eList:=["U", "2U", "A2"], k:=2));
    Add(ListRec, rec(eList:=["U", "2U", "A3"], k:=2));
    Add(ListRec, rec(eList:=["U", "2U", "A2", "A2"], k:=2));
    Add(ListRec, rec(eList:=["U", "U", "E7"], k:=2));
    Add(ListRec, rec(eList:=["U", "2U", "2E8"], k:=2)); # Enriques, should work as we did in the paper
    for eRec in ListRec
    do
        Print("Working with eRec=", eRec, "\n");
        eList:=eRec.eList;
        k:=eRec.k;
        result:=TestStab(eList);
        if result<>true then
            Print("Failed the TestStab\n");
            return result;
        fi;
        result:=TestEquiv(eList);
        if result<>true then
            Print("Failed the TestEquiv\n");
            return result;
        fi;
        result:=TestOrbitRepresentative(eList);
        if result<>true then
            Print("Failed the TestOrbitRepresentative\n");
            return result;
        fi;
        if TestOrbitIsotropic(eList, k)=false then
            return false;
        fi;
        if TestEquivNegative(eList)=false then
            return false;
        fi;
    od;
    # Definite forms: the vectors of norm Xnorm form nr_orbit orbits.
    ListRec:=[];
    Add(ListRec, rec(eList:=["A2"], Xnorm:=2, nr_orbit:=1));
    Add(ListRec, rec(eList:=["A1", "A1"], Xnorm:=2, nr_orbit:=1));
    Add(ListRec, rec(eList:=["A1", "A2"], Xnorm:=2, nr_orbit:=2));
    Add(ListRec, rec(eList:=["D4"], Xnorm:=2, nr_orbit:=1));
    Add(ListRec, rec(eList:=["E8"], Xnorm:=2, nr_orbit:=1));
    for eRec in ListRec
    do
        Print("Working with eRec=", eRec, "\n");
        if TestDefiniteOrbits(eRec)=false then
            return false;
        fi;
        result:=TestEquiv(eRec.eList);
        if result<>true then
            Print("Failed the TestEquiv\n");
            return result;
        fi;
        if TestEquivNegative(eRec.eList)=false then
            return false;
        fi;
    od;
    return true;
end;

result:=FullTest();
CI_Decision_Reset();
if result<>true then
    # Error case
    Print("Error case\n");
else
    # No error case
    Print("Normal case\n");
    CI_Write_Ok();
fi;
