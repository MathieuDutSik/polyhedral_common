Read("../common.g");
Print("Beginning of AllTests\n");

TestStab:=function(eList)
    local eProg, eGram, TmpDir, FileIn, FileOut, TheCommand, U, eP;
    eProg:="../../src_indefinite/INDEF_FORM_AutomorphismGroup";
    eGram:=GetGramMatrixFromList(eList);
    TmpDir:=DirectoryTemporary();
    FileIn:=Filename(TmpDir, "Test.in");
    FileOut:=Filename(TmpDir, "Test.out");
    WriteMatrixFile(FileIn, eGram);
    TheCommand:=Concatenation(eProg, " gmp ", FileIn, " GAP ", FileOut);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    U:=ReadAsFunction(FileOut)();
    RemoveFile(FileIn);
    RemoveFile(FileOut);
    for eP in U
    do
        if eP * eGram * TransposedMat(eP) <> eGram then
            Print("TestEquiv: Failed at U * eGram * TransposedMat(U) <> eGram\n");
            return false;
        fi;
    od;
    return true;
end;

TestEquiv:=function(eList)
    local TmpDir, eProg, eGram1, eP, eGram2, FileM1, FileM2, FileOut, TheCommand, U;
    TmpDir:=DirectoryTemporary();
    eProg:="../../src_indefinite/INDEF_FORM_TestEquivalence";
    eGram1:=GetGramMatrixFromList(eList);
    eP:=RandomIntegralUnimodularMatrix(Length(eGram1));
    eGram2:=eP * eGram1 * TransposedMat(eP);
    FileM1:=Filename(TmpDir, "Mat1.in");
    FileM2:=Filename(TmpDir, "Mat2.in");
    FileOut:=Filename(TmpDir, "Test.out");
    WriteMatrixFile(FileM1, eGram1);
    WriteMatrixFile(FileM2, eGram2);
    TheCommand:=Concatenation(eProg, " gmp ", FileM1, " ", FileM2, " GAP ", FileOut);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    U:=ReadAsFunction(FileOut)();
    RemoveFile(FileM1);
    RemoveFile(FileM2);
    RemoveFile(FileOut);
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
    local TmpDir, eProg, FileM, FileOut, TheCommand, U;
    TmpDir:=DirectoryTemporary();
    eProg:="../../src_indefinite/INDEF_FORM_GetOrbitRepresentative";
    FileM:=Filename(TmpDir, "Mat.in");
    FileOut:=Filename(TmpDir, "Test.out");
    WriteMatrixFile(FileM, eGram);
    TheCommand:=Concatenation(eProg, " gmp ", FileM, " ", String(Xnorm), " GAP ", FileOut);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return rec(is_correct:=false);
    fi;
    U:=ReadAsFunction(FileOut)();
    RemoveFile(FileM);
    RemoveFile(FileOut);
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
    local TmpDir, eProg, FileM, FileOut, TheCommand, U;
    TmpDir:=DirectoryTemporary();
    eProg:="../../src_indefinite/INDEF_FORM_GetOrbit_IsotropicKplane";
    FileM:=Filename(TmpDir, "Mat.in");
    FileOut:=Filename(TmpDir, "Test.out");
    WriteMatrixFile(FileM, eGram);
    TheCommand:=Concatenation(eProg, " gmp ", FileM, " ", String(kDim), " ", choice, " GAP ", FileOut);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return rec(is_correct:=false);
    fi;
    U:=ReadAsFunction(FileOut)();
    RemoveFile(FileM);
    RemoveFile(FileOut);
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
