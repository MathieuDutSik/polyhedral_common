Read("../common.g");
Print("Beginning of AllTests\n");

TestStab:=function(eList)
    local eProg, eGram, FileIn, FileOut, TheCommand, U, eP;
    eProg:="../../src_indefinite/INDEF_FORM_AutomorphismGroup";
    eGram:=GetGramMatrixFromList(eList);
    FileIn:=Filename(DirectoryTemporary(), "Test.in");
    FileOut:=Filename(DirectoryTemporary(), "Test.out");
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
            return false;
        fi;
    od;
    return true;
end;

TestEquiv:=function(eList)
    local eProg, eGram1, eP, eGram2, FileM1, FileM2, FileOut, TheCommand, U;
    eProg:="../../src_indefinite/INDEF_FORM_TestEquivalence";
    eGram1:=GetGramMatrixFromList(eList);
    eP:=RandomIntegralUnimodularMatrix(Length(eGram1));
    eGram2:=eP * eGram1 * TransposedMat(eP);
    FileM1:=Filename(DirectoryTemporary(), "Mat1.in");
    FileM2:=Filename(DirectoryTemporary(), "Mat2.in");
    FileOut:=Filename(DirectoryTemporary(), "Test.out");
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
        return false;
    fi;
    if U * eGram1 * TransposedMat(U) <> eGram2 then
        return false;
    fi;
    return true;
end;

GetNrOrbitRepresentative:=function(eGram, Xnorm)
    local eProg, FileM, FileOut, TheCommand, U;
    eProg:="../../src_indefinite/INDEF_FORM_GetOrbitRepresentative";
    FileM:=Filename(DirectoryTemporary(), "Mat.in");
    FileOut:=Filename(DirectoryTemporary(), "Test.out");
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
            return false;
        fi;
        nr_orbit1:=res1.nr_orbit;
        #
        res2:=GetNrOrbitRepresentative(eGram2, Xnorm);
        if res2.is_correct=false then
            return false;
        fi;
        nr_orbit2:=res2.nr_orbit;
        #
        if nr_orbit1<>nr_orbit2 then
            return false;
        fi;
    od;
    return true;
end;

GetNrOrbitIsotropic:=function(eGram, kDim, choice)
    local eProg, FileM, FileOut, TheCommand, U;
    eProg:="../../src_indefinite/INDEF_FORM_GetOrbit_IsotropicKplane";
    FileM:=Filename(DirectoryTemporary(), "Mat.in");
    FileOut:=Filename(DirectoryTemporary(), "Test.out");
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

TestOrbitIsotropic:=function(eList, k)
    local eGram1, eP, eGram2, Xnorm, choice, res1, res2, nr_orbit1, nr_orbit2;
    eGram1:=GetGramMatrixFromList(eList);
    eP:=RandomIntegralUnimodularMatrix(Length(eGram1));
    eGram2:=eP * eGram1 * TransposedMat(eP);
    Xnorm:=2;
    for choice in ["plane", "flag"]
    do
        res1:=GetNrOrbitIsotropic(eGram1, Xnorm);
        if res1.is_correct=false then
            return false;
        fi;
        nr_orbit1:=res1.nr_orbit;
        #
        res2:=GetNrOrbitIsotropic(eGram2, Xnorm);
        if res2.is_correct=false then
            return false;
        fi;
        nr_orbit2:=res2.nr_orbit;
        #
        if nr_orbit1<>nr_orbit2 then
            return false;
        fi;
    od;
    return true;
end;

FullTest:=function()
    local rec1, rec2, rec3, rec4, rec5, rec6, ListL, ListRk2, eList;
    ListL:=[];
    Add(ListL, ["U", "2U"]);
    Add(ListL, ["U", "2U", "A2"]);
    Add(ListL, ["U", "2U", "A3"]);
    Add(ListL, ["U", "2U", "A2", "A2"]);
    Add(ListL, ["U", "U", "E7"]);
    Add(ListL, ["U", "2U", "2E8"]); # Enriques, should work as we did in the paper
    for eList in ListL
    do
        if TestStab(eList)=false then
            return false;
        fi;
        if TestEquiv(eList)=false then
            return false;
        fi;
        if TestOrbitRepresentative(eList)=false then
            return false;
        fi;
    od;
    ListRk2:=[];
#    Add(ListL, ["U", "2U"]);
#    Add(ListL, ["U", "2U", "A2"]);
#    Add(ListL, ["U", "2U", "A3"]);
#    Add(ListL, ["U", "2U", "A2", "A2"]); # Too slow apparently.
    Add(ListL, ["U", "U", "E7"]);
#    Add(ListL, ["U", "2U", "2E8"]); # Enriques, should work as we did in the paper
    for eList in ListRk2
    do
        if TestOrbitIsotropic(eList)=false then
            return false;
        fi;
    od;
    return true;
end;

result:=FullTest();
CI_Decision_Reset();
if result=false then
    # Error case
    Print("Error case\n");
else
    # No error case
    Print("Normal case\n");
    CI_Write_Ok();
fi;
