Read("../common.g");
Read("../access_points.g");
Print("Beginning Test enumeration of Lorentzian indefinite equivalence\n");

Equi_SingleTest:=function(eRecEqui)
    local mat1, mat2, test, U;
    #
    mat1:=eRecEqui.mat1;
    mat2:=eRecEqui.mat2;
    test:=eRecEqui.test;
    if test<>fail then
        if mat2 <> test * mat1 * TransposedMat(test) then
            Error("erroneous input: test is not an equivalence. Debug the debug code");
        fi;
    fi;
    #
    U:=test_equivalent_lorentzian_matrices(mat1, mat2);
    if is_error(U) then
        return false;
    fi;
    if U=fail then
        if test<>fail then
            Print("Different equivalence conclusions\n");
            return false;
        fi;
    fi;
    if mat2 <> U * mat1 * TransposedMat(U) then
        Print("Different equivalence conclusions\n");
        return false;
    fi;
    #
    return true;
end;




Stab_SingleTest:=function(eRecStab)
    local mat, DirTemp, FileMat, FileOut, eProg, TheCommand, U, eMat;
    #
    mat:=eRecStab.mat;
    U:=stabilizer_lorentzian_matrix(mat);
    if is_error(U) then
        return false;
    fi;
    for eMat in U
    do
        if mat <> eMat * mat * TransposedMat(eMat) then
            Print("eMat does not stabilize the space\n");
            return false;
        fi;
    od;
    #
    # Would be nice to have some consistency checks here
    #
    return true;
end;



Equi_AllTests:=function()
    local TheDir, ListEquiFiles, n_error, iRec, eFile, FullFile, eRec, test;
    TheDir:="TestCasesEqui";
    ListEquiFiles:=ListFileDirectory(TheDir);
    n_error:=0;
    iRec:=0;
    for eFile in ListEquiFiles
    do
        FullFile:=Concatenation(TheDir, "/", eFile);
        eRec:=ReadAsFunction(FullFile)();
        Print("iRec=", iRec, " / ", Length(ListEquiFiles), "\n");
        test:=Equi_SingleTest(eRec);
        if test=false then
            n_error:=n_error+1;
            return n_error;
        fi;
        iRec:=iRec + 1;
    od;
    Print("Equi_AllTests, n_error=", n_error, "\n");
    return n_error;
end;

Stab_AllTests:=function()
    local TheDir, ListStabFiles, n_error, iRec, eFile, FullFile, eRec, test;
    TheDir:="TestCasesStab";
    ListStabFiles:=ListFileDirectory(TheDir);
    n_error:=0;
    iRec:=0;
    for eFile in ListStabFiles
    do
        FullFile:=Concatenation(TheDir, "/", eFile);
        eRec:=ReadAsFunction(FullFile)();
        Print("iRec=", iRec, " / ", Length(ListStabFiles), "\n");
        test:=Stab_SingleTest(eRec);
        if test=false then
            n_error:=n_error+1;
            return n_error;
        fi;
        iRec:=iRec + 1;
    od;
    Print("Stab_AllTests, n_error=", n_error, "\n");
    return n_error;
end;

AllTests:=function()
    return Equi_AllTests() + Stab_AllTests();
end;


n_error:=AllTests();
Print("n_error=", n_error, "\n");
CI_Decision_Reset();
if n_error > 0 then
    # Error case
    Print("Error case\n");
else
    # No error case
    Print("Normal case\n");
    CI_Write_Ok();
fi;

