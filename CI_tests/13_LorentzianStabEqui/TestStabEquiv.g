Read("../common.g");
Print("Beginning Test enumeration of Lorentzian indefinite equivalence\n");

Equi_SingleTest:=function(eRecEqui)
    local mat1, mat2, test, DirTemp, FileMat1, FileMat2, FileTest, eProg, TheCommand, U;
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
    DirTemp:=DirectoryTemporary();
    FileMat1:=Filename(DirTemp, "Indef_Equi.mat1");
    FileMat2:=Filename(DirTemp, "Indef_Equi.mat2");
    FileTest:=Filename(DirTemp, "Indef_Equi.test");
    #
    WriteMatrixFile(FileMat1, mat1);
    WriteMatrixFile(FileMat2, mat2);
    #
    eProg:="../../src_lorentzian/LORENTZ_PERF_Isomorphism";
    TheCommand:=Concatenation(eProg, " ", FileMat1, " ", FileMat2, " GAP ", FileTest);
    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    #
    if IsExistingFile(FileTest)=false then
        Print("The output file is not existing. That qualifies as a fail for Equi_SingleTest\n");
        return false;
    fi;
    U:=ReadAsFunction(FileTest)();
    RemoveFile(FileMat1);
    RemoveFile(FileMat2);
    RemoveFile(FileTest);
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
    #
    DirTemp:=DirectoryTemporary();
    FileMat:=Filename(DirTemp, "Indef_Stab.mat");
    FileOut:=Filename(DirTemp, "Indef_Stab.out");
    RemoveFileIfExist(FileMat);
    RemoveFileIfExist(FileOut);
    #
    WriteMatrixFile(FileMat, mat);
    #
    eProg:="../../src_lorentzian/LORENTZ_PERF_Automorphism";
    TheCommand:=Concatenation(eProg, " ", FileMat, " GAP ", FileOut);
    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    #
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail for Stab_SingleTest\n");
        return false;
    fi;
    U:=ReadAsFunction(FileOut)();
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

