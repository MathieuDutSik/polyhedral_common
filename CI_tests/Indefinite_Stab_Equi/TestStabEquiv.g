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
    RemoveFileIfExist(FileMat1);
    RemoveFileIfExist(FileMat2);
    RemoveFileIfExist(FileTest);
    #
    WriteMatrixFile(FileMat1, mat1);
    WriteMatrixFile(FileMat2, mat2);
    #
    eProg:="../../src_lorentzian/LORENTZ_PERF_Isomorphism";
    TheCommand:=Concatenation(eProg, " ", FileMat1, " ", FileMat2, " GAP ", FileTest);
    Exec(TheCommand);
    #
    if IsExistingFile(FileTest)=false then
        Print("The output file is not existing. That qualifies as a fail for Equi_SingleTest\n");
        return rec(is_correct:=false);
    fi;
    U:=ReadAsFunction(FileTest)();
    RemoveFileIfExist(FileMat1);
    RemoveFileIfExist(FileMat2);
    RemoveFileIfExist(FileTest);
    if U=fail then
        if test<>fail then
            Print("Different equivalence conclusions\n");
            return rec(is_correct:=false);
        fi;
    fi;
    if mat2 <> U * mat1 * TransposedMat(U) then
        Print("Different equivalence conclusions\n");
        return rec(is_correct:=false);
    fi;
    #
    return rec(is_correct:=true);
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
    Exec(TheCommand);
    #
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail for Stab_SingleTest\n");
        return rec(is_correct:=false);
    fi;
    U:=ReadAsFunction(FileOut)();
    for eMat in U
    do
        if mat <> eMat * mat * TransposedMat(eMat) then
            Print("eMat does not stabilize the space\n");
            return rec(is_correct:=false);
        fi;
    od;
    #
    # Would be nice to have some consistency checks here
    #
    return rec(is_correct:=true);
end;



Equi_AllTests:=function()
    local TheDir, ListEquiFiles, n_error, iRec, eFile, FullFile, eRec, RecReply;
    TheDir:="TestCasesEqui";
    ListEquiFiles:=ListFileDirectory(TheDir);
    n_error:=0;
    iRec:=0;
    for eFile in ListEquiFiles
    do
        FullFile:=Concatenation(TheDir, "/", eFile);
        eRec:=ReadAsFunction(FullFile)();
        Print("iRec=", iRec, " / ", Length(ListEquiFiles), "\n");
        RecReply:=Equi_SingleTest(eRec);
        if RecReply.is_correct=false then
            n_error:=n_error+1;
            return n_error;
        fi;
        iRec:=iRec + 1;
    od;
    return n_error;
end;

Stab_AllTests:=function()
    local TheDir, ListStabFiles, n_error, iRec, eFile, FullFile, eRec, RecReply;
    TheDir:="TestCasesStab";
    ListStabFiles:=ListFileDirectory(TheDir);
    n_error:=0;
    iRec:=0;
    for eFile in ListStabFiles
    do
        FullFile:=Concatenation(TheDir, "/", eFile);
        eRec:=ReadAsFunction(FullFile)();
        Print("iRec=", iRec, " / ", Length(ListStabFiles), "\n");
        RecReply:=Stab_SingleTest(eRec);
        if RecReply.is_correct=false then
            n_error:=n_error+1;
            return n_error;
        fi;
        iRec:=iRec + 1;
    od;
    return n_error;
end;

AllTests:=function()
    return Equi_AllTests() + Stab_AllTests();
end;


n_error:=AllTests();
if n_error > 0 then
    # Error case
    GAP_EXIT_CODE(1);
else
    # No error case
    GAP_EXIT_CODE(0);
fi;

