Read("../common.g");
Print("Beginning Test enumeration of Lorentzian indefinite equivalence\n");

Equi_SingleTest:=function(eRecEqui)
    local mat1, mat2, test, FileMat1, FileMat2, FileTest, eProg, TheCommand, U;
    #
    mat1:=eRecEqui.mat1;
    mat2:=eRecEqui.mat2;
    test:=eRecEqui.test;
    if test<>fail then
        if mat2 <> test * mat1 * TransposedMat(test) then
            Error("test is not an equivalence");
        fi;
    fi;
    #
    FileMat1:=Filename(DirectoryTemporary(), "LorentzEquiv.mat1");
    FileMat2:=Filename(DirectoryTemporary(), "LorentzEquiv.mat2");
    FileTest:=Filename(DirectoryTemporary(), "LorentzEquiv.test");
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
    if IsExistingFile(FileResult)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return rec(is_correct:=false);
    fi;
    U:=ReadAsFunction(FileTest)();
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




Stab_SingleTest:=function(eRecEqui)
    local mat1, mat2, test, FileMat1, FileMat2, FileTest, eProg, TheCommand, U;
    #
    mat1:=eRecEqui.mat1;
    mat2:=eRecEqui.mat2;
    test:=eRecEqui.test;
    if test<>fail then
        if mat2 <> test * mat1 * TransposedMat(test) then
            Error("test is not an equivalence");
        fi;
    fi;
    #
    FileMat1:=Filename(DirectoryTemporary(), "LorentzEquiv.mat1");
    FileMat2:=Filename(DirectoryTemporary(), "LorentzEquiv.mat2");
    FileTest:=Filename(DirectoryTemporary(), "LorentzEquiv.test");
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
    if IsExistingFile(FileResult)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return rec(is_correct:=false);
    fi;
    U:=ReadAsFunction(FileTest)();
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






Equi_AllTests:=function()
    local TheDir, ListEquiFiles, n_error, iRec, eRec, RecReply;
    TheDir:="TestCasesEqui";
    ListEquiFiles:=ListFileDirectory(TheDir);
    n_error:=0;
    iRec:=0;
    for eFile in ListEquiFiles
    do
        Print("iRec=", iRec, " / ", Length(ListRec), "\n");
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
    local TheDir, ListEquiFiles, n_error, iRec, eRec, RecReply;
    TheDir:="TestCasesStab";
    ListStabFiles:=ListFileDirectory(TheDir);
    n_error:=0;
    iRec:=0;
    for eFile in ListStabFiles
    do
        Print("iRec=", iRec, " / ", Length(ListRec), "\n");
        RecReply:=Stab_SingleTest(eRec);
        if RecReply.is_correct=false then
            n_error:=n_error+1;
            return n_error;
        fi;
        iRec:=iRec + 1;
    od;
    return n_error;
end;




n_error:=FullTest();
if n_error > 0 then
    # Error case
    GAP_EXIT_CODE(1);
else
    # No error case
    GAP_EXIT_CODE(0);
fi;

