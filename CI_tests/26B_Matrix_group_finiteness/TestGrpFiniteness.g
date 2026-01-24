Read("../common.g");
Print("Beginning Test for testing finiteness of matrix groups\n");

Finiteness_SingleTest:=function(eRecFin)
    local mat1, mat2, test, DirTemp, FileMat1, FileMat2, FileTest, eProg, TheCommand, U;
    #
    DirTemp:=DirectoryTemporary();
    FileMat:=Filename(DirTemp, "Finiteness_test.input");
    FileTest:=Filename(DirTemp, "Finiteness_test.result");
    #
    WriteMatrixFile(FileMat, ErecFin.GRPmatr);
    #
    eProg:="../../src_latt/GRP_TestFiniteness";
    TheCommand:=Concatenation(eProg, " gmp ", FileMat, " GAP ", FileTest);
    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    #
    if IsExistingFile(FileTest)=false then
        Print("The output file is not existing. That qualifies as a fail for Equi_SingleTest\n");
        return false;
    fi;
    U:=ReadAsFunction(FileTest)();
    RemoveFile(FileMat);
    RemoveFile(FileTest);
    if U.is_finite<>eRecFin.is_finite then
        Print("Different finiteness results\n");
        return false;
    fi;
    #
    return true;
end;

Finiteness_AllTests:=function()
    local TheDir, ListFinitenessFiles, n_error, iRec, eFile, FullFile, eRec, test;
    TheDir:="Finiteness";
    ListFinitenessFiles:=ListFileDirectory(TheDir);
    n_error:=0;
    iRec:=0;
    for eFile in ListFinitenessFiles
    do
        FullFile:=Concatenation(TheDir, "/", eFile);
        eRecFin:=ReadAsFunction(FullFile)();
        Print("iRec=", iRec, " / ", Length(ListStabFiles), "\n");
        test:=Finiteness_SingleTest(eRecFin);
        if test=false then
            n_error:=n_error+1;
            return n_error;
        fi;
        iRec:=iRec + 1;
    od;
    Print("Finiteness_AllTests, n_error=", n_error, "\n");
    return n_error;
end;

n_error:=Finiteness_AllTests();
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

