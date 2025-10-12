Read("../common.g");
Print("Beginning TestGroupSimplification\n");

eDir:="BigExamples";

TestSimplification:=function(eProg, eFile)
    local FileOut, fFile, fProg, TheCommand, ListMatrix;

    FileOut:=Filename(DirectoryTemporary(), "Test.out");
    #
    fFile:=Concatenation(eDir, "/", eFile);
    fProg:=Concatenation("../../src_group/", eProg);
    TheCommand:=Concatenation(fProg, " mpz_class ", fFile, " GAP ", FileOut);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    ListMatrix:=ReadAsFunction(FileOut)();
    RemoveFile(FileOut);
    return true;
end;

FullTest:=function(eDir)
    local ListProg, ListFiles, n_error, eFile, eProg, test;
    ListProg:=["GRP_MatrixGroupSimplification", "GRP_MatrixGroupSimplificationOnline", "GRP_MatrixGroupSimplificationOnlineOpt"];
    ListFiles:=ListFileDirectory(eDir);

    n_error:=0;
    for eFile in ListFiles
    do
        for eProg in ListProg
        do
            test:=TestSimplification(eProg, eFile);
            if test=false then
                n_error:=n_error + 1;
            fi;
        od;
    od;
    if n_error>0 then
        return false;
    else
        return true;
    fi;
end;

test:=FullTest(eDir);
Print("test=", test, "\n");

CI_Decision_Reset();
if test=false then
    # Error case
    Print("Error case\n");
else
    # No error case
    Print("Normal case\n");
    CI_Write_Ok();
fi;

