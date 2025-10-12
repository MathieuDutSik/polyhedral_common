Read("../common.g");
Print("Beginning TestGroupSimplification\n");

eDir:="BigExamples";

PrintComplexity:=function(ListMatrix)
    local n_matrix, max_coeff, sum_coeff, eMat, eLine, eVal, absVal;
    n_matrix:=Length(ListMatrix);
    max_coeff:=0;
    sum_coeff:=0;
    for eMat in ListMatrix
    do
        for eLine in eMat
        do
            for eVal in eLine
            do
                absVal:=eVal;
                if eVal < 0 then
                    absVal:=-eVal;
                fi;
                if absVal > max_coeff then
                    max_coeff:=absVal;
                fi;
                sum_coeff:=sum_coeff + absVal;
            od;
        od;
    od;
    Print("    n_matrix=", n_matrix, " max_coeff=", max_coeff, " sum_coeff=", sum_coeff, "\n");
end;



TestSimplification:=function(eProg, eFile)
    local TmpDir, FileOut, FileErr, fFile, fProg, TheCommand, ListMatrix;

    TmpDir:=DirectoryTemporary();
    FileOut:=Filename(TmpDir, "Test.out");
    FileErr:=Filename(TmpDir, "Test.err");
    #
    fFile:=Concatenation(eDir, "/", eFile);
    fProg:=Concatenation("../../src_group/", eProg);
    TheCommand:=Concatenation(fProg, " mpz_class ", fFile, " GAP ", FileOut, " 2> ", FileErr);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    ListMatrix:=ReadAsFunction(FileOut)();
    PrintComplexity(ListMatrix);
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
        Print("eFile=", eFile, "\n");
        for eProg in ListProg
        do
            Print("  eProg=", eProg, "\n");
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

