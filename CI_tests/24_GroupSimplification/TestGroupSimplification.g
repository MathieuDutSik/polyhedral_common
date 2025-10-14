Read("../common.g");
Print("Beginning TestGroupSimplification\n");

eDir:="BigExamples";

PrintComplexity:=function(ListMatrix, n_space)
    local n_matrix, max_coeff, sum_coeff, eMat, eLine, eVal, absVal, i;
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
    for i in [1..n_space]
    do
        Print(" ");
    od;
    Print("n_matrix=", n_matrix, " max_coeff=", max_coeff, " sum_coeff=", sum_coeff, "\n");
end;



TestSimplification:=function(eProg, eFile)
    local TmpDir, FileOut, fFile, fProg, TheCommand, ListMatrix;

    TmpDir:=DirectoryTemporary();
    FileOut:=Filename(TmpDir, "Test.out");
    #
    fProg:=Concatenation("../../src_group/", eProg);
    TheCommand:=Concatenation(fProg, " mpz_class ", eFile, " GAP ", FileOut);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return rec(result:=false);
    fi;
    ListMatrix:=ReadAsFunction(FileOut)();
    RemoveFile(FileOut);
    return rec(result:=true, ListMatrix:=ListMatrix);
end;

TestFamilySimplification:=function(ListProg, eFile)
    local TheInput, TmpDir, FileIn, TotalListMatrix, eProg, rec_result;
    TheInput:=ReadListMatrixFile(eFile);
    PrintComplexity(TheInput, 6);
    TmpDir:=DirectoryTemporary();
    FileIn:=Filename(TmpDir, "TotalListMatrix.in");
    TotalListMatrix:=[];
    for eProg in ListProg
    do
        Print("         ------------------\n");
        rec_result:=TestSimplification(eProg, eFile);
        if rec_result.result=false then
            return false;
        fi;
        PrintComplexity(rec_result.ListMatrix, 6);
        Append(TotalListMatrix, rec_result.ListMatrix);
    od;
    Print("         = - = - = - = - = - =\n");
    eProg:=ListProg[1];
    WriteListMatrixFile(FileIn, TotalListMatrix);
    rec_result:=TestSimplification(eProg, FileIn);
    if rec_result.result=false then
        return false;
    fi;
    PrintComplexity(rec_result.ListMatrix, 6);
    return true;
end;





FullTest:=function(eDir)
    local ListProg, ListFiles, n_error, eFile, fFile, TheInput, eProg, test;
    ListProg:=["GRP_MatrixGroupSimplification", "GRP_MatrixGroupSimplificationOnline", "GRP_MatrixGroupSimplificationOnlineOpt"];
    ListFiles:=ListFileDirectory(eDir);

    n_error:=0;
    for eFile in ListFiles
    do
        Print("================================\n");
        Print("eFile=", eFile, "\n");
        fFile:=Concatenation(eDir, "/", eFile);
        test:=TestFamilySimplification(ListProg, fFile);
        if test=false then
            n_error:=n_error + 1;
        fi;
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

