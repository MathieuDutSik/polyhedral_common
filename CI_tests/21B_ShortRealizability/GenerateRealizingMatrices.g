Read("../common.g");

ListSHV:=ReadAsFunction("ListSHV_n10_rnk10")();

get_mat:=function(SHV)
    local TmpDir, FileIn, FileOut, eProg, TheCommand, eRec;
    TmpDir:=DirectoryTemporary();
    FileIn:=Filename(TmpDir, "Test.in");
    FileOut:=Filename(TmpDir, "Test.out");
    #
    WriteMatrixFile(FileIn, SHV);
    eProg:="../../src_short/SHORT_TestRealizability";
    TheCommand:=Concatenation(eProg, " gmp ", FileIn, " GAP ", FileOut);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    #
    eRec:=ReadAsFunction(FileOut)();
    if eRec.realizable=false then
        return fail;
    fi;
    RemoveFile(FileIn);
    RemoveFile(FileOut);
    return eRec.matrix;
end;


generate_examples:=function()
    local ListMat, idx, eSHV, eMat, FileMat;
    ListMat:=[];
    idx:=0;
    for eSHV in ListSHV
    do
        idx:=idx+1;
        Print("Before get_mat at idx=", idx, "\n");
        Print("eSHV=\n");
        PrintArray(eSHV);
        eMat:=get_mat(eSHV);
        if eMat=fail then
            Print("eSHV=\n");
            PrintArray(eSHV);
            Print("eSHV is not realizable. Returning fail\n");
            return false;
        fi;
        Add(ListMat, eMat);
    od;
    FileMat:="ListGram_n10_rnk10";
    SaveDataToFile(FileMat, ListMat);
    return true;
end;

test:=generate_examples();


CI_Decision_Reset();
if test=false then
    # Error case
    Print("Error case\n");
else
    # No error case
    Print("Normal case\n");
    CI_Write_Ok();
fi;
