Read("../common.g");
Read("../access_points.g");

ListSHV:=ReadAsFunction("ListSHV_n10_rnk10")();

get_mat:=function(SHV)
    local TmpDir, FileIn, FileOut, eProg, TheCommand, eRec;
    eRec:=test_shortest_realizability(SHV);
    if is_error(eRec) then
        return false;
    fi;
    if eRec.realizable=false then
        return false;
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
        Print("Before get_mat at idx=", idx, " det=", DeterminantMat(eSHV), "\n");
        Print("eSHV=\n");
        PrintArray(eSHV);
        eMat:=test_shortest_realizability(eSHV);
        if eMat=false then
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
