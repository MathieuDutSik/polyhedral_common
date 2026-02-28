Read("../common.g");
Read("../access_points.g");
Print("Beginning Test for testing finiteness of matrix groups\n");

Finiteness_SingleTest:=function(eRecFin)
    local DirTemp, FileMat, FileTest, eProg, TheCommand, U;
    #
    U:=test_matrix_group_finiteness(eRecFin.GRPmatr);
    if is_error(U) then
        return false;
    fi;
    if U.is_finite<>eRecFin.is_finite then
        Print("Different finiteness results\n");
        return false;
    fi;
    #
    return true;
end;

Finiteness_AllTests:=function()
    local TheDir, ListFinitenessFiles, n_error, iRec, eFile, FullFile, eRecFin, test;
    TheDir:="Finiteness";
    ListFinitenessFiles:=ListFileDirectory(TheDir);
    n_error:=0;
    iRec:=0;
    for eFile in ListFinitenessFiles
    do
        FullFile:=Concatenation(TheDir, "/", eFile);
        eRecFin:=ReadAsFunction(FullFile)();
        Print("iRec=", iRec, " / ", Length(ListFinitenessFiles), " n_error=", n_error, " finite=", eRecFin.is_finite, "\n");
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

