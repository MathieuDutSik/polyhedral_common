Read("../common.g");
Read("../access_points.g");
Print("Beginning ComputeCoveringDensity\n");


TestComputation:=function(eRec)
    local TmpDir, FileIn, FileNml, FileOut, output, strOut, eProg, TheCommand, U, is_correct, UseMpi;
    result:=get_lattice_covering(eRec.eMat);
    if is_error(result) then
        return false;
    fi;
    Print("TheCov=", U.TheCov, "  det=", U.TheDet, "  CovDensity=", U.CovDensity, "\n");
    return true;
end;

eDir:="PerfectMatrices";
ListFiles:=ListFileDirectory(eDir);
ListRec:=[];
for eFile in ListFiles
do
    fFile:=Concatenation(eDir, "/", eFile);
    eMat:=ReadMatrixFile(fFile);
    Add(ListRec, rec(eMat:=eMat, eFile:=eFile));
od;


FullTest:=function()
    local n_error, iRec, eRec, RecReply;
    n_error:=0;
    iRec:=0;
    for eRec in ListRec
    do
        Print("iRec=", iRec, " / ", Length(ListRec), "\n");
        reply:=TestComputation(eRec);
        Print("reply(B)=", reply, "\n");
        if reply=false then
            n_error:=n_error+1;
            Print("1: n_error=", n_error, "\n");
            return n_error;
        fi;
        iRec:=iRec + 1;
    od;
    Print("FullTest: n_error=", n_error, "\n");
    return n_error;
end;

NestFunction:=function()
    local n_error;
    n_error:=FullTest();
    Print("2: n_error=", n_error, "\n");
    CI_Decision_Reset();
    if n_error > 0 then
        # Error case
        Print("Error case\n");
    else
        # No error case
        Print("Normal case\n");
        CI_Write_Ok();
    fi;
    CI_PrintExistConclusion();
end;

NestFunction();


