Read("../common.g");
Print("Beginning Test enumeration of iso-Delaunay domains\n");

TestEnumeration:=function(eRec)
    local FileLinSpa, FileNml, FileResult, output, eProg, TheCommand, U, is_correct;
    #
    FileLinSpa:="TSPACE_LinSpa";
    FileNml:="Enum_Tspaces_CI.nml";
    FileResult:="Result";
    Print("eRec=", eRec, "\n");
    TheCommand:=Concatenation("cp ", eRec.eFile, " ", FileLinSpa);
    Exec(TheCommand);
    #
    eProg:="../../src_latt/LATT_MPI_Lattice_IsoDelaunayDomain";
    TheCommand:=Concatenation(eProg, " ", FileNml);
    Exec(TheCommand);
    #
    if IsExistingFile(FileResult)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return rec(is_correct:=false);
    fi;
    U:=ReadAsFunction(FileResult)();
    RemoveFile(FileLinSpa);
    RemoveFile(FileResult);
    is_correct:=eRec.nb = U.nb;
    Print("eRec.nb=", eRec.nb, " U.nb=", U.nb, " is_correct=", is_correct, "\n");
    return rec(is_correct:=is_correct);
end;

ListRec:=[];
for eData in [rec(prefix:="TSPACES_Bravais", FileS:="ListCases_Bravais"),
              rec(prefix:="TSPACES_Coxeter", FileS:="ListCases_Coxeter")]
do
    Lst:=ReadAsFunction(eData.FileS)();
    for elst in Lst
    do
        eFile:=Concatenation(eData.prefix, "/", elst.eFile);
        eRec:=rec(eFile:=eFile, nb:=elst.nb);
        Add(ListRec, eRec);
    od;
od;

FullTest:=function()
    local n_error, iRec, eRec, RecReply;
    n_error:=0;
    iRec:=0;
    for eRec in ListRec
    do
        Print("iRec=", iRec, " / ", Length(ListRec), "\n");
        RecReply:=TestEnumeration(eRec);
        if RecReply.is_correct=false then
            n_error:=n_error+1;
            return n_error;
        fi;
        iRec:=iRec + 1;
    od;
    return n_error;
end;

n_error:=FullTest();
Print("n_error=", n_error, "\n");

CI_Decision_Reset();
if n_error > 0 then
    Print("Error case\n");
else
    Print("Normal case\n");
    CI_Write_Ok();
fi;
