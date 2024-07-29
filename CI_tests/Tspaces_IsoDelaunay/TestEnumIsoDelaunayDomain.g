Read("../common.g");
Print("Beginning Test enumeration of iso-Delaunay domains\n");

TestEnumeration:=function(eRec)
    local n, FileLinSpa, FileNml, FileResult, output, eProg, TheCommand, U, is_correct;
    #
    FileLinSpa:="TSPACE_LinSpa";
    FileNml:="Enum_Tspaces_CI.nml";
    FileResult:="Result";
    TheCommand:=Concatenation("cp ", eRec.FileLinSpa, " ", FileLinSpa);
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
    U:=ReadAsFunction(FileOut)();
    RemoveFile(FileLinSpa);
    RemoveFile(FileResult);
    is_correct:=eRec.n_domain = U.nb;
    Print("n=", n, " n_domain=", eRec.n_domain, " is_correct=", is_correct, "\n");
    return rec(is_correct:=is_correct);
end;

ListRec:=[];
for eFile in ["ListCases_Bravais", "ListCases_Coxeter"]
do
    Append(ListRec, ReadAsFunction(eFile)());
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
if n_error > 0 then
    # Error case
    Print("Error case\n");
    GAP_EXIT_CODE(1);
else
    # No error case
    Print("Normal case\n");
    GAP_EXIT_CODE(0);
fi;

