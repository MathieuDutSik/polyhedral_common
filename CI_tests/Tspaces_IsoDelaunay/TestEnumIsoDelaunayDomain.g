Read("../common.g");
Print("Beginning Test enumeration of iso-Delaunay domains\n");

get_nb_domains:=function(eFile)
    local FileLinSpa, FileNml, FileResult, output, eProg, TheCommand, U, is_correct;
    #
    FileLinSpa:="TSPACE_LinSpa";
    FileNml:="Enum_Tspaces_CI.nml";
    FileResult:="Result";
    Print("eFile=", eFile, "\n");
    TheCommand:=Concatenation("cp ", eFile, " ", FileLinSpa);
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
    return U.nb;
end;



TestEnumeration:=function(eRec)
    local comp_nb, is_correct;
    comp_nb:=get_nb_domains(eRec.eFile);
    is_correct:=eRec.nb = comp_nb;
    Print("eRec.nb=", eRec.nb, " comp_nb=", comp_nb, " is_correct=", is_correct, "\n");
    if is_correct=false then
        Print("  FOUND SOME ERROR\n");
    fi;
    return is_correct;
end;

eData1:=rec(prefix:="TSPACES_Bravais", FileS:="ListCases_Bravais");
eData2:=rec(prefix:="TSPACES_Coxeter", FileS:="ListCases_Coxeter");

ListRec:=[];
for eData in [eData1, eData2]
do
    Lst:=ReadAsFunction(eData.FileS)();
    for elst in Lst
    do
        eFile:=Concatenation(eData.prefix, "/", elst.eFile);
        eRec:=rec(eFile:=eFile, nb:=elst.nb);
        Add(ListRec, eRec);
    od;
od;

Recompute:=function(eData)
    local Lst, ListEntries, iter, elst, eFile, comp_nb, eEntry;
    Lst:=ReadAsFunction(eData.FileS)();
    ListEntries:=[];
    iter:=0;
    for elst in Lst
    do
        eFile:=Concatenation(eData.prefix, "/", elst.eFile);
        Print("Recompute iter=", iter, " / ", Length(Lst), " eFile=", elst.eFile, "\n");
        comp_nb:=get_nb_domains(eFile);
        eEntry:=rec(eFile:=elst.eFile, nb:=comp_nb);
        Add(ListEntries, eEntry);
        iter:=iter + 1;
    od;
    SaveDataToFile(eData.FileS, ListEntries);
end;
Recompute(eData1);
Recompute(eData2);




FullTest:=function()
    local n_error, iRec, eRec, RecReply;
    n_error:=0;
    iRec:=0;
    for eRec in ListRec
    do
        Print("iRec=", iRec, " / ", Length(ListRec), " n_error=", n_error, "\n");
        RecReply:=TestEnumeration(eRec);
        if RecReply=false then
            n_error:=n_error+1;
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
