Read("../common.g");
Print("Beginning Test enumeration of iso-Delaunay domains\n");

GetRecInfo:=function(FileLinSpa_input)
    local n, FileLinSpa, FileNml, FileResult, output, eProg, TheCommand, U, is_correct;
    #
    FileLinSpa:="TSPACE_LinSpa";
    FileNml:="Enum_Tspaces_GRPperm.nml";
    FileResult:="Result";
    RemoveFileIfExist(FileLinSpa);
    RemoveFileIfExist(FileResult);
    TheCommand:=Concatenation("cp ", FileLinSpa_input, " ", FileLinSpa);
    Exec(TheCommand);
    #
    eProg:="../../src_latt/LATT_MPI_Lattice_IsoDelaunayDomain";
    TheCommand:=Concatenation(eProg, " ", FileNml);
    Exec(TheCommand);
    #
    if IsExistingFile(FileResult)=false then
        Error("The output file is not existing. That qualifies as a fail");
    fi;
    U:=ReadAsFunction(FileResult)();
    RemoveFile(FileLinSpa);
    RemoveFile(FileResult);
    return U;
end;

DirName:="TSPACES_Bravais";
ListFile:=ListFileDirectory(DirName);

ListRec:=[];
for eFile in ListFile
do
    FullFile:=Concatenation(DirName, "/", eFile);
    Print("eFile=", eFile, "\n");
    U:=GetRecInfo(FullFile);
    nb:=Length(U);
    ListOrder:=List(U, x->Order(x.GRPperm));
    if Maximum(ListOrder) > 1 then
        Error("We found what were looking for, a non-trivial GRPperm");
    fi;
    Print("eFile=", eFile, " nb=", nb, "\n");
    eRec:=rec(eFile:=eFile, nb:=nb);
    Add(ListRec, eRec);
od;
