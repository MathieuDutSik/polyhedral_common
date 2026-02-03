Read("../common.g");
Print("Beginning Test enumeration of iso-Delaunay domains\n");

method:="serial";
#method:="mpi";


get_saturated_space:=function(ListMat)
    local TmpDir, FileI, FileO, FileE, eProg, TheCommand, ListMatRet;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Saturation.in");
    FileO:=Filename(TmpDir, "Saturation.out");
    FileE:=Filename(TmpDir, "Saturation.err");
    WriteListMatrixFile(FileI, ListMat);

    eProg:="../../src_latt/TSPACE_IntegralSaturation";
    TheCommand:=Concatenation(eProg, " mpq_class ", FileI, " GAP ", FileO);
    Exec(TheCommand);

    ListMatRet:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return ListMatRet;
end;

get_pairwise_scalar_inv:=function(ListMat, SuperMat)
    local ListB, pairwise_scalar, eB, eLine;
    ListB:=List(ListMat, x->x*Inverse(SuperMat));
    pairwise_scalar:=[];
    for eB in ListB
    do
        eLine:=List(ListB, x->Trace(x * eB));
        Add(pairwise_scalar, eLine);
    od;
    return Inverse(pairwise_scalar);
end;


# Reading the files from the storage
# It is fixed because we do not want to change them
read_existing_file:=function(eFile)
    local is, SuperMat, ListMat_pre, ListMat, ListComm, ListSubspaces, PtStabGens, l_spanning_elements, PairwiseScalarInv;
    is:=InputTextFile(eFile);
    SuperMat:=InputStreamMatrix(is);
    ListMat_pre:=InputStreamListMatrix(is);
    ListMat:=get_saturated_space(ListMat_pre);
    ListComm:=InputStreamListMatrix(is);
    ListSubspaces:=InputStreamListMatrix(is);
    PtStabGens:=InputStreamListMatrix(is);
    l_spanning_elements:=[];
    PairwiseScalarInv:=get_pairwise_scalar_inv(ListMat, SuperMat);
    return rec(SuperMat:=SuperMat, ListMat:=ListMat, ListComm:=ListComm, ListSubspaces:=ListSubspaces, PtStabGens:=PtStabGens, l_spanning_elements:=l_spanning_elements, PairwiseScalarInv:=PairwiseScalarInv);
end;


# Write to the file. That should follow the update to the C++ code.
write_linear_space_input:=function(eFile, RecLinSpa)
    local os;
    os:=OutputTextFile(eFile, true);
    OutputStreamMatrix(os, RecLinSpa.SuperMat);
    OutputStreamListMatrix(os, RecLinSpa.ListMat);
    OutputStreamListMatrix(os, RecLinSpa.ListComm);
    OutputStreamListMatrix(os, RecLinSpa.ListSubspaces);
    OutputStreamListMatrix(os, RecLinSpa.PtStabGens);
    OutputStreamListMatrix(os, RecLinSpa.l_spanning_elements);
    OutputStreamMatrix(os, RecLinSpa.PairwiseScalarInv);
    CloseStream(os);
end;





get_nb_domains:=function(eFile)
    local FileLinSpa, FileNml, FileResult, RecLinSpa, eProg, TheCommand, U, is_correct;
    #
    FileLinSpa:="TSPACE_LinSpa";
    FileNml:="Enum_Tspaces_CI.nml";
    FileResult:="Result";
    RecLinSpa:=read_existing_file(eFile);
    RemoveFileIfExist(FileLinSpa);
    write_linear_space_input(FileLinSpa, RecLinSpa);
    #
    if method = "serial" then
        eProg:="../../src_delaunay/LATT_SerialLattice_IsoDelaunayDomain";
    fi;
    if method = "mpi" then
        eProg:="../../src_delaunay/LATT_MPI_Lattice_IsoDelaunayDomain";
    fi;
    TheCommand:=Concatenation(eProg, " ", FileNml);
    Exec(TheCommand);
    #
    if IsExistingFile(FileResult)=false then
        Print("method=", method, "\n");
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

ListRec:=ListRec{[1..3]};

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
#Recompute(eData1);
#Recompute(eData2);




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
