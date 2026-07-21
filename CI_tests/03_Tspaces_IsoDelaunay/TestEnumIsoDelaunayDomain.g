Read("../common.g");
Read("../access_points.g");
Print("Beginning Test enumeration of iso-Delaunay / L-type domains\n");

method:="serial";
#method:="mpi";

# ==========================================================================
# Phase 1: iso-Delaunay domains of stored T-spaces (Bravais / Coxeter).
# The T-space is read from a fixed storage file, written back out in the
# LinSpa format, and LATT_SerialLattice_IsoDelaunayDomain is run to count the
# iso-Delaunay domains, which is compared against the stored expected number.
# ==========================================================================

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
    local is, SuperMat, ListMat_pre, ListMat, ListComm, ListSubspaces, PtStabGens, l_spanning_elements, l_outer_elements, PairwiseScalarInv;
    is:=InputTextFile(eFile);
    SuperMat:=InputStreamMatrix(is);
    ListMat_pre:=InputStreamListMatrix(is);
    ListMat:=get_saturated_space(ListMat_pre);
    ListComm:=InputStreamListMatrix(is);
    ListSubspaces:=InputStreamListMatrix(is);
    PtStabGens:=InputStreamListMatrix(is);
    l_spanning_elements:=[];
    l_outer_elements:=[];
    PairwiseScalarInv:=get_pairwise_scalar_inv(ListMat, SuperMat);
    return rec(SuperMat:=SuperMat, ListMat:=ListMat, ListComm:=ListComm, ListSubspaces:=ListSubspaces, PtStabGens:=PtStabGens, l_spanning_elements:=l_spanning_elements, l_outer_elements:=l_outer_elements, PairwiseScalarInv:=PairwiseScalarInv);
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
    OutputStreamListMatrix(os, RecLinSpa.l_outer_elements);
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
        eProg:=GetBinaryFilename("LATT_SerialLattice_IsoDelaunayDomain");
    fi;
    if method = "mpi" then
        eProg:=GetBinaryFilename("LATT_MPI_Lattice_IsoDelaunayDomain");
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




TestStoredTspaces:=function()
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


# ==========================================================================
# Phase 2: L-type domains for T-spaces from algebraic theory (formerly the
# separate CI entry 07B_L_domains_algebraic_rings). It relies on the same
# binary LATT_SerialLattice_IsoDelaunayDomain. The T-space is described by an
# imaginary quadratic field (discriminant d, dimension n); the namelist is
# built on the fly and the number of L-type domains is compared to the
# expected value.
#   comm_choice = "Use_realimag": group GL_n(Z[x]).
#   comm_choice = "Trivial":      group GL_n(Z[x]) with the conjugation.
# ==========================================================================

#keep_err:=false;
keep_err:=true;

get_rec_info:=function(fProg, d, n, comm_choice, keep_error)
    local strRun, FileNml, FileResult, FileErr, output, eProg, TheCommand, U, is_correct, info;
    #
    strRun:=Concatenation("_", String(n), "_", String(d));
    FileNml:=Concatenation("LtypeDomains", strRun , ".nml");
    FileResult:=Concatenation("Result", strRun);
    FileErr:=Concatenation("ERR_enumeration", strRun);
    RemoveFileIfExist(FileNml);
    RemoveFileIfExist(FileResult);
    #
    info:=GetFundamentalInfo(d);
    if info.IsCorrect=false then
        Print("Discriminant d=", d, " is not valid, skipping\n");
        return fail;
    fi;
    #
    # Create the namelist file
    output:=OutputTextFile(FileNml, true);
    AppendTo(output, "&SYSTEM\n");
    AppendTo(output, " max_runtime_second = 0\n");
    AppendTo(output, " ApplyStdUnitbuf = T\n");
    AppendTo(output, " Saving = F\n");
    AppendTo(output, " Prefix = \"/tmp/LtypeDomain/\"\n");
    AppendTo(output, " OutFile = \"", FileResult, "\"\n");
    AppendTo(output, " OutFormat = \"ObjectGAP\"\n");
    AppendTo(output, "/\n");
    AppendTo(output, "\n");
    AppendTo(output, "&DATA\n");
    AppendTo(output, " arithmetic = \"gmp\"\n");
    AppendTo(output, "/\n");
    AppendTo(output, "\n");
    AppendTo(output, "&TSPACE\n");
    AppendTo(output, " TypeTspace = \"", info.type_tspace, "\"\n");
    AppendTo(output, " FileLinSpa = \"unset.linspa\"\n");
    AppendTo(output, " SuperMatMethod = \"NotNeeded\"\n");
    AppendTo(output, " ListComm = \"", comm_choice, "\"\n");
    AppendTo(output, " PtGroupMethod = \"Trivial\"\n");
    AppendTo(output, " FileListSubspaces = \"unset\"\n");
    AppendTo(output, " RealImagDim = ", n, "\n");
    AppendTo(output, " RealImagSum = ", info.eSum, "\n");
    AppendTo(output, " RealImagProd = ", info.eProd, "\n");
    AppendTo(output, "/\n");
    CloseStream(output);
    #
    eProg:=GetBinaryFilename(fProg);
    TheCommand:=Concatenation(eProg, " ", FileNml);
    if keep_error then
        TheCommand:=Concatenation(TheCommand, " 2> ", FileErr);
    fi;
    Print("keep_error=", keep_error, " TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    #
    if IsExistingFile(FileResult)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return fail;
    fi;
    U:=ReadAsFunction(FileResult)();
    RemoveFile(FileErr);
    RemoveFile(FileNml);
    RemoveFile(FileResult);
    return U;
end;





ListCases:=[];

set_canonical_examples:=function()
    Add(ListCases, rec(n:=4, d:=-4, comm_choice:="Use_realimag"));
    Add(ListCases, rec(n:=3, d:=-3, comm_choice:="Use_realimag", n_domains:=23413));
    Add(ListCases, rec(n:=3, d:=-4, comm_choice:="Use_realimag", n_domains:=206));
end;

set_dim2_examples:=function()
    Add(ListCases, rec(n:=2, d:=-3, comm_choice:="Use_realimag", n_domains:=1));
    Add(ListCases, rec(n:=2, d:=-4, comm_choice:="Use_realimag", n_domains:=1));
    Add(ListCases, rec(n:=2, d:=-7, comm_choice:="Use_realimag", n_domains:=16));
    Add(ListCases, rec(n:=2, d:=-8, comm_choice:="Use_realimag", n_domains:=5));
    Add(ListCases, rec(n:=2, d:=-11, comm_choice:="Use_realimag", n_domains:=58));
    Add(ListCases, rec(n:=2, d:=-15, comm_choice:="Use_realimag", n_domains:=127));
    Add(ListCases, rec(n:=2, d:=-19, comm_choice:="Use_realimag", n_domains:=198));
    Add(ListCases, rec(n:=2, d:=-20, comm_choice:="Use_realimag", n_domains:=61));
    Add(ListCases, rec(n:=2, d:=-23, comm_choice:="Use_realimag", n_domains:=343));
    Add(ListCases, rec(n:=2, d:=-24, comm_choice:="Use_realimag", n_domains:=86));
end;

set_dim2_nocomm_examples:=function()
    Add(ListCases, rec(n:=2, d:=-3, comm_choice:="Trivial", n_domains:=1));
    Add(ListCases, rec(n:=2, d:=-4, comm_choice:="Trivial", n_domains:=1));
    Add(ListCases, rec(n:=2, d:=-7, comm_choice:="Trivial", n_domains:=11));
    Add(ListCases, rec(n:=2, d:=-8, comm_choice:="Trivial", n_domains:=5));
    Add(ListCases, rec(n:=2, d:=-11, comm_choice:="Trivial", n_domains:=35));
    Add(ListCases, rec(n:=2, d:=-15, comm_choice:="Trivial", n_domains:=74));
    Add(ListCases, rec(n:=2, d:=-19, comm_choice:="Trivial", n_domains:=112));
    Add(ListCases, rec(n:=2, d:=-20, comm_choice:="Trivial", n_domains:=43));
    Add(ListCases, rec(n:=2, d:=-23, comm_choice:="Trivial", n_domains:=190));
    Add(ListCases, rec(n:=2, d:=-24, comm_choice:="Trivial", n_domains:=61));
end;

#set_canonical_examples();
set_dim2_examples();
set_dim2_nocomm_examples();



#ListProg:=["LATT_SerialLattice_IsoDelaunayDomain", "LATT_MPI_Lattice_IsoDelaunayDomain"];
ListProg:=["LATT_SerialLattice_IsoDelaunayDomain"];

TestAlgebraicRings:=function()
    local n_error, n_case, i_case, eCase, fProg, eRec;
    n_error:=0;
    n_case:=Length(ListCases);
    for i_case in [1..n_case]
    do
        eCase:=ListCases[i_case];
        Print("----------------------------------------------------------------------\n");
        Print("i_case=", i_case, "/", n_case, " eCase=", eCase, "\n");
        for fProg in ListProg
        do
            eRec:=get_rec_info(fProg, eCase.d, eCase.n, eCase.comm_choice, keep_err);
            if eRec=fail then
                Print("Failing because eRec=fail\n");
                n_error:=n_error+1;
            elif IsBound(eCase.n_domains) then
                if Length(eRec)<>eCase.n_domains then
                    Print("Enumeration, |eRec|=", Length(eRec), " n_perf=", eCase.n_domains, "\n");
                    n_error:=n_error+1;
                fi;
            else
                Print("Number of L-type domains found=", Length(eRec), "\n");
            fi;
        od;
    od;
    return n_error;
end;


# ==========================================================================
# Combined decision: both phases must pass.
# ==========================================================================

n_error:=TestStoredTspaces() + TestAlgebraicRings();
Print("n_error=", n_error, "\n");

CI_Decision_Reset();
if n_error > 0 then
    Print("Error case\n");
else
    Print("Normal case\n");
    CI_Write_Ok();
fi;
