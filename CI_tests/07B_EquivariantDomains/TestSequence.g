Read("../common.g");
Print("Beginning Test enumeration of L-type domains in perfect forms in T-spaces\n");

#keep_err:=false;
keep_err:=true;

get_rec_info:=function(fProg, d, n, keep_error)
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
    AppendTo(output, " arithmetic_T = \"gmp_rational\"\n");
    AppendTo(output, " arithmetic_Tint = \"gmp_integer\"\n");
    AppendTo(output, "/\n");
    AppendTo(output, "\n");
    AppendTo(output, "&TSPACE\n");
    AppendTo(output, " TypeTspace = \"", info.type_tspace, "\"\n");
    AppendTo(output, " FileLinSpa = \"unset.linspa\"\n");
    AppendTo(output, " SuperMatMethod = \"NotNeeded\"\n");
    AppendTo(output, " ListComm = \"Use_realimag\"\n");
    AppendTo(output, " PtGroupMethod = \"Trivial\"\n");
    AppendTo(output, " FileListSubspaces = \"unset\"\n");
    AppendTo(output, " RealImagDim = ", n, "\n");
    AppendTo(output, " RealImagSum = ", info.eSum, "\n");
    AppendTo(output, " RealImagProd = ", info.eProd, "\n");
    AppendTo(output, "/\n");
    CloseStream(output);
    #
    eProg:=Concatenation("../../src_latt/", fProg);
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
    Add(ListCases, rec(n:=3, d:=-3));
    Add(ListCases, rec(n:=3, d:=-4));
    Add(ListCases, rec(n:=3, d:=-7));
    Add(ListCases, rec(n:=3, d:=-8));
    Add(ListCases, rec(n:=3, d:=-11));
    Add(ListCases, rec(n:=3, d:=-15));
    Add(ListCases, rec(n:=3, d:=-19));
    Add(ListCases, rec(n:=3, d:=-20));
    Add(ListCases, rec(n:=3, d:=-23));
    Add(ListCases, rec(n:=3, d:=-24));
    Add(ListCases, rec(n:=4, d:=-3));
    Add(ListCases, rec(n:=4, d:=-4));
end;

set_dim2_examples:=function()
    Add(ListCases, rec(n:=2, d:=-3));
    Add(ListCases, rec(n:=2, d:=-4));
    Add(ListCases, rec(n:=2, d:=-7));
    Add(ListCases, rec(n:=2, d:=-8));
    Add(ListCases, rec(n:=2, d:=-11));
    Add(ListCases, rec(n:=2, d:=-15));
    Add(ListCases, rec(n:=2, d:=-19));
    Add(ListCases, rec(n:=2, d:=-20));
    Add(ListCases, rec(n:=2, d:=-23));
    Add(ListCases, rec(n:=2, d:=-24));
end;

set_canonical_examples();
#set_dim2_examples();



#ListProg:=["LATT_SerialLattice_IsoDelaunayDomain", "LATT_MPI_Lattice_IsoDelaunayDomain"];
ListProg:=["LATT_SerialLattice_IsoDelaunayDomain"];

f_compute:=function()
    local n_case, i_case, eCase, fProg, eRec;
    n_case:=Length(ListCases);
    for i_case in [1..n_case]
    do
        eCase:=ListCases[i_case];
        Print("----------------------------------------------------------------------\n");
        Print("i_case=", i_case, "/", n_case, " eCase=", eCase, "\n");
        for fProg in ListProg
        do
            eRec:=get_rec_info(fProg, eCase.d, eCase.n, keep_err);
            if eRec=fail then
                Print("Failing because eRec=fail\n");
                return false;
            fi;
            if IsBound(eCase.n_domains) then
                if Length(eRec)<>eCase.n_domains then
                    Print("Enumeration, |eRec|=", Length(eRec), " n_perf=", eCase.n_domains, "\n");
                    return false;
                fi;
            else
                Print("Number of L-type domains found=", Length(eRec), "\n");
            fi;
        od;
    od;
    return true;
end;



test:=f_compute();
CI_Decision_Reset();
if test=false then
    Print("Error case\n");
else
    Print("Normal case\n");
    CI_Write_Ok();
fi;
