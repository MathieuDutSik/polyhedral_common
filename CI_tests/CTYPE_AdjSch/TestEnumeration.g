Read("../common.g");
Print("Beginning TestDelaunayEnumeration\n");

TestEnumeration:=function(eRec)
    local n, FileNml, FileOut, output, eProg, TheCommand, U, is_correct;
    FileNml:=Filename(DirectoryTemporary(), "DelaunayEnum.nml");
    FileOut:=Filename(DirectoryTemporary(), "DelaunayEnum.out");
    RemoveFileIfExist(FileNml);
    RemoveFileIfExist(FileOut);
    #
    output:=OutputTextFile(FileNml, true);
    AppendTo(output, "&DATA\n");
    AppendTo(output, " arithmetic_T = \"gmp_rational\"n");
    AppendTo(output, " arithmetic_Tint = \"gmp_integer\"n");
    AppendTo(output, " OutFormat = \"NumberGAP\"n");
    AppendTo(output, " OutFile = \"", FileOut, "\"n");
    AppendTo(output, " n = ", eRec.n, "n");
    AppendTo(output, " max_runtime_second = 0\n");
    AppendTo(output, " ApplyStdUnitbuf = T\n");
    AppendTo(output, " Saving = F\n");
    AppendTo(output, " Prefix = \"DATA/\"");
    AppendTo(output, "/\n");
    CloseStream(output);
    #
    eProg:="../../src_ctype_mpi/CTYP_MPI_AdjScheme";
    TheCommand:=Concatenation("mpirun -np 2 ", eProg, " ", FileNml);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return rec(is_correct:=false);
    fi;
    U:=ReadAsFunction(FileOut)();
    RemoveFile(FileNml);
    RemoveFile(FileOut);
    is_correct:=eRec.n_types = Length(U.nb);
    return rec(is_correct:=is_correct);
end;

ListRec:=[rec(n:=2, n_types:=1),
          rec(n:=3, n_types:=1),
          rec(n:=4, n_types:=3),
          rec(n:=5, n_types:=76)
         ];

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
        fi;
        iRec:=iRec + 1;
    od;
    return n_error;
end;

n_error:=FullTest();
if n_error > 0 then
    # Error case
    GAP_EXIT_CODE(1);
else
    # No error case
    GAP_EXIT_CODE(0);
fi;

