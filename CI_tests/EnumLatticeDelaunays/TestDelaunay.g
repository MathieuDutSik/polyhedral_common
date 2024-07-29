Read("../common.g");
Print("Beginning TestDelaunayEnumeration\n");

TestEnumeration:=function(eRec)
    local n, FileIn, FileNml, FileOut, output, eProg, TheCommand, U, is_correct;
    n:=Length(eRec.eG);
    FileIn:=Filename(DirectoryTemporary(), "DelaunayEnum.in");
    FileNml:=Filename(DirectoryTemporary(), "DelaunayEnum.nml");
    FileOut:=Filename(DirectoryTemporary(), "DelaunayEnum.out");
    RemoveFileIfExist(FileIn);
    RemoveFileIfExist(FileNml);
    RemoveFileIfExist(FileOut);
    #
    WriteMatrixFile(FileIn, eRec.eG);
    #
    output:=OutputTextFile(FileNml, true);
    AppendTo(output, "&DATA\n");
    AppendTo(output, " arithmetic_T = \"gmp_rational\"\n");
    AppendTo(output, " arithmetic_Tint = \"gmp_integer\"\n");
    AppendTo(output, " GRAMfile = \"", FileIn, "\"\n");
    AppendTo(output, " SVRfile = \"unset.svr\"\n");
    AppendTo(output, " OutFormat = \"GAP\"\n");
    AppendTo(output, " OutFile = \"", FileOut, "\"\n");
    AppendTo(output, " max_runtime_second = 10800\n");
    AppendTo(output, "/\n");
    AppendTo(output, "\n");
    AppendTo(output, "&STORAGE\n");
    AppendTo(output, " Saving = F\n");
    AppendTo(output, " Prefix = \"DATA/\"\n");
    AppendTo(output, "/\n");
    CloseStream(output);
    #
    eProg:="../../src_latt/LATT_MPI_ComputeDelaunay";
    TheCommand:=Concatenation("mpirun -np 2 ", eProg, " ", FileNml);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return rec(is_correct:=false);
    fi;
    U:=ReadAsFunction(FileOut)();
    RemoveFile(FileIn);
    RemoveFile(FileNml);
    RemoveFile(FileOut);
    is_correct:=eRec.n_del = Length(U);
    Print("n=", n, " n_del=", eRec.n_del, " is_correct=", is_correct, "\n");
    return rec(is_correct:=is_correct);
end;

ListRec:=ReadAsFunction("ListCases")();;


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

