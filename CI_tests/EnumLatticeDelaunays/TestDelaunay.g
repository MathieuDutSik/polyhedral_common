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
    AppendTo(output, " arithmetic_T = \"gmp_rational\"");
    AppendTo(output, " arithmetic_Tint = \"gmp_integer\"");
    AppendTo(output, " GRAMfile = \"", FileIn, "\"");
    AppendTo(output, " SVRfile = \"unset.svr\"");
    AppendTo(output, " OutFormat = \"GAP\"");
    AppendTo(output, " OutFile = \"", FileOut, "\"");
    AppendTo(output, " max_runtime_second = 10800\n");
    AppendTo(output, "/\n");
    AppendTo(output, "\n");
    AppendTo(output, "&STORAGE\n");
    AppendTo(output, " Saving = F\n");
    AppendTo(output, " Prefix = \"DATA/\"");
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
    local n_error, iRec, eRec, RecReply, ListIsotropicCases, eIsotropicCase;
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
    GAP_EXIT_CODE(1);
else
    # No error case
    GAP_EXIT_CODE(0);
fi;

