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
    strOut:="&DATA\n";
    strOut:=Concatenation(strOut, " arithmetic_T = \"gmp_rational\"\n");
    strOut:=Concatenation(strOut, " arithmetic_Tint = \"gmp_integer\"\n");
    strOut:=Concatenation(strOut, " GRAMfile = \"", FileIn, "\"\n");
    strOut:=Concatenation(strOut, " SVRfile = \"unset.svr\"\n");
    strOut:=Concatenation(strOut, " OutFormat = \"GAP\"\n");
    strOut:=Concatenation(strOut, " OutFile = \"", FileOut, "\"\n");
    strOut:=Concatenation(strOut, " max_runtime_second = 10800\n");
    strOut:=Concatenation(strOut, "/\n");
    strOut:=Concatenation(strOut, "\n");
    strOut:=Concatenation(strOut, "&STORAGE\n");
    strOut:=Concatenation(strOut, " Saving = F\n");
    strOut:=Concatenation(strOut, " Prefix = \"DATA/\"\n");
    strOut:=Concatenation(strOut, "/\n");
    #
    output:=OutputTextFile(FileNml, true);
    WriteAll(output, strOut);
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
        Print("RecReply=", RecReply, "\n");
        if RecReply.is_correct=false then
            n_error:=n_error+1;
            Print("1: n_error=", n_error, "\n");
            return n_error;
        fi;
        iRec:=iRec + 1;
    od;
    return n_error;
end;

n_error:=FullTest();
Print("2: n_error=", n_error, "\n");
for i in [1..10]
do
    Print("Loop i=", i, "\n");
od;

CI_Decision_Reset();
if n_error > 0 then
    # Error case
    Print("Error case\n");
else
    # No error case
    Print("Normal case\n");
    CI_Write_Ok();
fi;

