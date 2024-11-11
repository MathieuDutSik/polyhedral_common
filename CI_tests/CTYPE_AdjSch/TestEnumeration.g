Read("../common.g");
Print("Beginning Iso-Delaunay Enumeration\n");

TestEnumeration:=function(eRec)
    local n, FileNml, FileOut, output, strOut, eProg, TheCommand, U, is_correct;
    FileNml:=Filename(DirectoryTemporary(), "DelaunayEnum.nml");
    FileOut:=Filename(DirectoryTemporary(), "DelaunayEnum.out");
    RemoveFileIfExist(FileNml);
    RemoveFileIfExist(FileOut);
    #
    strOut:="&DATA\n";
    strOut:=Concatenation(strOut, " arithmetic_T = \"gmp_rational\"\n");
    strOut:=Concatenation(strOut, " arithmetic_Tint = \"gmp_integer\"\n");
    strOut:=Concatenation(strOut, " OutFormat = \"NumberGAP\"\n");
    strOut:=Concatenation(strOut, " OutFile = \"", FileOut, "\"\n");
    strOut:=Concatenation(strOut, " n = ", String(eRec.n), "\n");
    strOut:=Concatenation(strOut, " max_runtime_second = 0\n");
    strOut:=Concatenation(strOut, " ApplyStdUnitbuf = T\n");
    strOut:=Concatenation(strOut, " Saving = F\n");
    strOut:=Concatenation(strOut, " Prefix = \"DATA/\"\n");
    strOut:=Concatenation(strOut, "/\n");
    #
    output:=OutputTextFile(FileNml, true);
    WriteAll(output, strOut);
    CloseStream(output);
    #
    TheCommand:=Concatenation("cat ", FileNml);
    Exec(TheCommand);
    #
    eProg:="../../src_ctype/CTYP_MPI_AdjScheme";
    TheCommand:=Concatenation("mpirun -np 1 ", eProg, " ", FileNml);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return rec(is_correct:=false);
    fi;
    U:=ReadAsFunction(FileOut)();
    RemoveFile(FileNml);
    RemoveFile(FileOut);
    is_correct:=eRec.n_types = U.nb;
    return rec(is_correct:=is_correct);
end;

ListRec:=[rec(n:=2, n_types:=1),
          rec(n:=3, n_types:=1),
          rec(n:=4, n_types:=3),
          rec(n:=5, n_types:=77)
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
    Print("Error case\n");
    GAP_EXIT_CODE(1);
else
    # No error case
    Print("Normal case\n");
    GAP_EXIT_CODE(0);
fi;

