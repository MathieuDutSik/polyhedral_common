Print("Beginning TestReflectivity\n");

WriteMatrix:=function(TheFile, TheMat)
    local output, nRow, nCol, iRow, iCol;
    nRow:=Length(TheMat);
    nCol:=Length(TheMat[1]);
    output:=OutputTextFile(TheFile, true);
    AppendTo(output, nRow, " ", nCol, "\n");
    for iRow in [1..nRow]
    do
        for iCol in [1..nCol]
        do
            AppendTo(output, " ", TheMat[iRow][iCol]);
        od;
        AppendTo(output, "\n");
    od;
    CloseStream(output);
end;

TestReflectivity:=function(eRec)
    local n, FileIn, FileNml, FileOut, output, i, j, eProg, TheCommand, U;
    n:=Length(eRec.LorMat);
    FileIn:="Test.in";
    FileNml:="Test.nml";
    FileOut:="Test.out";
    RemoveFileIfExist(FileIn);
    RemoveFileIfExist(FileNml);
    RemoveFileIfExist(FileOut);
    #
    WriteMatrix(FileIn, eRec.LorMat);
    #
    output:=OutputTextFile(FileNml, true);
    AppendTo(output, "&PROC\n");
    AppendTo(output, " FileLorMat = \"", FileIn, "\"\n");
    AppendTo(output, " OptionInitialVertex = \"vinberg\"\n");
    AppendTo(output, " OutFormat = \"GAP\"\n");
    AppendTo(output, " FileOut = \"", FileOut, "\"\n");
    AppendTo(output, " OptionNorms = \"all\"\n");
    AppendTo(output, " EarlyTerminationIfNotReflective = T\n");
    AppendTo(output, " ComputeAllSimpleRoots = T\n");
    AppendTo(output, "/\n");
    CloseStream(output);
    #
    eProg:="../../src_lorentzian/LORENTZ_FundDomain_AllcockEdgewalk";
    TheCommand:=Concatenation(eProg, " ", FileNml);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    U:=ReadAsFunction(FileOut)();
    RemoveFile(FileIn);
    RemoveFile(FileNml);
    RemoveFile(FileOut);
    GRPmatr:=Group(eRec.ListIsomCox);
    Print("|GRPmatr|=", Order(GRPmatr), "\n");
    return eRec.n_simple = U.n_simple;
end;

ListRec:=ReadAsFunction("ListReflect")();;

n_error:=0;
for eRec in ListRec
do
    test:=TestReflectivity(eRec);
    if test=false then
        n_error:=n_error+1;
    fi;
od;
if n_error > 0 then
    # Error case
    GAP_EXIT_CODE(1);
else
    # No error case
    GAP_EXIT_CODE(0);
fi;

