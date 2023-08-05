
Print("Beginning TestReflectivity\n");


TestReflectivity:=function(eRec)
    local n, FileIn, FileNml, FileOut, output, i, j, eProg, TheCommand, U;
    n:=Length(eRec.LorMat);
    FileIn:="Test.in";
    FileNml:="Test.nml";
    FileOut:="Test.out";
    #
    output:=OutputTextFile(FileIn, true);
    AppendTo(output, n, " ", n, "\n");
    for i in [1..n]
    do
        for j in [1..n]
        do
            AppendTo(output, " ", eRec.LorMat[i][j]);
        od;
        AppendTo(output, "\n");
    od;
    CloseStream(output);
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
    return eRec.n_simple = U.n_simple;
end;




ListRec:=ReadAsFunction("ListReflect")();;


for iRec in [1..Length(ListRec)]
do
    eRec:=ListRec[iRec];
    test:=TestReflectivity(eRec);
    if test=false then
        # Error case
        GAP_EXIT_CODE(1);
    fi;
od;
# No error case
GAP_EXIT_CODE(0);
