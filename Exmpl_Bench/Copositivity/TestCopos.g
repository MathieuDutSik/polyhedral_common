
Print("Beginning TestCanonicalization\n");

Mat_Hoffman_Pereira:=[ [ 1, -1, 1, 0, 0, 1, -1 ], [ -1, 1, -1, 1, 0, 0, 1 ], [ 1, -1, 1, -1, 1, 0, 0 ], [ 0, 1, -1, 1, -1, 1, 0 ],
                       [ 0, 0, 1, -1, 1, -1, 1 ], [ 1, 0, 0, 1, -1, 1, -1 ], [ -1, 1, 0, 0, 1, -1, 1 ] ];


Mat_Horn:=[ [ 1, -1, 1, 1, -1 ], [ -1, 1, -1, 1, 1 ], [ 1, -1, 1, -1, 1 ], [ 1, 1, -1, 1, -1 ], [ -1, 1, 1, -1, 1 ] ];
Mat_Horn_Perturb:=[ [ 1, -1, 1, 1, -1 ], [ -1, 1, -1, 1, 1 ], [ 1, -1, 1, -1, 1 ], [ 1, 1, -1, 1, -1 ], [ -1, 1, 1, -1, 1 - 1/100 ] ];


Mat_test1:=[ [ 3, 4, 3, -3, -2 ], [ 4, 2, 0, 1, -2 ], [ 3, 0, 3, -1, -2 ], [ -3, 1, -1, 5, 3 ], [ -2, -2, -2, 3, 3 ] ];
Mat_test2:=[ [ 100, -72, -59, 120 ], [ -72, 100, -60, -46 ], [ -59, -60, 100, -60 ], [ 120, -46, -60, 100 ] ];




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
