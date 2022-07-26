
Print("Beginning Test copositivity\n");


case1:=rec(eMat:=[ [ 1, -1, 1, 0, 0, 1, -1 ], [ -1, 1, -1, 1, 0, 0, 1 ], [ 1, -1, 1, -1, 1, 0, 0 ], [ 0, 1, -1, 1, -1, 1, 0 ],
                   [ 0, 0, 1, -1, 1, -1, 1 ], [ 1, 0, 0, 1, -1, 1, -1 ], [ -1, 1, 0, 0, 1, -1, 1 ] ],
           name:="Hoffman_Pereira", reply:=true);

case2:=rec(eMat:=[ [ 1, -1, 1, 1, -1 ], [ -1, 1, -1, 1, 1 ], [ 1, -1, 1, -1, 1 ], [ 1, 1, -1, 1, -1 ], [ -1, 1, 1, -1, 1 ] ],
           name:="Horn", reply:=true);

case3:=rec(eMat:=[ [ 1, -1, 1, 1, -1 ], [ -1, 1, -1, 1, 1 ], [ 1, -1, 1, -1, 1 ], [ 1, 1, -1, 1, -1 ], [ -1, 1, 1, -1, 1 - 1/100 ] ],
           name:="Horn_perturb", reply:=false);

case4:=rec(eMat:=[ [ 3, 4, 3, -3, -2 ], [ 4, 2, 0, 1, -2 ], [ 3, 0, 3, -1, -2 ], [ -3, 1, -1, 5, 3 ], [ -2, -2, -2, 3, 3 ] ],
           name:="Dannenberg1", reply:=true);

case5:=rec(eMat:=[ [ 100, -72, -59, 120 ], [ -72, 100, -60, -46 ], [ -59, -60, 100, -60 ], [ 120, -46, -60, 100 ] ],
           name:="Dannenberg1", reply:=true);




TestReflectivity:=function(eCase)
    local n, FileIn, FileOut, output, i, j, eProg, TheCommand, U;
    n:=Length(eRec.eMat);
    FileIn:="Test.in";
    FileOut:="Test.out";
    #
    output:=OutputTextFile(FileIn, true);
    AppendTo(output, n, " ", n, "\n");
    for i in [1..n]
    do
        for j in [1..n]
        do
            AppendTo(output, " ", eRec.eMat[i][j]);
        od;
        AppendTo(output, "\n");
    od;
    CloseStream(output);
    #
    eProg:="../../src_copos/CP_TestCopositivity";
    TheCommand:=Concatenation(eProg, " ", FileIn, " GAP ", FileOut);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    U:=ReadAsFunction(FileOut)();
    RemoveFile(FileIn);
    RemoveFile(FileOut);
    return eCasec.reply = U.isCopositive;
end;


ListCase:=[case1, case2, case3, case4, case5];


for eCase in ListCase
do
    test:=TestReflectivity(eCase);
    if test=false then
        # Error case
        GAP_EXIT_CODE(1);
    fi;
od;
# No error case
GAP_EXIT_CODE(0);
