
Print("Beginning TestCompletePositivity\n");



case1:=rec(eMat:=[ [ 2, 1, 1, 1, 2 ], [ 1, 2, 2, 1, 1 ], [ 1, 2, 6, 5, 1 ], [ 1, 1, 5, 6, 2 ], [ 2, 1, 1, 2, 3 ] ],
           name:="Sec41_Sponsel_Dur", reply:=true);;
case2:=rec(eMat:=[ [ 1, 1, 0, 0, 1 ], [ 1, 2, 1, 0, 0 ], [ 0, 1, 3, 1, 0 ], [ 0, 0, 1, 4, 1 ], [ 1, 0, 0, 1, 5 ] ],
           name:="Berman_Example2_7", reply:=true);;
case3:=rec(eMat:=[ [ 6, 4, 1, 2, 2 ], [ 4, 6, 0, 1, 3 ], [ 1, 0, 3, 1, 2 ], [ 2, 1, 1, 2, 1 ], [ 2, 3, 2, 1, 5 ] ],
           name:="Example6_1_Nie", reply:=true);;
case4:=rec(eMat:=[ [ 1, 1, 0, 0, 1 ], [ 1, 2, 1, 0, 0 ], [ 0, 1, 2, 1, 0 ], [ 0, 0, 1, 2, 1 ], [ 1, 0, 0, 1, 6 ] ],
           name:="Example6_2_Nie", reply:=false);;
case5:=rec(eMat:=[ [ 2, 1, 0, 0, 1 ], [ 1, 2, 1, 0, 0 ], [ 0, 1, 2, 1, 0 ], [ 0, 0, 1, 2, 1 ], [ 1, 0, 0, 1, 2 ] ],
           name:="Berman_Example2_22", reply:=true);;
case6:=rec(eMat:=[ [ 6, 0, 2, 2 ], [ 0, 5, 4, 2 ], [ 2, 4, 6, 0 ], [ 2, 2, 0, 6 ] ],
           name:="Berman_Example2_28", reply:=true);;
case7:=rec(eMat:=[ [ 2, 1, 1, 0, 0 ], [ 1, 2, 0, 0, 1 ], [ 1, 0, 2, 1, 0 ], [ 0, 0, 1, 2, 1 ], [ 0, 1, 0, 1, 2 ] ],
           name:="Berman_Example2_16", reply:=true);;
case8:=rec(eMat:=[ [ 2, 0, 0, 1, 1 ], [ 0, 2, 0, 1, 1 ], [ 0, 0, 2, 1, 1 ], [ 1, 1, 1, 2, 0 ], [ 1, 1, 1, 0, 2 ] ],
           name:="Berman_p122_Example2_17", reply:=false);;
case9:=rec(eMat:=[ [ 6, 5, 3, 0 ], [ 5, 11, 4, 0 ], [ 3, 4, 2, 0 ], [ 0, 0, 0, 0 ] ],
           name:="Berman_p128_Exercise2_45", reply:=true);;
case10:=rec(eMat:=[ [ 4, 5, 4, 6, 4, 2 ], [ 5, 1, 4, 7, 4, 6 ], [ 4, 4, 4, 2, 5, 4 ], [ 6, 7, 2, 0, 3, 7 ], [ 4, 4, 5, 3, 1, 6 ], [ 2, 6, 4, 7, 6, 4 ] ],
            name:="Fan_Zhou_Example5_4", reply:=false);;
case11:=rec(eMat:=[ [ 2, 0, 0, 1, 1 ], [ 0, 2, 0, 1, 1 ], [ 0, 0, 2, 1, 1 ], [ 1, 1, 1, 3, 0 ], [ 1, 1, 1, 0, 3 ] ],
            name:="Prakash_etal_Example3_1_cprank4", reply:=true);;
case12:=rec(eMat:=[ [ 2, 1, 0, 0, 0, 0, 1 ], [ 1, 2, 1, 0, 0, 0, 0 ], [ 0, 1, 2, 1, 0, 0, 0 ], [ 0, 0, 1, 2, 1, 0, 0 ], [ 0, 0, 0, 1, 2, 1, 0 ], [ 0, 0, 0, 0, 1, 2, 1 ], [ 1, 0, 0, 0, 0, 1, 2 ] ],
            name:="Zhou_Fan_Example5_2", reply:=true);;



TestCompletePositivity:=function(eCase)
    local n, FileIn, FileOut, output, i, j, eProg, TheCommand, U;
    n:=Length(eCase.eMat);
    FileIn:="Test.in";
    FileOut:="Test.out";
    #
    output:=OutputTextFile(FileIn, true);
    AppendTo(output, n, " ", n, "\n");
    for i in [1..n]
    do
        for j in [1..n]
        do
            AppendTo(output, " ", eCase.eMat[i][j]);
        od;
        AppendTo(output, "\n");
    od;
    CloseStream(output);
    #
    eProg:="../../src_copos/CP_TestCompletePositivity";
    TheCommand:=Concatenation(eProg, " ", FileIn, " GAP ", FileOut);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    U:=ReadAsFunction(FileOut)();
    RemoveFile(FileIn);
    RemoveFile(FileOut);
    return U.result = eCase.reply;
end;


ListCase:=[case1, case2, case3, case4, case5, case6,
           case7, case8, case9, case10, case11, case12];;

n_error:=0;
for eCase in ListCase
do
    test:=TestCompletePositivity(eCase);
    if test=false then
        n_error:=n_error + 1;
    fi;
od;
Print("n_error=", n_error, "\n");
if n_error > 0 then
    # Error case
    GAP_EXIT_CODE(1);
else
    # No error case
    GAP_EXIT_CODE(0);
fi;

