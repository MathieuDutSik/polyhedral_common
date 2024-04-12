
Print("Beginning Test copositivity\n");


case1:=rec(eMat:=[ [ 1, -1, 1, 0, 0, 1, -1 ], [ -1, 1, -1, 1, 0, 0, 1 ], [ 1, -1, 1, -1, 1, 0, 0 ], [ 0, 1, -1, 1, -1, 1, 0 ],
                   [ 0, 0, 1, -1, 1, -1, 1 ], [ 1, 0, 0, 1, -1, 1, -1 ], [ -1, 1, 0, 0, 1, -1, 1 ] ],
           name:="Hoffman_Pereira", reply:=true);;

case2:=rec(eMat:=[ [ 1, -1, 1, 1, -1 ], [ -1, 1, -1, 1, 1 ], [ 1, -1, 1, -1, 1 ], [ 1, 1, -1, 1, -1 ], [ -1, 1, 1, -1, 1 ] ],
           name:="Horn", reply:=true);;

case3:=rec(eMat:=[ [ 1, -1, 1, 1, -1 ], [ -1, 1, -1, 1, 1 ], [ 1, -1, 1, -1, 1 ], [ 1, 1, -1, 1, -1 ], [ -1, 1, 1, -1, 1 - 1/100 ] ],
           name:="Horn_perturb", reply:=false);;

case4:=rec(eMat:=[ [ 3, 4, 3, -3, -2 ], [ 4, 2, 0, 1, -2 ], [ 3, 0, 3, -1, -2 ], [ -3, 1, -1, 5, 3 ], [ -2, -2, -2, 3, 3 ] ],
           name:="Dannenberg1", reply:=true);;

case5:=rec(eMat:=[ [ 100, -72, -59, 120 ], [ -72, 100, -60, -46 ], [ -59, -60, 100, -60 ], [ 120, -46, -60, 100 ] ],
           name:="Dannenberg2", reply:=true);;




TestCopositivity:=function(eCase)
    local n, FileIn, FileOut, output, i, j, eProg, TheCommand, U, test;
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
    test:=eCase.reply = U.isCopositive;
    if test=false then
        Print("We have eCase.reply=", eCase.reply, " but U.isCopositive=", U.isCopositive, "\n");
        Print("That is inconsistent\n");
    fi;
    return test;
end;


ListCase:=[case1, case2, case3, case4, case5];

n_error:=0;
ListCaseError:=[];
for eCase in ListCase
do
    test:=TestCopositivity(eCase);
    if test=false then
        n_error:=n_error + 1;
        Add(ListCaseError, eCase);
    fi;
od;
Print("n_case=", Length(ListCase), " n_error=", n_error, "\n");
if n_error > 0 then
    # Error case
    GAP_EXIT_CODE(1);
else
    # No error case
    GAP_EXIT_CODE(0);
fi;

