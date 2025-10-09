Read("../common.g");
Print("Beginning Test copositivity\n");


case1:=rec(eMat:=[ [ 1, -1, 1, 0, 0, 1, -1 ], [ -1, 1, -1, 1, 0, 0, 1 ], [ 1, -1, 1, -1, 1, 0, 0 ], [ 0, 1, -1, 1, -1, 1, 0 ],
                   [ 0, 0, 1, -1, 1, -1, 1 ], [ 1, 0, 0, 1, -1, 1, -1 ], [ -1, 1, 0, 0, 1, -1, 1 ] ],
           name:="Hoffman_Pereira", reply:=true);;

case2:=rec(eMat:=[ [ 1, -1, 1, 1, -1 ], [ -1, 1, -1, 1, 1 ], [ 1, -1, 1, -1, 1 ], [ 1, 1, -1, 1, -1 ], [ -1, 1, 1, -1, 1 ] ],
           name:="Horn", reply:=true);;

case3:=rec(eMat:=[ [ 1, -1, 1, 1, -1 ], [ -1, 1, -1, 1, 1 ], [ 1, -1, 1, -1, 1 ], [ 1, 1, -1, 1, -1 ], [ -1, 1, 1, -1, 1 - 1/100 ] ],
           name:="Horn_perturb", reply:=false);;

case4:=rec(eMat:=[ [ 3, 4, 3, -3, -2 ], [ 4, 2, 0, 1, -2 ], [ 3, 0, 3, -1, -2 ], [ -3, 1, -1, 5, 3 ], [ -2, -2, -2, 3, 3 ] ],
           name:="Dannenberg1", reply:=false);;

case5:=rec(eMat:=[ [ 100, -72, -59, 120 ], [ -72, 100, -60, -46 ], [ -59, -60, 100, -60 ], [ 120, -46, -60, 100 ] ],
           name:="Dannenberg2", reply:=false);;

case6:=rec(eMat:=[
[17, -91/5, 33/2, 38/3, -36/5],
[91/5, 59/3, -53/4, 8, 33/4],
[33/2, -53/4, 39/4, -13/2, 8],
[38/3, 8, -13/2, 16/3, -13/3],
[-36/5, 33/4, 8, -13/3, 1373628701/353935575]],
           name:="Strekelj_Zalar_Construction_of_exceptional_copositive_matrices", reply:=true);



ListCase_block1:=[case1, case2, case3, case4, case5, case6];
#ListCase_block1:=[case6];

TheDir:="Examples_from_Alexander_Oertel";


GetCases_90:=function()
    local TheCommand, TheFileLS, ListFiles, ListCases, eFile, fFile, eMat, eCase;
    TheFileLS:=Filename(DirectoryTemporary(), "LS_file");
    TheCommand:=Concatenation("(cd ", TheDir, " && find . -name \"cop*\" > ", TheFileLS, ")");
    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    ListFiles:=ReadTextFile(TheFileLS);
    #
    ListCases:=[];
    for eFile in ListFiles
    do
      	fFile:=Concatenation(TheDir, "/", eFile);
	Print("fFile=", fFile, "\n");
        eMat:=ReadMatrixFile(fFile);
        eCase:=rec(eMat:=eMat, reply:=true, name:=eFile);
        Add(ListCases, eCase);
    od;
    return ListCases;
end;

ListCase_block2:=GetCases_90();

ListCase:=Concatenation(ListCase_block1, ListCase_block2);
#ListCase:=ListCase_block1;


TestCopositivity:=function(eCase)
    local n, FileIn, FileOut, output, i, j, eProg, TheCommand, U, test;
    n:=Length(eCase.eMat);
    FileIn:=Filename(DirectoryTemporary(), "Test.in");
    FileOut:=Filename(DirectoryTemporary(), "Test.out");
    #
    WriteMatrixFile(FileIn, eCase.eMat);
    #
    eProg:="../../src_copos/CP_TestCopositivity";
    TheCommand:=Concatenation(eProg, " gmp ", FileIn, " GAP ", FileOut);
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


ProcessAllCases:=function()
    local n_error, ListCaseError, i_case, eCase, test;
    n_error:=0;
    ListCaseError:=[];
    i_case:=0;
    for eCase in ListCase
    do
        i_case:=i_case + 1;
        Print("TestCopositivity, i_case=", i_case, " name=", eCase.name, "\n");
        test:=TestCopositivity(eCase);
        if test=false then
            n_error:=n_error + 1;
            Add(ListCaseError, eCase);
        fi;
    od;
    Print("n_case=", Length(ListCase), " n_error=", n_error, "\n");
    return rec(n_error:=n_error, ListCaseError:=ListCaseError);
end;

RecResult:=ProcessAllCases();
Print("RecResult.n_error=", RecResult.n_error, "\n");

if RecResult.n_error > 0 then
    # Error case
    Print("Error case\n");
    GAP_EXIT_CODE(1);
else
    # No error case
    Print("Normal case\n");
    GAP_EXIT_CODE(0);
fi;

