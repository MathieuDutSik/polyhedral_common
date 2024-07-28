Read("../common.g");

TestEquiStabDatabase:=function(ListMat, n_equiv)
    local FullListMat, eMat, n, i, Pmat, NewMat, ListBlock, eBlock, pos, eProg, FileIn, FileOut, TheCommand, CompListBlock;
    FullListMat:=[];
    ListBlock:=[];
    pos:=0;
    for eMat in ListMat
    do
        n:=Length(eMat);
        eBlock:=[];
        for i in [1..n_equiv]
        do
            Pmat:=RandomIntegralUnimodularMatrix(n);
            NewMat:=Pmat * eMat * TransposedMat(Pmat);
            Add(FullListMat, NewMat);
            pos := pos+1;
            Add(ListBlock, pos);
        od;
        Add(ListBlock, eBlock);
    od;
    FileIn:=Filename(DirectoryTemporary(), "Test.in");
    FileOut:=Filename(DirectoryTemporary(), "Test.out");
    RemoveFileIfExist(FileIn);
    RemoveFileIfExist(FileOut);
    WriteListMatrixFile(FileIn, FullListMat);
    #
    eProg:="../../src_indefinite_models/TEST_EquiStabFamily";
    TheCommand:=Concatenation(eProg, " rational ", FileIn, " GAP ", FileOut);
    Exec(TheCommand);
    #
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    CompListBlock:=ReadAsFunction(FileOut)();
    if Set(ListBlock) <> CompListBlock then
        Print("The computed sets are not the same\n");
        return false;
    fi;
    return true;
end;


ListMat:=[ClassicalSporadicLattices("A4"), ClassicalSporadicLattices("D4")];
n_equiv:=10;

result:=TestEquiStabDatabase(ListMat, n_equiv);
if result = false then
    # Error case
    GAP_EXIT_CODE(1);
else
    # No error case
    GAP_EXIT_CODE(0);
fi;

