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
            Print("Pmat=\n");
            PrintArray(Pmat);
            NewMat:=Pmat * eMat * TransposedMat(Pmat);
            Add(FullListMat, NewMat);
            pos := pos+1;
            Add(eBlock, pos);
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
    Print("Before the Effective run\n");
    TheCommand:=Concatenation(eProg, " rational ", FileIn, " GAP ", FileOut);
    Exec(TheCommand);
    #
    if IsExistingFile(FileOut)=false then
        Error("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    CompListBlock:=ReadAsFunction(FileOut)();
    if Set(ListBlock) <> CompListBlock then
        Print("ListBlock=", ListBlock, "\n");
        Print("CompListBlock=", CompListBlock, "\n");
        Print("The computed sets are not the same\n");
        return false;
    fi;
    return true;
end;


ListMat:=[ClassicalSporadicLattices("A4"), ClassicalSporadicLattices("D4")];
n_equiv:=4;

result:=TestEquiStabDatabase(ListMat, n_equiv);
if result = false then
    # Error case
    Print("Error case\n")
    GAP_EXIT_CODE(1);
else
    # No error case
    Print("Normal case\n")
    GAP_EXIT_CODE(0);
fi;

