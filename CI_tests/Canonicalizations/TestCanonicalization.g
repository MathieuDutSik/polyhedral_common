
Print("Beginning TestCanonicalization\n");

G_A4:=[ [ 2, -1, 0, 0 ], [ -1, 2, -1, 0 ], [ 0, -1, 2, -1 ], [ 0, 0, -1, 2 ] ];
G_A5:=[ [ 2, -1, 0, 0, 0 ], [ -1, 2, -1, 0, 0 ], [ 0, -1, 2, -1, 0 ], [ 0, 0, -1, 2, -1 ], [ 0, 0, 0, -1, 2 ] ];
G_E6:=[ [ 2, 1, 1, 0, 1, 1 ], [ 1, 4, 1, 1, 1, 3 ], [ 1, 1, 2, 1, 1, 1 ], [ 0, 1, 1, 2, 1, 2 ], [ 1, 1, 1, 1, 2, 2 ], [ 1, 3, 1, 2, 2, 4 ] ];
G_D4:=[ [ 2, 1, 1, 0 ], [ 1, 2, 1, -1 ], [ 1, 1, 2, -1 ], [ 0, -1, -1, 2 ] ];
G_D5:=[ [ 2, 1, 1, 1, 0 ], [ 1, 2, 1, 1, -1 ], [ 1, 1, 2, 1, -1 ], [ 1, 1, 1, 2, -1 ], [ 0, -1, -1, -1, 2 ] ];




ListMat:=[G_A4, G_A5, G_E6, G_D4, G_D5];


GetCanonicalForm:=function(eMat)
    local n, FileIn, FileOut, output, i, j, eProg, TheCommand, U;
    n:=Length(eMat);
    FileIn:="Test.in";
    FileOut:="Test.out";
    output:=OutputTextFile(FileIn, true);
    AppendTo(output, n, " ", n, "\n");
    for i in [1..n]
    do
        for j in [1..n]
        do
            AppendTo(output, " ", eMat[i][j]);
        od;
        AppendTo(output, "\n");
    od;
    CloseStream(output);
    #
    eProg:="./src_latt/LATT_canonicalize";
    TheCommand:=Concatenation(eProg, " 2 ", FileIn, " ", FileOut);
    Exec(TheCommand);
    U:=ReadAsFunction(FileOut)();
    RemoveFile(FileIn);
    RemoveFile(FileOut);
    RemoveFile(FileOut);
    return U.eG;
end;




n_error:=0;
for iMat in [1..Length(ListMat)]
do
    eMat:=ListMat[iMat];
    TheCan:=GetCanonicalForm(eMat);
    n:=Length(eMat);
    GRP:=GeneralLinearGroup(n, Integers);
    LGen:=GeneratorsOfGroup(GRP);
    for iter in [1..10]
    do
        Print("iMat=", iMat, " |eMat|=", Length(eMat), " iter=", iter, "\n");
        len:=Random([1..2*n]);
        eP:=IdentityMat(n);
        for i in [1..len]
        do
            eP := eP * Random(LGen);
        od;
        eMat_B:=eP * eMat * TransposedMat(eP);
        TheCan_B:=GetCanonicalForm(eMat_B);
        if TheCan_B<>TheCan then
            n_error:=n_error + 1;
        fi;
    od;
od;
if n_error > 0 then
    # Error case
    Print("Error case\n");
    GAP_EXIT_CODE(1);
else
    # No error case
    Print("Normal case\n");
    GAP_EXIT_CODE(0);
fi;

