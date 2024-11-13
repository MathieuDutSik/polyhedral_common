Read("../common.g");
Print("Beginning TestCanonicalization\n");

G_A4:=[ [ 2, -1, 0, 0 ], [ -1, 2, -1, 0 ], [ 0, -1, 2, -1 ], [ 0, 0, -1, 2 ] ];
G_A5:=[ [ 2, -1, 0, 0, 0 ], [ -1, 2, -1, 0, 0 ], [ 0, -1, 2, -1, 0 ], [ 0, 0, -1, 2, -1 ], [ 0, 0, 0, -1, 2 ] ];
G_E6:=[ [ 2, 1, 1, 0, 1, 1 ], [ 1, 4, 1, 1, 1, 3 ], [ 1, 1, 2, 1, 1, 1 ], [ 0, 1, 1, 2, 1, 2 ], [ 1, 1, 1, 1, 2, 2 ], [ 1, 3, 1, 2, 2, 4 ] ];
G_D4:=[ [ 2, 1, 1, 0 ], [ 1, 2, 1, -1 ], [ 1, 1, 2, -1 ], [ 0, -1, -1, 2 ] ];
G_D5:=[ [ 2, 1, 1, 1, 0 ], [ 1, 2, 1, 1, -1 ], [ 1, 1, 2, 1, -1 ], [ 1, 1, 1, 2, -1 ], [ 0, -1, -1, -1, 2 ] ];

# Conway/Sloane lattice in dimension 11 with no basis
G_CS11:=[ [ 60, 5, 5, 5, 5, 5, -12, -12, -12, -12, -7 ],
  [ 5, 60, 5, 5, 5, 5, -12, -12, -12, -12, -7 ],
  [ 5, 5, 60, 5, 5, 5, -12, -12, -12, -12, -7 ],
  [ 5, 5, 5, 60, 5, 5, -12, -12, -12, -12, -7 ],
  [ 5, 5, 5, 5, 60, 5, -12, -12, -12, -12, -7 ],
  [ 5, 5, 5, 5, 5, 60, -12, -12, -12, -12, -7 ],
  [ -12, -12, -12, -12, -12, -12, 60, -1, -1, -1, -13 ],
  [ -12, -12, -12, -12, -12, -12, -1, 60, -1, -1, -13 ],
  [ -12, -12, -12, -12, -12, -12, -1, -1, 60, -1, -13 ],
  [ -12, -12, -12, -12, -12, -12, -1, -1, -1, 60, -13 ],
  [ -7, -7, -7, -7, -7, -7, -13, -13, -13, -13, 96 ] ];

ListMat:=[G_A4, G_A5, G_E6, G_D4, G_D5, G_CS11];


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
    eProg:="../../src_latt/LATT_Canonicalize";
    TheCommand:=Concatenation(eProg, " gmp ", FileIn, " GAP_full ", FileOut);
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
Print("n_error=", n_error, "\n");

CI_Decision_Reset();
if n_error > 0 then
    Print("Error case\n");
else
    Print("Normal case\n");
    CI_Write_Ok();
fi;

