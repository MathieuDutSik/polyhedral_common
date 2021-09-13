
for dim in [6..8]
do
    if dim=6 then
        nbCase:=1;
    fi;
    if dim=7 then
        nbCase:=2;
    fi;
    if dim=8 then
        nbCase:=27;
    fi;
    for iCase in [1..nbCase]
    do
        Print("dim=", dim, " iCase=", iCase, "\n");
        EXT:=ClassicalExtremeDelaunayPolytopes([dim, iCase]);
        FAC:=DualDescription(EXT);
        FileSave:=Concatenation("FACextreme_", String(dim), "_", String(iCase));
        SYMPOL_PrintMatrix(FileSave, FAC);
    od;
od;
