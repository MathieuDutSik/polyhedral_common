
nType:=76;
for i in [1..nType]
do
    eFile:=Concatenation("TheCtype_5_", String(i));
    eFileGRP:=Concatenation("TheCtype_5_", String(i), ".grp");
    eProg:="../../src_poly/GRP_LinPolytope_Automorphism";
    eCommand:=Concatenation(eProg, " rational ", eFile, " ", eFileGRP);
    Exec(eCommand);
    #
    eFileIrred1:=Concatenation("Irred_1_", String(i));
    fProg:="../../src_poly/POLY_redundancyClarkson";
    fCommand:=Concatenation(fProg, " GAP rational ", eFile, " ", eFileIrred1);
    Exec(fCommand);
    #
    eFileIrred2:=Concatenation("Irred_2_", String(i));
    gProg:="../../src_poly/POLY_redundancyClarksonGroup";
    gCommand:=Concatenation(gProg, " rational ", eFile, " ", eFileGRP, " GAP ", eFileIrred2);
    Exec(gCommand);
    #
    eFileIrred3:=Concatenation("Irred_3_", String(i));
    hProg:="../../src_poly/POLY_redundancy_Equivariant";
    hCommand:=Concatenation(hProg, " rational ", eFile, " ", eFileGRP, " GAP ", eFileIrred3);
    Exec(hCommand);
    #
    U1:=ReadAsFunction(eFileIrred1)();
    U2:=ReadAsFunction(eFileIrred2)();
    U3:=ReadAsFunction(eFileIrred3)();
    if U1<>U2 then
        Error("Inconsistency problem between method 1 and 2");
    fi;
    if U1<>U3 then
        Error("Inconsistency problem between method 1 and 3");
    fi;
od;
