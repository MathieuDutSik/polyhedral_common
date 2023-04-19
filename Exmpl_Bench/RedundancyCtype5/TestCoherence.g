
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
    U1:=ReadAsFunction(eFileIrred1)();
    U2:=ReadAsFunction(eFileIrred2)();
    if U1<>U2 then
        Error("Inconsistency problem between the two methods");
    fi;
od;
