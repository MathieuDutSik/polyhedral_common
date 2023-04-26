
DoTest2:=true;
DoTest3:=false;


nType:=76;
for i in [1..nType]
do
    eFile:=Concatenation("TheCtype_5_", String(i));
    eFileGRP:=Concatenation("TheCtype_5_", String(i), ".grp");
    eProg:="../../src_poly/GRP_LinPolytope_Automorphism";
    eCommand:=Concatenation(eProg, " rational ", eFile, " ", eFileGRP);
    Print("eCommand=", eCommand, "\n");
    Exec(eCommand);
    #
    eFileIrred1:=Concatenation("Irred_1_", String(i));
    fProg:="../../src_poly/POLY_redundancyClarkson";
    fCommand:=Concatenation(fProg, " GAP rational ", eFile, " ", eFileIrred1);
    Print("fCommand=", fCommand, "\n");
    Exec(fCommand);
    U1:=ReadAsFunction(eFileIrred1)();
    #
    if DoTest2 then
        eFileIrred2:=Concatenation("Irred_2_", String(i));
        gProg:="../../src_poly/POLY_redundancyClarksonGroup";
        gCommand:=Concatenation(gProg, " rational ", eFile, " ", eFileGRP, " GAP ", eFileIrred2);
        Print("gCommand=", gCommand, "\n");
        Exec(gCommand);
        U2:=ReadAsFunction(eFileIrred2)();
        #
        if U1<>U2 then
            Error("Inconsistency problem between method 1 and 2");
        fi;
    fi;
    #
    if DoTest3 then
        eFileIrred3:=Concatenation("Irred_3_", String(i));
        hProg:="../../src_poly/POLY_redundancy_Equivariant";
        hCommand:=Concatenation(hProg, " rational ", eFile, " ", eFileGRP, " GAP ", eFileIrred3);
        Print("hCommand=", hCommand, "\n");
        Exec(hCommand);
        U3:=ReadAsFunction(eFileIrred3)();
        #
        if U1<>U3 then
            Error("Inconsistency problem between method 1 and 3");
        fi;
    fi;
od;
