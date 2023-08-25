
DoTest2_Clarkson:=true;
DoTest3_Equivariant:=true;
DoTest4_ClarksonBlock:=true;


ListIdx:=[1..76];
#ListIdx:=[7];
for i in ListIdx
do
    eFile:=Concatenation("TheCtype_5_", String(i));
    eFileGRP:=Concatenation("TheCtype_5_", String(i), ".grp");
    eProg:="../../src_poly/GRP_LinPolytope_Automorphism";
    eCommand:=Concatenation(eProg, " rational ", eFile, " ", eFileGRP);
    Print("eCommand=", eCommand, "\n");
    Exec(eCommand);
    #
    eFileIrred1:=Concatenation("Irred_1_", String(i));
    fProg:="../../src_poly/POLY_redundancy";
    fCommand:=Concatenation(fProg, " HitAndRun rational ", eFile, " GAP ", eFileIrred1);
    Print("fCommand=", fCommand, "\n");
    Exec(fCommand);
    U1:=ReadAsFunction(eFileIrred1)();
    #
    if DoTest2_Clarkson then
        eFileIrred2:=Concatenation("Irred_2_", String(i));
        gProg:="../../src_poly/POLY_redundancy";
        gCommand:=Concatenation(gProg, " Clarkson rational ", eFile, " GAP ", eFileIrred2);
        Print("gCommand=", gCommand, "\n");
        Exec(gCommand);
        U2:=ReadAsFunction(eFileIrred2)();
        #
        if U1<>U2 then
            Error("Inconsistency problem between method 1 and 2");
        fi;
    fi;
    #
    if DoTest3_Equivariant then
        eFileIrred3:=Concatenation("Irred_3_", String(i));
        gProg:="../../src_poly/POLY_redundancyGroup";
        gCommand:=Concatenation(gProg, " Equivariant rational ", eFile, " ", eFileGRP, " GAP ", eFileIrred2);
        Print("gCommand=", gCommand, "\n");
        Exec(gCommand);
        U3:=ReadAsFunction(eFileIrred3)();
        #
        if U1<>U3 then
            Error("Inconsistency problem between method 1 and 3");
        fi;
    fi;
    #
    if DoTest4_ClarksonBlock then
        eFileIrred3:=Concatenation("Irred_3_", String(i));
        hProg:="../../src_poly/POLY_redundancyGroup";
        hCommand:=Concatenation(hProg, " ClarksonBlock rational ", eFile, " ", eFileGRP, " GAP ", eFileIrred3);
        Print("hCommand=", hCommand, "\n");
        Exec(hCommand);
        U3:=ReadAsFunction(eFileIrred3)();
        #
        if U1<>U4 then
            Error("Inconsistency problem between method 1 and 4");
        fi;
    fi;
od;
