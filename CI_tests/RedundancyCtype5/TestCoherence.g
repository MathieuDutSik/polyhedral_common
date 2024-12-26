Read("../common.g");
Print("Beginning TestCoherencen");

DoTest2_Clarkson:=true;
DoTest3_Equivariant:=true;
DoTest4_ClarksonBlock:=false; # This requires more work


n_error:=0;
ListIdx:=[1..76];
#ListIdx:=[7];
for i in ListIdx
do
    eFile:=Concatenation("TheCtype_5_", String(i));
    eFileGRP:=Concatenation("TheCtype_5_", String(i), ".grp");
    eProg:="../../src_group/GRP_LinPolytope_Automorphism";
    eCommand:=Concatenation(eProg, " rational ", eFile, " GAP ", eFileGRP);
    Print("eCommand=", eCommand, "\n");
    Exec(eCommand);
    if IsExistingFile(eFileGRP)=false then
        Print("Missing file eFileGRP=", eFileGRP, "\n");
        Error("Correct the creation file");
    fi;
    #
    eFileIrred1:=Concatenation("Irred_1_", String(i));
    fProg:="../../src_poly/POLY_redundancy";
    fCommand:=Concatenation(fProg, " HitAndRun rational ", eFile, " GAP ", eFileIrred1);
    Print("fCommand=", fCommand, "\n");
    Exec(fCommand);
    U1:=ReadAsFunction(eFileIrred1)();
    RemoveFile(eFileIrred1);
    #
    if DoTest2_Clarkson then
        eFileIrred2:=Concatenation("Irred_2_", String(i));
        gProg:="../../src_poly/POLY_redundancy";
        gCommand:=Concatenation(gProg, " Clarkson rational ", eFile, " GAP ", eFileIrred2);
        Print("gCommand=", gCommand, "\n");
        Exec(gCommand);
        U2:=ReadAsFunction(eFileIrred2)();
        RemoveFile(eFileIrred2);
        #
        if U1<>U2 then
            Print("Inconsistency problem between method 1 and 2 at i=", i, "\n");
            n_error:=n_error + 1;
        fi;
    fi;
    #
    if DoTest3_Equivariant then
        eFileIrred3:=Concatenation("Irred_3_", String(i));
        gProg:="../../src_poly/POLY_redundancyGroup";
        gCommand:=Concatenation(gProg, " Equivariant rational ", eFile, " ", eFileGRP, " GAP ", eFileIrred3);
        Exec(gCommand);
        Print("gCommand=", gCommand, "\n");
        U3:=ReadAsFunction(eFileIrred3)();
        RemoveFile(eFileIrred3);
        #
        if U1<>U3 then
            Print("Inconsistency problem between method 1 and 3 at i=", i, "\n");
            n_error:=n_error + 1;
        fi;
    fi;
    #
    if DoTest4_ClarksonBlock then
        eFileIrred4:=Concatenation("Irred_4_", String(i));
        hProg:="../../src_poly/POLY_redundancyGroup";
        hCommand:=Concatenation(hProg, " ClarksonBlock rational ", eFile, " ", eFileGRP, " GAP ", eFileIrred4);
        Print("hCommand=", hCommand, "\n");
        Exec(hCommand);
        U4:=ReadAsFunction(eFileIrred4)();
        RemoveFile(eFileIrred4);
        #
        if U1<>U4 then
            Print("Inconsistency problem between method 1 and 4 at i=", i, "\n");
            n_error:=n_error + 1;
        fi;
    fi;
od;
Print("n_error=", n_error, "\n");

CI_Decision_Reset();
if n_error > 0 then
    Print("Error case\n");
else
    Print("Normal case\n");
    CI_Write_Ok();
fi;


