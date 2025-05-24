Read("../common.g");
Print("Beginning TestCoherencen");

DoTest2_Clarkson:=true;
DoTest3_Equivariant:=true;
DoTest4_ClarksonBlock:=false; # This requires more work



TestIdx:=function(i)
    local eFile, eFileGRP, eProg, eCommand, eFileIrred1, fProg, fCommand, U1, eFileIrred2, gProg2, gCommand2, U2, eFileIrred3, gProg3, gCommand3, U3, eFileIrred4, gProg4, gCommand4, U4;
    eFile:=Concatenation("TheCtype_5_", String(i));
    eFileGRP:=Concatenation("TheCtype_5_", String(i), ".grp");
    eProg:="../../src_group/GRP_LinPolytope_Automorphism";
    eCommand:=Concatenation(eProg, " rational ", eFile, " GAP ", eFileGRP);
    Print("eCommand=", eCommand, "\n");
    Exec(eCommand);
    if IsExistingFile(eFileGRP)=false then
        Print("Missing file eFileGRP=", eFileGRP, "\n");
        return false;
    fi;
    #
    eFileIrred1:=Concatenation("Irred_1_", String(i));
    fProg:="../../src_poly/POLY_redundancy";
    fCommand:=Concatenation(fProg, " HitAndRun rational ", eFile, " GAP ", eFileIrred1);
    Print("fCommand=", fCommand, "\n");
    Exec(fCommand);
    if IsExistingFile(eFileIrred1)=false then
        Print("Missing file eFileIrred1\n");
        return false;
    fi;
    U1:=ReadAsFunction(eFileIrred1)();
    RemoveFile(eFileIrred1);
    #
    if DoTest2_Clarkson then
        eFileIrred2:=Concatenation("Irred_2_", String(i));
        gProg2:="../../src_poly/POLY_redundancy";
        gCommand2:=Concatenation(gProg2, " Clarkson rational ", eFile, " GAP ", eFileIrred2);
        Print("gCommand2=", gCommand2, "\n");
        Exec(gCommand2);
        if IsExistingFile(eFileIrred2)=false then
            Print("Missing file eFileIrred2\n");
            return false;
        fi;
        U2:=ReadAsFunction(eFileIrred2)();
        RemoveFile(eFileIrred2);
        #
        if U1<>U2 then
            Print("Inconsistency problem between method 1 and 2 at i=", i, "\n");
            return false;
        fi;
    fi;
    #
    if DoTest3_Equivariant then
        eFileIrred3:=Concatenation("Irred_3_", String(i));
        gProg3:="../../src_poly/POLY_redundancyGroup";
        gCommand3:=Concatenation(gProg3, " Equivariant rational ", eFile, " ", eFileGRP, " GAP ", eFileIrred3);
        Exec(gCommand3);
        Print("gCommand3=", gCommand3, "\n");
        if IsExistingFile(eFileIrred3)=false then
            Print("Missing file eFileIrred3\n");
            return false;
        fi;
        U3:=ReadAsFunction(eFileIrred3)();
        RemoveFile(eFileIrred3);
        #
        if U1<>U3 then
            Print("Inconsistency problem between method 1 and 3 at i=", i, "\n");
            return false;
        fi;
    fi;
    #
    if DoTest4_ClarksonBlock then
        eFileIrred4:=Concatenation("Irred_4_", String(i));
        gProg4:="../../src_poly/POLY_redundancyGroup";
        gCommand4:=Concatenation(gProg4, " ClarksonBlock rational ", eFile, " ", eFileGRP, " GAP ", eFileIrred4);
        Print("gCommand4=", gCommand4, "\n");
        Exec(gCommand4);
        if IsExistingFile(eFileIrred4)=false then
            Print("Missing file eFileIrred4\n");
            return false;
        fi;
        U4:=ReadAsFunction(eFileIrred4)();
        RemoveFile(eFileIrred4);
        #
        if U1<>U4 then
            Print("Inconsistency problem between method 1 and 4 at i=", i, "\n");
            return false;
        fi;
    fi;
    return true;
end;



n_error:=0;
ListIdx:=[1..76];
#ListIdx:=[7];

for i in ListIdx
do
    test:=TestIdx(i);
    if test=false then
        n_error:=n_error + 1;
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


