Read("../common.g");
Print("Beginning TestCoherence (redundancy)\n");

# ==========================================================================
# Phase 1: C-type-5 polytopes, cross-method consistency of the irredundant
# facet computation. For each stored polytope TheCtype_5_<i> we compute the
# automorphism group (GRP_LinPolytope_Automorphism) and the set of irredundant
# facets by several methods, checking they all agree:
#   1. POLY_redundancy       HitAndRun
#   2. POLY_redundancy       Clarkson
#   3. POLY_redundancyGroup  Equivariant   (uses the automorphism group)
#   4. POLY_redundancyGroup  ClarksonBlock (disabled, needs more work)
# ==========================================================================

DoTest2_Clarkson:=true;
DoTest3_Equivariant:=true;
DoTest4_ClarksonBlock:=true;



TestIdx:=function(i)
    local eFile, eFileGRP, eProg, eCommand, eFileIrred1, fProg, fCommand, U1, eFileIrred2, gProg2, gCommand2, U2, eFileIrred3, gProg3, gCommand3, U3, eFileIrred4, gProg4, gCommand4, U4;
    eFile:=Concatenation("TheCtype_5_", String(i));
    eFileGRP:=Concatenation("TheCtype_5_", String(i), ".grp");
    eProg:=GetBinaryFilename("GRP_LinPolytope_Automorphism");
    # The CPP output format is the native WriteGroup format, which is the
    # one that POLY_redundancyGroup can read back (GAP format cannot).
    eCommand:=Concatenation(eProg, " rational ", eFile, " CPP ", eFileGRP);
    Print("eCommand=", eCommand, "\n");
    Exec(eCommand);
    if IsExistingFile(eFileGRP)=false then
        Print("Missing file eFileGRP=", eFileGRP, "\n");
        return false;
    fi;
    #
    eFileIrred1:=Concatenation("Irred_1_", String(i));
    fProg:=GetBinaryFilename("POLY_redundancy");
    fCommand:=Concatenation(fProg, " HitAndRun mpq_class ", eFile, " GAP ", eFileIrred1);
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
        gProg2:=GetBinaryFilename("POLY_redundancy");
        gCommand2:=Concatenation(gProg2, " Clarkson mpq_class ", eFile, " GAP ", eFileIrred2);
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
        gProg3:=GetBinaryFilename("POLY_redundancyGroup");
        gCommand3:=Concatenation(gProg3, " Equivariant mpq_class ", eFile, " ", eFileGRP, " GAP ", eFileIrred3);
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
        gProg4:=GetBinaryFilename("POLY_redundancyGroup");
        gCommand4:=Concatenation(gProg4, " ClarksonBlock mpq_class ", eFile, " ", eFileGRP, " GAP ", eFileIrred4);
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

TestCtype5:=function()
    local n_error, i;
    n_error:=0;
    for i in [1..76]
    do
        if TestIdx(i)=false then
            n_error:=n_error + 1;
        fi;
    od;
    return n_error;
end;


# ==========================================================================
# Phase 2: equivariant redundancy on the larger "walls" polytope (formerly the
# separate CI entry 25_RedundantEquiv). It relies on the same binaries; we
# compute the automorphism group and the equivariant irredundant facets and
# check the number of irredundant facets.
# ==========================================================================

TestRedundancy:=function(eRec)
    local eFileGRP, eFileIrred, eProg, eCommand, fProg, fCommand, U;
    eFileGRP:=Filename(DirectoryTemporary(), "Test.grp");
    eFileIrred:=Filename(DirectoryTemporary(), "Test.irred");
    RemoveFileIfExist(eFileGRP);
    RemoveFileIfExist(eFileIrred);
    #
    eProg:=GetBinaryFilename("GRP_LinPolytope_Automorphism");
    eCommand:=Concatenation(eProg, " rational ", eRec.eFile, " CPP ", eFileGRP);
    Print("eCommand=", eCommand, "\n");
    Exec(eCommand);
    if IsExistingFile(eFileGRP)=false then
        Print("Missing file eFileGRP=", eFileGRP, "\n");
        return false;
    fi;
    Print("We have eFileGRP=", eFileGRP, "\n");
    #
    fProg:=GetBinaryFilename("POLY_redundancyGroup");
    fCommand:=Concatenation(fProg, " Equivariant mpq_class ", eRec.eFile, " ", eFileGRP, " GAP ", eFileIrred);
    Print("fCommand=", fCommand, "\n");
    Exec(fCommand);
    if IsExistingFile(eFileIrred)=false then
        Print("Missing file eFileIrred=", eFileIrred, "\n");
        return false;
    fi;
    Print("We have eFileIrred=", eFileIrred, "\n");
    U:=ReadAsFunction(eFileIrred)();
    RemoveFileIfExist(eFileGRP);
    RemoveFileIfExist(eFileIrred);
    if Length(U)<>eRec.n_irred then
        Print("Wrong number of entries");
        return false;
    fi;
    return true;
end;

TestWalls:=function()
    local n_error, ListRec, eRec;
    n_error:=0;
    # The value 2576 is confirmed by the Clarkson, ClarksonBlock and
    # Equivariant methods which all agree on the walls polytope.
    ListRec:=[rec(eFile:="walls", n_irred:=2576)];
    for eRec in ListRec
    do
        if TestRedundancy(eRec)=false then
            n_error:=n_error + 1;
        fi;
    od;
    return n_error;
end;


# ==========================================================================
# Combined decision: both phases must pass.
# ==========================================================================

n_error:=TestCtype5() + TestWalls();
Print("n_error=", n_error, "\n");

CI_Decision_Reset();
if n_error > 0 then
    Print("Error case\n");
else
    Print("Normal case\n");
    CI_Write_Ok();
fi;
