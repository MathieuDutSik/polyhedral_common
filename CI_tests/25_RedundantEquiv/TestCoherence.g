Read("../common.g");
Print("Beginning TestCoherence\n");

TestRedundancy:=function(eRec)
    local eFileGRP, eFileIrred, eProg, eCommand, fProg, fCommand, U;
    eFileGRP:=Filename(DirectoryTemporary(), "Test.grp");
    eFileIrred:=Filename(DirectoryTemporary(), "Test.irred");
    RemoveFileIfExist(eFileGRP);
    RemoveFileIfExist(eFileIrred);
    #
    eProg:="../../src_group/GRP_LinPolytope_Automorphism";
    eCommand:=Concatenation(eProg, " rational ", eRec.eFile, " GAP ", eFileGRP);
    Print("eCommand=", eCommand, "\n");
    Exec(eCommand);
    if IsExistingFile(eFileGRP)=false then
        Print("Missing file eFileGRP=", eFileGRP, "\n");
        return false;
    fi;
    Print("We have eFileGRP=", eFileGRP, "\n");
    #
    fProg:="../../src_poly/POLY_redundancyGroup";
    fCommand:=Concatenation(fProg, " Equivariant rational ", eRec.eFile, " ", eFileGRP, " GAP ", eFileIrred);
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

eRec1:=rec(eFile:="walls", n_irred:=2400); #Need to put the correct value
ListRec:=[eRec1];

FullTest:=function()
    local iRec, eRec, test;
    iRec:=0;
    for eRec in ListRec
    do
        Print("iRec=", iRec, " / ", Length(ListRec), "\n");
        test:=TestRedundancy(eRec);
        if test=false then
            return false;
        fi;
        iRec:=iRec + 1;
    od;
    return true;
end;

test:=FullTest();
Print("test=", test, "\n");

CI_Decision_Reset();
if test=false then
    # Error case
    Print("Error case\n");
else
    # No error case
    Print("Normal case\n");
    CI_Write_Ok();
fi;
