Read("../common.g");


TestOrbitShortest:=function(eRec)
    local GramMat, FileIn, FileOut, eProg, TheCommand, U, V, eNorm, eNormSqr;
    GramMat:=GetGramMatrixFromString(eRec.name);
    FileIn:=Filename(DirectoryTemporary(), "Test.in");
    FileOut:=Filename(DirectoryTemporary(), "Test.out");
    WriteMatrixFile(FileIn, GramMat);
    #
    eProg:="../../src_latt/LATT_ComputeShortestOrbits";
    TheCommand:=Concatenation(eProg, " gmp ", FileIn, " GAP ", FileOut);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    U:=ReadAsFunction(FileOut)();
    RemoveFile(FileIn);
    RemoveFile(FileOut);
    Print("|U.vf|=", Length(U.vf), " n_orbit=", eRec.n_orbit, "\n");
    if Length(U.vf)<>eRec.n_orbit then
        Print("The number of orbits is incorrect\n");
        return false;
    fi;
    return true;
end;

ListRec:=[];
Add(ListRec, rec(name:="A2", n_orbit:=1));
Add(ListRec, rec(name:="A3", n_orbit:=2));

FullTest:=function()
    local iRec, eRec, test;
    iRec:=0;
    for eRec in ListRec
    do
        Print("iRec=", iRec, " / ", Length(ListRec), "\n");
        test:=TestOrbitShortest(eRec);
        if test=false then
            return false;
        fi;
        iRec:=iRec + 1;
    od;
    return true;
end;

result:=FullTest();
Print("result=", result, "\n");

CI_Decision_Reset();
if result = false then
    Print("Error case\n");
else
    Print("Normal case\n");
    CI_Write_Ok();
fi;

