Read("../common.g");


TestOrbitShortest:=function(eRec)
    local GramMat, FileIn, FileOut, eProg, TheCommand, U, V, eNorm, eNormSqr;
    GramMat:=GetGramMatrixFromString(eRec.name);
    FileIn:=Filename(DirectoryTemporary(), "Test.in");
    FileOut:=Filename(DirectoryTemporary(), "Test.out");
    RemoveFileIfExist(FileIn);
    RemoveFileIfExist(FileOut);
    #
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
    if Length(U.vf)<>eRec.n_orbit then
        Print("The number of orbits is incorrect\n");
        Print("|U.vf|=", Length(U.vf), " n_orbit=", eRec.n_orbit, "\n");
        return false;
    fi;
    return true;
end;


eRec1:=rec(name:="A2", n_orbit:=1);
eRec2:=rec(name:="A3", n_orbit:=1);
ListRec:=[eRec1, eRec2];

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
if result = false then
    # Error case
    Print("Error case\n");
    GAP_EXIT_CODE(1);
else
    # No error case
    Print("Normal case\n");
    GAP_EXIT_CODE(0);
fi;

