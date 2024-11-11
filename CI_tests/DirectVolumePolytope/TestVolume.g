TestVolume:=function(eRec)
    local FileOut, eProg, TheCommand, the_volume;
    FileOut:=Filename(DirectoryTemporary(), "Test.out");
    #
    eProg:="../../src_poly/POLY_lrs_volume";
    TheCommand:=Concatenation(eProg, " rational ", eRec.FileIn, " GAP ", FileOut);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    the_volume:=ReadAsFunction(FileOut)();
    RemoveFile(FileOut);
    if the_volume<>eRec.the_volume then
        Print("Thevolume is incorrect\n");
        return false;
    fi;
    return true;
end;

eRec1:=rec(FileIn:="G6.ext", the_volume:=1/2);
eRec2:=rec(FileIn:="24cell_poly.ext", the_volume:=8);
eRec3:=rec(FileIn:="H3.ext", the_volume:=1);
ListRec:=[eRec1, eRec2, eRec3];

FullTest:=function()
    local iRec, eRec, test;
    iRec:=0;
    for eRec in ListRec
    do
        Print("iRec=", iRec, " / ", Length(ListRec), "\n");
        test:=TestVolume(eRec);
        if test=false then
            return false;
        fi;
        iRec:=iRec + 1;
    od;
    return true;
end;

test:=FullTest();
if test=false then
    # Error case
    Print("Error case\n");
    GAP_EXIT_CODE(1);
else
    # No error case
    Print("Normal case\n");
    GAP_EXIT_CODE(0);
fi;

