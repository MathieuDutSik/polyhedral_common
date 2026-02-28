Read("../common.g");
Read("../access_points.g");
Print("Beginning of TestVolume\n");

TestVolume:=function(eRec)
    local the_volume;
    the_volume:=get_ext_volume(eRec.EXT);
    if is_error(the_volume) then
        return false;
    fi;
    if the_volume<>eRec.the_volume then
        Print("Thevolume is incorrect\n");
        return false;
    fi;
    return true;
end;

eRec1:=rec(FileIn:="G6.ext", the_volume:=1/2);
eRec2:=rec(FileIn:="24cell_poly.ext", the_volume:=8);
eRec3:=rec(FileIn:="H3.ext", the_volume:=1);
ListRec:=[];
for eRec in [eRec1, eRec2, eRec3]
do
    EXT:=ReadMatrixFile(eRec.FileIn);
    fRec:=rec(EXT:=EXT, the_volume:=eRec.the_volume);
    Add(ListRec, fRec);
od;

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

CI_Decision_Reset();
if test=false then
    Print("Error case\n");
else
    # No error case
    Print("Normal case\n");
    CI_Write_Ok();
fi;

