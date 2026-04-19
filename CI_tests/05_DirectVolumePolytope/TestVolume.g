Read("../common.g");
Read("../access_points.g");
Print("Beginning of TestVolume\n");

Eulerian:=function(n, k)
    if n = 0 then
        if k = 0 then
            return 1;
        else
            return 0;
        fi;
    fi;

    if k < 0 or k >= n then
        return 0;
    fi;

    return (n - k) * Eulerian(n - 1, k - 1) + (k + 1) * Eulerian(n - 1, k);
end;


# It is well known that the volume of the hypersimplex is related to the Eulerian numbers
# See e.g.
# 1: Thomas Lam and Alexander Postnikov, "Alcoved Polytopes I", https://arxiv.org/pdf/math/0501246
# 2: Gaku Liu, "Mixed Volumes of Hypersimplices", https://www.combinatorics.org/ojs/index.php/eljc/article/view/v23i3p19
#
# Original reference seems to be from Laplace:
# M. de Laplace, Oeuvres completes, Vol. 7, reedite par Gauthier-Villars, Paris, 1886
get_hypersimplex_case:=function(n,k)
    local the_volume, EXT, eSet, eEXT, pos, name;
    the_volume:=Eulerian(n-1,k-1) / Factorial(n-1);
    EXT:=[];
    for eSet in Combinations([1..n], k)
    do
        eEXT:=ListWithIdenticalEntries(n, 0);
        for pos in eSet
        do
            eEXT[pos]:=1;
        od;
        eEXT[1]:=1;
        Add(EXT, eEXT);
    od;
    name:=Concatenation("Hypersimplex(", String(n), ",", String(k), ")");
    return rec(EXT:=EXT, the_volume:=the_volume, name:=name);
end;



TestVolume:=function(eRec)
    local the_volume;
    the_volume:=get_ext_volume(eRec.EXT);
    if is_error(the_volume) then
        return false;
    fi;
    if the_volume<>eRec.the_volume then
        Print("the_colume=", the_volume, " eRec.the_volume=", eRec.the_volume, "\n");
        Print("Incohency in the computation\n");
        return false;
    fi;
    return true;
end;

eRec1:=rec(FileIn:="G6.ext", the_volume:=1/2, name:="G6");
eRec2:=rec(FileIn:="24cell_poly.ext", the_volume:=8, name:="24cell");
eRec3:=rec(FileIn:="H3.ext", the_volume:=1, name:="H3");
ListRec:=[];
for eRec in [eRec1, eRec2, eRec3]
do
    EXT:=ReadMatrixFile(eRec.FileIn);
    fRec:=rec(EXT:=EXT, the_volume:=eRec.the_volume, name:=eRec.name);
    Add(ListRec, fRec);
od;

# Add the hypersimplex cases.
Add(ListRec, get_hypersimplex_case(6, 2));
Add(ListRec, get_hypersimplex_case(6, 3));
Add(ListRec, get_hypersimplex_case(6, 3));





FullTest:=function()
    local iRec, eRec, test;
    iRec:=0;
    for eRec in ListRec
    do
        Print("iRec=", iRec, " / ", Length(ListRec), " name=", eRec.name, "\n");
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

