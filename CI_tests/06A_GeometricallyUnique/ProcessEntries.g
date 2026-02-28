Read("../common.g");
Read("../access_points.g");

TheDir:="AllExamples";
ListFiles:=ListFileDirectory(TheDir);

Print("ListFiles=", ListFiles);

TestOneMat:=function(FAC)
    local n_fac, dim, GRP_fac, eCent, i, ePerm, eMatr, eMatr_cgr, FACprod, NewFAC, eCentProd, eCentMap;
    n_fac:=Length(FAC);
    dim:=Length(FAC[1]);
    GRP_fac:=SymmetricGroup(n_fac);
    eCent:=get_interior_point(FAC);
    if is_error(eCent) then
        return false;
    fi;
    for i in [1..20]
    do
        ePerm:=Random(GRP_fac);
        eMatr:=RandomIntegralUnimodularMatrix(dim);
        eMatr_cgr:=Inverse(TransposedMat(eMatr));
        FACprod:=FAC * eMatr;
        NewFAC:=Permuted(FACprod, ePerm);
        eCentProd:=get_interior_point(NewFAC);
        if is_error(eCentProd) then
            return false;
        fi;
        eCentMap:=eCent * eMatr_cgr;
        if eCentMap<>eCentProd then
            return false;
        fi;
    od;
    return true;
end;



iFile:=0;
n_error:=0;
n_file:=Length(ListFiles);
for eFile in ListFiles
do
    eFileFac:=Concatenation(TheDir, "/", eFile);
    FAC:=ReadMatrixFile(eFileFac);
    result:=TestOneMat(FAC);
    if result=false then
        n_error:=n_error + 1;
    fi;
    iFile:=iFile+1;
    Print("iFile=", iFile, "/", n_file, " n_error=", n_error, "\n");
od;
CI_Decision_Reset();
if n_error > 0 then
    # Error case
    Print("Error case\n");
else
    # No error case
    Print("Normal case\n");
    CI_Write_Ok();
fi;
