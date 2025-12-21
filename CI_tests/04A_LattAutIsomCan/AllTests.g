Read("../common.g");
Print("Beginning TestCanonicalization\n");

ListMat:=[];

ClassicMatrices:=true;
ConwayExample:=true;
WellRoundedDim10:=true;

if ClassicMatrices then
    Add(ListMat, ClassicalSporadicLattices("A4"));
    Add(ListMat, ClassicalSporadicLattices("A5"));
    Add(ListMat, ClassicalSporadicLattices("E6"));
    Add(ListMat, ClassicalSporadicLattices("D4"));
    Add(ListMat, ClassicalSporadicLattices("D5"));
fi;

if ConwayExample then
    # Conway/Sloane lattice in dimension 11 with no basis
    Add(ListMat, ClassicalSporadicLattices("ConwaySloane11"));
fi;

if WellRoundedDim10 then
    FileSave:="../21B_ShortRealizability/ListGram_n10_rnk10";
    if IsExistingFile(FileSave) then
        ListGram:=ReadAsFunction(FileSave)();
        Append(ListMat, ListGram);
    fi;
fi;




get_canonical_form:=function(eMat)
    local n, FileIn, FileOut, output, i, j, eProg, TheCommand, U;
    n:=Length(eMat);
    FileIn:="Test.in";
    FileOut:="Test.out";
    output:=OutputTextFile(FileIn, true);
    AppendTo(output, n, " ", n, "\n");
    for i in [1..n]
    do
        for j in [1..n]
        do
            AppendTo(output, " ", eMat[i][j]);
        od;
        AppendTo(output, "\n");
    od;
    CloseStream(output);
    #
    eProg:="../../src_latt/LATT_Canonicalize";
    TheCommand:=Concatenation(eProg, " gmp ", FileIn, " GAP_full ", FileOut);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        return fail;
    fi;
    U:=ReadAsFunction(FileOut)();
    RemoveFile(FileIn);
    RemoveFile(FileOut);
    RemoveFile(FileOut);
    return U.eG;
end;


test_canonicalization_function:=function(eMat)
    local TheCan, n, GRP, iter, len, eP, eMat_B, TheCan_B, LGen, i;
    TheCan:=get_canonical_form(eMat);
    if TheCan=fail then
        return false;
    fi;
    n:=Length(eMat);
    GRP:=GeneralLinearGroup(n, Integers);
    LGen:=GeneratorsOfGroup(GRP);
    for iter in [1..10]
    do
        Print("  |eMat|=", Length(eMat), " iter=", iter, "\n");
        len:=Random([1..2*n]);
        eP:=IdentityMat(n);
        for i in [1..len]
        do
            eP := eP * Random(LGen);
        od;
        eMat_B:=eP * eMat * TransposedMat(eP);
        TheCan_B:=get_canonical_form(eMat_B);
        if TheCan_B=fail then
            return false;
        fi;
        if TheCan_B<>TheCan then
            return false;
        fi;
    od;
    return true;
end;




n_error:=0;
nMat:=Length(ListMat);
for iMat in [1..nMat]
do
    eMat:=ListMat[iMat];
    Print("         iMat=", iMat, "/", nMat, " |eMat|=", Length(eMat), "\n");
    test:=test_canonicalization_function(eMat);
    if test=false then
        n_error:=n_error+1;
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

