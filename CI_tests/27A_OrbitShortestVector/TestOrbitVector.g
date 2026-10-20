Read("../common.g");
Read("../access_points.g");


TestOrbitShortest:=function(eRec)
    local GramMat, U;
    GramMat:=GetGramMatrixFromString(eRec.name);
    U:=get_orbit_shortest(GramMat);
    if is_error(U) then
        return false;
    fi;
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

# The certified sieve of LATT_ShvecSieve against the ShortestVectors of GAP:
# no vector may have a norm below the minimum returned, and the minimal
# vectors must be exactly the ones of GAP, taken with both signs.
TestShvecSieve:=function(name, GramMat)
    local U, d, sv, SetGAP;
    U:=get_shvec_sieve(GramMat);
    if is_error(U) then
        Print("TestShvecSieve, ", name, ": program failure\n");
        return false;
    fi;
    # ShortestVectors wants an integral bound: the matrix and the bound are
    # scaled by the common denominator, which keeps the minimal vectors.
    d:=Lcm(List(Concatenation(GramMat), DenominatorRat));
    sv:=ShortestVectors(d * GramMat, d * U.min);
    Print("  ", name, " min=", U.min, " |SHV|=", Length(U.SHV), "\n");
    if Minimum(sv.norms) <> d * U.min then
        Print("TestShvecSieve, ", name, ": the minimum is wrong\n");
        return false;
    fi;
    SetGAP:=Set(Concatenation(sv.vectors, -sv.vectors));
    if Set(U.SHV) <> SetGAP then
        Print("TestShvecSieve, ", name, ": the minimal vectors are wrong\n");
        return false;
    fi;
    return true;
end;

get_random_gram:=function(n)
    local B;
    while true
    do
        B:=List([1..n], i->List([1..n], j->Random([-3..3])));
        if DeterminantMat(B) <> 0 then
            return B * TransposedMat(B);
        fi;
    od;
end;

FullTestShvecSieve:=function()
    local ListGram, eRec, P, iter;
    ListGram:=[];
    for eRec in ["A2", "A4", "D4", "D5", "E6", "E7", "E8"]
    do
        Add(ListGram, rec(name:=eRec, G:=ClassicalSporadicLattices(eRec)));
    od;
    Add(ListGram, rec(name:="Z5", G:=IdentityMat(5)));
    for iter in [1..4]
    do
        Add(ListGram, rec(name:=Concatenation("random", String(iter)), G:=get_random_gram(4 + iter)));
    od;
    Add(ListGram, rec(name:="random_half", G:=get_random_gram(5) / 2));
    for eRec in ListGram
    do
        if TestShvecSieve(eRec.name, eRec.G)=false then
            return false;
        fi;
        P:=RandomIntegralUnimodularMatrix(Length(eRec.G));
        if TestShvecSieve(Concatenation(eRec.name, "_conj"), P * eRec.G * TransposedMat(P))=false then
            return false;
        fi;
    od;
    return true;
end;

result:=FullTest() and FullTestShvecSieve();
Print("result=", result, "\n");

CI_Decision_Reset();
if result = false then
    Print("Error case\n");
else
    Print("Normal case\n");
    CI_Write_Ok();
fi;

