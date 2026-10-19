Read("../common.g");
Read("../access_points.g");
Print("Beginning Test of the presentations by generators and relators\n");

# The presentation of GL_n(Z) from the perfect form complex is checked by
# --- The relators holding for the matrices of the generators.
# --- The abelian invariants.
# --- The orders of the quotients by the normal closure of E_{12}^N. For n=2,
#     it is GL_2(Z) / Gamma(N) for N <= 5, of order 2 |SL_2(Z/N)| for N >= 3.
#     For n >= 3, it is GL_n(Z) / Gamma(N) by the congruence subgroup property
#     (Bass, Milnor, Serre), of order 2 |SL_n(Z/N)|.
# --- For n=2, the numbers of subgroups of small index, compared with the ones
#     of GL_2(Z) = D_8 *_{C_2 x C_2} D_12 (Serre, Trees, I.4.2).

# A word in the generators for the matrix target, by breadth first search.
FindWord:=function(R, target)
    local F, gens, lett, n, seen, frontier, nf, x, i, y;
    F:=FreeGroupOfFpGroup(R.fp);
    gens:=Concatenation(R.ListGen, List(R.ListGen, Inverse));
    lett:=Concatenation(GeneratorsOfGroup(F), List(GeneratorsOfGroup(F), Inverse));
    n:=Length(R.ListGen[1]);
    seen:=NewDictionary(IdentityMat(n), true);
    AddDictionary(seen, IdentityMat(n), One(F));
    frontier:=[IdentityMat(n)];
    while LookupDictionary(seen, target)=fail
    do
        nf:=[];
        for x in frontier
        do
            for i in [1..Length(gens)]
            do
                y:=x*gens[i];
                if LookupDictionary(seen, y)=fail then
                    AddDictionary(seen, y, LookupDictionary(seen, x)*lett[i]);
                    Add(nf, y);
                fi;
            od;
        od;
        frontier:=nf;
    od;
    return LookupDictionary(seen, target);
end;

CongruenceQuotientOrder:=function(R, N)
    local F, n, E, w, G;
    F:=FreeGroupOfFpGroup(R.fp);
    n:=Length(R.ListGen[1]);
    E:=IdentityMat(n);
    E[1][2]:=1;
    w:=FindWord(R, E);
    G:=F/Concatenation(RelatorsOfFpGroup(R.fp), [w^N]);
    return Size(G);
end;

StandardGL2Z:=function()
    local F, a, b, c;
    F:=FreeGroup("a", "b", "c");
    a:=F.1;
    b:=F.2;
    c:=F.3;
    return F/[a^2, b^4, (a*b)^2, c^6, (a*c)^2, b^2*c^-3];
end;

ListCases:=[rec(n:=2, abelian_invariants:=[2,2], quotient_orders:=[[2,12], [3,48], [4,96], [5,240]], max_index:=8),
            rec(n:=3, abelian_invariants:=[2], quotient_orders:=[[2,336], [3,11232]])];

test_case:=function(eCase)
    local desc, R, hom, inv, pair, order, G_std, k, n_ours, n_std;
    desc:=GenerateTspaceDescription_classic(eCase.n, true);
    R:=PERFCOMP_group_generator_relation(desc);
    Print("n=", eCase.n, " |gens|=", Length(R.ListGen), " |rels|=", Length(RelatorsOfFpGroup(R.fp)), "\n");
    hom:=GroupHomomorphismByImages(R.fp, Group(R.ListGen), GeneratorsOfGroup(R.fp), R.ListGen);
    if hom=fail then
        Print("The relators do not hold for the matrices of the generators\n");
        return false;
    fi;
    inv:=AbelianInvariants(R.fp);
    Print("AbelianInvariants=", inv, "\n");
    if inv<>eCase.abelian_invariants then
        Print("The abelian invariants should be ", eCase.abelian_invariants, "\n");
        return false;
    fi;
    for pair in eCase.quotient_orders
    do
        order:=CongruenceQuotientOrder(R, pair[1]);
        Print("N=", pair[1], " order=", order, "\n");
        if order<>pair[2] then
            Print("The order of the quotient should be ", pair[2], "\n");
            return false;
        fi;
    od;
    if IsBound(eCase.max_index) then
        G_std:=StandardGL2Z();
        for k in [1..eCase.max_index]
        do
            n_ours:=Length(LowIndexSubgroupsFpGroup(R.fp, k));
            n_std:=Length(LowIndexSubgroupsFpGroup(G_std, k));
            Print("index<=", k, " n_ours=", n_ours, " n_std=", n_std, "\n");
            if n_ours<>n_std then
                Print("The numbers of subgroups of small index should be equal\n");
                return false;
            fi;
        od;
    fi;
    return true;
end;

f_compute:=function()
    local eCase;
    for eCase in ListCases
    do
        Print("----------------------------------------------------------------------\n");
        if test_case(eCase)=false then
            return false;
        fi;
    od;
    return true;
end;

test:=f_compute();
CI_Decision_Reset();
if test=false then
    Print("Error case\n");
else
    Print("Normal case\n");
    CI_Write_Ok();
fi;
