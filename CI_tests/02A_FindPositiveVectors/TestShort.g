Read("../common.g");
Read("../access_points.g");
Print("Beginning TestShort\n");

# The programs of src_short on the root lattices of ListPositiveDefinite,
# scrambled by random unimodular matrices, each result checked in GAP:
# * SHORT_GetShortestVector: the minimum 2 and the shortest vectors of
#   GAP, up to sign;
# * SHORT_AutomorphismGroup: generators integral, unimodular, preserving
#   the configuration of roots, and generating a group of the order of the
#   automorphism group of the root system (RootLatticeAutomorphismOrders);
# * SHORT_ReduceVectorFamilyGAP: transformed copies of a configuration,
#   rows permuted and signs changed, reduced to one per configuration;
# * SHORT_SplitVectorFamily: the configurations of index i mod N.

Reset(GlobalMersenneTwister, 1);;
Reset(GlobalRandomSource, 1);;

ListGram:=ReadAsFunction("ListPositiveDefinite")();;
ListOrder:=ReadAsFunction("RootLatticeAutomorphismOrders")();;

Scramble:=function(M)
    local P;
    P:=RandomIntegralUnimodularMatrix(Length(M));
    return P * M * TransposedMat(P);
end;

# The set of the vectors +-v for the rows v of SHV.
SymmetricSet:=SHV -> Set(Concatenation(SHV, -SHV));

# A copy of the configuration SHV under a random element of GL_n(Z), its
# rows permuted and some of their signs changed.
TransformConfiguration:=function(SHV)
    local P, M, perm;
    P:=RandomIntegralUnimodularMatrix(Length(SHV[1]));
    M:=SHV * P;
    perm:=Random(SymmetricGroup(Length(M)));
    M:=Permuted(M, perm);
    return List(M, v -> Random([1, -1]) * v);
end;

TestShortestVector:=function(G)
    local res, ref;
    res:=short_get_shortest_vector(G);
    if is_error(res) then
        Print("  GetShortestVector: ", res, "\n");
        return fail;
    fi;
    if res.min <> 2 then
        Print("  GetShortestVector: min=", res.min, " instead of 2\n");
        return fail;
    fi;
    ref:=ShortestVectors(G, 2).vectors;
    if SymmetricSet(res.SHV) <> SymmetricSet(ref) then
        Print("  GetShortestVector: the vectors differ from those of GAP\n");
        return fail;
    fi;
    return res.SHV;
end;

TestAutomorphismGroup:=function(SHV, expected_order)
    local ListGen, eSet, g, order;
    ListGen:=short_automorphism_group(SHV);
    if is_error(ListGen) then
        Print("  AutomorphismGroup: ", ListGen, "\n");
        return false;
    fi;
    eSet:=SymmetricSet(SHV);
    for g in ListGen do
        if ForAll(Flat(g), IsInt)=false or AbsInt(DeterminantMat(g)) <> 1 then
            Print("  AutomorphismGroup: a generator is not integral unimodular\n");
            return false;
        fi;
        if SymmetricSet(SHV * g) <> eSet then
            Print("  AutomorphismGroup: a generator does not preserve the configuration\n");
            return false;
        fi;
    od;
    if Length(ListGen)=0 then
        order:=1;
    else
        order:=Size(Group(ListGen));
    fi;
    if order <> expected_order then
        Print("  AutomorphismGroup: order ", order, " instead of ", expected_order, "\n");
        return false;
    fi;
    return true;
end;

FullTest:=function()
    local n_error, ListBase, i, G, SHV, ListSHV, eSHV, ListRed, N, r, ListSplit;
    n_error:=0;
    ListBase:=[];
    for i in [1..Length(ListGram)] do
        G:=Scramble(ListGram[i]);
        Print("Root lattice ", i, " n=", Length(G), " det=", DeterminantMat(G), "\n");
        SHV:=TestShortestVector(G);
        if SHV=fail then
            n_error:=n_error + 1;
        else
            Add(ListBase, SHV);
            if TestAutomorphismGroup(SHV, ListOrder[i])=false then
                n_error:=n_error + 1;
            fi;
        fi;
    od;
    # Three transformed copies of the configurations of A3, A4, D4 and D5:
    # four classes.
    ListSHV:=[];
    for eSHV in ListBase{[2, 3, 6, 7]} do
        Append(ListSHV, [eSHV, TransformConfiguration(eSHV), TransformConfiguration(eSHV)]);
    od;
    Print("ReduceVectorFamily on ", Length(ListSHV), " configurations\n");
    ListRed:=short_reduce_vector_family(ListSHV);
    if is_error(ListRed) then
        Print("  ReduceVectorFamily: ", ListRed, "\n");
        n_error:=n_error + 1;
    elif Length(ListRed) <> 4 then
        Print("  ReduceVectorFamily: ", Length(ListRed), " classes instead of 4\n");
        n_error:=n_error + 1;
    fi;
    N:=3;
    for r in [0..N-1] do
        ListSplit:=short_split_vector_family(ListSHV, N, r);
        if is_error(ListSplit) or ListSplit <> ListSHV{Filtered([1..Length(ListSHV)], i -> (i-1) mod N = r)} then
            Print("  SplitVectorFamily: wrong block for residue ", r, "\n");
            n_error:=n_error + 1;
        fi;
    od;
    return n_error;
end;

n_error:=FullTest();;
Print("n_error=", n_error, "\n");

CI_Decision_Reset();
if n_error > 0 then
    # Error case
    Print("Error case\n");
else
    # No error case
    Print("Normal case\n");
    CI_Write_Ok();
fi;
