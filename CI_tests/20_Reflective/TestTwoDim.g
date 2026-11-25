Read("../common.g");
Read("../access_points.g");
Print("Beginning TestTwoDim\n");

# The algorithms of Section 8 of Allcock, "An alternative to Vinberg's
# algorithm", for the 2-dimensional Lorentzian forms a x^2 + 2 b x y + c y^2,
# by the program LORENTZ_TwoDimAnisotropic_Allcock. Each result is checked
# by enumerating the lattice vectors of a bounded region.
#
# TwoDimAnisotropicCases (a b c M, anisotropic forms): the start vector r0,
# the generator g of SO^+(L) and the list R of Algorithm 8.7. g preserves
# the form and is not +-1, and R is exactly the set of the primitive vectors
# of norm <= M after r0 and before g(R[1]), in ascending order. If there is
# no result, no vector of norm in ]0,M] lies in a box.
#
# TwoDimFirstNextCases (a b c x y N): the first vector v of norm N after
# r0 = (x,y) in its sector, for isotropic and anisotropic forms, multiples
# included. No vector of norm N lies strictly between r0 and v. If there is
# no result, none lies after r0 in a box.

GetGram:=c -> [[c[1], c[2]], [c[2], c[3]]];
Det2:=function(u, v) return u[1]*v[2] - u[2]*v[1]; end;
NoneBox:=60;

# A box containing the vectors of norm <= M between the rays of u and w:
# they lie in the triangle of 0 and the points of norm M of the two rays.
GetBox:=function(G, u, w, M)
    local B, v;
    B:=0;
    for v in [u, w] do
        B:=Maximum(B, Maximum(AbsInt(v[1]), AbsInt(v[2])) * (RootInt(Int(M / (v*G*v)) + 1) + 1));
    od;
    return B;
end;

# The vectors x of the box with f(x) true.
BoxVectors:=function(B, f)
    local L, x1, x2;
    L:=[];
    for x1 in [-B..B] do
        for x2 in [-B..B] do
            if f([x1, x2]) then
                Add(L, [x1, x2]);
            fi;
        od;
    od;
    return L;
end;

TestAnisotropic:=function(eCase, res)
    local G, M, r0, g, R, n, i, nextR, B, ListRef;
    G:=GetGram(eCase);
    M:=eCase[4];
    r0:=res.r0;
    if r0*G*r0 <= 0 then
        Print("  r0 is not of positive norm\n");
        return false;
    fi;
    if not IsBound(res.g) then
        ListRef:=BoxVectors(NoneBox, x -> x*G*x > 0 and x*G*x <= M);
        if Length(ListRef) > 0 then
            Print("  no result, but ", ListRef[1], " has norm ", ListRef[1]*G*ListRef[1], " <= M\n");
            return false;
        fi;
        return true;
    fi;
    g:=res.g;
    R:=res.R;
    if ForAll(Flat(g), IsInt)=false or DeterminantMat(g) <> 1 or g*G*TransposedMat(g) <> G then
        Print("  g is not an element of SO(L)\n");
        return false;
    fi;
    if g = IdentityMat(2) or g = -IdentityMat(2) then
        Print("  g is +-1\n");
        return false;
    fi;
    n:=Length(R);
    for i in [1..n] do
        if Gcd(R[i]) <> 1 or R[i]*G*R[i] <= 0 or R[i]*G*R[i] > M or R[i]*G*r0 <= 0 or Det2(r0, R[i]) <= 0 then
            Print("  R[", i, "]=", R[i], " is not a primitive vector of S after r0 of norm <= M\n");
            return false;
        fi;
        if i < n and Det2(R[i], R[i+1]) <= 0 then
            Print("  R is not in ascending order\n");
            return false;
        fi;
    od;
    nextR:=R[1]*g;
    if Det2(R[n], nextR) <= 0 then
        Print("  g(R[1]) does not come after R\n");
        return false;
    fi;
    B:=GetBox(G, r0, nextR, M);
    ListRef:=BoxVectors(B, x -> Gcd(x)=1 and x*G*x > 0 and x*G*x <= M and x*G*r0 > 0 and Det2(r0, x) > 0 and Det2(x, nextR) > 0);
    if Set(ListRef) <> Set(R) then
        Print("  R is not the set of the vectors before g(R[1]): ", Length(ListRef), " found by enumeration, ", n, " in R\n");
        return false;
    fi;
    return true;
end;

TestFirstNext:=function(eCase, res)
    local G, r0, N, v, B, ListRef;
    G:=GetGram(eCase);
    r0:=[eCase[4], eCase[5]];
    N:=eCase[6];
    if not IsBound(res.v) then
        ListRef:=BoxVectors(NoneBox, x -> x*G*x = N and x*G*r0 > 0 and Det2(r0, x) > 0);
        if Length(ListRef) > 0 then
            Print("  no result, but ", ListRef[1], " is a vector of norm N after r0\n");
            return false;
        fi;
        return true;
    fi;
    v:=res.v;
    if v*G*v <> N or v*G*r0 <= 0 or Det2(r0, v) <= 0 then
        Print("  v=", v, " is not a vector of norm N after r0 in its sector\n");
        return false;
    fi;
    B:=GetBox(G, r0, v, N);
    ListRef:=BoxVectors(B, x -> x*G*x = N and x*G*r0 > 0 and Det2(r0, x) > 0 and Det2(x, v) > 0);
    if Length(ListRef) > 0 then
        Print("  ", ListRef[1], " of norm N comes before v=", v, "\n");
        return false;
    fi;
    return true;
end;

ReadCases:=function(FileCases)
    local L;
    L:=List(SplitString(StringFile(FileCases), "\n"), x -> List(SplitString(x, " "), Int));
    L:=Filtered(L, x -> Length(x) > 1);
    return L;
end;

FullTest:=function()
    local n_error, ListCase, ListRes, i;
    n_error:=0;
    ListCase:=ReadCases("TwoDimAnisotropicCases");
    ListRes:=two_dim_lorentzian_cases("anisotropic", "TwoDimAnisotropicCases");
    if is_error(ListRes) then
        Print(ListRes, "\n");
        return 1;
    fi;
    for i in [1..Length(ListCase)] do
        if TestAnisotropic(ListCase[i], ListRes[i])=false then
            Print("anisotropic case ", ListCase[i], " is wrong\n");
            n_error:=n_error + 1;
        fi;
    od;
    Print("anisotropic: ", Length(ListCase), " cases\n");
    ListCase:=ReadCases("TwoDimFirstNextCases");
    ListRes:=two_dim_lorentzian_cases("first_next", "TwoDimFirstNextCases");
    if is_error(ListRes) then
        Print(ListRes, "\n");
        return n_error + 1;
    fi;
    for i in [1..Length(ListCase)] do
        if TestFirstNext(ListCase[i], ListRes[i])=false then
            Print("first_next case ", ListCase[i], " is wrong\n");
            n_error:=n_error + 1;
        fi;
    od;
    Print("first_next: ", Length(ListCase), " cases\n");
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
