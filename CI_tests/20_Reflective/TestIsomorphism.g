Read("../common.g");
Read("../access_points.g");

# The edgewalk programs are built with ENABLE_FLINT_SUPPORT=1 and run on
# flint, the faster arithmetic for these computations.
ReflectiveOptions:=rec(arith:="flint");;
Print("Beginning TestIsomorphism\n");

# The isometry test of LORENTZ_FundDomain_AllcockEdgewalk_Isomorphism on
# the reflective lattices of ListReflect:
# * a lattice and a scrambled copy: an isometry is expected, and checked
#   to be integral, unimodular and to map the first form to the second;
# * two lattices of the same dimension and determinant, one even and one
#   odd, and two lattices of different determinants: fail is expected,
#   the parity or the determinant certifying that they are not isometric.

Reset(GlobalMersenneTwister, 1);;
Reset(GlobalRandomSource, 1);;
NbPerDim:=5;;
ListDim:=[3 .. 8];;

ListRec:=ReadAsFunction("ListReflect")();;

IsEven:=G -> ForAll([1..Length(G)], i -> IsEvenInt(G[i][i]));;

TestIsometric:=function(G1)
    local P, G2, U;
    P:=RandomIntegralUnimodularMatrix(Length(G1));
    G2:=P * G1 * TransposedMat(P);
    U:=lorentzian_isomorphism(G1, G2, ReflectiveOptions);
    if is_error(U) then
        Print("  ", U, "\n");
        return false;
    fi;
    if U = fail then
        Print("  no isometry found between isometric lattices G1=", G1, "\n");
        return false;
    fi;
    if ForAll(Flat(U), IsInt)=false or AbsInt(DeterminantMat(U)) <> 1 then
        Print("  the isometry is not integral unimodular\n");
        return false;
    fi;
    if U * G1 * TransposedMat(U) <> G2 then
        Print("  U G1 U^T <> G2 for G1=", G1, "\n");
        return false;
    fi;
    return true;
end;

TestNonIsometric:=function(G1, G2)
    local U;
    U:=lorentzian_isomorphism(G1, G2, ReflectiveOptions);
    if is_error(U) then
        Print("  ", U, "\n");
        return false;
    fi;
    if U <> fail then
        Print("  an isometry returned for non isometric lattices G1=", G1, " G2=", G2, "\n");
        return false;
    fi;
    return true;
end;

FullTest:=function()
    local n_error, dim, ListG, G, i, j, pair_found, Ga, Gb;
    n_error:=0;
    for dim in ListDim do
        ListG:=List(Filtered(ListRec, x -> Length(x.LorMat)=dim), x -> x.LorMat);
        for G in ListG{[1..Minimum(NbPerDim, Length(ListG))]} do
            Print("dim=", dim, " isometric pair, det=", DeterminantMat(G), "\n");
            if TestIsometric(G)=false then
                n_error:=n_error + 1;
            fi;
        od;
        # Same determinant, different parity.
        pair_found:=false;
        for i in [1..Length(ListG)] do
            for j in [i+1..Length(ListG)] do
                if pair_found=false and DeterminantMat(ListG[i])=DeterminantMat(ListG[j]) and IsEven(ListG[i]) <> IsEven(ListG[j]) then
                    pair_found:=true;
                    Print("dim=", dim, " non isometric pair of det=", DeterminantMat(ListG[i]), " (parity)\n");
                    if TestNonIsometric(ListG[i], ListG[j])=false then
                        n_error:=n_error + 1;
                    fi;
                fi;
            od;
        od;
        # Different determinants.
        Ga:=ListG[1];
        Gb:=First(ListG, x -> AbsInt(DeterminantMat(x)) <> AbsInt(DeterminantMat(Ga)));
        if Gb <> fail then
            Print("dim=", dim, " non isometric pair of dets ", DeterminantMat(Ga), " and ", DeterminantMat(Gb), "\n");
            if TestNonIsometric(Ga, Gb)=false then
                n_error:=n_error + 1;
            fi;
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
