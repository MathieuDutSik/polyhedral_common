Read("../common.g");
Read("../access_points.g");
Print("Beginning TestReductions\n");

# The reductions of src_isotropy, each result being checked in GAP:
# the transformation is unimodular (or rational for the determinant
# minimization), it maps the input to the returned form, and the
# condition the method claims holds.
#
# The inputs are the indefinite Gram matrices of LGramReflect (one in
# Stride), the root lattices of ListPositiveDefinite, and forms built
# from them, all scrambled by random unimodular matrices.

Reset(GlobalMersenneTwister, 1);;
Reset(GlobalRandomSource, 1);;
Stride:=8;;

GetSqrMat:=function(eV)
    local dim;
    dim:=Sqrt(Length(eV));
    return List([1..dim], i -> eV{[(i-1)*dim+1..i*dim]});
end;

Scramble:=function(M)
    local P;
    P:=RandomIntegralUnimodularMatrix(Length(M));
    return P * M * TransposedMat(P);
end;

IsUnimodular:=function(B)
    return ForAll(Flat(B), IsInt) and AbsInt(DeterminantMat(B))=1;
end;

ListIndefinite:=List(ReadAsFunction("LGramReflect")(), GetSqrMat);;
ListIndefinite:=List(ListIndefinite{Filtered([1..Length(ListIndefinite)], i -> i mod Stride = 1)}, Scramble);;
ListPositiveDefinite:=List(ReadAsFunction("ListPositiveDefinite")(), Scramble);;

# The Gram-Schmidt data of the standard basis for the form M: the norms
# of the b*_i and the coefficients mu. fail if a b*_i is isotropic.
GramSchmidtForm:=function(M)
    local n, Bstar, Norms, mu, i, j, v;
    n:=Length(M);
    Bstar:=[];
    Norms:=[];
    mu:=NullMat(n, n);
    for i in [1..n] do
        v:=IdentityMat(n)[i];
        for j in [1..i-1] do
            mu[i][j]:=(IdentityMat(n)[i] * M * Bstar[j]) / Norms[j];
            v:=v - mu[i][j] * Bstar[j];
        od;
        if v * M * v = 0 then
            return fail;
        fi;
        Add(Bstar, v);
        Add(Norms, v * M * v);
    od;
    return rec(Norms:=Norms, mu:=mu);
end;

# The LLL conditions on the form Mred, with absolute values of the norms
# so that they apply to Simon's indefinite LLL as well.
IsLLLReduced:=function(Mred, c)
    local GS, n, i, j;
    GS:=GramSchmidtForm(Mred);
    if GS=fail then
        Print("  isotropic Gram-Schmidt vector\n");
        return false;
    fi;
    n:=Length(Mred);
    for i in [2..n] do
        for j in [1..i-1] do
            if AbsoluteValue(GS.mu[i][j]) > 1/2 then
                Print("  not size reduced: mu[", i, "][", j, "]=", GS.mu[i][j], "\n");
                return false;
            fi;
        od;
        if AbsoluteValue(GS.Norms[i] + GS.mu[i][i-1]^2 * GS.Norms[i-1]) < c * AbsoluteValue(GS.Norms[i-1]) then
            Print("  Lovasz condition fails at k=", i, "\n");
            return false;
        fi;
    od;
    return true;
end;

CheckTransform:=function(name, M, res)
    if is_error(res) then
        Print(name, ": ", res, "\n");
        return false;
    fi;
    if IsUnimodular(res.B)=false then
        Print(name, ": B is not unimodular\n");
        return false;
    fi;
    if res.B * M * TransposedMat(res.B) <> res.Mred then
        Print(name, ": Mred <> B M B^T\n");
        return false;
    fi;
    return true;
end;

TestIndefiniteLLL:=function(M)
    local res, X;
    res:=indefinite_lll(M);
    if CheckTransform("IndefiniteLLL", M, res)=false then
        return false;
    fi;
    if IsBound(res.Xisotrop) then
        X:=res.Xisotrop;
        if IsZero(X) or X * M * X <> 0 then
            Print("IndefiniteLLL: Xisotrop is not a nonzero isotropic vector\n");
            return false;
        fi;
        return true;
    fi;
    return IsLLLReduced(res.Mred, 7/8);
end;

TestIndefiniteReduction:=function(M)
    return CheckTransform("IndefiniteReduction", M, indefinite_reduction(M));
end;

TestLLLreduceBasis:=function(M)
    local res;
    res:=lll_reduce_basis(M);
    if CheckTransform("LLLreduceBasis", M, res)=false then
        return false;
    fi;
    return IsLLLReduced(res.Mred, 99/100);
end;

# A determinant minimization of Q: Mred integral, P Q P^T a nonzero
# multiple of Mred, and |det Mred| <= |det Q| when require_decrease. The
# equivariant reduction of LATT_AutomorphyReduction can increase the
# determinant (for example -49 -> -343), which is only reported.
CheckDetMin:=function(name, Q, res, require_decrease)
    local Prod, pos, lambda;
    if ForAll(Flat(res.Mred), IsInt)=false then
        Print(name, ": Mred is not integral\n");
        return false;
    fi;
    if res.det1 <> DeterminantMat(Q) or res.det2 <> DeterminantMat(res.Mred) then
        Print(name, ": inconsistent det1/det2\n");
        return false;
    fi;
    if AbsoluteValue(res.det2) > AbsoluteValue(res.det1) then
        if require_decrease then
            Print(name, ": the determinant increased\n");
            return false;
        fi;
        Print(name, ": note, the determinant increased from ", res.det1, " to ", res.det2, "\n");
    fi;
    Prod:=res.P * Q * TransposedMat(res.P);
    pos:=PositionProperty(Flat(res.Mred), x -> x <> 0);
    lambda:=Flat(Prod)[pos] / Flat(res.Mred)[pos];
    if lambda = 0 or Prod <> lambda * res.Mred then
        Print(name, ": P Q P^T is not a multiple of Mred\n");
        return false;
    fi;
    return true;
end;

TestDetMinimization:=function(M)
    local res;
    res:=determinant_minimization(M);
    if is_error(res) then
        Print("DetMinimization: ", res, "\n");
        return false;
    fi;
    return CheckDetMin("DetMinimization", M, res, true);
end;

TestAutomorphyReduction:=function(M)
    local list_rec, Minv;
    list_rec:=automorphy_reduction(M);
    if is_error(list_rec) then
        Print("AutomorphyReduction: ", list_rec, "\n");
        return false;
    fi;
    if Length(list_rec) <> 2 then
        Print("AutomorphyReduction: ", Length(list_rec), " records instead of 2\n");
        return false;
    fi;
    Minv:=Inverse(M);
    Minv:=Minv * Lcm(List(Flat(Minv), DenominatorRat));
    Minv:=Minv / Gcd(Flat(Minv));
    return CheckDetMin("AutomorphyReduction(M)", M, list_rec[1], false)
       and CheckDetMin("AutomorphyReduction(M^-1)", Minv, list_rec[2], false);
end;

TestPositiveSemiDefinite:=function(M, expected)
    local res;
    res:=check_positive_semidefinite(M);
    if res <> expected then
        Print("CheckPositiveSemiDefinite: got ", res, " expected ", expected, "\n");
        return false;
    fi;
    return true;
end;

# A basis of a sublattice, given by generators: the reduced basis spans
# the same lattice and has no larger L1 norm.
TestSublatticeBasisReduction:=function(M)
    local res, L1;
    res:=sublattice_basis_reduction(M);
    if is_error(res) then
        Print("SublatticeBasisReduction: ", res, "\n");
        return false;
    fi;
    if HermiteNormalFormIntegerMat(res.Mred) <> HermiteNormalFormIntegerMat(M) then
        Print("SublatticeBasisReduction: the lattice changed\n");
        return false;
    fi;
    L1:=X -> Sum(List(Flat(X), AbsoluteValue));
    if L1(res.Mred) > L1(M) then
        Print("SublatticeBasisReduction: the L1 norm increased\n");
        return false;
    fi;
    return true;
end;

# The guarantee of Lenstra, Lenstra, Lovasz (Proposition 1.39):
# an integer q >= 1 with |q alpha_i - p_i| <= epsilon for all i. The
# program works in double precision, hence the 1/10^6 of tolerance.
TestDiophantine:=function(ListStr, epsilon_str)
    local res, alpha, eps, q, i;
    res:=diophantine_approximation(ListStr, epsilon_str);
    if is_error(res) then
        Print("DiophantApprox: ", res, "\n");
        return false;
    fi;
    if Length(Set(res.Denominators)) <> 1 then
        Print("DiophantApprox: the denominators differ\n");
        return false;
    fi;
    q:=res.Denominators[1];
    alpha:=List(ListStr, Rat);
    eps:=Rat(epsilon_str);
    if q < 1 then
        Print("DiophantApprox: q=", q, "\n");
        return false;
    fi;
    for i in [1..Length(alpha)] do
        if AbsoluteValue(q * alpha[i] - res.Numerators[i]) > eps + 1/10^6 then
            Print("DiophantApprox: |q alpha - p| > epsilon for ", ListStr, " epsilon=", epsilon_str, " q=", q, "\n");
            return false;
        fi;
    od;
    return true;
end;

FullTest:=function()
    local n_error, M, n, k, R, Mpsd, Z, eCase;
    n_error:=0;
    for M in ListIndefinite do
        Print("Indefinite n=", Length(M), " det=", DeterminantMat(M), "\n");
        if TestIndefiniteLLL(M)=false then n_error:=n_error+1; fi;
        if TestIndefiniteReduction(M)=false then n_error:=n_error+1; fi;
        if TestDetMinimization(M)=false then n_error:=n_error+1; fi;
        if TestAutomorphyReduction(M)=false then n_error:=n_error+1; fi;
        if TestPositiveSemiDefinite(M, "not")=false then n_error:=n_error+1; fi;
    od;
    for M in ListPositiveDefinite do
        n:=Length(M);
        Print("Positive definite n=", n, " det=", DeterminantMat(M), "\n");
        if TestLLLreduceBasis(M)=false then n_error:=n_error+1; fi;
        if TestIndefiniteLLL(M)=false then n_error:=n_error+1; fi;
        if TestPositiveSemiDefinite(M, "definite")=false then n_error:=n_error+1; fi;
        # A positive semidefinite form: M plus a zero block, scrambled.
        Z:=NullMat(n+1, n+1);
        Z{[1..n]}{[1..n]}:=M;
        Mpsd:=Scramble(Z);
        if TestPositiveSemiDefinite(Mpsd, "semidefinite")=false then n_error:=n_error+1; fi;
        # A sublattice of rank k < n given by a skewed basis.
        k:=n - 1;
        R:=RandomIntegralUnimodularMatrix(k) * RandomIntegralUnimodularMatrix(n){[1..k]};
        if TestSublatticeBasisReduction(R)=false then n_error:=n_error+1; fi;
    od;
    # LATT_DiophantApprox works in double precision: the term
    # eps^(2n+2) / 2^(n(n+1)/2) of its Gram matrix must stay above the
    # rounding of the sum of the alpha_i^2, which limits eps from below.
    for eCase in [ [["1.41421356237", "1.73205080757"], "0.01"],
                   [["1.41421356237", "1.73205080757", "3.14159265359"], "0.1"],
                   [["0.61803398875", "2.71828182846", "0.57721566490", "1.20205690316"], "0.1"] ] do
        Print("Diophantine approximation of ", eCase[1], " epsilon=", eCase[2], "\n");
        if TestDiophantine(eCase[1], eCase[2])=false then n_error:=n_error+1; fi;
    od;
    Print("n_error=", n_error, "\n");
    return n_error;
end;

n_error:=FullTest();;
CI_Decision_Reset();
if n_error > 0 then
    # Error case
    Print("Error case\n");
else
    # No error case
    Print("Normal case\n");
    CI_Write_Ok();
fi;
