Read("../common.g");
Read("../access_points.g");
Print("Beginning TestReductions\n");

# The reductions of LATT_Reduction are run on scrambled presentations of lattices
# whose good basis is known. For every run two things are checked here, in
# GAP and without sharing any code with the C++ side:
#   * the transformation is integral, unimodular and produces the returned
#     form;
#   * the returned form satisfies the condition the method claims: the Lovasz
#     condition, the deep condition, the BKZ block condition, the two families
#     of slide reduction, Minkowski's condition, the local minimality of
#     Seysen's measure, or for "best" that it is no worse than its input.

delta:=99/100;

#
# The lattices.
#

get_gram_Zn:=function(n)
    return IdentityMat(n);
end;

# G = R R^T with R integral and close to the identity: a form of low symmetry,
# so that the test does not run only on lattices with large automorphism
# groups.
get_gram_random:=function(rs, n)
    local R, i, j;
    while true
    do
        R:=IdentityMat(n);
        for i in [1..n]
        do
            for j in [1..n]
            do
                if i<>j then
                    R[i][j]:=Random(rs, [-1..1]);
                fi;
            od;
        od;
        if DeterminantMat(R)<>0 then
            return R * TransposedMat(R);
        fi;
    od;
end;

# G = R R^T with R integral with entries in [-spread, spread]: far from any
# good basis, so that the projected blocks of BKZ and slide reduction carry
# integers of dozens of digits. On such a form of dimension 16 the
# enumeration used to overflow a double in the seeding of its bounds and
# then never returned.
get_gram_random_wide:=function(rs, n, spread)
    local R;
    while true
    do
        R:=List([1..n], i->List([1..n], j->Random(rs, [-spread..spread])));
        if DeterminantMat(R)<>0 then
            return R * TransposedMat(R);
        fi;
    od;
end;

# A random element of GL_n(Z), as a product of transvections, swaps and sign
# changes. n_ops says how far the good presentation is pushed away.
get_random_unimodular:=function(rs, n, n_ops)
    local U, i_op, kind, i, j, coef;
    U:=IdentityMat(n);
    if n=1 then
        return U;
    fi;
    for i_op in [1..n_ops]
    do
        kind:=Random(rs, [1..10]);
        i:=Random(rs, [1..n]);
        j:=Random(rs, [1..n]);
        if kind <= 8 then
            if i<>j then
                coef:=Random(rs, [-2,-1,1,2]);
                U[i]:=U[i] + coef * U[j];
            fi;
        elif kind=9 then
            if i<>j then
                U{[i,j]}:=U{[j,i]};
            fi;
        else
            U[i]:=-U[i];
        fi;
    od;
    return U;
end;

#
# Gram-Schmidt data of a Gram matrix over the rationals: Bn[i] = |b_i^*|^2
# and mu[i][j] for j < i.
#
get_gram_schmidt:=function(G)
    local n, mu, Bn, i, j, k, val;
    n:=Length(G);
    mu:=NullMat(n, n);
    Bn:=[];
    for i in [1..n]
    do
        for j in [1..i-1]
        do
            val:=G[i][j];
            for k in [1..j-1]
            do
                val:=val - mu[j][k] * mu[i][k] * Bn[k];
            od;
            mu[i][j]:=val / Bn[j];
        od;
        val:=G[i][i];
        for k in [1..i-1]
        do
            val:=val - mu[i][k]^2 * Bn[k];
        od;
        Bn[i]:=val;
    od;
    return rec(mu:=mu, Bn:=Bn);
end;

# |pi_i(b_k)|^2, the norm of b_k projected orthogonally to b_1, ..., b_{i-1}.
get_projected_norm:=function(G, gs, k, i)
    local val, l;
    val:=G[k][k];
    for l in [1..i-1]
    do
        val:=val - gs.mu[k][l]^2 * gs.Bn[l];
    od;
    return val;
end;

# The Gram matrix of pi_j(b_j), ..., pi_j(b_kend): a Schur complement.
get_projected_block:=function(G, j, kend)
    local P, J;
    J:=[j..kend];
    if j=1 then
        return G{J}{J};
    fi;
    P:=[1..j-1];
    return G{J}{J} - G{J}{P} * Inverse(G{P}{P}) * G{P}{J};
end;

# The minimum of a positive definite form, rational entries allowed.
get_lattice_minimum:=function(G)
    local scal, Gs, bound, rec_shv;
    scal:=Lcm(List(Flat(G), DenominatorRat));
    Gs:=scal * G;
    bound:=Minimum(List([1..Length(Gs)], i->Gs[i][i]));
    rec_shv:=ShortestVectors(Gs, bound);
    return Minimum(rec_shv.norms) / scal;
end;

get_orth_defect_sq:=function(G)
    return Product(List([1..Length(G)], i->G[i][i])) / DeterminantMat(G);
end;

get_seysen_measure:=function(G)
    local Ginv;
    Ginv:=Inverse(G);
    return Sum(List([1..Length(G)], i->G[i][i] * Ginv[i][i]));
end;

#
# The conditions. Each returns true or a string describing the failure.
#

check_size_reduced:=function(G)
    local gs, n, k, j;
    gs:=get_gram_schmidt(G);
    n:=Length(G);
    for k in [2..n]
    do
        for j in [1..k-1]
        do
            if AbsInt(gs.mu[k][j]) > 1/2 then
                return Concatenation("not size reduced at k=", String(k), " j=", String(j));
            fi;
        od;
    od;
    return true;
end;

check_lll:=function(G)
    local gs, n, k, test;
    test:=check_size_reduced(G);
    if test<>true then
        return test;
    fi;
    gs:=get_gram_schmidt(G);
    n:=Length(G);
    for k in [2..n]
    do
        if gs.Bn[k] < (delta - gs.mu[k][k-1]^2) * gs.Bn[k-1] then
            return Concatenation("the Lovasz condition fails at k=", String(k));
        fi;
    od;
    return true;
end;

# depth = 0 is unrestricted deep insertion. With depth > 0 the positions
# tried for index k are the first depth ones and the last depth ones before
# k, stated here with 0-based indices as in the C++ code.
check_deep:=function(G, depth)
    local gs, n, k, i, test, is_admissible;
    test:=check_lll(G);
    if test<>true then
        return test;
    fi;
    gs:=get_gram_schmidt(G);
    n:=Length(G);
    is_admissible:=function(i0, k0)
        if depth <= 0 then
            return true;
        fi;
        return i0 < depth or i0 >= k0 - depth;
    end;
    for k in [2..n]
    do
        for i in [1..k-1]
        do
            if is_admissible(i-1, k-1) then
                if get_projected_norm(G, gs, k, i) < delta * gs.Bn[i] then
                    return Concatenation("the deep condition fails at k=", String(k), " i=", String(i));
                fi;
            fi;
        od;
    od;
    return true;
end;

# b_j^* is, up to delta, a shortest vector of the projected block [j, kend].
check_block_condition:=function(G, gs, j, kend)
    local Gblock;
    Gblock:=get_projected_block(G, j, kend);
    if get_lattice_minimum(Gblock) < delta * gs.Bn[j] then
        return Concatenation("the block [", String(j), ",", String(kend), "] has a vector shorter than b_j^*");
    fi;
    return true;
end;

check_bkz:=function(G, beta)
    local gs, n, j, test;
    test:=check_lll(G);
    if test<>true then
        return test;
    fi;
    gs:=get_gram_schmidt(G);
    n:=Length(G);
    for j in [1..n-1]
    do
        test:=check_block_condition(G, gs, j, Minimum(j + beta - 1, n));
        if test<>true then
            return test;
        fi;
    od;
    return true;
end;

# Slide reduction is defined for a block size k >= 2 dividing n; a form of
# dimension at most one is reduced whatever k.
is_slide_applicable:=function(n, k)
    return n <= 1 or (k >= 2 and n mod k = 0);
end;

# Primal: every block [ik+1, (i+1)k] is HKZ reduced. Dual: every shifted
# block [ik+2, (i+1)k+1] has its last Gram-Schmidt norm maximal up to delta,
# that is the dual of the block has no vector shorter, up to delta, than the
# last vector of the dual basis, whose norm is 1/|b_last^*|^2.
check_slide:=function(G, k)
    local gs, n, p, i, j, test, Gblock, Gdual, m;
    test:=check_size_reduced(G);
    if test<>true then
        return test;
    fi;
    gs:=get_gram_schmidt(G);
    n:=Length(G);
    if n <= 1 then
        return true;
    fi;
    p:=n / k;
    for i in [0..p-1]
    do
        for j in [i*k+1..(i+1)*k-1]
        do
            test:=check_block_condition(G, gs, j, (i+1)*k);
            if test<>true then
                return Concatenation("primal: ", test);
            fi;
        od;
    od;
    for i in [0..p-2]
    do
        Gblock:=get_projected_block(G, i*k+2, (i+1)*k+1);
        Gdual:=Inverse(Gblock);
        m:=Length(Gdual);
        if get_lattice_minimum(Gdual) < delta * Gdual[m][m] then
            return Concatenation("dual: the condition fails at boundary ", String(i+1));
        fi;
    od;
    return true;
end;

# At every index i, no vector v with gcd(v_i, ..., v_n) = 1, the condition
# for (b_1, ..., b_{i-1}, v) to extend to a basis, is shorter than b_i.
check_minkowski:=function(G)
    local n, i, rec_shv, iV, eV;
    n:=Length(G);
    for i in [1..n]
    do
        rec_shv:=ShortestVectors(G, G[i][i]);
        for iV in [1..Length(rec_shv.vectors)]
        do
            eV:=rec_shv.vectors[iV];
            if rec_shv.norms[iV] < G[i][i] and Gcd(eV{[i..n]})=1 then
                return Concatenation("an admissible vector is shorter than b_", String(i));
            fi;
        od;
    od;
    return true;
end;

# Seysen's descent stops when no transvection b_i <- b_i + lambda b_j lowers
# the measure. The change of the measure is a convex quadratic in lambda that
# vanishes at 0, so it suffices to test lambda = 1 and lambda = -1.
check_seysen:=function(G)
    local n, Ginv, i, j, lambda, gain, den, num;
    n:=Length(G);
    Ginv:=Inverse(G);
    for i in [1..n]
    do
        for j in [1..n]
        do
            if i<>j then
                den:=2 * G[j][j] * Ginv[i][i];
                num:=G[j][j] * Ginv[i][j] - G[i][j] * Ginv[i][i];
                for lambda in [-1,1]
                do
                    gain:=lambda * (lambda * den - 2 * num);
                    if gain < 0 then
                        return Concatenation("the move i=", String(i), " j=", String(j), " lambda=", String(lambda), " lowers the Seysen measure");
                    fi;
                od;
            fi;
        od;
    od;
    return true;
end;

#
# The methods and the condition each one claims.
#
get_method_check:=function(method)
    local prefix, rest;
    if method="direct" then
        return check_lll;
    fi;
    if method="seysen" or method="seysen_best" then
        return check_seysen;
    fi;
    if method="deep" then
        return G->check_deep(G, 0);
    fi;
    if method="minkowski" then
        return check_minkowski;
    fi;
    for prefix in ["deep-", "bkz-", "slide-"]
    do
        rest:=starts_with(method, prefix);
        if rest<>fail then
            if prefix="deep-" then
                return G->check_deep(G, Int(rest));
            fi;
            if prefix="bkz-" then
                return G->check_bkz(G, Int(rest));
            fi;
            return G->check_slide(G, Int(rest));
        fi;
    od;
    # dual, seysen_lll and best claim no condition of their own; what they
    # guarantee is checked in test_reduction.
    return G->true;
end;

ListMethod:=["direct", "dual", "seysen", "seysen_best", "seysen_lll",
             "deep", "deep-3", "bkz-2", "bkz-4", "bkz-8",
             "slide-2", "slide-3", "slide-4", "minkowski", "best"];

# The block size of a slide method, or fail for another method.
get_slide_parameter:=function(method)
    local rest;
    rest:=starts_with(method, "slide-");
    if rest=fail then
        return fail;
    fi;
    return Int(rest);
end;

test_reduction:=function(eMat, method)
    local res, Pmat, test, input_defect, direct, k;
    # A slide method whose block size does not divide the dimension is not
    # defined, and must be rejected rather than replaced by something else.
    k:=get_slide_parameter(method);
    if k<>fail and is_slide_applicable(Length(eMat), k)=false then
        res:=get_lattice_reduction(eMat, method);
        if is_error(res)=false then
            Print("  method=", method, " should have been rejected in dimension ", Length(eMat), "\n");
            return false;
        fi;
        return true;
    fi;
    res:=get_lattice_reduction(eMat, method);
    if is_error(res) then
        return false;
    fi;
    Pmat:=res.Pmat;
    if ForAll(Flat(Pmat), IsInt)=false then
        Print("  method=", method, " the transformation is not integral\n");
        return false;
    fi;
    if AbsInt(DeterminantMat(Pmat))<>1 then
        Print("  method=", method, " the transformation is not unimodular\n");
        return false;
    fi;
    if Pmat * eMat * TransposedMat(Pmat)<>res.GramMat then
        Print("  method=", method, " the transformation does not produce the returned form\n");
        return false;
    fi;
    test:=get_method_check(method)(res.GramMat);
    if test<>true then
        Print("  method=", method, " ", test, "\n");
        return false;
    fi;
    if method="seysen_lll" then
        if get_seysen_measure(res.GramMat) > get_seysen_measure(eMat) then
            Print("  method=seysen_lll the Seysen measure went up\n");
            return false;
        fi;
    fi;
    if method="best" then
        # The input and every candidate are in the search, so the result is
        # no worse than either by the criterion of the search.
        input_defect:=get_orth_defect_sq(eMat);
        if get_orth_defect_sq(res.GramMat) > input_defect then
            Print("  method=best the defect went up\n");
            return false;
        fi;
        direct:=get_lattice_reduction(eMat, "direct");
        if is_error(direct) then
            return false;
        fi;
        if get_orth_defect_sq(res.GramMat) > get_orth_defect_sq(direct.GramMat) then
            Print("  method=best is worse than direct\n");
            return false;
        fi;
        if res.method<>"none" and Position(ListMethod, res.method)=fail then
            Print("  method=best returned the unknown method ", res.method, "\n");
            return false;
        fi;
        k:=get_slide_parameter(res.method);
        if k<>fail and is_slide_applicable(Length(eMat), k)=false then
            Print("  method=best returned ", res.method, " which does not apply in dimension ", Length(eMat), "\n");
            return false;
        fi;
    fi;
    return true;
end;

#
# The names the parser must reject: no hyphen, a parameter below the
# meaningful minimum, a missing or malformed parameter, an unknown name.
#
ListRejected:=["bkz8", "deep5", "deep-0", "bkz-1", "slide-1", "bkz-", "bkz-4x", "bkz--4", "foo"];

test_rejected:=function(method)
    local res;
    res:=get_lattice_reduction(ClassicalSporadicLattices("A4"), method);
    if is_error(res)=false then
        Print("  the method name ", method, " should have been rejected\n");
        return false;
    fi;
    return true;
end;

#
# The instances.
#
rs:=RandomSource(IsMersenneTwister, 16);

ListCase:=[];
Add(ListCase, rec(name:="dim1", G:=[[5]]));
Add(ListCase, rec(name:="Z5", G:=get_gram_Zn(5)));
# A11 has no block size below 11 dividing its dimension, so every slide
# method of the list must be rejected on it.
for eName in ["A4", "A7", "D4", "D6", "E6", "E7", "E8", "A11"]
do
    Add(ListCase, rec(name:=eName, G:=ClassicalSporadicLattices(eName)));
od;
Add(ListCase, rec(name:="random6", G:=get_gram_random(rs, 6)));
Add(ListCase, rec(name:="random9", G:=get_gram_random(rs, 9)));
Add(ListCase, rec(name:="wide16", G:=get_gram_random_wide(rs, 16, 10)));

ListInstance:=[];
for eCase in ListCase
do
    n:=Length(eCase.G);
    for n_ops in [20, 80]
    do
        U:=get_random_unimodular(rs, n, n_ops);
        Add(ListInstance, rec(name:=Concatenation(eCase.name, "/", String(n_ops)), G:=U * eCase.G * TransposedMat(U)));
    od;
od;

FullTest:=function()
    local n_error, eInst, method, reply, iInst;
    n_error:=0;
    iInst:=0;
    for eInst in ListInstance
    do
        iInst:=iInst + 1;
        Print("iInst=", iInst, " / ", Length(ListInstance), " name=", eInst.name, " n_error=", n_error, "\n");
        for method in ListMethod
        do
            reply:=test_reduction(eInst.G, method);
            if reply=false then
                Print("  FAILURE name=", eInst.name, " method=", method, "\n");
                n_error:=n_error + 1;
            fi;
        od;
    od;
    for method in ListRejected
    do
        reply:=test_rejected(method);
        if reply=false then
            n_error:=n_error + 1;
        fi;
    od;
    return n_error;
end;

n_error:=FullTest();
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
