Read("../common.g");
Read("../access_points.g");
Print("Beginning TestLocalDensity\n");

# The local density code of src_genera (GENUS_LocalDensity) is checked here
# piece by piece, each against something computed independently of it: GAP's
# own Bernoulli numbers, the definition of the generalized Bernoulli numbers
# evaluated in GAP, class numbers, the orders of the orthogonal groups over
# F_p, and reference values computed with Sage (sage/make_density_references.py)
# from Sage's own special values.

# The arithmetics the binary was built with, read off its usage text: flint
# is there only in a build with ENABLE_FLINT_SUPPORT.
available_arith:=function()
    local TmpDir, FileE, eProg, text, ListArith;
    TmpDir:=DirectoryTemporary();
    FileE:=Filename(TmpDir, "usage.err");
    eProg:=GetBinaryFilename("GENUS_LocalDensity");
    run_binary(eProg, [], FileE);
    text:=StringFile(FileE);
    ListArith:=["gmp"];
    if PositionSublist(text, "FLINT") <> fail then
        Add(ListArith, "flint");
    fi;
    RemoveFile(FileE);
    return ListArith;
end;

ListArith:=available_arith();
Print("Testing the arithmetics ", ListArith, "\n");

#
# Independent GAP versions of the quantities.
#

# The Kronecker symbol (a / n) for n > 0.
kronecker_symbol:=function(a, n)
    local res, r;
    res:=1;
    while n mod 2 = 0
    do
        n:=n / 2;
        if a mod 2 = 0 then
            return 0;
        fi;
        r:=a mod 8;
        if r = 3 or r = 5 then
            res:=-res;
        fi;
    od;
    if n = 1 then
        return res;
    fi;
    return res * Jacobi(a, n);
end;

is_fundamental_discriminant:=function(D)
    local d;
    if D = 1 or D = 0 then
        return false;
    fi;
    if D mod 4 = 1 then
        return ForAll(Collected(FactorsInt(AbsInt(D))), x->x[2]=1);
    fi;
    if D mod 4 <> 0 then
        return false;
    fi;
    d:=D / 4;
    if d mod 4 <> 2 and d mod 4 <> 3 then
        return false;
    fi;
    return ForAll(Collected(FactorsInt(AbsInt(d))), x->x[2]=1);
end;

bernoulli_polynomial_value:=function(k, x)
    return Sum([0..k], j->Binomial(k, j) * Bernoulli(j) * x^(k-j));
end;

# B_{k, chi_D} = f^{k-1} sum_{a=1}^{f} chi(a) B_k(a/f), straight from the
# definition.
generalized_bernoulli:=function(k, D)
    local f;
    f:=AbsInt(D);
    return f^(k-1) * Sum([1..f], a->kronecker_symbol(D, a) * bernoulli_polynomial_value(k, a/f));
end;

# |O(V)| for a non-degenerate quadratic space of dimension dim and
# determinant det over F_p, p odd.
orthogonal_group_order:=function(dim, det, p)
    local r, eps;
    if dim = 0 then
        return 1;
    fi;
    r:=QuoInt(dim, 2);
    if dim mod 2 = 1 then
        return 2 * p^(r*r) * Product([1..r], i->p^(2*i)-1);
    fi;
    eps:=Legendre((-1)^r * det, p);
    return 2 * p^(r*(r-1)) * (p^r - eps) * Product([1..r-1], i->p^(2*i)-1);
end;

# The density of the representations of a form S (rank m, det detS) by a
# form T (rank n, det detT) at a prime p not dividing 2 detS detT: by Witt,
# the F_p embeddings are |O(V_T)| / |O(W)|, W the orthogonal complement,
# and Hensel lifts them.
good_prime_density:=function(n, m, detS, detT, p)
    return orthogonal_group_order(n, detT, p) / orthogonal_group_order(n-m, detS*detT, p) / p^(n*m - m*(m+1)/2);
end;

#
# The checks.
#

test_bernoulli:=function(arith)
    local ListK, ListVal, n_error, i;
    ListK:=[0..40];
    ListVal:=get_genus_local_density(arith, "bernoulli", List(ListK, k->[k]));
    if is_error(ListVal) then
        return 1;
    fi;
    n_error:=0;
    for i in [1..Length(ListK)]
    do
        if ListVal[i] <> Bernoulli(ListK[i]) then
            Print("  bernoulli k=", ListK[i], ": got ", ListVal[i], " expected ", Bernoulli(ListK[i]), "\n");
            n_error:=n_error + 1;
        fi;
    od;
    return n_error;
end;

test_generalized_bernoulli:=function(arith)
    local ListRec, ListVal, n_error, D, k, i, expected, ListClass;
    ListRec:=[];
    for D in [-200..200]
    do
        if is_fundamental_discriminant(D) then
            for k in [1..6]
            do
                Add(ListRec, [k, D]);
            od;
        fi;
    od;
    ListVal:=get_genus_local_density(arith, "generalized_bernoulli", ListRec);
    if is_error(ListVal) then
        return 1;
    fi;
    n_error:=0;
    for i in [1..Length(ListRec)]
    do
        expected:=generalized_bernoulli(ListRec[i][1], ListRec[i][2]);
        if ListVal[i] <> expected then
            Print("  generalized_bernoulli ", ListRec[i], ": got ", ListVal[i], " expected ", expected, "\n");
            n_error:=n_error + 1;
        fi;
    od;
    # B_{1, chi_D} = -2 h(D) / w(D) for D < 0: class numbers as a check
    # that shares nothing with the definition above.
    ListClass:=[[-3, 1, 6], [-4, 1, 4], [-7, 1, 2], [-8, 1, 2], [-15, 2, 2],
                [-23, 3, 2], [-47, 5, 2], [-71, 7, 2], [-104, 6, 2], [-163, 1, 2]];
    ListVal:=get_genus_local_density(arith, "generalized_bernoulli", List(ListClass, x->[1, x[1]]));
    if is_error(ListVal) then
        return n_error + 1;
    fi;
    for i in [1..Length(ListClass)]
    do
        expected:=-2 * ListClass[i][2] / ListClass[i][3];
        if ListVal[i] <> expected then
            Print("  class number check D=", ListClass[i][1], ": got ", ListVal[i], " expected ", expected, "\n");
            n_error:=n_error + 1;
        fi;
    od;
    return n_error;
end;

test_euler_value:=function(arith)
    local ListRec, n, m, detS, detT, p, ListVal, n_error, i, expected;
    ListRec:=[];
    for n in [8, 13, 14, 16]
    do
        for m in [1..n-3]
        do
            for detS in [1, 2, 3, 4, 6]
            do
                for detT in [1, 3, 180, 243, 507]
                do
                    for p in [5, 7, 11, 17, 19, 23]
                    do
                        if (2 * detS * detT) mod p <> 0 then
                            Add(ListRec, [n, m, detS, detT, p]);
                        fi;
                    od;
                od;
            od;
        od;
    od;
    ListVal:=get_genus_local_density(arith, "euler_value", ListRec);
    if is_error(ListVal) then
        return 1;
    fi;
    n_error:=0;
    for i in [1..Length(ListRec)]
    do
        expected:=CallFuncList(good_prime_density, ListRec[i]);
        if ListVal[i] <> expected then
            Print("  euler_value ", ListRec[i], ": got ", ListVal[i], " expected ", expected, "\n");
            n_error:=n_error + 1;
        fi;
    od;
    Print("  euler_value: ", Length(ListRec), " cases\n");
    return n_error;
end;

test_generic_part:=function(arith)
    local ListRef, ListVal, n_error, i;
    ListRef:=ReadAsFunction("density_generic_part.g")();
    ListVal:=get_genus_local_density(arith, "generic_part", List(ListRef, x->x{[1..4]}));
    if is_error(ListVal) then
        return 1;
    fi;
    n_error:=0;
    for i in [1..Length(ListRef)]
    do
        if ListVal[i] <> ListRef[i][5] then
            Print("  generic_part ", ListRef[i]{[1..4]}, ": got ", ListVal[i], " expected ", ListRef[i][5], "\n");
            n_error:=n_error + 1;
        fi;
    od;
    Print("  generic_part: ", Length(ListRef), " cases\n");
    return n_error;
end;

FullTest:=function()
    local n_error, arith, eTest, n_err;
    n_error:=0;
    for arith in ListArith
    do
        for eTest in [["bernoulli", test_bernoulli],
                      ["generalized_bernoulli", test_generalized_bernoulli],
                      ["euler_value", test_euler_value],
                      ["generic_part", test_generic_part]]
        do
            n_err:=eTest[2](arith);
            Print("arith=", arith, " test=", eTest[1], " n_error=", n_err, "\n");
            n_error:=n_error + n_err;
        od;
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
