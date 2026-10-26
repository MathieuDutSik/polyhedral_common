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


#
# Katsurada's Siegel series and the densities at a unimodular prime.
#

test_siegel_series:=function(arith)
    local ListRef, ListRec, ListVal, n_error, i;
    ListRef:=ReadAsFunction("density_siegel_series.g")();
    ListRec:=List(ListRef, x->Concatenation([x[1], x[2]], x[3]));
    ListVal:=get_genus_local_density(arith, "siegel_series", ListRec);
    if is_error(ListVal) then
        return 1;
    fi;
    n_error:=0;
    for i in [1..Length(ListRef)]
    do
        if ListVal[i] <> ListRef[i][4] then
            Print("  siegel_series p=", ListRef[i][1], " S=", ListRef[i][3], ": got ", ListVal[i], " expected ", ListRef[i][4], "\n");
            n_error:=n_error + 1;
        fi;
    od;
    Print("  siegel_series: ", Length(ListRef), " cases\n");
    return n_error;
end;

test_alpha_unimodular:=function(arith)
    local ListRef, ListRec, ListVal, n_error, i;
    ListRef:=ReadAsFunction("density_alpha_unimodular.g")();
    ListRec:=List(ListRef, x->Concatenation([x[1], x[2], x[3], x[4]], x[5]));
    ListVal:=get_genus_local_density(arith, "alpha_unimodular", ListRec);
    if is_error(ListVal) then
        return 1;
    fi;
    n_error:=0;
    for i in [1..Length(ListRef)]
    do
        if ListVal[i] <> ListRef[i][6] then
            Print("  alpha_unimodular ", ListRef[i]{[1..4]}, " S=", ListRef[i][5], ": got ", ListVal[i], " expected ", ListRef[i][6], "\n");
            n_error:=n_error + 1;
        fi;
    od;
    Print("  alpha_unimodular: ", Length(ListRef), " cases\n");
    return n_error;
end;

#
# Representation numbers, counted in GAP, against the whole chain
#     M_S / mass = 2 * (generic part) * prod_{p | 2 det S} alpha_p(S, T)
# for the genus of E8 (one class) and of the even unimodular lattices of
# rank 16 (two classes, E8 + E8 and D16+).
#

cartan_A:=function(n)
    return List([1..n], i->List([1..n], function(j)
        if i=j then return 2; fi;
        if AbsInt(i-j)=1 then return -1; fi;
        return 0; end));
end;

cartan_D:=function(n)
    local M;
    M:=cartan_A(n);
    M[n-1][n]:=0; M[n][n-1]:=0;
    M[n-2][n]:=-1; M[n][n-2]:=-1;
    return M;
end;

# E8 with the Bourbaki labelling: the chain 1-3-4-5-6-7-8 and 2 on 4.
cartan_E8:=function()
    local M, edges, e;
    M:=2 * IdentityMat(8);
    edges:=[[1,3],[3,4],[4,5],[5,6],[6,7],[7,8],[2,4]];
    for e in edges
    do
        M[e[1]][e[2]]:=-1;
        M[e[2]][e[1]]:=-1;
    od;
    return M;
end;

block_diagonal:=function(ListM)
    local n, M, pos, A, i, j;
    n:=Sum(List(ListM, Length));
    M:=NullMat(n, n);
    pos:=0;
    for A in ListM
    do
        for i in [1..Length(A)]
        do
            for j in [1..Length(A)]
            do
                M[pos+i][pos+j]:=A[i][j];
            od;
        od;
        pos:=pos + Length(A);
    od;
    return M;
end;

# The Gram matrix of D16+, from a Z-basis of D16 + Z (1/2, ..., 1/2).
gram_D16plus:=function()
    local gens, i, v, B;
    gens:=[];
    for i in [1..15]
    do
        v:=ListWithIdenticalEntries(16, 0);
        v[i]:=2; v[i+1]:=-2;
        Add(gens, v);
    od;
    v:=ListWithIdenticalEntries(16, 0);
    v[15]:=2; v[16]:=2;
    Add(gens, v);
    Add(gens, ListWithIdenticalEntries(16, 1));
    B:=BaseIntMat(gens) / 2;
    return B * TransposedMat(B) ;
end;

# The roots of G (both signs) and their inner products.
roots_and_products:=function(G)
    local sv, roots;
    sv:=ShortestVectors(G, 2).vectors;
    roots:=Concatenation(sv, -sv);
    return rec(roots:=roots, IP:=List(roots, x->List(roots, y->x * G * y)));
end;

# The number of ordered tuples of roots with Gram matrix S (diagonal 2).
count_representations:=function(RP, S)
    local m, n_root, rec_count;
    m:=Length(S);
    n_root:=Length(RP.roots);
    rec_count:=function(chosen)
        local k, total, j;
        k:=Length(chosen) + 1;
        total:=0;
        for j in [1..n_root]
        do
            if ForAll([1..k-1], i->RP.IP[chosen[i]][j] = S[i][k]) then
                if k = m then
                    total:=total + 1;
                else
                    total:=total + rec_count(Concatenation(chosen, [j]));
                fi;
            fi;
        od;
        return total;
    end;
    return rec_count([]);
end;

# 2 * (generic part) * prod_{p | 2 det S} alpha_p(S, T), T of rank n = 2k
# and determinant 1, from GENUS_LocalDensity.
predicted_mass_ratio:=function(arith, S, k)
    local m, detS, gp, ListP, ListRec, ListAlpha;
    m:=Length(S);
    detS:=DeterminantMat(S);
    gp:=get_genus_local_density(arith, "generic_part", [[2*k, m, detS, 1]]);
    if is_error(gp) then
        return fail;
    fi;
    ListP:=Set(FactorsInt(2 * detS));
    ListRec:=List(ListP, p->Concatenation([p, k, kronecker_symbol((-1)^k, p), m], Concatenation(S)));
    ListAlpha:=get_genus_local_density(arith, "alpha_unimodular", ListRec);
    if is_error(ListAlpha) then
        return fail;
    fi;
    return 2 * gp[1] * Product(ListAlpha);
end;

test_representation_numbers:=function(arith)
    local n_error, E8, RP8, G16a, G16b, RP16a, RP16b, aut_a, aut_b, mass16, ListS8, ListS16, S, truth, pred;
    n_error:=0;
    E8:=cartan_E8();
    RP8:=roots_and_products(E8);
    ListS8:=[cartan_A(1), cartan_A(2), block_diagonal([cartan_A(1), cartan_A(1)]),
             cartan_A(3), block_diagonal([cartan_A(1), cartan_A(2)]),
             block_diagonal([cartan_A(1), cartan_A(1), cartan_A(1)]),
             cartan_A(4), cartan_D(4), block_diagonal([cartan_A(2), cartan_A(2)])];
    for S in ListS8
    do
        truth:=count_representations(RP8, S);
        pred:=predicted_mass_ratio(arith, S, 4);
        if pred <> truth then
            Print("  E8, S=", S, ": counted ", truth, " predicted ", pred, "\n");
            n_error:=n_error + 1;
        fi;
    od;
    G16a:=block_diagonal([E8, E8]);
    G16b:=gram_D16plus();
    RP16a:=roots_and_products(G16a);
    RP16b:=roots_and_products(G16b);
    aut_a:=2 * 696729600^2;
    aut_b:=2^15 * Factorial(16);
    mass16:=1/aut_a + 1/aut_b;
    ListS16:=[cartan_A(1), cartan_A(2), block_diagonal([cartan_A(1), cartan_A(1)]), cartan_A(3)];
    for S in ListS16
    do
        truth:=(count_representations(RP16a, S) / aut_a + count_representations(RP16b, S) / aut_b) / mass16;
        pred:=predicted_mass_ratio(arith, S, 8);
        if pred <> truth then
            Print("  rank 16, S=", S, ": counted ", truth, " predicted ", pred, "\n");
            n_error:=n_error + 1;
        fi;
    od;
    Print("  representation numbers: ", Length(ListS8), " forms in E8, ", Length(ListS16), " in rank 16\n");
    return n_error;
end;


#
# Densities at the primes where T is not unimodular.
#

# A reference table whose entries are [inputs..., value]; value may be fail,
# and the C++ program must then answer fail as well.
compare_table:=function(arith, task, FileRef, f_record, n_input)
    local ListRef, ListVal, n_error, i, expected;
    ListRef:=ReadAsFunction(FileRef)();
    ListVal:=get_genus_local_density(arith, task, List(ListRef, f_record));
    if is_error(ListVal) then
        return 1;
    fi;
    n_error:=0;
    for i in [1..Length(ListRef)]
    do
        expected:=ListRef[i][n_input + 1];
        if ListVal[i] <> expected then
            Print("  ", task, " record ", i, ": got ", ListVal[i], " expected ", expected, "\n");
            n_error:=n_error + 1;
        fi;
    od;
    Print("  ", task, ": ", Length(ListRef), " cases\n");
    return n_error;
end;

test_rank_one_density:=function(arith)
    return compare_table(arith, "rank_one_density", "density_rank_one_density.g",
                         x->Concatenation([x[1], x[2]], x[3], [x[4]]), 4);
end;

test_unit_chain:=function(arith)
    return compare_table(arith, "unit_chain", "density_unit_chain.g",
                         x->Concatenation([x[1], x[2]], x[3], [x[4]], x[5]), 5);
end;

# Certified values M_S / mass of completely enumerated genera: a value the
# program computes must be exact; fail (outside the unit chain) is allowed.
test_certified_mass_ratio:=function(arith)
    local ListRef, ListVal, n_error, n_covered, i;
    ListRef:=ReadAsFunction("density_certified_mass_ratio.g")();
    ListVal:=get_genus_local_density(arith, "mass_ratio",
        List(ListRef, x->Concatenation([x[1]], x[2], [x[3]], x[4])));
    if is_error(ListVal) then
        return 1;
    fi;
    n_error:=0;
    n_covered:=0;
    for i in [1..Length(ListRef)]
    do
        if ListVal[i] <> fail then
            n_covered:=n_covered + 1;
            if ListVal[i] <> ListRef[i][5] then
                Print("  certified mass ratio ", i, ": got ", ListVal[i], " expected ", ListRef[i][5], "\n");
                n_error:=n_error + 1;
            fi;
        fi;
    od;
    Print("  certified mass ratio: ", Length(ListRef), " cases, ", n_covered, " covered\n");
    if n_covered < 80 then
        Print("  certified mass ratio: coverage dropped to ", n_covered, "\n");
        n_error:=n_error + 1;
    fi;
    return n_error;
end;

# The classes of a genus from GENUS_Enumerate, the representations counted in
# GAP class by class: M_S / mass against the densities, for genera that are
# not unimodular at the primes dividing their determinant.
enumerate_genus:=function(name)
    local TmpDir, FileO, FileE, eProg, list_arg;
    TmpDir:=DirectoryTemporary();
    FileO:=Filename(TmpDir, "genus.out");
    FileE:=Filename(TmpDir, "genus.err");
    eProg:=GetBinaryFilename("GENUS_Enumerate");
    list_arg:=["gmp", Concatenation(name, "_genus.txt"), Concatenation(name, "_lattice.txt"),
               Concatenation(name, "_mass.txt"), "GAP", FileO];
    run_binary_checked(eProg, list_arg, FileE);
    return ReadAsFunction(FileO)();
end;

test_enumerated_genera:=function(arith)
    local n_error, name, R, ListRP, ListS, S, truth, ListVal, i, n_checked;
    n_error:=0;
    n_checked:=0;
    ListS:=[cartan_A(1), cartan_A(2), block_diagonal([cartan_A(1), cartan_A(1)]),
            cartan_A(3), block_diagonal([cartan_A(1), cartan_A(2)]),
            block_diagonal([cartan_A(1), cartan_A(1), cartan_A(1)]), cartan_D(4),
            cartan_A(4), block_diagonal([cartan_A(2), cartan_A(2)])];
    for name in ["rank8_det9", "rank8_det25"]
    do
        R:=enumerate_genus(name);
        if R.complete <> true then
            Print("  ", name, ": the enumeration is not complete\n");
            return 1;
        fi;
        ListRP:=List(R.ListGram, roots_and_products);
        ListVal:=get_genus_local_density(arith, "mass_ratio",
            List(ListS, S->Concatenation([8], Concatenation(R.ListGram[1]), [Length(S)], Concatenation(S))));
        if is_error(ListVal) then
            return n_error + 1;
        fi;
        for i in [1..Length(ListS)]
        do
            if ListVal[i] = fail then
                continue;
            fi;
            S:=ListS[i];
            truth:=Sum([1..Length(R.ListGram)], j->count_representations(ListRP[j], S) / R.ListAutOrder[j]) / R.accumulated_mass;
            n_checked:=n_checked + 1;
            if ListVal[i] <> truth then
                Print("  ", name, ", S=", S, ": counted ", truth, " predicted ", ListVal[i], "\n");
                n_error:=n_error + 1;
            fi;
        od;
    od;
    Print("  enumerated genera: ", n_checked, " values checked against counted representations\n");
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
                      ["generic_part", test_generic_part],
                      ["siegel_series", test_siegel_series],
                      ["alpha_unimodular", test_alpha_unimodular],
                      ["representation_numbers", test_representation_numbers],
                      ["rank_one_density", test_rank_one_density],
                      ["unit_chain", test_unit_chain],
                      ["certified_mass_ratio", test_certified_mass_ratio],
                      ["enumerated_genera", test_enumerated_genera]]
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
