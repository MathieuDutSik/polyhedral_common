Read("../common.g");
Read("../access_points.g");
Print("Beginning Hecke operator tests\n");

# The Hecke matrices are compared with the known eigenvalues of the weight 2
# modular forms: for Gamma_0(11) the cusp form 11a has a_2=-2, a_3=-1,
# a_5=1, a_7=-2 and the Eisenstein class has eigenvalue 1+p. For Gamma_1(13)
# the two cusp forms have a_2 root of x^2+3x+3.

# The level of the complex whose homology carries the weight 2 forms is
# index 1 (the edges of the Voronoi tessellation).
get_hecke_level:=function(TheResult, index)
    local eLevel;
    for eLevel in TheResult.ListLevel
    do
        if eLevel.index=index then
            return eLevel;
        fi;
    od;
    Error("Failed to find the level");
end;

test_gamma0_11:=function()
    local desc, p, a_p, ListPrime, ListAp, i, x, TheResult, eLevel, M, dim, P1, P2;
    desc:=GenerateTspaceDescription_classic(2, false);
    ListPrime:=[2, 3, 5, 7];
    ListAp:=[-2, -1, 1, -2];
    for i in [1..Length(ListPrime)]
    do
        p:=ListPrime[i];
        a_p:=ListAp[i];
        x:=[[1,0],[0,p]];
        TheResult:=PERFCOMP_hecke_operators(desc, x, "Gamma0", 11, true);
        Print("Gamma0(11) p=", p, " n_coset_gamma=", TheResult.n_coset_gamma, " n_hecke_coset=", TheResult.n_hecke_coset, "\n");
        if TheResult.n_coset_gamma<>12 then
            Print("The index of Gamma0(11) should be 12\n");
            return false;
        fi;
        if TheResult.n_hecke_coset<>p+1 then
            Print("The number of Hecke cosets should be p+1\n");
            return false;
        fi;
        eLevel:=get_hecke_level(TheResult, 1);
        M:=eLevel.HeckeMatrix;
        dim:=eLevel.dim_homology;
        Print("  dim_homology=", dim, " HeckeMatrix=", M, "\n");
        if dim<>2 then
            Print("The homology should be of dimension 2\n");
            return false;
        fi;
        P1:=M - (1+p)*IdentityMat(dim);
        P2:=M - a_p*IdentityMat(dim);
        if P1*P2<>NullMat(dim, dim) then
            Print("The eigenvalues should be 1+p and a_p\n");
            return false;
        fi;
    od;
    return true;
end;

test_gamma1_13:=function()
    local desc, x, TheResult, eLevel, M, dim, P, nullity;
    desc:=GenerateTspaceDescription_classic(2, false);
    x:=[[1,0],[0,2]];
    TheResult:=PERFCOMP_hecke_operators(desc, x, "Gamma1", 13, true);
    Print("Gamma1(13) p=2 n_coset_gamma=", TheResult.n_coset_gamma, " n_hecke_coset=", TheResult.n_hecke_coset, "\n");
    eLevel:=get_hecke_level(TheResult, 1);
    M:=eLevel.HeckeMatrix;
    dim:=eLevel.dim_homology;
    Print("  dim_homology=", dim, "\n");
    if dim<>13 then
        Print("The homology should be of dimension 13\n");
        return false;
    fi;
    # The two cusp forms give the factor x^2+3x+3 of the characteristic
    # polynomial, that is a kernel of dimension 2 for M^2+3M+3.
    P:=M*M + 3*M + 3*IdentityMat(dim);
    nullity:=dim - RankMat(P);
    Print("  nullity of M^2+3M+3 = ", nullity, "\n");
    if nullity<>2 then
        Print("The cusp forms of Gamma1(13) should give the factor x^2+3x+3\n");
        return false;
    fi;
    return true;
end;

test_gl3_full:=function()
    local desc, x, TheResult, eLevel, M, dim;
    desc:=GenerateTspaceDescription_classic(3, false);
    x:=[[1,0,0],[0,1,0],[0,0,2]];
    TheResult:=PERFCOMP_hecke_operators(desc, x, "Full", 1, true);
    Print("GL3 full p=2 n_hecke_coset=", TheResult.n_hecke_coset, "\n");
    if TheResult.n_hecke_coset<>7 then
        Print("The number of Hecke cosets should be 7\n");
        return false;
    fi;
    eLevel:=get_hecke_level(TheResult, 0);
    M:=eLevel.HeckeMatrix;
    dim:=eLevel.dim_homology;
    Print("  dim_homology=", dim, " HeckeMatrix=", M, "\n");
    if dim<>1 or M<>[[7]] then
        Print("The trivial class should have eigenvalue 7\n");
        return false;
    fi;
    return true;
end;

# The sharbly method computes the Hecke operator at the symbol level (the
# top degree), by the Ash-Rudolph reduction to unimodular symbols, and at
# the level below it by filling. It computes the transposed side operator
# relative to the Dual chain map, so the matrices differ but the
# characteristic polynomials coincide.
test_sharbly:=function()
    local desc, p, a_p, ListPrime, ListAp, i, x, ResS, ResD, eLevS, eLevD, M, dim, P1, P2, N, ListN, ListChar, index, FileCache;
    desc:=GenerateTspaceDescription_classic(2, false);
    ListPrime:=[2, 3];
    ListAp:=[-2, -1];
    for i in [1..Length(ListPrime)]
    do
        p:=ListPrime[i];
        a_p:=ListAp[i];
        x:=[[1,0],[0,p]];
        ResS:=PERFCOMP_hecke_operators(desc, x, "Gamma0", 11, true, "Sharbly");
        ResD:=PERFCOMP_hecke_operators(desc, x, "Gamma0", 11, true, "Dual");
        eLevS:=get_hecke_level(ResS, 1);
        eLevD:=get_hecke_level(ResD, 1);
        M:=eLevS.HeckeMatrix;
        dim:=eLevS.dim_homology;
        Print("Sharbly Gamma0(11) p=", p, " dim=", dim, " M=", M, "\n");
        if dim<>2 then
            Print("The homology should be of dimension 2\n");
            return false;
        fi;
        P1:=M - (1+p)*IdentityMat(dim);
        P2:=M - a_p*IdentityMat(dim);
        if P1*P2<>NullMat(dim, dim) then
            Print("The eigenvalues should be 1+p and a_p\n");
            return false;
        fi;
        if CharacteristicPolynomial(M)<>CharacteristicPolynomial(eLevD.HeckeMatrix) then
            Print("The Sharbly and Dual characteristic polynomials should coincide\n");
            return false;
        fi;
    od;
    # GL2, Principal(3) and Principal(4): the level below the symbol level
    # (index 0) carries the trivial class, with eigenvalue 1+p, and both
    # levels are compared with the Dual method.
    for N in [3, 4]
    do
        p:=N-1;
        x:=[[1,0],[0,p]];
        ResS:=PERFCOMP_hecke_operators(desc, x, "Principal", N, true, "Sharbly");
        ResD:=PERFCOMP_hecke_operators(desc, x, "Principal", N, true, "Dual");
        eLevS:=get_hecke_level(ResS, 0);
        Print("Sharbly Principal(", N, ") p=", p, " index=0 M=", eLevS.HeckeMatrix, "\n");
        if eLevS.HeckeMatrix<>[[1+p]] then
            Print("The trivial class should have eigenvalue 1+p\n");
            return false;
        fi;
        for index in [0, 1]
        do
            eLevS:=get_hecke_level(ResS, index);
            eLevD:=get_hecke_level(ResD, index);
            if CharacteristicPolynomial(eLevS.HeckeMatrix)<>CharacteristicPolynomial(eLevD.HeckeMatrix) then
                Print("The Sharbly and Dual characteristic polynomials should coincide\n");
                return false;
            fi;
        od;
    od;
    # GL2, Gamma0(11) intersected with SL2: the level below (index 0) has
    # the trivial class and the symbol level (index 1) has dimension 3 with
    # T2 acting by (x-3)(x+2)^2 (Eisenstein class and the form 11a).
    x:=[[1,0],[0,2]];
    ResS:=PERFCOMP_hecke_operators(desc, x, "Gamma0SL", 11, true, "Sharbly");
    ResD:=PERFCOMP_hecke_operators(desc, x, "Gamma0SL", 11, true, "Dual");
    ListChar:=[[-3, 1], [-12, -8, 1, 1]];
    for index in [0, 1]
    do
        eLevS:=get_hecke_level(ResS, index);
        eLevD:=get_hecke_level(ResD, index);
        Print("Sharbly Gamma0SL(11) index=", index, " M=", eLevS.HeckeMatrix, "\n");
        if CoefficientsOfUnivariatePolynomial(CharacteristicPolynomial(eLevS.HeckeMatrix))<>ListChar[index+1] then
            Print("The characteristic polynomial does not match the expected one\n");
            return false;
        fi;
        if CharacteristicPolynomial(eLevS.HeckeMatrix)<>CharacteristicPolynomial(eLevD.HeckeMatrix) then
            Print("The Sharbly and Dual characteristic polynomials should coincide\n");
            return false;
        fi;
    od;
    # The chain map does not depend on the subgroup: computed once by the
    # Sharbly method for Gamma0(11), it is read from the cache for
    # Gamma0SL(11). The second run asks for the Dual method, so getting the
    # Sharbly method and matrices back shows that the cache was used.
    FileCache:=Filename(DirectoryTemporary(), "hecke_chain_map_cache");
    ResS:=PERFCOMP_hecke_operators(desc, x, "Gamma0", 11, true, "Sharbly", FileCache);
    ResD:=PERFCOMP_hecke_operators(desc, x, "Gamma0SL", 11, true, "Dual", FileCache);
    Print("Cache: method of the second run=", ResD.method, "\n");
    if ResS.method<>"Sharbly" or ResD.method<>"Sharbly" then
        Print("The chain map should have been read from the cache\n");
        return false;
    fi;
    ResS:=PERFCOMP_hecke_operators(desc, x, "Gamma0SL", 11, true, "Sharbly");
    for index in [0, 1]
    do
        if get_hecke_level(ResS, index).HeckeMatrix<>get_hecke_level(ResD, index).HeckeMatrix then
            Print("The Hecke matrices with and without the cache should be equal\n");
            return false;
        fi;
    od;
    RemoveFile(FileCache);
    # GL3, Gamma0(9) and Gamma0(11): nonzero homology at the symbol level
    # (index 3), compared between the Sharbly and the Dual methods.
    desc:=GenerateTspaceDescription_classic(3, false);
    x:=[[1,0,0],[0,1,0],[0,0,2]];
    ListN:=[9, 11];
    ListChar:=[[-7, 1], [-6, 1, 1]];
    for i in [1..Length(ListN)]
    do
        N:=ListN[i];
        ResS:=PERFCOMP_hecke_operators(desc, x, "Gamma0", N, true, "Sharbly");
        ResD:=PERFCOMP_hecke_operators(desc, x, "Gamma0", N, true, "Dual");
        eLevS:=get_hecke_level(ResS, 3);
        eLevD:=get_hecke_level(ResD, 3);
        Print("Sharbly GL3 Gamma0(", N, ") dim=", eLevS.dim_homology, " M=", eLevS.HeckeMatrix, "\n");
        if CharacteristicPolynomial(eLevS.HeckeMatrix)<>CharacteristicPolynomial(eLevD.HeckeMatrix) then
            Print("The Sharbly and Dual characteristic polynomials should coincide\n");
            return false;
        fi;
        if CoefficientsOfUnivariatePolynomial(CharacteristicPolynomial(eLevS.HeckeMatrix))<>ListChar[i] then
            Print("The characteristic polynomial does not match the stored one\n");
            return false;
        fi;
    od;
    return true;
end;

# The homology modulo a prime p by sparse linear algebra: the Hecke
# matrices are given by their residues, and their characteristic
# polynomials must be those of the rational Hecke matrices reduced mod p,
# for the Sharbly and the Dual (transposed) chain maps.
# The residues modulo p of a rational matrix, as elements of GF(p).
matrix_mod_p:=function(M, p)
    return List(M, row -> List(row, x -> (NumeratorRat(x) * PowerMod(DenominatorRat(x), -1, p)) mod p) * One(GF(p)));
end;

test_homology_mod:=function()
    local p, desc, x, ListCase, eCase, ResQ, ResP, index, eLevQ, eLevP, CharQ, CharP;
    p:=2147483647;
    ListCase:=[[2, [[1,0],[0,2]], "Gamma0SL", 11, "Sharbly", [0, 1]],
               [2, [[1,0],[0,2]], "Gamma0SL", 11, "Dual", [0, 1]],
               [3, [[1,0,0],[0,1,0],[0,0,2]], "Gamma0", 11, "Sharbly", [2, 3]]];
    for eCase in ListCase
    do
        desc:=GenerateTspaceDescription_classic(eCase[1], false);
        ResQ:=PERFCOMP_hecke_operators(desc, eCase[2], eCase[3], eCase[4], true, eCase[5], "none", 0);
        ResP:=PERFCOMP_hecke_operators(desc, eCase[2], eCase[3], eCase[4], true, eCase[5], "none", p);
        if ResP.prime<>p then
            Print("The modular result should record the prime\n");
            return false;
        fi;
        for index in eCase[6]
        do
            eLevQ:=get_hecke_level(ResQ, index);
            eLevP:=get_hecke_level(ResP, index);
            if eLevQ.dim_homology<>eLevP.dim_homology then
                Print("The dimensions of the homology differ, ", eCase[5], " index=", index, "\n");
                return false;
            fi;
            if eLevQ.dim_homology > 0 then
                CharQ:=CharacteristicPolynomial(matrix_mod_p(eLevQ.HeckeMatrix, p));
                CharP:=CharacteristicPolynomial(matrix_mod_p(eLevP.HeckeMatrix, p));
                Print("Mod p ", eCase[3], "(", eCase[4], ") ", eCase[5], " index=", index,
                      " dim=", eLevP.dim_homology, " same=", CharQ=CharP, "\n");
                if CharQ<>CharP then
                    Print("The characteristic polynomials mod p differ\n");
                    return false;
                fi;
            fi;
        od;
    od;
    return true;
end;

# GL3 over the Gaussian integers, Gamma0(11 O_K) in SL3(Z[i]), the Sharbly
# method: the homology at index 5 (the degree 5 of the cohomology) is of
# dimension 2, carrying the Eisenstein classes of the Bianchi form 121.1-a
# (the base change of 11a): T_pi has the eigenvalues a + q^2 and 1 + q a,
# with q = N(pi) and a = a_pi, which is a_2 = -2 for pi = 1+i and a_5 = 1
# for pi = 2+i. The index of Gamma is 59052: the homology is computed
# modulo p (the rational one uses dense matrices).
test_sharbly_gaussian:=function()
    local p, desc, ListX, ListChar, i, Res, eLev, CharP;
    p:=2147483647;
    desc:=GenerateTspaceDescription_imag_quad(3, -4, true);
    # The realified matrices diag(1, 1, pi), in the layout
    # (x_1, x_2, x_3, y_1, y_2, y_3) of z_j = x_j + i y_j.
    ListX:=[[[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,1],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,-1,0,0,1]],
            [[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,2,0,0,1],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,-1,0,0,2]]];
    ListChar:=[[-6, 1, 1], [156, -32, 1]];
    for i in [1..Length(ListX)]
    do
        Res:=PERFCOMP_hecke_operators(desc, ListX[i], "Gamma0SL", 11, true, "Sharbly", "none", p);
        eLev:=get_hecke_level(Res, 5);
        Print("Sharbly SL3(Z[i]) Gamma0(11) i=", i, " dim=", eLev.dim_homology, " M=", eLev.HeckeMatrix, "\n");
        if eLev.dim_homology<>2 then
            Print("The homology should be of dimension 2\n");
            return false;
        fi;
        CharP:=CharacteristicPolynomial(matrix_mod_p(eLev.HeckeMatrix, p));
        if CoefficientsOfUnivariatePolynomial(CharP)<>List(ListChar[i], c->c*One(GF(p))) then
            Print("The characteristic polynomial does not match the expected one\n");
            return false;
        fi;
    od;
    return true;
end;

f_compute:=function()
    local n_error;
    n_error:=0;
    if not test_gamma0_11() then
        n_error:=n_error+1;
    fi;
    if not test_gamma1_13() then
        n_error:=n_error+1;
    fi;
    if not test_gl3_full() then
        n_error:=n_error+1;
    fi;
    if not test_sharbly() then
        n_error:=n_error+1;
    fi;
    if not test_homology_mod() then
        n_error:=n_error+1;
    fi;
    if not test_sharbly_gaussian() then
        n_error:=n_error+1;
    fi;
    return n_error=0;
end;

test:=f_compute();
CI_Decision_Reset();
if test=false then
    Print("Error case\n");
else
    Print("Normal case\n");
    CI_Write_Ok();
fi;
