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
