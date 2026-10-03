Read("../common.g");
Read("../access_points.g");
Print("Beginning Test of the stabilizers of non-symmetric matrices\n");

# Every case has a known group order, and every generator returned is
# checked here to preserve the matrix, independently of the C++ code.

#
# The programs
#

get_direct_matrix_stabilizer:=function(M)
    local TmpDir, FileI, FileO, FileE, eProg, TheCommand, GRP;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, M);
    eProg:=GetBinaryFilename("GRP_DirectMatrix_Stabilizer");
    TheCommand:=Concatenation(eProg, " rational ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: GRP_DirectMatrix_Stabilizer did not return anything, likely crash";
    fi;
    GRP:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return GRP;
end;

get_grammat_automorphism:=function(EXT, GramMat)
    local TmpDir, FileI, FileG, FileO, FileE, eProg, TheCommand, GRP;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.ext");
    FileG:=Filename(TmpDir, "Test.gram");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, EXT);
    WriteMatrixFile(FileG, GramMat);
    eProg:=GetBinaryFilename("GRP_LinPolytope_Automorphism_GramMat");
    TheCommand:=Concatenation(eProg, " rational ", FileI, " ", FileG, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: GRP_LinPolytope_Automorphism_GramMat did not return anything, likely crash";
    fi;
    GRP:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileG);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return GRP;
end;

#
# The test matrices
#

# The Paley tournament of a prime p = 3 mod 4: M[i][j] is the Legendre symbol
# of j - i. Its automorphism group is {x -> a x + b, a a nonzero square}, of
# order p (p-1) / 2. M is antisymmetric: ignoring the orientation would give
# the complete graph and the symmetric group.
PaleyTournament:=function(p)
    local M, i, j;
    M:=NullMat(p, p);
    for i in [1..p] do
        for j in [1..p] do
            if i<>j then
                M[i][j]:=Legendre(j-i, p);
            fi;
        od;
    od;
    return M;
end;

# M[i][j] = (j - i) mod n, a circulant with n distinct values: only the
# n translations preserve it.
DistinctCirculant:=function(n)
    return List([1..n], i->List([1..n], j->(j-i) mod n));
end;

# All entries distinct: trivial group, and n^2 distinct weights, more than an
# 8-bit weight index can hold when n >= 16.
AllDistinct:=function(n)
    return List([1..n], i->List([1..n], j->n*(i-1) + (j-1)));
end;

# The adjacency matrix of the Petersen graph: symmetric, group of order 120.
PetersenMatrix:=function()
    local LV;
    LV:=Combinations([1..5], 2);
    return List(LV, x->List(LV, function(y)
        if Intersection(x, y)=[] then
            return 1;
        fi;
        return 0;
    end));
end;

# The integer points of [-16,16]^2 but the origin, 1088 vectors.
BoxVectors:=function()
    local EXT, a, b;
    EXT:=[];
    for a in [-16..16] do
        for b in [-16..16] do
            if [a,b]<>[0,0] then
                Add(EXT, [a,b]);
            fi;
        od;
    od;
    return EXT;
end;

#
# The checks
#

IsMatrixAutomorphism:=function(M, g)
    local n, i, j;
    n:=Length(M);
    for i in [1..n] do
        for j in [1..n] do
            if M[i^g][j^g]<>M[i][j] then
                return false;
            fi;
        od;
    od;
    return true;
end;

check_group:=function(name, GRP, ScalMat, expected_order)
    local eGen;
    if is_error(GRP) then
        return false;
    fi;
    for eGen in GeneratorsOfGroup(GRP) do
        if IsMatrixAutomorphism(ScalMat, eGen)=false then
            Print("  ", name, ": a generator does not preserve the matrix\n");
            return false;
        fi;
    od;
    if Order(GRP)<>expected_order then
        Print("  ", name, ": |GRP|=", Order(GRP), " but expected ", expected_order, "\n");
        return false;
    fi;
    Print("  ", name, ": |GRP|=", Order(GRP), " correct\n");
    return true;
end;

DirectMatrixCase:=function(name, M, expected_order)
    return check_group(name, get_direct_matrix_stabilizer(M), M, expected_order);
end;

GramMatCase:=function(name, EXT, GramMat, expected_order)
    local ScalMat;
    # The weight of the pair (i,j) for the automorphisms of the configuration:
    # a permutation preserves EXT GramMat EXT^T exactly when it preserves its
    # transpose, so the convention of the C++ code does not matter.
    ScalMat:=EXT * GramMat * TransposedMat(EXT);
    return check_group(name, get_grammat_automorphism(EXT, GramMat), ScalMat, expected_order);
end;

#
# The tests
#

NonSymmetric_AllTests:=function()
    local n_error, f_case, p, n, Box, Qrot, Qid, Q4;
    n_error:=0;
    # The arguments are evaluated first, so the details of the case are
    # printed before its conclusion.
    f_case:=function(name, test)
        if test=false then
            Print("Case ", name, ": FAILED\n");
            n_error:=n_error+1;
        else
            Print("Case ", name, ": ok\n");
        fi;
    end;
    # The stabilizer of a matrix given directly.
    for p in [7, 11, 19, 23] do
        f_case(Concatenation("Paley(", String(p), ")"),
               DirectMatrixCase("Paley", PaleyTournament(p), p*(p-1)/2));
    od;
    f_case("Paley(7) symmetrized",
           DirectMatrixCase("Paley sym", List(PaleyTournament(7), x->List(x, AbsInt)), Factorial(7)));
    for n in [10, 17] do
        f_case(Concatenation("DistinctCirculant(", String(n), ")"),
               DirectMatrixCase("Circulant", DistinctCirculant(n), n));
    od;
    for n in [16, 20, 22] do
        f_case(Concatenation("AllDistinct(", String(n), ")"),
               DirectMatrixCase("AllDistinct", AllDistinct(n), 1));
    od;
    f_case("Petersen", DirectMatrixCase("Petersen", PetersenMatrix(), 120));
    # The automorphisms of a configuration of vectors for a non-symmetric
    # Gram matrix.
    Box:=BoxVectors();
    # A rotation-invariant non-symmetric form: the rotations by pi/2 preserve
    # it and the reflections exchange it with its transpose.
    Qrot:=[[1,1],[-1,1]];
    Qid:=[[1,0],[0,1]];
    Q4:=[[3,2,-1,0],[-2,5,1,1],[1,-1,4,-2],[0,-1,2,6]];
    # 1088 vectors: above 1000 rows the heuristic scheme is used.
    f_case("Box, rotation-invariant form", GramMatCase("Box rot", Box, Qrot, 4));
    f_case("Box, identity form", GramMatCase("Box id", Box, Qid, 8));
    # 550 random antipodal pairs, only -Id preserves them: 1100 vectors.
    f_case("PlusMinus_1100, rotation-invariant form",
           GramMatCase("PlusMinus_1100", ReadMatrixFile("PlusMinus_1100.ext"), Qrot, 2));
    # 20 generic vectors with 316 distinct weights: more than an 8-bit
    # weight index can hold.
    f_case("Generic_20, non-symmetric form",
           GramMatCase("Generic_20", ReadMatrixFile("Generic_20.ext"), Q4, 1));
    return n_error;
end;

n_error:=NonSymmetric_AllTests();
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
