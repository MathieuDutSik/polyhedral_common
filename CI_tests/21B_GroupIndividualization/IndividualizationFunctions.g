Read("../common.g");

# The automorphism group of large point sets with a small group, computed
# by GRP_LinPolytope_Automorphism. Above THRESHOLD_INDIVIDUALIZATION_STAB
# points with no small spanning subset, it goes through the orbit-stabilizer
# computation by individualization (GetStabilizerWeightMatrix_Individualization)
# and passes the order of the group to its construction. The groups are
# compared with the action of the known group on the points.

# The hyperoctahedral group B_n as signed permutation matrices.
GroupB:=function(n)
    local ListGen, i, D;
    ListGen:=[];
    for i in [1..n-1]
    do
        Add(ListGen, PermutationMat((i, i+1), n));
    od;
    D:=IdentityMat(n);
    D[1][1]:=-1;
    Add(ListGen, D);
    return Group(ListGen);
end;

# The matrix group G acting on R^n, extended by the identity on R^k.
ExtendGroup:=function(G, k)
    local n, ListGen, eGen, M;
    n:=DimensionOfMatrixGroup(G);
    ListGen:=[];
    for eGen in GeneratorsOfGroup(G)
    do
        M:=IdentityMat(n + k);
        M{[1..n]}{[1..n]}:=eGen;
        Add(ListGen, M);
    od;
    return Group(ListGen);
end;

# A test case from points: the polytope of the vectors [1, pt] and the
# action of G on the points as reference.
CaseFromPoints:=function(name, G, ListPts)
    return rec(name:=name, EXT:=List(ListPts, x->Concatenation([1], x)),
               arith:="rational", GRPref:=Action(G, ListPts, OnRight));
end;

# The orbit of a point under B_n.
CaseOrbit:=function(name, n, pt)
    local G;
    G:=GroupB(n);
    return CaseFromPoints(name, G, Set(Orbit(G, pt, OnRight)));
end;

# B_n acting on R^n x R: the 2n points [+-c e_i, 0], whose span does not
# contain the other ones, and the regular orbit of [pt, 1]. The vertices
# outside the large block are few but not spanning, so that the subset
# method needs the large block and the individualization uses the small
# subset made of them and one vertex of the large block.
CaseSmallSubset:=function(n, pt, c)
    local G, GB, ListPts, i, v;
    GB:=GroupB(n);
    G:=ExtendGroup(GB, 1);
    ListPts:=Set(Orbit(G, Concatenation(pt, [1]), OnRight));
    for i in [1..n]
    do
        v:=ListWithIdenticalEntries(n + 1, 0);
        v[i]:=c;
        Add(ListPts, v);
        Add(ListPts, -v);
    od;
    return CaseFromPoints(Concatenation("B", String(n), " small subset"), G,
                          ListPts);
end;

# The WythoffH4 polytope of the test 28A, with its group. Its coordinates
# are in Q(sqrt(5)), given as the real algebraic field of
# 28A_WythoffH4/FileDescSqrt5.
ReadGroupFile:=function(eFile)
    local ListLines, eHead, n, nGen, ListGen, i, eLine;
    ListLines:=Filtered(SplitString(StringFile(eFile), "\n"),
                        x->Length(NormalizedWhitespace(x)) > 0);
    eHead:=List(SplitString(NormalizedWhitespace(ListLines[1]), " "), Int);
    n:=eHead[1];
    nGen:=eHead[2];
    ListGen:=[];
    for i in [1..nGen]
    do
        eLine:=List(SplitString(NormalizedWhitespace(ListLines[i + 1]), " "),
                    Int);
        Add(ListGen, PermList(eLine + 1));
    od;
    return Group(ListGen);
end;

CaseWythoffH4:=function()
    return rec(name:="WythoffH4", FileEXT:="../28A_WythoffH4/WythoffH4.ext",
               arith:="RealAlgebraic=../28A_WythoffH4/FileDescSqrt5",
               GRPref:=ReadGroupFile("../28A_WythoffH4/WythoffH4.grp"));
end;

# Runs GRP_LinPolytope_Automorphism and compares the group with the
# reference one.
TestCase:=function(eRec)
    local TmpDir, FileI, FileO, FileE, eProg, TheCommand, GRP;
    TmpDir:=DirectoryTemporary();
    if IsBound(eRec.FileEXT) then
        FileI:=eRec.FileEXT;
    else
        FileI:=Filename(TmpDir, "Test.ext");
        RemoveFileIfExist(FileI);
        WriteMatrixFile(FileI, eRec.EXT);
    fi;
    FileO:=Filename(TmpDir, "Test.grp");
    FileE:=Filename(TmpDir, "Test.err");
    RemoveFileIfExist(FileO);
    eProg:=GetBinaryFilename("GRP_LinPolytope_Automorphism");
    TheCommand:=Concatenation(eProg, " ", eRec.arith, " ", FileI, " GAP ",
                              FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO) = false then
        Print(eRec.name, ": no output file\n");
        return false;
    fi;
    GRP:=ReadAsFunction(FileO)();
    Print(eRec.name, ": |GRP|=", Order(GRP), " |GRPref|=", Order(eRec.GRPref),
          "\n");
    if GRP <> eRec.GRPref then
        Print(eRec.name, ": the groups are different\n");
        return false;
    fi;
    return true;
end;

RunIndividualizationTests:=function(ListRec)
    local n_error, eRec;
    n_error:=0;
    for eRec in ListRec
    do
        if TestCase(eRec) = false then
            n_error:=n_error + 1;
        fi;
    od;
    Print("n_error=", n_error, "\n");
    CI_Decision_Reset();
    if n_error > 0 then
        Print("Error case\n");
    else
        Print("Normal case\n");
        CI_Write_Ok();
    fi;
end;
