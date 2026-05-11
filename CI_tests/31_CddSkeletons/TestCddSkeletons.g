Read("../common.g");
Read("../access_points.g");
Print("Beginning TestCddSkeletons\n");

# Given a vertex count and an edge list (1-based pairs), return the
# multiset of vertex degrees as Collected(degs) = [[deg, count], ...]
# sorted by deg.  This is a much stronger combinatorial fingerprint of
# the graph than just nbVert.
DegSequence:=function(nbVert, ListEdges)
    local degs, eEdge;
    degs:=ListWithIdenticalEntries(nbVert, 0);
    for eEdge in ListEdges do
        degs[eEdge[1]]:=degs[eEdge[1]] + 1;
        degs[eEdge[2]]:=degs[eEdge[2]] + 1;
    od;
    return Collected(degs);
end;

TestSkeletons:=function(eFile, expected_FAC_dims, expected_SkelDegSeq,
                        expected_RidgeDegSeq)
    local EXT, eRec, ok, SkelDeg, RidgeDeg;
    EXT:=ReadAsFunction(eFile)();
    if EXT=fail then
        Print("Could not parse input file ", eFile, "\n");
        return false;
    fi;
    eRec:=get_cdd_skeletons(EXT);
    if IsString(eRec) then
        Print("eFile=", eFile, ": ", eRec, "\n");
        return false;
    fi;
    ok:=true;
    if Length(eRec.FAC) <> expected_FAC_dims[1]
       or Length(eRec.FAC[1]) <> expected_FAC_dims[2] then
        Print("eFile=", eFile, ": FAC is ", Length(eRec.FAC), "x",
              Length(eRec.FAC[1]), ", expected ", expected_FAC_dims, "\n");
        ok:=false;
    fi;
    SkelDeg:=DegSequence(eRec.nbVertSkel, eRec.SkelEdges);
    if SkelDeg <> expected_SkelDegSeq then
        Print("eFile=", eFile, ": SkelDegSeq=", SkelDeg, ", expected ",
              expected_SkelDegSeq, "\n");
        ok:=false;
    fi;
    RidgeDeg:=DegSequence(eRec.nbVertRidge, eRec.RidgeEdges);
    if RidgeDeg <> expected_RidgeDegSeq then
        Print("eFile=", eFile, ": RidgeDegSeq=", RidgeDeg, ", expected ",
              expected_RidgeDegSeq, "\n");
        ok:=false;
    fi;
    if ok then
        Print("eFile=", eFile, ": FAC=", Length(eRec.FAC), "x",
              Length(eRec.FAC[1]), " SkelDegSeq=", SkelDeg,
              " RidgeDegSeq=", RidgeDeg, " OK\n");
    fi;
    return ok;
end;

n_error:=0;

# 3-cube:
#  - FAC: 6 facets in dim 4 (homogeneous);
#  - skeleton: 8 vertices, each of degree 3 (3-regular);
#  - ridge: 6 facet-nodes, each of degree 4 (each facet shares a ridge
#           with the 4 facets that are not itself or its antipode).
if TestSkeletons("Example1_cube3.g",
                 [6, 4],
                 [[3, 8]],
                 [[4, 6]])=false then
    n_error:=n_error+1;
fi;

# 3-cross-polytope:
#  - FAC: 8 facets in dim 4;
#  - skeleton: 6 vertices (the +/-e_i), every pair adjacent except
#              antipodal pairs, so degree 4 (4-regular);
#  - ridge: 8 facet-nodes (octants), each octant shares a ridge with
#           the 3 octants differing in exactly one sign, so 3-regular.
if TestSkeletons("Example2_cross3.g",
                 [8, 4],
                 [[4, 6]],
                 [[3, 8]])=false then
    n_error:=n_error+1;
fi;

Print("n_error=", n_error, "\n");

CI_Decision_Reset();
if n_error > 0 then
    Print("Error case\n");
else
    Print("Normal case\n");
    CI_Write_Ok();
fi;
