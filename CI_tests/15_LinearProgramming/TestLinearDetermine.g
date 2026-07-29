Read("../common.g");
Read("../access_points.g");
Print("Beginning TestLinearDetermine\n");

# POLY_LinearDetermineByInequalities returns a matrix whose rows span
# the affine span of the polytope defined by the H-rep FAC.  When the
# polytope is full-dimensional this is a square identity; when an
# implicit linearity x.v = 0 is present, the returned matrix has fewer
# rows and every row v satisfies <v, eqv> = 0 for each forced equality.
TestLinearDetermine:=function(eFile, expected_rows, expected_cols,
                              expected_rank)
    local FAC, Mat, ok, r;
    FAC:=ReadAsFunction(eFile)();
    if FAC=fail then
        Print("Could not parse input file ", eFile, "\n");
        return false;
    fi;
    Mat:=get_linear_determine_by_inequalities(FAC);
    if IsString(Mat) then
        Print("eFile=", eFile, ": ", Mat, "\n");
        return false;
    fi;
    ok:=true;
    if Length(Mat) <> expected_rows then
        Print("eFile=", eFile, ": got ", Length(Mat), " rows, expected ",
              expected_rows, "\n");
        ok:=false;
    elif Length(Mat[1]) <> expected_cols then
        Print("eFile=", eFile, ": got ", Length(Mat[1]), " cols, expected ",
              expected_cols, "\n");
        ok:=false;
    else
        r:=RankMat(Mat);
        if r <> expected_rank then
            Print("eFile=", eFile, ": RankMat=", r, ", expected ",
                  expected_rank, "\n");
            ok:=false;
        fi;
    fi;
    if ok then
        Print("eFile=", eFile, ": matrix is ", Length(Mat), "x",
              Length(Mat[1]), " rank=", RankMat(Mat), " OK\n");
    fi;
    return ok;
end;

n_error:=0;

# Full-dimensional 3-cross-polytope, no implicit linearities.  The
# kernel function returns the 4x4 identity (rank 4 in 4-col ambient).
if TestLinearDetermine("Example1_cross3_full.g", 4, 4, 4)=false then
    n_error:=n_error+1;
fi;

# Same polytope intersected with z=0 (two extra inequalities forcing
# the equality).  The implicit linearity drops one dimension: 3
# generators of the {z=0} hyperplane in the 4-col ambient.
if TestLinearDetermine("Example2_slab.g", 3, 4, 3)=false then
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
