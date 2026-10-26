Read("ErdahlTestFunctions.g");
Print("Beginning TestErdahl_zero_set\n");

# Runs ERDAHL_zero_set on the function f(x) = F[(1, x)] and compares the
# result with the expected one. eRec has the fields name, F and
# nonnegative, and when nonnegative is true min, n_ext (the number of
# representatives of the zero set) and rank_L (the rank of its lattice).
# The points of the zero set have to be zeros of f, and the points given
# when f takes negative values have to be negative values.
TestZeroSet:=function(eRec)
    local FileF, FileOut, eProg, TheCommand, answer, eEXT;
    Print("Case ", eRec.name, "\n");
    FileF:=Filename(DirectoryTemporary(), "Test.F");
    FileOut:=Filename(DirectoryTemporary(), "Test.out");
    WriteMatrixFile(FileF, eRec.F);
    eProg:=GetBinaryFilename("ERDAHL_zero_set");
    TheCommand:=Concatenation(eProg, " ", FileF, " ", FileOut);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    answer:=ReadAsFunction(FileOut)();
    RemoveFile(FileOut);
    RemoveFile(FileF);
    if answer.nonnegative<>eRec.nonnegative then
        Print("nonnegative=", answer.nonnegative, " expected=", eRec.nonnegative, "\n");
        return false;
    fi;
    if answer.nonnegative=false then
        if Length(answer.points)=0 then
            Print("No point of negative value is given\n");
            return false;
        fi;
        for eEXT in answer.points
        do
            if eEXT*eRec.F*eEXT>=0 then
                Print("The point ", eEXT, " is not of negative value\n");
                return false;
            fi;
        od;
        return true;
    fi;
    for eEXT in answer.EXT
    do
        if eEXT*eRec.F*eEXT<>0 then
            Print("The point ", eEXT, " is not a zero of the function\n");
            return false;
        fi;
    od;
    if answer.min<>eRec.min or Length(answer.EXT)<>eRec.n_ext or Length(answer.L)<>eRec.rank_L then
        Print("min=", answer.min, " |EXT|=", Length(answer.EXT), " rank(L)=", Length(answer.L),
              " expected ", eRec.min, " ", eRec.n_ext, " ", eRec.rank_L, "\n");
        return false;
    fi;
    return true;
end;

# The function sum_i x_i (x_i - 1), whose zero set is the cube {0,1}^n.
CubeFunction:=function(n)
    local F, i;
    F:=NullMat(n+1, n+1);
    for i in [1..n]
    do
        F[i+1][i+1]:=1;
        F[1][i+1]:=-1/2;
        F[i+1][1]:=-1/2;
    od;
    return F;
end;

ListCase:=[
  # The indefinite form x1^2 - x2^2.
  rec(name:="indefinite", F:=[[0,0,0],[0,1,0],[0,0,-1]], nonnegative:=false),
  # x1^2 + x2: the linear part does not vanish on the kernel.
  rec(name:="linear_on_kernel", F:=[[0,0,1/2],[0,1,0],[1/2,0,0]], nonnegative:=false),
  # (x1 - 1/2)^2 + x2^2 - 1: a negative minimum.
  rec(name:="negative_minimum", F:=[[-3/4,-1/2,0],[-1/2,1,0],[0,0,1]], nonnegative:=false),
  # x1^2 + x2^2 + 1: a positive minimum, the zero set is empty.
  rec(name:="positive_minimum", F:=[[1,0,0],[0,1,0],[0,0,1]], nonnegative:=true,
      min:=1, n_ext:=0, rank_L:=0),
  # The zero function: the zero set is Z^3.
  rec(name:="zero_function", F:=NullMat(4,4), nonnegative:=true,
      min:=0, n_ext:=1, rank_L:=3),
  # The strip {0,1} x Z.
  rec(name:="strip", F:=[[0,-1/2,0],[-1/2,1,0],[0,0,0]], nonnegative:=true,
      min:=0, n_ext:=2, rank_L:=1),
  # The square and the cube.
  rec(name:="square", F:=CubeFunction(2), nonnegative:=true,
      min:=0, n_ext:=4, rank_L:=0),
  rec(name:="cube", F:=CubeFunction(3), nonnegative:=true,
      min:=0, n_ext:=8, rank_L:=0)
];

# All the cases are run, so that every failure is reported.
test:=ForAll(List(ListCase, TestZeroSet), x->x);
ConcludeErdahl(test);
