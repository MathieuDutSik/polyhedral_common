WriteMatrix:=function(output, EXT)
    local eEXT, eVal;
    AppendTo(output, Length(EXT), " ", Length(EXT[1]), "\n");
    for eEXT in EXT
    do
        for eVal in eEXT
        do
            AppendTo(output, " ", eVal);
        od;
        AppendTo(output, "\n");
    od;
end;

WriteMatrixFile:=function(eFile, EXT)
    local output, eEXT, eVal;
    output:=OutputTextFile(eFile, true);
    WriteMatrix(output, EXT);
    CloseStream(output);
end;

WriteListMatrixFile:=function(eFile, ListMat)
    local output, eMat;
    output:=OutputTextFile(eFile, true);
    AppendTo(output, Length(ListMat), "\n");
    for eMat in ListMat
    do
        WriteMatrix(output, eMat);
    od;
    CloseStream(output);
end;

WriteGroupFile:=function(eFile, n, GRP)
  local ListGen, output, eGen, i, j;
  ListGen:=GeneratorsOfGroup(GRP);
  output:=OutputTextFile(eFile, true);
  AppendTo(output, n, " ", Length(ListGen), "\n");
  for eGen in ListGen
  do
    for i in [1..n]
    do
      j:=OnPoints(i, eGen);
      if j>n then
          Error("We have j=", j, " but n=", n);
      fi;
      AppendTo(output, " ", j-1);
    od;
    AppendTo(output, "\n");
  od;
  CloseStream(output);
end;


WriteVector:=function(output, V)
    local eVal;
    AppendTo(output, Length(V), "\n");
    for eVal in V
    do
        AppendTo(output, " ", eVal);
    od;
    AppendTo(output, "\n");
end;

IsIntegralMat:=function(eMat)
  local eLine, eVal;
  for eLine in eMat
  do
    for eVal in eLine
    do
      if IsInt(eVal)=false then
        return false;
      fi;
    od;
  od;
  return true;
end;

RemoveFileIfExist:=function(FileName)
    if IsExistingFile(FileName) then
        RemoveFile(FileName);
    fi;
end;

SaveDataToFile:=function(FileName, OBJ)
  local output;
  Exec("rm -f ", FileName,"\n");
  output:=OutputTextFile(FileName, true);;
  AppendTo(output, "return ", OBJ, ";\n");
  CloseStream(output);
end;

ListFileDirectory:=function(TheDir)
    local FileOUT, TheCommand, ListFiles, file, TheRead, len, TheReadRed;
    FileOUT:=Filename(DirectoryTemporary(), "LSfile");
    TheCommand:=Concatenation("ls ", TheDir, " > ", FileOUT);
    Exec(TheCommand);
    ListFiles:=[];
    file := InputTextFile(FileOUT);
    while(true)
    do
        TheRead:=ReadLine(file);
        if TheRead=fail then
            break;
        fi;
        len:=Length(TheRead);
        TheReadRed:=TheRead{[1..len-1]};
        Add(ListFiles, TheReadRed);
    od;
    CloseStream(file);
    RemoveFileIfExist(FileOUT);
    return ListFiles;
end;



ClassicalSporadicLattices:=function(TheName)
  local ListNames, ListGram, len, i;
  ListNames:=[];
  ListGram:=[];
  #
  # root lattices
  #
  Add(ListNames, "E6");
  Add(ListGram,
[[2,1,1,0,1,1],[1,4,1,1,1,3],[1,1,2,1,1,1],
 [0,1,1,2,1,2],[1,1,1,1,2,2],[1,3,1,2,2,4]]);
  #
  Add(ListNames, "E7");
  Add(ListGram,
[[4,3,2,3,1,1,1],[3,4,3,3,1,1,1],
[2,3,4,3,2,2,1],[3,3,3,4,2,1,2],[1,1,2,2,2,1,1],
[1,1,2,1,1,2,0],[1,1,1,2,1,0,2]]);
  #
  Add(ListNames, "A2");
  Add(ListGram,
[ [ 2, -1 ], [ -1, 2 ] ]);
  #
  Add(ListNames, "E8");
  Add(ListGram, [[2,0,0,0,0,0,1,0],[0,2,-1,0,-1,1,0,1],
                 [0,-1,2,1,0,0,0,-1],[0,0,1,2,0,0,1,-1],
                 [0,-1,0,0,2,-1,1,-1],[0,1,0,0,-1,2,0,0],
                 [1,0,0,1,1,0,2,-1],[0,1,-1,-1,-1,0,-1,2]]);
  #
  Add(ListNames, "A4");
  Add(ListGram, [ [ 2, -1, 0, 0 ], [ -1, 2, -1, 0 ], [ 0, -1, 2, -1 ], [ 0, 0, -1, 2 ] ]);
  #
  Add(ListNames, "A5");
  Add(ListGram, [ [ 2, -1, 0, 0, 0 ], [ -1, 2, -1, 0, 0 ], [ 0, -1, 2, -1, 0 ], [ 0, 0, -1, 2, -1 ], [ 0, 0, 0, -1, 2 ] ]);
  #
  Add(ListNames, "A6");
  Add(ListGram, [ [ 2, -1, 0, 0, 0, 0 ], [ -1, 2, -1, 0, 0, 0 ], [ 0, -1, 2, -1, 0, 0 ], [ 0, 0, -1, 2, -1, 0 ], [ 0, 0, 0, -1, 2, -1 ], [ 0, 0, 0, 0, -1, 2 ] ]);
  #
  Add(ListNames, "D4");
  Add(ListGram, [ [ 2, -1, 0, 0 ], [ -1, 2, -1, -1 ], [ 0, -1, 2, 0 ], [ 0, -1, 0, 2 ] ]);
  #
  Add(ListNames, "D5");
  Add(ListGram,  [ [ 2, -1, 0, 0, 0 ], [ -1, 2, -1, 0, -1 ], [ 0, -1, 2, -1, 0 ], [ 0, 0, -1, 2, 0 ], [ 0, -1, 0, 0, 2 ] ]);
  #
  Add(ListNames, "D6");
  Add(ListGram, [ [ 2, -1, 0, 0, 0, 0 ], [ -1, 2, -1, 0, 0, -1 ], [ 0, -1, 2, -1, 0, 0 ], [ 0, 0, -1, 2, -1, 0 ], [ 0, 0, 0, -1, 2, 0 ],[ 0, -1, 0, 0, 0, 2 ] ]);
  #
  len:=Length(ListNames);
  for i in [1..len]
  do
      if ListNames[i] = TheName then
          return ListGram[i];
      fi;
  od;
  Error("Failed to find the entry in the database");
end;


RandomIntegralUnimodularMatrix:=function(n)
    local GetRandomTriangular, GRP, ePerm, eMat1, eMat2, i, j, eMat3;
    GetRandomTriangular:=function()
        local TheMat, i, j;
        TheMat:=NullMat(n, n);
        for i in [1..n]
        do
            TheMat[i][i]:=Random([-1,1]);
            for j in [i+1..n]
            do
                TheMat[i][j]:=Random([-2..2]);
            od;
        od;
        return TheMat;
    end;
    GRP:=SymmetricGroup(n);
    ePerm:=Random(GRP);
    eMat1:=GetRandomTriangular();
    eMat2:=NullMat(n,n);
    for i in [1..n]
    do
        j:=OnPoints(i, ePerm);
        eMat2[i][j]:=1;
    od;
    eMat3:=TransposedMat(GetRandomTriangular());
    return eMat1 * eMat2 * eMat3;
end;
