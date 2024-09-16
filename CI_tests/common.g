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

STRING_Split:=function(eStr, ePat)
  local lenStr, lenPat, lenTest, ListMatch, i, eStrPart, ListStrInter, iVal, eVal, FirstVal, nbMatch, eStrInter;
  lenStr:=Length(eStr);
  lenPat:=Length(ePat);
  ListMatch:=[];
  lenTest:=lenStr - (lenPat-1);
  for i in [1..lenTest]
  do
    eStrPart:=eStr{[i..i+lenPat-1]};
    if eStrPart=ePat then
      Add(ListMatch, i);
    fi;
  od;
  nbMatch:=Length(ListMatch);
  if nbMatch=0 then
    return rec(ListStrInter:=[eStr], ListMatch:=[]);
  fi;
  Print("nbMatch=", nbMatch, "\n");
  ListStrInter:=[];
  for iVal in [1..nbMatch]
  do
    eVal:=ListMatch[iVal];
    if eVal>1 then
      if iVal>1 then
        FirstVal:=ListMatch[iVal-1]+lenPat;
      else
        FirstVal:=1;
      fi;
      eStrInter:=eStr{[FirstVal..eVal-1]};
      Add(ListStrInter, eStrInter);
    fi;
    if iVal=nbMatch then
      if eVal+lenPat<lenStr then
        FirstVal:=eVal+lenPat;
        eStrInter:=eStr{[FirstVal..lenStr]};
        Add(ListStrInter, eStrInter);
      fi;
    fi;
  od;
  return rec(ListStrInter:=ListStrInter, ListMatch:=ListMatch);
end;

ReadMatrixFile:=function(eFile)
    local file, line, eRec, nbRow, nbCol, TheMat, iRow, eLine, eStr, val;
    file:=InputTextFile(eFile);
    line:=ReadLine(file);
    eRec:=SplitString(line, " ");
    nbRow:=Int(eRec[1]);
    nbCol:=Int(eRec[2]);
    TheMat:=[];
    for iRow in [1..nbRow]
    do
        line:=ReadLine(file);
        eRec:=SplitString(line, " ");
        eLine:=[];
        for eStr in eRec.ListStrInter
        do
            if Length(eStr) > 0 then
                val:=Rat(eStr);
                Add(eLine, val);
            fi;
        od;
        if Length(eLine)<>nbCol then
            Error("Inconsistent length of vector");
        fi;
        Add(TheMat, eLine);
    od;
    return TheMat;
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
    local TheNameRed, dim, TheMat, ListPair, ePair, ListNames, ListGram, len, i;
    if TheName[1]='A' then
        TheNameRed:=TheName{[2..Length(TheName)]};
        dim:=Int(TheNameRed);
        if dim<>fail then
            TheMat:=NullMat(dim,dim);
            for i in [1..dim]
            do
                TheMat[i][i]:=2;
            od;
            for i in [1..dim-1]
            do
                TheMat[i][i+1]:=-1;
                TheMat[i+1][i]:=-1;
            od;
            return TheMat;
        fi;
    fi;
    if TheName[1]='D' then
        TheNameRed:=TheName{[2..Length(TheName)]};
        dim:=Int(TheNameRed);
        if dim<>fail then
            TheMat:=NullMat(dim,dim);
            for i in [1..dim]
            do
                TheMat[i][i]:=2;
            od;
            ListPair:=[[1,3],[2,3]];
            for i in [4..dim]
            do
                Add(ListPair, [i-1,i]);
            od;
            for ePair in ListPair
            do
                TheMat[ePair[1]][ePair[2]]:=-1;
                TheMat[ePair[2]][ePair[1]]:=-1;
            od;
            return TheMat;
        fi;
    fi;
    #
    ListNames:=[];
    ListGram:=[];
    #
    # Hyperbolic plane
    #
    Add(ListNames, "U");
    Add(ListGram, [[0,1],[1,0]]);
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
    Add(ListNames, "E8");
    Add(ListGram, [[2,0,0,0,0,0,1,0],[0,2,-1,0,-1,1,0,1],
                   [0,-1,2,1,0,0,0,-1],[0,0,1,2,0,0,1,-1],
                   [0,-1,0,0,2,-1,1,-1],[0,1,0,0,-1,2,0,0],
                   [1,0,0,1,1,0,2,-1],[0,1,-1,-1,-1,0,-1,2]]);
    #
    len:=Length(ListNames);
    for i in [1..len]
    do
        if ListNames[i] = TheName then
            return ListGram[i];
        fi;
    od;
    Print("TheName=", TheName, "\n");
    Error("Failed to find the entry in the database");
end;

GetGramMatrixFromString:=function(eStr)
    local ListChar, i, IsNumber, scal, eStrRed;
    ListChar:=[];
    for i in [0..9]
    do
        Add(ListChar, String(i));
    od;
    IsNumber:=function(eChar)
        local fChar;
        for fChar in ListChar
        do
            if fChar=eChar then
                return true;
            fi;
        od;
        return false;
    end;
    if IsNumber(eStr[1])=false then
        return ClassicalSporadicLattices(eStr);
    fi;
    scal:=Number(eStr[1]);
    eStrRed:=eStr{[2..Length(eStr)]};
    return scal * ClassicalSporadicLattices(eStr);
end;



GetGramMatrixFromList:=function(eList)
    local ListGram, entry, dim, GramMat, shift, eGram, locdim, u, v;
    ListGram:=[];
    for entry in eList
    do
        if IsString(entry) then
            Add(ListGram, GetGramMatrixFromString(entry));
        else
            Add(ListGram, entry);
        fi;
    od;
    dim:=Sum(List(ListGram, Length));
    GramMat:=NullMat(dim, dim);
    shift:=0;
    for eGram in ListGram
    do
        locdim:=Length(eGram);
        for u in [1..locdim]
        do
            for v in [1..locdim]
            do
                GramMat[shift+u][shift+v]:=eGram[u][v];
            od;
        od;
        shift:=shift+locdim;
    od;
    return GramMat;
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
