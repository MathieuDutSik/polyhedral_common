CI_Decision_Reset:=function()
    local FileName;
    FileName:="CI_CONCLUSION";
    if IsExistingFile(FileName) then
        RemoveFile(FileName);
    fi;
end;

CI_Write_Ok:=function()
    local FileName, TheCommand;
    FileName:="CI_CONCLUSION";
    TheCommand:=Concatenation("touch ", FileName);
    Exec(TheCommand);
end;

CI_PrintExistConclusion:=function()
    local FileName, test;
    FileName:="CI_CONCLUSION";
    test:=IsExistingFile(FileName);
    Print("CI_PrintExistConclusion, FileName=", FileName, " IsExist=", test, "\n");
end;

WriteStringFile:=function(eFile, strOut)
    local output;
    output:=OutputTextFile(eFile, true);
    WriteAll(output, strOut);
    CloseStream(output);
end;


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

ReadLineRed:=function(file)
    local line, n_char, lineRed;
    line:=ReadLine(file);
    n_char:=Length(line) - 1;
    lineRed:=line{[1..n_char]};
    return lineRed;
end;



StreamMatrixFile:=function(file)
    local line, LStr, nbRow, nbCol, TheMat, iRow, eLine, eStr, val;
    line:=ReadLineRed(file);
    LStr:=SplitString(line, " ");
    nbRow:=Int(LStr[1]);
    nbCol:=Int(LStr[2]);
    TheMat:=[];
    for iRow in [1..nbRow]
    do
        line:=ReadLineRed(file);
        LStr:=SplitString(line, " ");
        eLine:=[];
        for eStr in LStr
        do
            if Length(eStr) > 0 then
                val:=Rat(eStr);
                Add(eLine, val);
            fi;
        od;
        if Length(eLine)<>nbCol then
            Print("|eLine|=", Length(eLine), " nbCol=", nbCol, "\n");
            Error("Inconsistent length of vector");
        fi;
        Add(TheMat, eLine);
    od;
    return TheMat;
end;

ReadMatrixFile:=function(eFile)
    local file, TheMat;
    file:=InputTextFile(eFile);
    TheMat:=StreamMatrixFile(file);
    CloseStream(file);
    return TheMat;
end;


ReadListMatrixFile:=function(eFile)
    local file, line, n_matrix, ListMat, iMat, TheMat;
    file:=InputTextFile(eFile);
    line:=ReadLineRed(file);
    n_matrix:=Int(line);
    ListMat:=[];
    for iMat in [1..n_matrix]
    do
        TheMat:=StreamMatrixFile(file);
        Add(ListMat, TheMat);
    od;
    CloseStream(file);
    return ListMat;
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

ReadTextFile:=function(FileIn)
    local ListLines, file, TheRead, len, TheReadRed;
    ListLines:=[];
    if IsExistingFile(FileIn)= false then
        Print("ReadTextFile error, the file FileIn=", FileIn, " is missing\n");
        Error("Correct your stuff");
    fi;
    file := InputTextFile(FileIn);
    while(true)
    do
        TheRead:=ReadLine(file);
        if TheRead=fail then
            break;
        fi;
        len:=Length(TheRead);
        TheReadRed:=TheRead{[1..len-1]};
        Add(ListLines, TheReadRed);
    od;
    CloseStream(file);
    return ListLines;
end;

GetFreeFile:=function(prefix)
    local iFile, eFile;
    iFile:=1;
    while(true)
    do
        eFile:=Concatenation(prefix, String(iFile));
        if IsExistingFile(eFile)=false then
            return eFile;
        fi;
        iFile:=iFile + 1;
    od;
end;

ListFileDirectory:=function(TheDir)
    local FileOUT, TheCommand, ListFiles;
    FileOUT:=Filename(DirectoryTemporary(), "LSfile");
    TheCommand:=Concatenation("ls ", TheDir, " > ", FileOUT);
    Exec(TheCommand);
    ListFiles:=ReadTextFile(FileOUT);
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
    local ListChar, ListVal, i, GetNumber, eChar, scal, eStrRed;
    ListChar:=[];
    ListVal:=[];
    for i in [0..9]
    do
        Add(ListChar, String(i));
        Add(ListVal, i);
    od;
    GetNumber:=function(eChar)
        local pos;
        pos:=Position(ListChar, eChar);
        if pos=fail then
            return fail;
        fi;
        return ListVal[pos];
    end;
    eChar:=eStr{[1]};
    scal:=GetNumber(eChar);
    if scal=fail then
        return ClassicalSporadicLattices(eStr);
    fi;
    eStrRed:=eStr{[2..Length(eStr)]};
    return scal * ClassicalSporadicLattices(eStrRed);
end;





GetGramMatrixFromListGram:=function(ListGram)
    local dim, GramMat, shift, eGram, locdim, u, v;
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
    return GetGramMatrixFromListGram(ListGram);
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

GetFundamentalInfo:=function(d)
  local res, IsCorrect, eSum, eProd, Dval, eQuot, type_tspace;
  res:=d mod 4;
  IsCorrect:=false;
  eSum:=0;
  eProd:=0;
  if res=0 then
    eSum:=0;
    eProd:=-d/4;
    Dval:=-eProd;
    IsCorrect:=true;
    if IsInt((Dval-1)/4)=true or IsInt(Dval/4)=true then
      IsCorrect:=false;
    fi;
  fi;
  if res=1 then
    eQuot:=(1-d)/4;
    eSum:=1;
    eProd:=eQuot;
    IsCorrect:=true;
  fi;
  if d < 0 then
      type_tspace:="ImagQuad";
  else
      type_tspace:="RealQuad";
  fi;
  return rec(eSum:=eSum, eProd:=eProd, IsCorrect:=IsCorrect, type_tspace:=type_tspace);
end;

CorrectnessRealQuadratic:=function(eSum, eProd)
  local TheDiscriminant, ListMult, pos, TheRes;
  if eSum=0 then
    TheDiscriminant:=-eProd;
  else
    TheDiscriminant:=eSum*eSum -4*eProd;
  fi;
  TheRes:=TheDiscriminant mod 4;
  if TheDiscriminant<=0 then
    Print("eSum=", eSum, " eProd=", eProd, " The Field is not real quadratic, impossible to work\n");
    return false;
  fi;
  ListMult:=List(Collected(FactorsInt(TheDiscriminant)), x->x[2]);
  pos:=First(ListMult, x->x mod 2=0);
  if pos<>fail then
    Print("TheDiscriminant=", TheDiscriminant, "\n");
    Print("The Discriminant has a square factor, illegal\n");
    return false;
  fi;
  return true;
end;

