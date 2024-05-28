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

WriteGroupFile:=function(eFile, n, GRP)
  local ListGen, eGen, i, j;
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
