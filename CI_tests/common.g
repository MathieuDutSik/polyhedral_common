WriteMatrixFile:=function(eFile, EXT)
    local output, eEXT, eVal;
    output:=OutputTextFile(eFile, true);
    AppendTo(output, Length(EXT), " ", Length(EXT[1]), "\n");
    for eEXT in EXT
    do
        for eVal in eEXT
        do
            AppendTo(output, " ", eVal);
        od;
        AppendTo(output, "\n");
    od;
    CloseStream(output);
end;



RemoveFileIfExist:=function(FileName)
    if IsExistingFile(FileName) then
        RemoveFile(FileName);
    fi;
end;


