


WriteEntry:=function(FileOut, EXT)
    local output;
    RemoveFileIfExist(FileOut);
    output:=OutputTextFile(FileOut, true);
    CPP_WriteMatrix(output, EXT);
    CloseStream(output);
end;

WriteEntry("Schlafli", ClassicalExtremeDelaunayPolytopes("G6"));
WriteEntry("Gosset", ClassicalExtremeDelaunayPolytopes("G7"));
