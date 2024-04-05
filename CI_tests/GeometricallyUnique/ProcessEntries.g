Read("../common.g");

TheDir:="AllExamples";
ListFiles:=ListFileDirectory(TheDir);

Print("ListFiles=", ListFiles);

eProg:="../../src_poly/TEST_GeometricallyUniquePoint";
iFile:=0;
n_error:=0;
for eFile in ListFiles
do
    eFileExt:=Concatenation(TheDir, "/", eFile);
    strFile:=Concatenation("result_", String(iFile));
    FileResult:=Filename(DirectoryTemporary(), strFile);
    TheCommand:=Concatenation(eProg, " ", eFileExt, " ", FileResult);
    Exec(TheCommand);
    TheResult:=ReadAsFunction(FileResult)();
    if TheResult=false then
        n_error:=n_error + 1;
    fi;
    RemoveFileIfExist(FileResult);
    iFile:=iFile + 1;
    Print("iFile=", iFile, " n_error=", n_error, "\n");
od;
if n_error > 0 then
    # Error case
    GAP_EXIT_CODE(1);
else
    # No error case
    GAP_EXIT_CODE(0);
fi;
