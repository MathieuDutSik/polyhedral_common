Read("../common.g");

TheDir:="AllExamples";
ListFiles:=ListFileDirectory(TheDir);

Print("ListFiles=", ListFiles);

eProg:="../../src_poly/TEST_GeometricallyUniquePoint";
iFile:=0;
n_error:=0;
n_file:=Length(ListFiles);
for eFile in ListFiles
do
    eFileExt:=Concatenation(TheDir, "/", eFile);
    strFile:=Concatenation("result_", String(iFile));
    FileResult:=Filename(DirectoryTemporary(), strFile);
    TheCommand:=Concatenation(eProg, " ", eFileExt, " ", FileResult);
    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    if IsExistingFile(FileResult)=false then
        Print("The FileResult does not exist\n");
        n_error:=n_error + 1;
    else
        TheResult:=ReadAsFunction(FileResult)();
        if TheResult=false then
            Print("The result is not correct\n");
            n_error:=n_error + 1;
        fi;
        RemoveFileIfExist(FileResult);
    fi;
    iFile:=iFile + 1;
    Print("iFile=", iFile, "/", n_file, " n_error=", n_error, "\n");
od;
CI_Decision_Reset();
if n_error > 0 then
    # Error case
    Print("Error case\n");
else
    # No error case
    Print("Normal case\n");
    CI_Write_Ok();
fi;
