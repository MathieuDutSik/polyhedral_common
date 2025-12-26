Read("../common.g");


AppendDelaunaySimplices:=true;
AppendReflectiveDim45:=false;

ListEXT:=[];

if AppendDelaunaySimplices then
    ListFiles:=[];
    for n in [5..7]
    do
        FileName:=Concatenation("../DATA/ClassificationSimplices", String(n));
        OneBlock:=ReadAsFunction(FileName)();
        Append(ListEXT, OneBlock);
    od;
fi;

if AppendReflectiveDim45 then
    OneBlock:=ReadAsFunction("ListSimpleRootSystem_4_56_X_5_47")();
    Append(ListEXT, OneBlock);
fi;


get_automorphism_result:=function(EXT)
    local TmpDir, FileI, FileO, eProg, TheCommand, rec_gap;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Can.in");
    FileO:=Filename(TmpDir, "Can.out");
    WriteMatrixFile(FileI, EXT);
    eProg:="../../src_group/GRP_LinPolytope_Automorphism";
    TheCommand:=Concatenation(eProg, " rational ", FileI, " RecGAP ", FileO);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        Print("FileO does not exist, qualifies as a fail\n");
        return false;
    fi;
    rec_gap:=ReadAsFunction(FileO)();
    RemoveFileIfExist(FileI);
    RemoveFileIfExist(FileO);
    return rec_gap;
end;


get_canonic_form:=function(EXT)
    local TmpDir, FileI, FileO, eProg, TheCommand, EXT_can;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Can.in");
    FileO:=Filename(TmpDir, "Can.out");
    WriteMatrixFile(FileI, EXT);
    eProg:="../../src_group/GRP_LinPolytope_Canonic";
    TheCommand:=Concatenation(eProg, " mpq_class ", FileI, " GAP ", FileO);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        Print("FileO does not exist, qualifies as a fail\n");
        return false;
    fi;
    EXT_can:=ReadAsFunction(FileO)();
    RemoveFileIfExist(FileI);
    RemoveFileIfExist(FileO);
    return EXT_can;
end;

get_isomorphism_result:=function(EXT1, EXT2)
    local TmpDir, File1, File2, FileO, eProg, TheCommand, result;
    TmpDir:=DirectoryTemporary();
    File1:=Filename(TmpDir, "Test1.in");
    File2:=Filename(TmpDir, "Test2.in");
    FileO:=Filename(TmpDir, "Test.out");
    WriteMatrixFile(File1, EXT1);
    WriteMatrixFile(File2, EXT2);
    eProg:="../../src_group/GRP_LinPolytope_Isomorphism";
    TheCommand:=Concatenation(eProg, " rational ", File1, " ", File2, " GAP ", FileO);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        Print("FileO does not exist, qualifies as a fail\n");
        return false;
    fi;
    result:=ReadAsFunction(FileO)();
    RemoveFileIfExist(File1);
    RemoveFileIfExist(File2);
    RemoveFileIfExist(FileO);
    return result;
end;


scramble_ext:=function(EXT)
    local n, n_vert, P;
    n:=Length(EXT[1]);
    n_vert:=Length(EXT);
    P:=RandomIntegralUnimodularMatrix(n);
    return Permuted(EXT * P, Random(SymmetricGroup(n_vert)));
end;

TestCase_Automorphy:=function(EXT)
    local EXT_img, rec_gap1, rec_gap2;
    EXT_img:=scramble_ext(EXT);
    rec_gap1:=get_automorphism_result(EXT);
    if rec_gap1=false then
        return false;
    fi;
    rec_gap2:=get_automorphism_result(EXT_img);
    if rec_gap2=false then
        return false;
    fi;
    if Order(rec_gap1.GAPperm)<>Order(rec_gap2.GAPperm) then
        return false;
    fi;
    return true;
end;

TestCase_Canonic:=function(EXT)
    local EXT_img, EXT_can1, EXT_can2;
    EXT_img:=scramble_ext(EXT);
    EXT_can1:=get_canonic_form(EXT);
    if EXT_can1=false then
        return false;
    fi;
    EXT_can2:=get_canonic_form(EXT_img);
    if EXT_can2=false then
        Print("No result, returning false\n");
        return false;
    fi;
    if EXT_can1<>EXT_can2 then
        Print("Different canonical form\n");
        return false;
    fi;
    return true;
end;




TestCase_Isomorphy:=function(EXT)
    local EXT_img, result;
    EXT_img:=scramble_ext(EXT);
    result:=get_isomorphism_result(EXT, EXT_img);
    if result=false then
        # No result obtained
        return false;
    fi;
    if result=fail then
        # Found to be non-isomorphic
        return false;
    fi;
    return true;
end;


#

n_error:=0;
nEXT:=Length(ListEXT);
for iEXT in [1..nEXT]
do
    EXT:=ListEXT[iEXT];
    Print("---------------------------------------- ", iEXT, "/", nEXT, " ", n_error, " ----------------------------------------\n");
    test:=TestCase_Automorphy(EXT);
    if test=false then
        n_error:=n_error+1;
    fi;
    Print("----\n");
    #
    test:=TestCase_Canonic(EXT);
    if test=false then
        n_error:=n_error+1;
    fi;
    Print("----\n");
    #
    test:=TestCase_Isomorphy(EXT);
    if test=false then
        n_error:=n_error+1;
    fi;
    Print("----\n");
od;
Print("-------------------------------------------------------\n");
Print("n_error=", n_error, "\n");
CI_Decision_Reset();
if n_error > 0 then
    Print("Error case\n");
else
    Print("Normal case\n");
    CI_Write_Ok();
fi;
