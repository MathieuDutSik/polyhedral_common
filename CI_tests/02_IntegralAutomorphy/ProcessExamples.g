Read("../common.g");


AppendDelaunaySimplices:=true;
AppendReflectiveDim45:=true;

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

GeneratorsPreservePolytope:=function(TheGRP,EXT)
    local gens, BasisExt, i, j, ImageExt, Mb, Mi, A;
    BasisEXT:=[EXT[1]];
    BasisIndices:=[1];
    n:=Length(EXT[1]);
    i:=1;
    for i in [2..n] do
        for j in [1..Length(EXT)] do
            M:=Concatenation(BasisEXT,[EXT[j]]);
            if RankMat(M)=Length(BasisEXT)+1 then
                Add(BasisEXT, EXT[j]);
                Add(BasisIndices, j);
                break;
            fi;
        od;
    od;
    gens:=GeneratorsOfGroup(TheGRP);
    gensCheck:=[];
    for g in gens do
        ImageIndices:=List(BasisIndices, i -> i^g);
        ImageEXT := List(ImageIndices, i -> EXT[i]);
        Mb:=TransposedMat(BasisEXT);
        Mi:=TransposedMat(ImageEXT);
        A:=Mi*Inverse(Mb);
        if ForAll(Flat(A), x -> x=Int(x)) then
            continue;
        else
            return false;
        fi;
    od;
    return true;
end;

TestCase_Automorphy:=function(EXT)
    local TmpDir, FileI, FileO, arith, eProg, TheCommand, TheGRP;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    WriteMatrixFile(FileI, EXT);
    eProg:="../../build/GRP_LinPolytopeIntegral_Automorphism";
    TheCommand:=Concatenation(eProg, " mpz_class ", FileI, " GAP ", FileO);
    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        RemoveFileIfExist(FileI);
        Print("FileO does not exist\n");
        return false;
    fi;
    TheGRP:=ReadAsFunction(FileO)();
    RemoveFileIfExist(FileI);
    RemoveFileIfExist(FileO);
    if not GeneratorsPreservePolytope(TheGRP,EXT) then
        return false;
    fi;
    Print("|TheGRP|=", Order(TheGRP), "\n");
    return true;
end;

TestCase_Isomorphy:=function(EXT)
    local n, n_vert, P, EXT_img, TmpDir, File1, File2, FileO, arith, eProg, TheCommand, TheGRP;
    n:=Length(EXT[1]);
    n_vert:=Length(EXT);
    P:=RandomIntegralUnimodularMatrix(n);
    EXT_img:=Permuted(EXT * P, Random(SymmetricGroup(n_vert)));
    TmpDir:=DirectoryTemporary();
    File1:=Filename(TmpDir, "Test1.in");
    File2:=Filename(TmpDir, "Test2.in");
    FileO:=Filename(TmpDir, "Test.out");
    WriteMatrixFile(File1, EXT);
    WriteMatrixFile(File2, EXT_img);
    eProg:="../../build/GRP_LinPolytopeIntegral_Isomorphism";
    TheCommand:=Concatenation(eProg, " mpz_class ", File1, " ", File2, " GAP ", FileO);
    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        RemoveFileIfExist(File1);
        RemoveFileIfExist(File2);
        Print("FileO does not exist\n");
        return false;
    fi;
    result:=ReadAsFunction(FileO)();
    RemoveFileIfExist(File1);
    RemoveFileIfExist(File2);
    RemoveFileIfExist(FileO);
    if result=fail then
        return false;
    fi;
    return true;
end;

TestCase_Automorphy_RightCoset:=function(EXT)
    local TmpDir, FileI, FileO, arith, OutFormat, eProg, TheCommand, TheGRP;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    WriteMatrixFile(FileI, EXT);
    eProg:="../../build/GRP_LinPolytopeIntegral_Automorphism_RightCoset";
    TheCommand:=Concatenation(eProg, " mpz_class ", FileI, " GAP ", FileO);
    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        RemoveFileIfExist(FileI);
        Print("FileO does not exist\n");
        return false;
    fi;
    TheGRP:=ReadAsFunction(FileO)();
    RemoveFileIfExist(FileI);
    RemoveFileIfExist(FileO);
    Print("|TheGRP|=", Order(TheGRP), "\n");
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
    test:=TestCase_Isomorphy(EXT);
    if test=false then
        n_error:=n_error+1;
    fi;
    Print("----\n");
    #
    test:=TestCase_Automorphy_RightCoset(EXT);
    if test=false then
        n_error:=n_error+1;
    fi;
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
