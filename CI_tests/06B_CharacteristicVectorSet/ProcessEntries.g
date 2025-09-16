Read("../common.g");
Print("Beginning TestDelaunayEnumeration\n");


TestGeneration:=function(matrix, method)
    local n, TmpDir, FileMat, FileOut, eProg, TheCommand, U;
    TmpDir:=DirectoryTemporary();
    FileMat:=Filename(TmpDir, "GenerateVectorSet.mat");
    FileOut:=Filename(TmpDir, "GenerateVectorSet.out");
    #
    WriteMatrixFile(FileMat, matrix);
    #
    eProg:="../../src_latt/LATT_GenerateCharacteristicVectorSet";
    TheCommand:=Concatenation(eProg, " mpq_class mpz_class ", method, " ", FileMat, " GAP ", FileOut);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    U:=ReadAsFunction(FileOut)();
    RemoveFile(FileMat);
    RemoveFile(FileOut);
    return true;
end;

ListRec:=ReadAsFunction("ListCases")();;
ListMethod:=["shortest", "relevant_voronoi", "filtered_relevant_voronoi", "fullrank", "spanning"];


FullTest:=function()
    local n_error, iRec, eRec, RecReply, method;
    iRec:=0;
    for eRec in ListRec
    do
        Print("----------------------------------------------------------------------------\n");
        Print("iRec=", iRec, "/", Length(ListRec), " Treating lattice named ", eRec.name, "\n");
        n_meth:=Length(ListMethod);
        for i_meth in [1..n_meth]
        do
            method:=ListMethod[i_meth];
            Print("    i_meth=", i_meth, "/", n_meth, " method=", method, "\n");
            result:=TestGeneration(eRec.eG, method);
            if result=false then
                return false;
            fi;
        od;
        iRec:=iRec+1;
    od;
    return true;
end;

result:=FullTest();
Print("2: result=", result, "\n");
CI_Decision_Reset();
if result=false then
    # Error case
    Print("Error case\n");
else
    # No error case
    Print("Normal case\n");
    CI_Write_Ok();
fi;
CI_PrintExistConclusion();

