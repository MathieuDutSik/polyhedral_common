Read("../common.g");
Print("Beginning Test enumeration of perfect forms in T-spaces\n");

GetFundamentalInfo:=function(d)
  local res, IsCorrect, eSum, eProd, Dval, eQuot;
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
  return rec(eSum:=eSum, eProd:=eProd, IsCorrect:=IsCorrect);
end;

GetRecInfo:=function(fProg, d, n)
    local FileNml, FileResult, output, eProg, TheCommand, U, is_correct, info;
    #
    FileNml:="PerfectForm.nml";
    FileResult:="Result";
    RemoveFileIfExist(FileNml);
    RemoveFileIfExist(FileResult);
    #
    info:=GetFundamentalInfo(d);
    if info.IsCorrect=false then
        Print("Discriminant d=", d, " is not valid, skipping\n");
        return fail;
    fi;
    #
    # Create the namelist file
    output:=OutputTextFile(FileNml, true);
    AppendTo(output, "&DATA\n");
    AppendTo(output, " arithmetic_T=\"gmp_rational\"\n");
    AppendTo(output, " arithmetic_Tint=\"gmp_integer\"\n");
    AppendTo(output, " OutFormat=\"ObjectGAP\"\n");
    AppendTo(output, " OutFile=\"", FileResult, "\"\n");
    AppendTo(output, " max_runtime_second=0\n");
    AppendTo(output, " ApplyStdUnitbuf=T\n");
    AppendTo(output, "/\n");
    AppendTo(output, "\n");
    AppendTo(output, "&STORAGE\n");
    AppendTo(output, " Saving=F\n");
    AppendTo(output, " Prefix=\"/tmp/PerfectForm/\"\n");
    AppendTo(output, "/\n");
    AppendTo(output, "\n");
    AppendTo(output, "&TSPACE\n");
    AppendTo(output, " TypeTspace=\"ImagQuad\"\n");
    AppendTo(output, " FileLinSpa=\"unset.linspa\"\n");
    AppendTo(output, " SuperMatMethod=\"NotNeeded\"\n");
    AppendTo(output, " ListComm=\"Trivial\"\n");
    AppendTo(output, " PtGroupMethod=\"Trivial\"\n");
    AppendTo(output, " FileListSubspaces=\"unset\"\n");
    AppendTo(output, " RealImagDim=", n, "\n");
    AppendTo(output, " RealImagSum=", info.eSum, "\n");
    AppendTo(output, " RealImagProd=", info.eProd, "\n");
    AppendTo(output, "/\n");
    CloseStream(output);
    #
    eProg:=Concatenation("../../src_perfect/", fProg);
    TheCommand:=Concatenation(eProg, " ", FileNml);
    Exec(TheCommand);
    #
    if IsExistingFile(FileResult)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
#        Print(NullMat(5));
        return fail;
    fi;
    U:=ReadAsFunction(FileResult)();
    RemoveFile(FileNml);
    RemoveFile(FileResult);
    return U;
end;

ListCases:=[];

Add(ListCases, rec(n:=3, d:=-3, n_perf:=1));
Add(ListCases, rec(n:=3, d:=-4, n_perf:=2));
Add(ListCases, rec(n:=3, d:=-7, n_perf:=3));
Add(ListCases, rec(n:=3, d:=-8, n_perf:=5));
Add(ListCases, rec(n:=3, d:=-11, n_perf:=8));
Add(ListCases, rec(n:=3, d:=-15, n_perf:=34));
Add(ListCases, rec(n:=3, d:=-19, n_perf:=43));
Add(ListCases, rec(n:=3, d:=-20, n_perf:=69));
Add(ListCases, rec(n:=3, d:=-23, n_perf:=204));
Add(ListCases, rec(n:=3, d:=-24, n_perf:=158));
Add(ListCases, rec(n:=4, d:=-3, n_perf:=2));
Add(ListCases, rec(n:=4, d:=-4, n_perf:=4));


ListProg:=["PERF_SerialEnumeratePerfectCones", "PERF_MPI_EnumeratePerfectCones"];


f_compute:=function()
    local eCase, fProg, eRec;
    for eCase in ListCases
    do
        for fProg in ListProg
        do
            eRec:=GetRecInfo(fProg, eCase.d, eCase.n);
            if eRec=fail then
                return false;
            fi;
            if eRec.nb<>eCase.n_perf then
                return false;
            fi;
        od;
    od;
    return true;
end;



test:=f_compute();
CI_Decision_Reset();
if test=false then
    Print("Error case\n");
else
    Print("Normal case\n");
    CI_Write_Ok();
fi;
