Read("../common.g");
Read("../access_points.g");
Print("Beginning TestIsotropic\n");

TestIsotropic:=function(eRec)
    local test_result, find_result, eNorm, eNormSqr;
    test_result:=quadratic_form_is_isotropic(eRec.M);
    if is_error(test_result) then
        return false;
    fi;
    find_result:=quadratic_form_isotropic_vector(eRec.M);
    if is_error(find_result) then
        return false;
    fi;
    if eRec.has_isotropic then
        if test_result<>true then
            Print("We should have found an isotropic vector A");
            return false;
        fi;
        if find_result=fail then
            Print("We should have found an isotropic vector B");
            return false;
        fi;
        eNorm:=find_result * eRec.M * find_result;
        if eNorm<>0 then
            Print("The vector is not isotropic\n");
            return false;
        fi;
        eNormSqr:=find_result * find_result;
        if eNormSqr=0 then
            Print("The vector is zero\n");
            return false;
        fi;
    else
        if test_result<>false then
            Print("We should not have found an isotropic vector A");
            return false;
        fi;
        if find_result<>fail then
            Print("We should not have found an isotropic vector B");
            return false;
        fi;
    fi;
    return true;
end;

ListRec:=ReadAsFunction("../DATA/IsotropicCases")();;
ListDim3:=Filtered(ListRec, x->Length(x.M) = 3);
ListDim4:=Filtered(ListRec, x->Length(x.M) = 4);
ListDim5:=Filtered(ListRec, x->Length(x.M) = 5);
ListDim5gt:=Filtered(ListRec, x->Length(x.M) > 5);
#ListRec:=Concatenation(ListDim3_B, ListDim4_B, ListDim5high);
#ListRec:=Concatenation(ListDim3_B, ListDim5high);


#ListRec:=ListDim3;
#ListRec:=ListDim4;
#ListRec:=ListDim5;
#ListRec:=ListDim5gt;
#ListRec:=ListDim3_A{[3271..3272]};
#ListRec:=ListDim4_A;
#ListRec:=ListDim5;


FullTest:=function()
    local n_error, iRec, eRec, reply;
    n_error:=0;
    iRec:=0;
    for eRec in ListRec
    do
        Print("iRec=", iRec, " / ", Length(ListRec), " n_error=", n_error, "\n");
        reply:=TestIsotropic(eRec);
        if reply=false then
            n_error:=n_error+1;
        fi;
        iRec:=iRec + 1;
    od;
    return n_error;
end;

n_error:=FullTest();
Print("n_error=", n_error, "\n");
CI_Decision_Reset();
if n_error > 0 then
    # Error case
    Print("Error case\n");
else
    # No error case
    Print("Normal case\n");
    CI_Write_Ok();
fi;

