Read("../common.g");
Read("../access_points.g");
Print("Beginning contracting homotopy tests\n");

random_grp_element:=function(l_gens)
    local n_gen, n_prod, eMat, i;
    n_gen:=Length(l_gens);
    n_prod:=2 * n_gen + Random([1..n_gen]);

    eMat:=Random(l_gens);

    for i in [2..n_prod]
    do
        eMat:=eMat * Random(l_gens);
    od;
    return eMat;
end;




test_series:=function(desc)
    local l_gens, GRPmatr, index, list_cells, eElt, EXT, EXT2, e_equiv, rec_dim, is_face, is_well_rounded, ffs, list_lower_cell1, list_lower_cell2, eStab;
    #
    l_gens:=PERFCOMP_group_generators(desc);
    Print("|l_gens|=", Length(l_gens), "\n");
    #
    index:=0;
    list_cells:=PERFCOMP_get_cells(desc, index);
    Print("list_cells=", list_cells, "\n");
    EXT:=list_cells[1].EXT;
    eElt:=random_grp_element(l_gens);
    EXT2:=EXT * eElt;
    e_equiv:=PERFCOMP_test_equivalence(desc, EXT, EXT2);
    if e_equiv=fail then
        Print("They should be equivalent\n");
        return false;
    fi;
    Print("e_equiv=", e_equiv, "\n");
    rec_dim:=PERFCOMP_dimension(desc, EXT);
    Print("rec_dim=", rec_dim, "\n");
    if rec_dim.index <> index then
        Print("The indices should be equal\n");
        return false;
    fi;
    is_face:=PERFCOMP_is_face(desc, EXT);
    Print("is_face=", is_face, "\n");
    if is_face=false then
        Print("The face should exist\n");
        return false;
    fi;
    is_well_rounded:=PERFCOMP_is_well_rounded(desc, EXT);
    Print("is_well_rounded=", is_well_rounded, "\n");
    if is_well_rounded=false then
        Print("The face should exist\n");
        return false;
    fi;
    ffs:=PERFCOMP_face_search(desc, EXT2);
    Print("ffs=", ffs, "\n");
    eStab:=PERFCOMP_stabilizer(desc, EXT2);
    Print("eStab=", eStab, "\n");
    list_lower_cell1:=PERFCOMP_lower_boundary_cell(desc, EXT);
    list_lower_cell2:=PERFCOMP_lower_boundary_cell(desc, EXT2);
    Print("|list_lower_cell1|=", Length(list_lower_cell1), " |list_lower_cell2|=", Length(list_lower_cell2), "\n");
    if Length(list_lower_cell1)<>Length(list_lower_cell2) then
        Print("The list of cells should be of equal length\n");
        return false;
    fi;
    EXT_index1:=list_lower_cell2[1];
    list_upper_cell:=PERFCOMP_upper_boundary_cell(desc, EXT_index1);
    Print("|list_upper_cell|=", Length(list_upper_cell), "\n");
    l_ffs:=List(list_upper_cell, x->PERFCOMP_face_search(desc, x));
    Print("l_ffs=", l_ffs, "\n");
    if Length(list_upper_cell)<>2 then
        Print("We should have a length equal to 2\n");
        return false;
    fi;


    #    Print("l_gens=", l_gens, "\n");
    return true;
end;




case1:=GenerateTspaceDescription_classic(4, false);
case2:=GenerateTspaceDescription_imag_quad(3,-7, false);

#ListCases:=[case1, case2];
ListCases:=[case1];


f_compute:=function()
    local n_error, case, test;
    n_error:=0;
    for case in ListCases
    do
        test:=test_series(case);
        if not test then
            n_error:=n_error+1;
        fi;
    od;
    if n_error > 0 then
        return false;
    fi;
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
