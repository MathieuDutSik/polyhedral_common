


Rank2_Un:=function(n)
    return [[0, n], [n, 0]];
end;

Rank2_Unk:=function(n, k)
    return [[0, n], [n, 2*k]];
end;


Rank3_Theorem4_3:=function()
    local ListM, A2, n, U;
    ListM:=[];
    A2:=ClassicalSporadicLattices("A2");
    for n in [5,8,10,11,14,15,20,24,26,30,42,48,60,90]
    do
        U:=[[2*n]];
        eM:=GetGramMatrixFromListGram([U, A2]);
        Add(ListM, eM);
    od;
    #
    A1:=ClassicalSporadicLattices("A1");
    for n in [6,7,12,15,24,36]
    do
        U:=[[2*n]];
        eM:=GetGramMatrixFromListGram([U, A1, A1]);
        Add(ListM, eM);
    od;
    #
    for l in [11,13,15,21,23,33]
    do
        eM:=[[2*l, 1, 0], [1, -2, 1], [0, 1, -2]];
        Add(ListM, eM);
    od;
    #
    for l in [5,12]
    do
        eM:=[[2*l, 1, 0], [1, -2, 0], [0, 0, -2]];
        Add(ListM, eM);
    od;
    #
    for l in [3,7]
    do
        eM:=[[2*l, 1, 1], [1, -2, 0], [1, 0, -2]];
        Add(ListM, eM);
    od;
    #
    for n in [11,14,15,20,24]
    do
        A:=Rank2_Un(n);
        B:=[[-2]];
        eM:=GetGramMatrixFromListGram([A, B]);
        Add(ListM, eM);
    od;
    #
    for p in [[8,2],[8,3],[8,6],[9,3],[12,3],[12,6],[12,8],[16,8]]
    do
        n:=p[1];
        k:=p[2];
        A:=Rank2_Unk(n, k);
        B:=[[-2]];
        eM:=GetGramMatrixFromListGram([A, B]);
        Add(ListM, eM);
    od;
    #
    eM:=[[-2, 3, 3], [3, -2, 2], [3, 2, -2]];
    Add(ListM, eM);
    eM:=[[-2, 3, 7], [3, -2, 2], [7, 2, -2]];
    Add(ListM, eM);
    #
    return ListM;
end;


Rank4_Theorem4_5:=function()
    A:=Rank2_Un(6);
    A1:=ClassicalSporadicLattices("A1");
    eM1:=GetGramMatrixFromListGram([A, A1, A1]);
    #
    eM2:=[[0, 4, 0, 0], [4, -6, 1, 1], [0, 1, -2, 0], [0, 1, 0, -2]];
    eM3:=[[0, 5, 0, 0], [5, -6, 2, 2], [0, 2, -2, 1], [0, 2, 1, -2]];
    return [eM1, eM2, eM3];
end;

