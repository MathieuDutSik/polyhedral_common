ePerm1:=(1,2,3,4,5,6,7,8,9,10);
ePerm2:=(1,2);

eG:=Group([ePerm1, ePerm2]);
wS:=MinimalStabChain(eG);

#u:=Random(eG);

