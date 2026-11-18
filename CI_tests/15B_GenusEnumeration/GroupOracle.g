# The oracle of the tests of GENUS_GlueIsometry factors for a finite group
# G, abelian or not: for a genus given by its classes Lambda and generators
# of the O(Lambda), the sum over the classes and over the faithful
# homomorphisms rho : G -> O(Lambda) of 1 / |O(Lambda)|, that is the sum
# over the classes of the pairs (Lambda, rho) of 1 / |C(rho)|, the total
# mass of the genera of pairs (L, G) of each character, by the
# multiplicities (m_1, ..., m_k) of the irreducible characters Irr(G) (in
# the order of GAP) in the character of rho; output lines "m1,...,mk mass",
# after a comment line with the degrees of Irr(G).
#
# The homomorphisms with image conjugate to a subgroup H of P = O(Lambda)
# isomorphic to G are the e o alpha, e : G -> H one isomorphism and alpha in
# Aut(G), each in |P : N_P(H)| conjugates: the faithful homomorphisms of
# character chi weigh sum_H #{alpha : e o alpha of character chi} /
# |N_P(H)|, H over the classes of subgroups isomorphic to G.
#
# Input (InFile): return rec(ListGram := [...], ListListGen := [...]), the
# generators acting on the rows, g A g^T = A, or on the columns (as for
# CharPolyOracle.g); GroupGens, the generators of G as permutations, in
# the order of the generators of the spec files.
#
# Usage: gap -q with InFile, OutFile and GroupGens bound first, e.g.
#   echo 'InFile:="charpoly/genus_E8.g";; OutFile:="out.txt";;
#         GroupGens:=[(1,2,3),(1,2)];; Read("GroupOracle.g");' | gap -q
G := Group(GroupGens);;
irr := Irr(G);;
ccG := ConjugacyClasses(G);;
autG := AsList(AutomorphismGroup(G));;
R := ReadAsFunction(InFile)();;
masses := rec();;
AddMass := function(s, m)
  if not IsBound(masses.(s)) then
    masses.(s) := 0;
  fi;
  masses.(s) := masses.(s) + m;
end;;
for i in [1 .. Length(R.ListGram)] do
  A := R.ListGram[i];
  gens := R.ListListGen[i];
  if ForAll(gens, g -> TransposedMat(g) * A * g = A) then
    gens := List(gens, TransposedMat);
  fi;
  if not ForAll(gens, g -> g * A * TransposedMat(g) = A) then
    Error("the generators are not isometries");
  fi;
  M := Group(gens);
  # a faithful permutation action: the vectors of norm at most the largest
  # diagonal entry, which span
  sv := ShortestVectors(A, Maximum(DiagonalOfMat(A))).vectors;
  sv := Concatenation(sv, -sv);
  Sort(sv);
  hom := ActionHomomorphism(M, sv, OnRight);
  P := Image(hom);
  for e in IsomorphicSubgroups(P, G) do
    H := Image(e);
    w := 1 / Size(Normalizer(P, H));
    for alpha in autG do
      # the character of e o alpha, by the traces on the classes of G
      tr := List(ccG, c -> TraceMat(PreImagesRepresentative(hom,
                Image(e, Image(alpha, Representative(c))))));
      chi := ClassFunction(G, tr);
      mult := List(irr, x -> ScalarProduct(chi, x));
      AddMass(JoinStringsWithSeparator(List(mult, String), ","), w);
    od;
  od;
od;
out := OutputTextFile(OutFile, false);;
SetPrintFormattingStatus(out, false);;
PrintTo(out, "# degrees ", JoinStringsWithSeparator(
              List(irr, x -> String(x[1])), ","), "\n");
for s in Set(RecNames(masses)) do
  PrintTo(out, s, " ", masses.(s), "\n");
od;
CloseStream(out);
QUIT;
