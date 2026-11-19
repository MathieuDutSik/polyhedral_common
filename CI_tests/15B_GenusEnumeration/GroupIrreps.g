# The irreducible characters of a finite group G as factors of the spec
# files of GENUS_GlueIsometry factors (MakeGroupTests.py): for each chi in
# Irr(G), in the order of GAP, a line of OutFile
#   orth s_1 ... s_k       a linear character of values s_j = +-1 on the
#                          generators;
#   herm d e_1 ... e_k     a linear character of values zeta_d^e_j on the
#                          generators, d >= 3, the first of its Galois
#                          orbit;
#   conj i                 a Galois conjugate of the character i;
#   irrep FILE m           a rational character of higher degree, of Schur
#                          index m (1, or 2 for the indicator -1), FILE
#                          (Prefix_i.txt) holding a G-lattice M in the
#                          rational irreducible representation, of
#                          character m chi: its Gram matrix, then the
#                          actions of the generators on the columns;
#   unsupported REASON     otherwise (no permutation character of
#                          multiplicity m).
# A character of degree > 1 and field K = Q(chi) != Q: its Galois orbit
# makes the rational irreducible representation, of character m sum_sigma
# chi^sigma, carried by the first of the orbit ("irrep FILE m"), the others
# "conj i".
# M is the projection by the central idempotent e_chi of the permutation
# lattice Z[G/H], H a subgroup with <chi, 1_H^G> = m, so that the image is
# irreducible of degree m chi(1); the permutation form restricted to it is
# G-invariant, scaled to be integral and primitive.
#
# Usage: gap -q with GroupGens, OutFile and Prefix bound first, e.g.
#   GroupGens:=[(1,2,3,4),(1,2)];; OutFile:="group/group_S4.txt";;
#   Prefix:="group/irrep_S4";; Read("GroupIrreps.g");
G := Group(GroupGens);;
irr := Irr(G);;
ccG := ConjugacyClasses(G);;
subs := List(ConjugacyClassesSubgroups(G), Representative);;
out := OutputTextFile(OutFile, false);;
SetPrintFormattingStatus(out, false);;
WriteMat := function(f, M)
  local row;
  AppendTo(f, Length(M), " ", Length(M[1]), "\n");
  for row in M do
    AppendTo(f, JoinStringsWithSeparator(List(row, String), " "), "\n");
  od;
end;;
first := [];;
for i in [1 .. Length(irr)] do
  chi := irr[i];
  if chi[1] = 1 then
    vals := List(GroupGens, g -> g^chi);
    if ForAll(vals, v -> v = 1 or v = -1) then
      PrintTo(out, "orth ", JoinStringsWithSeparator(List(vals, String), " "),
              "\n");
      continue;
    fi;
    # a Galois conjugate of an earlier linear character
    cond := Conductor(chi);
    j := First([1 .. i - 1], k -> irr[k][1] = 1 and
               ForAny(Filtered([1 .. cond], a -> Gcd(a, cond) = 1),
                      a -> List(GroupGens, g -> GaloisCyc(g^irr[k], a)) =
                           vals));
    if j <> fail then
      PrintTo(out, "conj ", j, "\n");
      continue;
    fi;
    d := Lcm(List(vals, v -> Order(v)));
    es := List(vals, v -> First([0 .. d - 1], e -> E(d)^e = v));
    PrintTo(out, "herm ", d, " ", JoinStringsWithSeparator(List(es, String),
                                                         " "), "\n");
    continue;
  fi;
  # the Galois orbit of chi: the conjugates of an earlier character are
  # marked, the first one carries the rational representation
  cond := Conductor(chi);
  conjs := Set(List(Filtered([1 .. cond], a -> Gcd(a, cond) = 1),
                    a -> GaloisCyc(chi, a)));;
  j := First([1 .. i - 1], k -> irr[k] in conjs);
  if j <> fail then
    PrintTo(out, "conj ", j, "\n");
    continue;
  fi;
  m := 1;
  if Indicator(CharacterTable(G), [chi], 2)[1] = -1 then
    m := 2;
  fi;
  H := First(subs, h -> ScalarProduct(chi, PermutationCharacter(G, h)) = m);
  if H = fail then
    PrintTo(out, "unsupported no permutation character of multiplicity ",
            m, "\n");
    continue;
  fi;
  hom := ActionHomomorphism(G, RightCosets(G, H), OnRight);;
  n := Index(G, H);
  pm := g -> PermutationMat(Image(hom, g), n);
  # the central idempotent of the rational representation (the sum over
  # the Galois orbit), times |G| / chi(1): integral
  Echi := Sum(ccG, c -> Sum(conjs, x -> (Representative(c)^-1)^x) *
                        Sum(AsList(c), pm));;
  B := BaseIntMat(Echi);;
  # the action on the rows of B (x -> x P), and the induced form
  F := B * TransposedMat(B);;
  F := F / Gcd(Flat(F));;
  if ForAny(DiagonalOfMat(F), x -> x mod 2 <> 0) then
    F := 2 * F;
  fi;
  file := Concatenation(Prefix, "_", String(i), ".txt");
  PrintTo(file, Length(GroupGens) + 1, "\n");
  WriteMat(file, F);
  for g in GroupGens do
    P := List(B, row -> SolutionMat(B, row * pm(g)));
    # on the columns: g^T F g = F with g = P^T
    WriteMat(file, TransposedMat(P));
  od;
  PrintTo(out, "irrep ", file, " ", m, "\n");
od;
CloseStream(out);
QUIT;
