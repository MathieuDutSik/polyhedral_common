// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Arithmetic_dispatch_types.h"
#include "Group.h"
#include "Permutation.h"
#include "IsometryGlue.h"
#include "IsometryGlueFactors.h"
#include "GenusSymbol.h"
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>
// clang-format on

/*
  The seeds of the enumeration of pairs (L, g) by gluing two pairs, see
  IsometryGlue.h.

  glue      : the distinct sigma-stable primitive extensions of (L_1, g_1)
              (+) (L_2, g_2), of the given determinant or of any.
  involutions: every genus of pairs (L, sigma) with L in the genus of the
              given lattice and sigma an involution of +1-eigenspace of
              rank r: the genera of the pieces L_+ and L_- (even, of ranks r
              and n - r, determinants d_+ d_- = d 4^a for a glue (Z/2)^a,
              their odd local symbols adding up to those of L, as the glue
              is 2-primary), glued with +1 and -1 and kept in the genus of
              L, then as enumerate.
  pieces    : the same for sigma = 1 (+) sigma_2, the +1-piece of rank r and
              (L_2, sigma_2) among the pairs of a file, one per genus of
              pairs (e.g. the hermitian genera over Z[zeta_d] of
              HermitianPieces.jl): the factors "file" and "orth 1 r".
  factors   : the same for sigma acting on several factors of pairwise
              coprime minimal polynomials, given by a spec file, the pieces
              glued one factor at a time (IsometryGlueFactors.h).
  enumerate : for each of them not yet met, the mass of its genus of pairs
              (GENUS_IsometryMass) and the enumeration of that genus
              (GENUS_EnumerateIsometry); the result is every class of pairs
              whose pieces lie in the genera of the (L_i, g_i).
              Exit status as GENUS_EnumerateIsometry: 2 when a genus did
              not close or a glued pair has no proved mass, 3 when the
              classes of a genus weigh more than its mass.
 */

template <typename T>
std::optional<T> ReadDeterminantArg(std::string const &s) {
  if (s == "any") {
    return {};
  }
  return T(s);
}

template <typename T>
void WriteGluedPairGAP(std::ostream &os, GluedPair<T> const &gp) {
  os << "rec(det:=" << gp.det << ", glue_order:=" << gp.glue_order
     << ", gram:=";
  WriteMatrixGAP(os, gp.gram);
  os << ", isom:=";
  WriteMatrixGAP(os, gp.isom);
  os << ")";
}

/*
  The odd local symbols of g1 + g2 (the components of equal scale merged)
  against those of s. A symbol lists only the primes of its determinant:
  at another odd p the lattice is unimodular, a component (0, rank, (det /
  p)).
 */
inline std::vector<std::array<int, 3>> OddLocalComponents(GenusSymbol const &g,
                                                          int64_t p) {
  auto iter = g.sym_odd.find(p);
  if (iter != g.sym_odd.end()) {
    return iter->second;
  }
  mpz_class d = GenusSymbolDeterminant<mpz_class>(g);
  mpz_class r = d % mpz_class(static_cast<long>(p));
  int sign = kronecker_symbol<int64_t>(r.get_si(), p);
  return {{0, g.rank, sign}};
}

inline bool OddSymbolsAddUp(GenusSymbol const &g1, GenusSymbol const &g2,
                            GenusSymbol const &s) {
  std::set<int64_t> primes;
  for (auto const *g : {&g1, &g2, &s}) {
    for (auto const &kv : g->sym_odd) {
      primes.insert(kv.first);
    }
  }
  for (int64_t p : primes) {
    std::map<int, std::pair<int, int>> merged;
    for (auto const *g : {&g1, &g2}) {
      for (auto const &c : OddLocalComponents(*g, p)) {
        auto &m = merged[c[0]];
        if (m.first == 0) {
          m.second = 1;
        }
        m.first += c[1];
        m.second *= c[2];
      }
    }
    std::vector<std::array<int, 3>> sum;
    for (auto const &[k, rs] : merged) {
      if (rs.first > 0) {
        sum.push_back({k, rs.first, rs.second});
      }
    }
    std::vector<std::array<int, 3>> target;
    for (auto const &c : OddLocalComponents(s, p)) {
      if (c[1] > 0) {
        target.push_back(c);
      }
    }
    std::sort(target.begin(), target.end());
    if (sum != target) {
      return false;
    }
  }
  return true;
}

/*
  The seeds of every genus of pairs (L, sigma), L in the genus of A and
  sigma an involution with a +1-eigenspace of rank r.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<GluedPair<T>> InvolutionSeeds(MyMatrix<T> const &A, int r,
                                          std::ostream &os) {
  int n = A.rows();
  int rm = n - r;
  if (r < 0 || rm < 0) {
    std::cerr << "GENUS_GlueIsometry: the rank " << r << " is not in [0, "
              << n << "]\n";
    throw TerminalException{1};
  }
  T det = DeterminantMat(A);
  std::vector<GluedPair<T>> res;
  if (r == 0 || rm == 0) {
    MyMatrix<T> g = IdentityMat<T>(n);
    if (r == 0) {
      g = -g;
    }
    res.push_back({A, g, det, 1});
    return res;
  }
  GenusSymbol S = GenusSymbolOfGram<T, Tint>(A);
  Tint d = UniversalScalarConversion<Tint, T>(det);
  std::map<std::string, MyMatrix<T>> reps;
  auto rep_of = [&](GenusSymbol const &g) -> MyMatrix<T> {
    std::string key = GenusSymbolString(g);
    auto iter = reps.find(key);
    if (iter != reps.end()) {
      return iter->second;
    }
    MyMatrix<T> G = GenusSymbolRepresentative<T, Tint>(g, os);
    reps[key] = G;
    return G;
  };
  IsometryGlueBounds bounds;
  size_t n_cand = 0;
  for (int a = 0; a <= std::min(r, rm); a++) {
    Tint N = d * MyPow<Tint>(Tint(2), 2 * a);
    Tint two_a = MyPow<Tint>(Tint(2), a);
    std::vector<Tint> divisors{Tint(1)};
    for (auto const &[q, mult] : FactorsIntMap(N)) {
      std::vector<Tint> next;
      for (auto const &x : divisors) {
        Tint y = x;
        for (size_t i = 0; i <= mult; i++) {
          next.push_back(y);
          y *= q;
        }
      }
      divisors = next;
    }
    for (Tint const &dp : divisors) {
      Tint dm = N / dp;
      if (ResInt(dp, two_a) != 0 || ResInt(dm, two_a) != 0) {
        continue;
      }
      std::vector<GenusSymbol> Lp = EnumerateEvenGenusSymbols<Tint>(r, dp);
      std::vector<GenusSymbol> Lm = EnumerateEvenGenusSymbols<Tint>(rm, dm);
      for (auto const &gp : Lp) {
        for (auto const &gm : Lm) {
          if (!OddSymbolsAddUp(gp, gm, S)) {
            continue;
          }
          n_cand++;
          std::cerr << "GENUS_GlueIsometry: pieces " << GenusSymbolString(gp)
                    << " and " << GenusSymbolString(gm) << "\n";
          MyMatrix<T> Ap = rep_of(gp), Am = rep_of(gm);
          std::vector<GluedPair<T>> glued = GlueIsometryPairs<T, Tint, Tgroup>(
              Ap, IdentityMat<T>(r), Am, MyMatrix<T>(-IdentityMat<T>(rm)), det,
              bounds, os);
          for (auto const &gl : glued) {
            if (GenusSymbolOfGram<T, Tint>(gl.gram) == S) {
              res.push_back(gl);
            }
          }
        }
      }
    }
  }
  std::cerr << "GENUS_GlueIsometry: " << n_cand << " pairs of genera of the "
            << "pieces, " << res.size() << " glued pairs in the genus of L\n";
  return res;
}

template <typename T, typename Tint, typename Tgroup>
void process(std::string const &mode, std::string const &FileGram1,
             std::string const &FileIsom1, std::string const &FileGram2,
             std::string const &FileIsom2, std::string const &det_arg,
             std::string const &OutFormat, QuadLValueMethod lvalue,
             std::ostream &os) {
  std::vector<GluedPair<T>> ListGlued;
  // an odd lattice (factors): the pairs of its even sublattice, then the
  // odd overlattices (IsometryGlueFactors.h)
  std::optional<MyMatrix<T>> odd_gram;
  if (mode == "factors") {
    // FileGram1 the lattice, FileIsom1 the spec file of the factors (the
    // pieces of the even sublattice when the lattice is odd)
    MyMatrix<T> A = ReadMatrixFile<T>(FileGram1);
    bool is_even = true;
    for (int i = 0; i < A.rows(); i++) {
      is_even = is_even && IsInteger(A(i, i) / T(2));
    }
    if (!is_even) {
      odd_gram = A;
      MyMatrix<T> B0 = EvenSublatticeBasis(A);
      A = B0 * A * B0.transpose();
    }
    ListGlued = FactorSeeds<T, Tint, Tgroup>(
        A, ReadIsometryFactors<T>(FileIsom1), std::cerr);
  } else if (mode == "pieces") {
    // FileGram1 the lattice, FileIsom1 the rank r of the +1-piece, FileGram2
    // the list of the other pieces: the factors "file FileGram2" and
    // "orth 1 r"
    std::vector<IsometryFactor<T>> factors{
        ParseIsometryFactor<T>("file " + FileGram2)};
    if (std::stoi(FileIsom1) > 0) {
      factors.push_back(ParseIsometryFactor<T>("orth 1 " + FileIsom1));
    }
    ListGlued = FactorSeeds<T, Tint, Tgroup>(ReadMatrixFile<T>(FileGram1),
                                             factors, std::cerr);
  } else if (mode == "involutions") {
    // FileGram1 the lattice, FileIsom1 the rank of the +1-eigenspace
    ListGlued = InvolutionSeeds<T, Tint, Tgroup>(
        ReadMatrixFile<T>(FileGram1), std::stoi(FileIsom1), std::cerr);
  } else {
    MyMatrix<T> A1 = ReadMatrixFile<T>(FileGram1);
    MyMatrix<T> g1 = ReadMatrixFile<T>(FileIsom1);
    MyMatrix<T> A2 = ReadMatrixFile<T>(FileGram2);
    MyMatrix<T> g2 = ReadMatrixFile<T>(FileIsom2);
    std::optional<T> det = ReadDeterminantArg<T>(det_arg);
    IsometryGlueBounds bounds;
    ListGlued = GlueIsometryPairs<T, Tint, Tgroup>(A1, g1, A2, g2, det, bounds,
                                                   std::cerr);
  }
  if (mode == "glue") {
    if (OutFormat == "Summary") {
      os << "glued pairs = " << ListGlued.size() << "\n";
      for (auto const &gp : ListGlued) {
        os << "det = " << gp.det << " glue order = " << gp.glue_order << "\n";
      }
      return;
    }
    if (OutFormat == "GAP") {
      os << "return [";
      for (size_t i = 0; i < ListGlued.size(); i++) {
        os << (i > 0 ? ",\n" : "\n");
        WriteGluedPairGAP(os, ListGlued[i]);
      }
      os << "];\n";
      return;
    }
    if (OutFormat == "ListMatrix") {
      // the Gram matrices, then the isometries, for GENUS_EnumerateIsometry
      std::vector<MyMatrix<T>> ListGram, ListIsom;
      for (auto const &gp : ListGlued) {
        ListGram.push_back(gp.gram);
        ListIsom.push_back(gp.isom);
      }
      WriteListMatrix(os, ListGram);
      WriteListMatrix(os, ListIsom);
      return;
    }
    std::cerr << "GENUS_GlueIsometry: unknown OutFormat=" << OutFormat
              << " (Summary, GAP or ListMatrix)\n";
    throw TerminalException{1};
  }
  if (mode != "enumerate" && mode != "involutions" && mode != "pieces" &&
      mode != "factors") {
    std::cerr << "GENUS_GlueIsometry: the mode must be glue, enumerate, "
                 "involutions, pieces or factors\n";
    throw TerminalException{1};
  }
  using Tresult = GenusEnumerationGroupResult<T, Tgroup>;
  using Tmass = typename Tresult::Tmass;
  PairMethods methods = DefaultPairMethods();
  std::set<std::vector<MyMatrix<T>>> known;
  auto key_of = [&](MyMatrix<T> const &A,
                    MyMatrix<T> const &g) -> std::vector<MyMatrix<T>> {
    return PairCanonicalForms<T, Tint, Tgroup>(
        A, std::vector<MyMatrix<T>>{g.transpose()}, methods.canonical,
        std::cerr);
  };
  struct GenusRecord {
    GluedPair<T> seed;
    T mass;
    Tresult result;
  };
  std::vector<GenusRecord> ListGenus;
  size_t n_without_mass = 0;
  for (auto const &gp : ListGlued) {
    if (known.count(key_of(gp.gram, gp.isom)) == 1) {
      continue;
    }
    IsometryMassBounds mass_bounds;
    mass_bounds.lvalue = lvalue;
    std::map<int, T> herm_mass;
    IsometryMassResult<T> mres = IsometryGenusMass<T, Tint>(
        gp.gram, gp.isom, herm_mass, mass_bounds, std::cerr);
    if (!mres.mass) {
      // no component mass, or local factors missing
      // (mass_without_local_factors, only a lower bound)
      std::cerr << "GENUS_GlueIsometry: a seed has no proved mass, its genus "
                   "is skipped\n";
      n_without_mass++;
      continue;
    }
    Tmass TotalMass = MassToGroupField<Tmass>(*mres.mass);
    std::vector<std::pair<MyMatrix<T>, MyMatrix<T>>> ListSeed{
        {gp.gram, gp.isom}};
    Tresult result = GenusEnumerationGroup<T, Tint, Tgroup>(
        ListSeed, TotalMass, 0, 6, methods, std::cerr);
    for (size_t i = 0; i < result.ListGram.size(); i++) {
      known.insert(key_of(result.ListGram[i], result.ListListIsom[i][0]));
    }
    // the seed itself, in case the enumeration stopped short
    known.insert(key_of(gp.gram, gp.isom));
    ListGenus.push_back({gp, *mres.mass, result});
  }
  if (odd_gram) {
    // the classes of the odd pairs (L, sigma), L in the genus of odd_gram:
    // the orbits of the odd overlattices of the classes of the even pairs
    GenusSymbol S_odd = GenusSymbolOfGram<T, Tint>(*odd_gram);
    Tmass total_odd(0);
    size_t n_odd = 0;
    bool all_complete = n_without_mass == 0;
    for (auto const &rec : ListGenus) {
      all_complete = all_complete && rec.result.complete;
      for (size_t i = 0; i < rec.result.ListGram.size(); i++) {
        MyMatrix<T> const &A0 = rec.result.ListGram[i];
        MyMatrix<T> P = rec.result.ListListIsom[i][0].transpose();
        std::vector<MyMatrix<Tint>> ListGen =
            PairCentralizer<T, Tint, Tgroup>(A0, std::vector<MyMatrix<T>>{P},
                                             methods.automorphism, std::cerr)
                .ListGen;
        for (auto const &ov :
             OddOverlatticeOrbits<T, Tint>(A0, P, ListGen)) {
          MyMatrix<T> AL = ov.basis * A0 * ov.basis.transpose();
          if (GenusSymbolOfGram<T, Tint>(AL) != S_odd) {
            continue;
          }
          n_odd++;
          total_odd += Tmass(static_cast<long>(ov.orbit_size)) /
                       Tmass(rec.result.ListCentralizerOrder[i]);
        }
      }
    }
    if (OutFormat != "Summary") {
      std::cerr << "GENUS_GlueIsometry: an odd lattice has only the Summary "
                   "output\n";
      throw TerminalException{1};
    }
    os << "genera of pairs of the even sublattice = " << ListGenus.size()
       << "\n";
    os << "classes of odd pairs = " << n_odd << "\n";
    os << "total mass = " << total_odd << "\n";
    os << "complete = " << (all_complete ? "true" : "false") << "\n";
    for (auto const &rec : ListGenus) {
      if (rec.result.mass_exceeded()) {
        throw GenusMassExceededException{};
      }
    }
    for (auto const &rec : ListGenus) {
      CheckGenusEnumerationGroupResult(rec.result);
    }
    if (n_without_mass > 0) {
      throw GenusMassNotClosedException{};
    }
    return;
  }
  // after the output: a genus whose classes weigh more than its mass first,
  // then the genera that did not close or have no proved mass
  auto check_genera = [&]() -> void {
    for (auto const &rec : ListGenus) {
      if (rec.result.mass_exceeded()) {
        throw GenusMassExceededException{};
      }
    }
    for (auto const &rec : ListGenus) {
      CheckGenusEnumerationGroupResult(rec.result);
    }
    if (n_without_mass > 0) {
      throw GenusMassNotClosedException{};
    }
  };
  if (OutFormat == "Summary") {
    os << "glued pairs = " << ListGlued.size() << "\n";
    os << "genera of pairs = " << ListGenus.size() << "\n";
    size_t n_class = 0;
    bool all_complete = true;
    for (auto const &rec : ListGenus) {
      os << "det = " << rec.seed.det << " glue order = " << rec.seed.glue_order
         << " mass = " << rec.mass
         << " classes = " << rec.result.ListGram.size()
         << " complete = " << (rec.result.complete ? "true" : "false") << "\n";
      n_class += rec.result.ListGram.size();
      all_complete = all_complete && rec.result.complete;
    }
    os << "classes of pairs = " << n_class << "\n";
    T total(0);
    for (auto const &rec : ListGenus) {
      total += rec.mass;
    }
    os << "total mass = " << total << "\n";
    os << "skipped without a mass = " << n_without_mass << "\n";
    os << "complete = "
       << (all_complete && n_without_mass == 0 ? "true" : "false") << "\n";
    check_genera();
    return;
  }
  if (OutFormat == "GAP") {
    os << "return [";
    for (size_t i = 0; i < ListGenus.size(); i++) {
      auto const &rec = ListGenus[i];
      os << (i > 0 ? ",\n" : "\n") << "rec(seed:=";
      WriteGluedPairGAP(os, rec.seed);
      os << ",\n mass:=" << rec.mass
         << ",\n complete:=" << (rec.result.complete ? "true" : "false")
         << ",\n ListGram:=[";
      for (size_t j = 0; j < rec.result.ListGram.size(); j++) {
        os << (j > 0 ? ",\n" : "");
        WriteMatrixGAP(os, rec.result.ListGram[j]);
      }
      os << "],\n ListIsom:=[";
      for (size_t j = 0; j < rec.result.ListListIsom.size(); j++) {
        os << (j > 0 ? ",\n" : "");
        WriteMatrixGAP(os, rec.result.ListListIsom[j][0]);
      }
      os << "],\n ListCentralizerOrder:=[";
      for (size_t j = 0; j < rec.result.ListCentralizerOrder.size(); j++) {
        os << (j > 0 ? ", " : "") << rec.result.ListCentralizerOrder[j];
      }
      os << "])";
    }
    os << "];\n";
    check_genera();
    return;
  }
  std::cerr << "GENUS_GlueIsometry: unknown OutFormat=" << OutFormat
            << " (Summary or GAP)\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    // the group orders are in TintGroup, whatever the arithmetic of the
    // coordinates T / Tint chosen by arithmetic=
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using TintGroup = mpz_class;
    using Tgroup = permutalib::Group<Telt, TintGroup>;
    std::string arith = ExtractArithmeticOption(argc, argv, "gmp");
    ArithmeticDispatch dispatch;
    QuadLValueMethod lvalue = ExtractQuadLValueOption(argc, argv);
    // involutions and factors: [Gram] [r or FileSpec] [OutFormat] [OutFile]
    bool short_args = (argc >= 2 && (std::string(argv[1]) == "involutions" ||
                                     std::string(argv[1]) == "factors"));
    bool pieces_mode = (argc >= 2 && std::string(argv[1]) == "pieces");
    if (short_args ? (argc < 4 || argc > 8)
                   : pieces_mode ? (argc < 5 || argc > 9)
                                 : (argc < 7 || argc > 10)) {
      std::cerr << "GENUS_GlueIsometry [mode] [Gram1] [Isom1] [Gram2] [Isom2] "
                   "[det] [OutFormat] [OutFile]\n";
      std::cerr << "GENUS_GlueIsometry involutions [Gram] [r] [OutFormat] "
                   "[OutFile]\n";
      std::cerr << "GENUS_GlueIsometry pieces [Gram] [r] [FilePieces] "
                   "[OutFormat] [OutFile]\n";
      std::cerr << "GENUS_GlueIsometry factors [Gram] [FileSpec] [OutFormat] "
                   "[OutFile]\n";
      std::cerr << "\n";
      std::cerr << "mode             : glue or enumerate\n";
      std::cerr
          << "Gram1, Isom1     : the pair (L_1, g_1), g_1^T A_1 g_1 = A_1\n";
      std::cerr << "Gram2, Isom2     : the pair (L_2, g_2)\n";
      std::cerr
          << "det              : the determinant of the glued lattices,\n";
      std::cerr << "                   or any\n";
      std::cerr << "OutFormat        : Summary (default), GAP, or for glue\n";
      std::cerr
          << "                   ListMatrix (the Gram matrices then the\n";
      std::cerr << "                   isometries)\n";
      std::cerr << "OutFile          : stderr (default), stdout or a file\n";
      std::cerr << "\n";
      std::cerr << "glue      : the classes of the g-stable primitive\n";
      std::cerr << "            extensions of L_1 (+) L_2, g = g_1 (+) g_2\n";
      std::cerr << "enumerate : the classes of all the genera of pairs met,\n";
      std::cerr << "            each by its mass\n";
      std::cerr << "involutions: every genus of pairs (L, sigma), L in the\n";
      std::cerr << "            genus of Gram, sigma an involution with a\n";
      std::cerr << "            +1-eigenspace of rank r, by the genera of\n";
      std::cerr << "            the pieces, glued, then as enumerate\n";
      std::cerr << "pieces    : the same with sigma = 1 (+) sigma_2, the\n";
      std::cerr << "            +1-piece of rank r and the other one among\n";
      std::cerr << "            the pairs of FilePieces (a list of Gram\n";
      std::cerr << "            matrices and isometries, one per genus)\n";
      std::cerr << "factors   : the same for sigma acting on several\n";
      std::cerr << "            factors, one per line of FileSpec:\n";
      std::cerr << "            orth 1 r, orth -1 r (the even genera with\n";
      std::cerr << "            +1 or -1) or file FilePieces, the pieces\n";
      std::cerr << "            glued one factor at a time\n";
      std::cerr << "\n";
      std::cerr << "lvalue=double / lvalue=mpfr / lvalue=arb (anywhere,\n";
      std::cerr << "        optional): the L-values of the quadratic masses\n";
      std::cerr << "        with E not abelian, in double precision (when\n";
      std::cerr << "        it suffices), MPFR (default) or certified with\n";
      std::cerr << "        Arb (needs flint, about twice the time)\n";
      std::cerr << "arithmetic=NAME (anywhere, optional, default gmp): the\n";
      std::cerr << "        arithmetic of the computation, among\n";
      dispatch.print_descriptions(std::cerr);
      return -1;
    }
    std::string mode = argv[1];
    std::string FileGram1 = argv[2];
    std::string FileIsom1 = argv[3];
    std::string FileGram2, FileIsom2, det_arg = "any", OutFormat, OutFile;
    if (short_args) {
      OutFormat = (argc >= 5) ? argv[4] : "Summary";
      OutFile = (argc >= 6) ? argv[5] : "stderr";
    } else if (pieces_mode) {
      FileGram2 = argv[4];
      OutFormat = (argc >= 6) ? argv[5] : "Summary";
      OutFile = (argc >= 7) ? argv[6] : "stderr";
    } else {
      FileGram2 = argv[4];
      FileIsom2 = argv[5];
      det_arg = argv[6];
      OutFormat = (argc >= 8) ? argv[7] : "Summary";
      OutFile = (argc >= 9) ? argv[8] : "stderr";
    }
    add_rational_integer_arithmetics_by_library(
        dispatch, [&]<typename T, typename Tint>() -> void {
          auto f = [&](std::ostream &os) -> void {
            process<T, Tint, Tgroup>(mode, FileGram1, FileIsom1, FileGram2, FileIsom2,
                                     det_arg, OutFormat, lvalue, os);
          };
          FILE_PrintStderrStdoutFile(OutFile, f);
        });
    dispatch.run(arith);
    std::cerr << "Normal termination of GENUS_GlueIsometry\n";
  } catch (GenusMassExceededException const &) {
    std::cerr << "GENUS_GlueIsometry: ERROR: the classes of a genus of pairs "
                 "weigh more than its mass\n";
    exit(GENUS_EXIT_MASS_EXCEEDED);
  } catch (GenusMassNotClosedException const &) {
    std::cerr << "GENUS_GlueIsometry: INCOMPLETE: a genus of pairs did not "
                 "close its mass or has no proved mass\n";
    exit(GENUS_EXIT_MASS_NOT_CLOSED);
  } catch (TerminalException const &e) {
    std::cerr << "Error in GENUS_GlueIsometry\n";
    exit(e.eVal);
  }
  runtime(time);
}
