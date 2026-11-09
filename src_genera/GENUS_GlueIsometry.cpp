// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Group.h"
#include "Permutation.h"
#include "IsometryGlue.h"
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
  enumerate : for each of them not yet met, the mass of its genus of pairs
              (GENUS_IsometryMass) and the enumeration of that genus
              (GENUS_EnumerateIsometry); the result is every class of pairs
              whose pieces lie in the genera of the (L_i, g_i).
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

template <typename T, typename Tint, typename Tgroup>
void process(std::string const &mode, std::string const &FileGram1,
             std::string const &FileIsom1, std::string const &FileGram2,
             std::string const &FileIsom2, std::string const &det_arg,
             std::string const &OutFormat, QuadLValueMethod lvalue,
             std::ostream &os) {
  MyMatrix<T> A1 = ReadMatrixFile<T>(FileGram1);
  MyMatrix<T> g1 = ReadMatrixFile<T>(FileIsom1);
  MyMatrix<T> A2 = ReadMatrixFile<T>(FileGram2);
  MyMatrix<T> g2 = ReadMatrixFile<T>(FileIsom2);
  std::optional<T> det = ReadDeterminantArg<T>(det_arg);
  IsometryGlueBounds bounds;
  std::vector<GluedPair<T>> ListGlued = GlueIsometryPairs<T, Tint, Tgroup>(
      A1, g1, A2, g2, det, bounds, std::cerr);
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
  if (mode != "enumerate") {
    std::cerr << "GENUS_GlueIsometry: the mode must be glue or enumerate\n";
    throw TerminalException{1};
  }
  using Tresult = GenusEnumerationGroupResult<T, Tgroup>;
  using Tmass = typename Tresult::Tmass;
  std::set<std::pair<int, std::vector<MyMatrix<T>>>> known;
  auto key_of =
      [&](MyMatrix<T> const &A,
          MyMatrix<T> const &g) -> std::pair<int, std::vector<MyMatrix<T>>> {
    CanonicalPairData<T> can = CanonicalPair<T, Tint, Tgroup>(
        A, std::vector<MyMatrix<T>>{g.transpose()}, std::cerr);
    return {can.method, PairForms(can.Acan, can.ListPcan)};
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
      std::cerr << "GENUS_GlueIsometry: a seed has a component without a "
                   "mass, its genus is skipped\n";
      n_without_mass++;
      continue;
    }
    Tmass TotalMass = UniversalScalarConversion<Tmass, T>(*mres.mass);
    std::vector<std::pair<MyMatrix<T>, MyMatrix<T>>> ListSeed{
        {gp.gram, gp.isom}};
    Tresult result = GenusEnumerationGroup<T, Tint, Tgroup>(ListSeed, TotalMass,
                                                            0, 6, std::cerr);
    for (size_t i = 0; i < result.ListGram.size(); i++) {
      known.insert(key_of(result.ListGram[i], result.ListListIsom[i][0]));
    }
    // the seed itself, in case the enumeration stopped short
    known.insert(key_of(gp.gram, gp.isom));
    ListGenus.push_back({gp, *mres.mass, result});
  }
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
    os << "skipped without a mass = " << n_without_mass << "\n";
    os << "complete = "
       << (all_complete && n_without_mass == 0 ? "true" : "false") << "\n";
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
    using T = mpq_class;
    using Tint = mpz_class;
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using TintGroup = mpz_class;
    using Tgroup = permutalib::Group<Telt, TintGroup>;
    QuadLValueMethod lvalue = ExtractQuadLValueOption(argc, argv);
    if (argc < 7 || argc > 10) {
      std::cerr << "GENUS_GlueIsometry [mode] [Gram1] [Isom1] [Gram2] [Isom2] "
                   "[det] [OutFormat] [OutFile]\n";
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
      std::cerr << "\n";
      std::cerr << "lvalue=mpfr / lvalue=arb (anywhere, optional): the\n";
      std::cerr << "        L-values of the quadratic masses with E not\n";
      std::cerr << "        abelian, MPFR (default) or certified with Arb\n";
      std::cerr << "        (needs flint, about twice the time)\n";
      return -1;
    }
    std::string mode = argv[1];
    std::string FileGram1 = argv[2];
    std::string FileIsom1 = argv[3];
    std::string FileGram2 = argv[4];
    std::string FileIsom2 = argv[5];
    std::string det_arg = argv[6];
    std::string OutFormat = (argc >= 8) ? argv[7] : "Summary";
    std::string OutFile = (argc >= 9) ? argv[8] : "stderr";
    auto f = [&](std::ostream &os) -> void {
      process<T, Tint, Tgroup>(mode, FileGram1, FileIsom1, FileGram2, FileIsom2,
                               det_arg, OutFormat, lvalue, os);
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    std::cerr << "Normal termination of GENUS_GlueIsometry\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GENUS_GlueIsometry\n";
    exit(e.eVal);
  }
  runtime(time);
}
