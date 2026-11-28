// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
# include "NumberTheoryFlint.h"
#endif
#include "HermitianEnumeration.h"
#include "HermitianGenera.h"
#include "Namelist.h"
#include <fstream>
#include <set>
#include <string>
#include <vector>
// clang-format on

/*
  GENUS_HermitianGenera: all the genera of Hermitian lattices over a CM field
  E of a given rank, determinant and signatures (HermitianGenera.h), as
  Hecke.jl's hermitian_genera, written as the genus symbols that
  GENUS_HermitianEnumerate takes.

  With CompareFile (one symbol per line, e.g. from Hecke.jl), the genera
  found are matched with those of the file up to the equivalence of the
  local symbols; the exit status is 2 when the two lists differ, after the
  output is written.
 */

namespace {

FullNamelist NAMELIST_GetStandard_HERMITIAN_GENERA() {
  std::map<std::string, SingleBlock> ListBlock;
  {
    std::map<std::string, std::string> ListStringValues_doc;
    std::map<std::string, std::string> ListIntValues_doc;
    std::map<std::string, std::string> ListListIntValues_doc;
    std::map<std::string, std::string> ListListStringValues_doc;
    ListListIntValues_doc["FieldPolynomial"] = "Default: 0\n\
The coefficients c_0, ..., c_d of the minimal polynomial of theta, monic\n\
(c_d = 1), by increasing degree: E = Q(theta), a CM field";
    ListListStringValues_doc["FieldConjugate"] = "Default: \"0\"\n\
The complex conjugate of theta on 1, theta, ..., theta^{d-1}, rationals";
    ListIntValues_doc["AbelianConductor"] = "Default: 0\n\
n > 0 for the abelian CM field fixed by the subgroup H of (Z/n)^* generated\n\
by AbelianGenerators, in place of FieldPolynomial and FieldConjugate";
    ListListIntValues_doc["AbelianGenerators"] = "Default: 1\n\
The generators of H (1 for H trivial, E = Q(zeta_n))";
    ListIntValues_doc["Rank"] = "Default: 1\n\
The rank m over E";
    ListStringValues_doc["Determinant"] = "Default: 1\n\
The determinant (volume) ideal D, invariant under the conjugation: an\n\
integer n for n O_E, or a product of factors (ell,k)^v, (ell,k) the k-th\n\
prime of E above ell (one prime of a conjugate pair is enough). Its\n\
valuation at P is the sum of s r over the Jordan blocks (s, r) at P";
    ListStringValues_doc["MinScale"] = "Default: 1\n\
The lower bound of the scales, as an ideal written as Determinant (1: the\n\
scales are at least 0)";
    ListStringValues_doc["MaxScale"] = "Default: unset\n\
The upper bound of the scales (unset: the determinant)";
    ListListIntValues_doc["Signatures"] = "Default: 0\n\
For each real place of the real subfield, the number of negative entries\n\
(a single 0: all the places definite). Only their parity enters";
    SingleBlock BlockDATA;
    BlockDATA.setListStringValues_doc(ListStringValues_doc);
    BlockDATA.setListIntValues_doc(ListIntValues_doc);
    BlockDATA.setListListIntValues_doc(ListListIntValues_doc);
    BlockDATA.setListListStringValues_doc(ListListStringValues_doc);
    ListBlock["DATA"] = BlockDATA;
  }
  {
    std::map<std::string, std::string> ListBoolValues_doc;
    ListBoolValues_doc["CheckSeed"] = "Default: F\n\
Build a lattice of every genus (HermSeedFromSymbol): a check that the\n\
genus exists";
    std::map<std::string, std::string> ListStringValues_doc;
    ListStringValues_doc["CompareFile"] = "Default: unset\n\
A file of genus symbols, one per line (from another program, Hecke.jl):\n\
each is matched with the genera found, prime by prime up to equivalence,\n\
and the counts of the matched and the unmatched ones are printed";
    SingleBlock BlockMETHOD;
    BlockMETHOD.setListBoolValues_doc(ListBoolValues_doc);
    BlockMETHOD.setListStringValues_doc(ListStringValues_doc);
    ListBlock["METHOD"] = BlockMETHOD;
  }
  {
    std::map<std::string, std::string> ListStringValues_doc;
    ListStringValues_doc["OutFormat"] = "Default: GAP\n\
GAP: the list of the genus symbols. TXT: one symbol per line";
    ListStringValues_doc["OutFile"] = "Default: stdout\n\
The output file, or stdout / stderr";
    SingleBlock BlockSYSTEM;
    BlockSYSTEM.setListStringValues_doc(ListStringValues_doc);
    ListBlock["SYSTEM"] = BlockSYSTEM;
  }
  return FullNamelist(ListBlock);
}

} // namespace

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    using T = mpq_class;
    using Tint = mpz_class;
    FullNamelist eFull = NAMELIST_GetStandard_HERMITIAN_GENERA();
    if (argc != 2) {
      std::cerr << "GENUS_HermitianGenera [file.nml]\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    NAMELIST_ReadNamelistFile(argv[1], eFull);
    SingleBlock const &BlockDATA = eFull.get_block("DATA");
    SingleBlock const &BlockMETHOD = eFull.get_block("METHOD");
    SingleBlock const &BlockSYSTEM = eFull.get_block("SYSTEM");
    std::string error;
    std::optional<HermFieldInput<T>> opt_field =
        HermFieldFromBlock<T>(BlockDATA, error);
    if (!opt_field) {
      std::cerr << "GENUS_HermitianGenera: " << error << "\n";
      throw TerminalException{1};
    }
    CMOrder<T, Tint> ord = CMOrderOfField<T, Tint>(opt_field->K, true, std::cerr);
    int m = BlockDATA.get_int("Rank");
    if (m < 1) {
      std::cerr << "GENUS_HermitianGenera: the rank must be positive\n";
      throw TerminalException{1};
    }
    // the ideals: factors (ell, k, v)
    auto parse = [&](std::string const &name, std::string const &str)
        -> std::vector<HermIdealFactor> {
      std::optional<std::vector<HermIdealFactor>> opt =
          HermParseIdealFactors(str, error);
      if (!opt) {
        std::cerr << "GENUS_HermitianGenera: " << name << ": " << error << "\n";
        throw TerminalException{1};
      }
      return *opt;
    };
    // an integer n: the ideal n O_E, v_P(n O_E) = e(P) v_ell(n)
    std::string str_det = BlockDATA.get_string("Determinant");
    std::vector<HermIdealFactor> l_det;
    if (str_det.find('(') == std::string::npos) {
      Tint n_det = ParseScalar<Tint>(str_det);
      if (n_det < 1) {
        std::cerr << "GENUS_HermitianGenera: the determinant n must be "
                     "positive\n";
        throw TerminalException{1};
      }
      for (auto const &f : FactorsInt(n_det)) {
        int64_t ell = UniversalScalarConversion<int64_t, Tint>(f);
        bool seen = false;
        for (auto const &g : l_det) {
          seen = seen || g.ell == ell;
        }
        if (seen) {
          continue;
        }
        int a = 0;
        Tint rest = n_det;
        while (ResInt(rest, f) == 0) {
          rest /= f;
          a++;
        }
        std::vector<HermIdealPrime<Tint>> l_prime = HermIdealPrimesAbove(ord, ell);
        for (size_t k = 0; k < l_prime.size(); k++) {
          l_det.push_back({ell, k, a * l_prime[k].e});
        }
      }
    } else {
      l_det = parse("Determinant", str_det);
    }
    std::vector<HermIdealFactor> l_min =
        parse("MinScale", BlockDATA.get_string("MinScale"));
    std::string str_max = BlockDATA.get_string("MaxScale");
    bool max_given = str_max != "unset";
    std::vector<HermIdealFactor> l_max =
        max_given ? parse("MaxScale", str_max) : std::vector<HermIdealFactor>{};
    std::vector<int> sig = BlockDATA.get_list_int("Signatures");
    int n_real = ord.deg / 2;
    if (sig.size() == 1 && sig[0] == 0) {
      sig.assign(n_real, 0);
    }
    if (static_cast<int>(sig.size()) != n_real) {
      std::cerr << "GENUS_HermitianGenera: Signatures needs " << n_real
                << " entries, one per real place of the real subfield\n";
      throw TerminalException{1};
    }
    int n_odd = 0;
    for (int s : sig) {
      if (s < 0 || s > m) {
        std::cerr << "GENUS_HermitianGenera: a signature is not in [0, " << m
                  << "]\n";
        throw TerminalException{1};
      }
      n_odd += s % 2;
    }
    // the relevant rational primes: 2, those under the ramified primes (the
    // divisors of the discriminant), those of the ideals
    std::set<int64_t> l_ell{2};
    {
      Tint disc = T_abs(DeterminantMat(ord.trace_norm_form()));
      for (auto const &f : FactorsInt(disc)) {
        l_ell.insert(UniversalScalarConversion<int64_t, Tint>(f));
      }
    }
    for (auto const *l : {&l_det, &l_min, &l_max}) {
      for (auto const &f : *l) {
        l_ell.insert(f.ell);
      }
    }
    // the bounds by rational prime and representative index of the prime of
    // K_0
    std::vector<std::pair<int64_t, std::map<size_t, HermGeneraBounds>>> l_bound;
    for (int64_t ell : l_ell) {
      std::vector<HermIdealPrime<Tint>> l_prime = HermIdealPrimesAbove(ord, ell);
      auto val_of = [&](std::vector<HermIdealFactor> const &l, size_t i_P,
                        int def) -> int {
        std::optional<int> v;
        for (auto const &f : l) {
          if (f.ell != ell) {
            continue;
          }
          if (f.k >= l_prime.size()) {
            std::cerr << "GENUS_HermitianGenera: there is no prime (" << ell
                      << "," << f.k << ")\n";
            throw TerminalException{1};
          }
          if (HermK0Representative(l_prime, f.k) == i_P) {
            if (v && *v != f.v) {
              std::cerr << "GENUS_HermitianGenera: an ideal is not invariant "
                           "under the conjugation at " << l_prime[f.k].name
                        << "\n";
              throw TerminalException{1};
            }
            v = f.v;
          }
        }
        return v ? *v : def;
      };
      std::map<size_t, HermGeneraBounds> bmap;
      for (auto const &p : HermK0PrimesAbove(ord, l_prime)) {
        int dv = val_of(l_det, p.i_P, 0);
        if (p.type == HermIdealType::Ramified && dv % 2 != 0) {
          // the determinant is the ideal of an element of K_0 up to norms:
          // an even valuation at a ramified prime (as Hecke.jl requires)
          std::cerr << "GENUS_HermitianGenera: the determinant has an odd "
                       "valuation at the ramified prime "
                    << l_prime[p.i_P].name << ": there is no genus\n";
          throw TerminalException{1};
        }
        int smin = val_of(l_min, p.i_P, 0);
        int smax = val_of(l_max, p.i_P, max_given ? 0 : dv);
        if (!max_given) {
          smax = dv;
        }
        bmap[p.i_P] = HermGeneraBounds{dv, smin, smax};
      }
      l_bound.push_back({ell, bmap});
    }
    std::vector<std::string> l_genus =
        HermitianGenera<T, Tint>(ord, m, l_bound, n_odd, std::cerr);
    std::cerr << "GENUS_HermitianGenera: " << l_genus.size() << " genera\n";
    if (BlockMETHOD.get_bool("CheckSeed")) {
      for (auto const &sym : l_genus) {
        HermSeed<T, Tint> seed = HermSeedFromSymbol<T, Tint>(ord, sym, std::cerr);
        (void)seed;
      }
      std::cerr << "GENUS_HermitianGenera: a lattice of every genus built\n";
    }
    bool compare_ok = true;
    std::string CompareFile = BlockMETHOD.get_string("CompareFile");
    if (CompareFile != "unset") {
      std::ifstream is(CompareFile);
      if (!is.good()) {
        std::cerr << "GENUS_HermitianGenera: cannot read " << CompareFile
                  << "\n";
        throw TerminalException{1};
      }
      HermPrimeCache<T, Tint> cache;
      cache.ord = &ord;
      std::vector<std::map<HermK0Ref, HermLocalSymbol>> l_ours, l_other;
      std::vector<std::string> l_other_str;
      for (auto const &sym : l_genus) {
        int mm;
        l_ours.push_back(HermGenusLocalMap(cache, sym, mm));
      }
      std::string line;
      while (std::getline(is, line)) {
        if (line.find("m=") == std::string::npos) {
          continue;
        }
        int mm;
        l_other.push_back(HermGenusLocalMap(cache, line, mm));
        l_other_str.push_back(line);
        if (mm != m) {
          std::cerr << "GENUS_HermitianGenera: " << line << " is not of rank "
                    << m << "\n";
          throw TerminalException{1};
        }
      }
      std::vector<int> ours_hit(l_ours.size(), 0);
      size_t n_matched = 0, n_only_other = 0;
      for (size_t j = 0; j < l_other.size(); j++) {
        int n_hit = 0;
        for (size_t i = 0; i < l_ours.size(); i++) {
          if (HermGenusEquivalent(cache, m, l_ours[i], l_other[j])) {
            ours_hit[i]++;
            n_hit++;
          }
        }
        if (n_hit == 1) {
          n_matched++;
        } else if (n_hit == 0) {
          n_only_other++;
          std::cerr << "COMPARE: only in " << CompareFile << ": "
                    << l_other_str[j] << "\n";
        } else {
          std::cerr << "COMPARE: " << l_other_str[j] << " matches " << n_hit
                    << " genera found\n";
        }
      }
      size_t n_only_ours = 0;
      for (size_t i = 0; i < l_ours.size(); i++) {
        if (ours_hit[i] == 0) {
          n_only_ours++;
          std::cerr << "COMPARE: only found here: " << l_genus[i] << "\n";
        }
      }
      compare_ok = n_matched == l_other.size() && n_only_ours == 0 &&
                   l_ours.size() == l_other.size();
      std::cerr << "COMPARE: found " << l_ours.size() << ", given "
                << l_other.size() << ", matched " << n_matched
                << ", only found here " << n_only_ours << ", only given "
                << n_only_other << "\n";
    }
    std::string OutFormat = BlockSYSTEM.get_string("OutFormat");
    if (OutFormat != "GAP" && OutFormat != "TXT") {
      std::cerr << "GENUS_HermitianGenera: OutFormat is GAP or TXT\n";
      throw TerminalException{1};
    }
    auto f = [&](std::ostream &os) -> void {
      if (OutFormat == "TXT") {
        for (auto const &sym : l_genus) {
          os << sym << "\n";
        }
        return;
      }
      os << "return [";
      for (size_t i = 0; i < l_genus.size(); i++) {
        os << (i > 0 ? ",\n" : "\n") << "\"" << l_genus[i] << "\"";
      }
      os << "];\n";
    };
    FILE_PrintStderrStdoutFile(BlockSYSTEM.get_string("OutFile"), f);
    if (!compare_ok) {
      std::cerr << "GENUS_HermitianGenera: the genera differ from those of "
                   "CompareFile\n";
      std::cerr << "Normal termination of GENUS_HermitianGenera\n";
      runtime(time);
      return 2;
    }
    std::cerr << "Normal termination of GENUS_HermitianGenera\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GENUS_HermitianGenera\n";
    exit(e.eVal);
  }
  runtime(time);
}
