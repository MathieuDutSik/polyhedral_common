// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Group.h"
#include "Permutation.h"
#include "HermitianEnumeration.h"
#include "Namelist.h"
#include <string>
#include <vector>
// clang-format on

/*
  GENUS_HermitianEnumerate: the classes of a genus of positive definite
  Hermitian lattices over a CM field, from its genus symbol
  (HermitianGenusSymbol.h), by HermitianEnumeration.h. The field is given
  by the minimal polynomial of theta and the conjugate of theta; the mass
  of the genus is given, and the enumeration is complete when the classes
  found reach it.
 */

namespace {

FullNamelist NAMELIST_GetStandard_HERMITIAN_ENUMERATE() {
  std::map<std::string, SingleBlock> ListBlock;
  {
    std::map<std::string, std::string> ListStringValues_doc;
    std::map<std::string, std::string> ListListIntValues_doc;
    std::map<std::string, std::string> ListListStringValues_doc;
    ListListIntValues_doc["FieldPolynomial"] = "Default: 0\n\
The coefficients c_0, ..., c_d of the minimal polynomial of theta, monic\n\
(c_d = 1), by increasing degree: E = Q(theta), a CM field";
    ListListStringValues_doc["FieldConjugate"] = "Default: \"0\"\n\
The complex conjugate of theta on 1, theta, ..., theta^{d-1}, rationals";
    ListStringValues_doc["GenusSymbol"] = "Default: unset (required)\n\
The genus, m=<rank>; <prime>: (s,r,+)(s,r,-)...; ... (see\n\
HermitianGenusSymbol.h), the blocks with their norm valuation (s,r,+,n) at\n\
the dyadic primes ramified in E / K_0";
    ListStringValues_doc["TotalMass"] = "Default: unset\n\
The mass of the genus, a rational. Over an abelian field (AbelianConductor)\n\
it is computed exactly, and checked against this value when one is given";
    std::map<std::string, std::string> ListIntValues_doc;
    ListIntValues_doc["AbelianConductor"] = "Default: 0\n\
n > 0 for the abelian CM field fixed by the subgroup H of (Z/n)^* generated\n\
by AbelianGenerators, in place of FieldPolynomial and FieldConjugate";
    ListListIntValues_doc["AbelianGenerators"] = "Default: 1\n\
The generators of H (1 for H trivial, E = Q(zeta_n))";
    SingleBlock BlockDATA;
    BlockDATA.setListStringValues_doc(ListStringValues_doc);
    BlockDATA.setListIntValues_doc(ListIntValues_doc);
    BlockDATA.setListListIntValues_doc(ListListIntValues_doc);
    BlockDATA.setListListStringValues_doc(ListListStringValues_doc);
    ListBlock["DATA"] = BlockDATA;
  }
  {
    std::map<std::string, std::string> ListIntValues_doc;
    ListIntValues_doc["max_prime"] = "Default: 100\n\
The largest rational prime under the neighbor primes";
    SingleBlock BlockMETHOD;
    BlockMETHOD.setListIntValues_doc(ListIntValues_doc);
    ListBlock["METHOD"] = BlockMETHOD;
  }
  {
    std::map<std::string, std::string> ListStringValues_doc;
    ListStringValues_doc["OutFormat"] = "Default: GAP\n\
GAP: a record with the field, the form h of the ambient space (h = sum_k\n\
theta^k h_k), the bases of the classes in it, their unitary group orders\n\
and the mass";
    ListStringValues_doc["OutFile"] = "Default: stdout\n\
The output file, or stdout / stderr";
    SingleBlock BlockSYSTEM;
    BlockSYSTEM.setListStringValues_doc(ListStringValues_doc);
    ListBlock["SYSTEM"] = BlockSYSTEM;
  }
  return FullNamelist(ListBlock);
}

template <typename T, typename Tint>
void WriteHermEnumResult(std::ostream &os, NumberField<T> const &K,
                         HermSeed<T, Tint> const &seed,
                         std::vector<MyMatrix<T>> const &h,
                         HermEnumResult<T, Tint> const &res) {
  os << "return rec(FieldPolynomial:=[";
  for (size_t k = 0; k < K.Phi.size(); k++) {
    os << (k > 0 ? ", " : "") << K.Phi[k];
  }
  os << "], FieldConjugate:=[";
  for (int k = 0; k < K.deg; k++) {
    os << (k > 0 ? ", " : "") << K.conj[k];
  }
  os << "], rank:=" << seed.ring.m << ",\nForm:=[";
  for (size_t k = 0; k < h.size(); k++) {
    os << (k > 0 ? ", " : "");
    WriteMatrixGAP(os, h[k]);
  }
  os << "],\nListBasis:=[";
  for (size_t i = 0; i < res.l_class.size(); i++) {
    os << (i > 0 ? ",\n" : "");
    WriteMatrixGAP(os, res.l_class[i].basis);
  }
  os << "],\nListAutOrder:=[";
  for (size_t i = 0; i < res.l_class.size(); i++) {
    os << (i > 0 ? ", " : "") << res.l_class[i].aut_order;
  }
  os << "],\nmass:=" << res.mass << ", complete:="
     << (res.complete ? "true" : "false") << ");\n";
}

} // namespace

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    using T = mpq_class;
    using Tint = mpz_class;
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt, mpz_class>;
    FullNamelist eFull = NAMELIST_GetStandard_HERMITIAN_ENUMERATE();
    if (argc != 2) {
      std::cerr << "GENUS_HermitianEnumerate [file.nml]\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    NAMELIST_ReadNamelistFile(argv[1], eFull);
    SingleBlock const &BlockDATA = eFull.get_block("DATA");
    SingleBlock const &BlockMETHOD = eFull.get_block("METHOD");
    SingleBlock const &BlockSYSTEM = eFull.get_block("SYSTEM");
    NumberField<T> K;
    std::optional<AbelianField> opt_F;
    int conductor = BlockDATA.get_int("AbelianConductor");
    if (conductor > 0) {
      std::vector<int64_t> gens;
      for (int g : BlockDATA.get_list_int("AbelianGenerators")) {
        gens.push_back(g);
      }
      opt_F = MakeAbelianCMField(conductor, gens);
      K = NumberFieldOfAbelianCM<T>(*opt_F);
    } else {
      K.d = 0;
      for (int c : BlockDATA.get_list_int("FieldPolynomial")) {
        K.Phi.push_back(c);
      }
      K.deg = static_cast<int>(K.Phi.size()) - 1;
      std::vector<std::string> l_conj =
          BlockDATA.get_list_string("FieldConjugate");
      if (K.deg < 1 || static_cast<int>(l_conj.size()) != K.deg) {
        std::cerr << "GENUS_HermitianEnumerate: the field is given by "
                     "AbelianConductor or by FieldPolynomial and "
                     "FieldConjugate (" << K.deg << " coefficients)\n";
        throw TerminalException{1};
      }
      K.conj = ElemE<T>(K.deg, T(0));
      for (int k = 0; k < K.deg; k++) {
        K.conj[k] = ParseScalar<T>(l_conj[k]);
      }
    }
    std::string symbol = BlockDATA.get_string("GenusSymbol");
    std::string str_mass = BlockDATA.get_string("TotalMass");
    std::optional<T> mass_input;
    if (str_mass != "unset") {
      mass_input = ParseScalar<T>(str_mass);
    }
    int max_prime = BlockMETHOD.get_int("max_prime");
    CMOrder<T, Tint> ord = CMOrderOfField<T, Tint>(K, true, std::cerr);
    HermSeed<T, Tint> seed = HermSeedFromSymbol<T, Tint>(ord, symbol, std::cerr);
    std::cerr << "GENUS_HermitianEnumerate: seed of " << symbol << " built\n";
    std::optional<T> mass = mass_input;
    if (opt_F) {
      HermLatticeData<T, Tint> L0 = HermLatticeOfBasis(ord, seed.ring, seed.basis);
      std::vector<MyMatrix<T>> Omega;
      for (auto const &O : L0.Omega) {
        Omega.push_back(UniversalMatrixConversion<T, Tint>(O));
      }
      std::optional<T> mass_exact = HermitianMassAbelianLattice<T>(
          UniversalMatrixConversion<T, Tint>(L0.G), Omega, *opt_F);
      if (mass_exact) {
        std::cerr << "GENUS_HermitianEnumerate: mass " << *mass_exact
                  << " (exact, abelian field)\n";
        if (mass_input && *mass_input != *mass_exact) {
          std::cerr << "GENUS_HermitianEnumerate: the given mass "
                    << *mass_input << " differs from the exact one\n";
          throw TerminalException{1};
        }
        mass = mass_exact;
      }
    }
    if (!mass) {
      std::cerr << "GENUS_HermitianEnumerate: the mass is not available for "
                   "this field and genus, TotalMass must be given\n";
      throw TerminalException{1};
    }
    HermEnumResult<T, Tint> res = HermitianEnumeration<T, Tint, Tgroup>(
        ord, seed, *mass, max_prime, std::cerr);
    std::cerr << "GENUS_HermitianEnumerate: " << res.l_class.size()
              << " classes, mass " << res.mass << " of " << *mass
              << (res.complete ? ", complete" : ", NOT complete") << "\n";
    std::vector<MyMatrix<T>> const &h = seed.h;
    auto f = [&](std::ostream &os_out) -> void {
      WriteHermEnumResult(os_out, K, seed, h, res);
    };
    FILE_PrintStderrStdoutFile(BlockSYSTEM.get_string("OutFile"), f);
    std::cerr << "Normal termination of GENUS_HermitianEnumerate\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GENUS_HermitianEnumerate\n";
    exit(e.eVal);
  }
  runtime(time);
}
