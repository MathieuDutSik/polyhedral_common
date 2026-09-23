// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Group.h"
#include "Permutation.h"
#include "LatticeAutomorphismVinberg.h"
#include "Positivity.h"
// clang-format on

/*
  Validates the Vinberg automorphism computation against the direct graph
  method: the order must match, and every returned generator must
  preserve the form. Times both to show the speedup.
 */
template <typename T,
          typename Tint = typename underlying_z_ring<T>::ring_type>
void process(std::string const &FileListGram, std::ostream &os) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  std::vector<MyMatrix<T>> ListGram = ReadListMatrixFile<T>(FileListGram);
  int n_error = 0;
  double tvin = 0, tdir = 0;
  for (size_t i = 0; i < ListGram.size(); i++) {
    MyMatrix<T> const &eMat = ListGram[i];
    MyMatrix<Tint> G_int = UniversalMatrixConversion<Tint, T>(eMat);
    // Vinberg
    MicrosecondTime tv;
    VinbergAutom<Tint> vin =
        ComputeAutomorphismVinberg<T, Tint, Tgroup>(eMat, os);
    double dv = static_cast<double>(tv.const_eval_int64());
    tvin += dv;
    // direct: order off the full family
    MicrosecondTime td;
    std::vector<MyMatrix<Tint>> gens =
        ArithmeticAutomorphismGroup<T, Tint, Tgroup>(eMat, os);
    CanonicVectorFamily<Tint> fam = GetCanonicVectorFamily<T, Tint>(eMat, os);
    mpz_class aut = OrderFromGens<Tint, Tgroup>(fam.get_full(), gens);
    double dd = static_cast<double>(td.const_eval_int64());
    tdir += dd;
    bool order_ok = (vin.order == aut);
    bool gens_ok = true;
    for (auto &g : vin.ListGen) {
      if (g * G_int * g.transpose() != G_int) {
        gens_ok = false;
        break;
      }
    }
    os << "class " << (i + 1) << ": vin_order=" << vin.order
       << " (|W|=" << vin.weyl_order << " x res=" << vin.residual_order
       << ") direct=" << aut << " match=" << (order_ok ? "yes" : "NO")
       << " gens_ok=" << (gens_ok ? "yes" : "NO")
       << " | vin=" << (dv / 1000) << "ms dir=" << (dd / 1000) << "ms\n";
    if (!order_ok || !gens_ok) {
      n_error++;
    }
  }
  os << "TOTAL: vinberg=" << (tvin / 1000) << "ms direct=" << (tdir / 1000)
     << "ms ratio=" << (tvin / tdir) << "\n";
  if (n_error > 0) {
    std::cerr << "AUTVINBERG TEST: " << n_error << " errors\n";
    throw TerminalException{1};
  }
  os << "AUTVINBERG TEST: all " << ListGram.size() << " checks passed\n";
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 2) {
      std::cerr << "TEST_AutVinberg [FileListGram]\n";
      return -1;
    }
    std::string FileListGram = argv[1];
    using T = mpq_class;
    process<T>(FileListGram, std::cerr);
    std::cerr << "Normal termination of TEST_AutVinberg\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TEST_AutVinberg\n";
    exit(e.eVal);
  }
  runtime(time);
}
