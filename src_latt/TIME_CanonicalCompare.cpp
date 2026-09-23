// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Time the monolithic vs the root-decomposed canonical form.
// clang-format off
#include "NumberTheory.h"
#include "Group.h"
#include "Permutation.h"
#include "LatticeRootDecomposition.h"
// clang-format on

template <typename T>
void process(std::string const &FileListGram, std::ostream &os) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  std::vector<MyMatrix<T>> ListGram = ReadListMatrixFile<T>(FileListGram);
  double tmono = 0, tdec = 0;
  for (size_t i = 0; i < ListGram.size(); i++) {
    MyMatrix<T> const &eMat = ListGram[i];
    MicrosecondTime tm;
    (void)ComputeCanonicalForm<T, Tint, Tgroup>(eMat, os);
    double dm = static_cast<double>(tm.const_eval_int64());
    MicrosecondTime td;
    (void)ComputeCanonicalFormRootDecomposed<T, Tint, Tgroup>(eMat, os);
    double dd = static_cast<double>(td.const_eval_int64());
    tmono += dm;
    tdec += dd;
    os << "class " << (i + 1) << ": monolithic=" << (dm / 1000)
       << "ms  decomposed=" << (dd / 1000) << "ms\n";
  }
  os << "TOTAL: monolithic=" << (tmono / 1000) << "ms  decomposed="
     << (tdec / 1000) << "ms  ratio=" << (tdec / tmono) << "\n";
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 2) {
      std::cerr << "TIME_CanonicalCompare [FileListGram]\n";
      return -1;
    }
    std::string FileListGram = argv[1];
    using T = mpq_class;
    process<T>(FileListGram, std::cerr);
    std::cerr << "Normal termination of TIME_CanonicalCompare\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TIME_CanonicalCompare\n";
    exit(e.eVal);
  }
  runtime(time);
}
