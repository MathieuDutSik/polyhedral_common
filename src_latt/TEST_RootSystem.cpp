// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Group.h"
#include "Permutation.h"
#include "LatticeRootSystem.h"
#include "LatticeStabEquiCan.h"
#include "Positivity.h"
// clang-format on

/*
  Validates the root-system recognition. For each Gram matrix it prints
  the Dynkin type and the analytic Weyl order, and checks that |W(R)|
  divides the true automorphism group order |Aut(L)| (computed by the
  graph method): the Weyl group is a subgroup of Aut(L), so this must
  hold, and it is a strong test of the recognition.
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
  for (size_t i = 0; i < ListGram.size(); i++) {
    MyMatrix<T> const &eMat = ListGram[i];
    RootSystemData<Tint> rs = ComputeRootSystem<T, Tint>(eMat, os);
    mpz_class w = WeylGroupOrder(rs);
    std::string type;
    for (auto &c : rs.components) {
      type += " ";
      type += c.letter;
      type += std::to_string(c.rank);
    }
    // |Aut(L)| via the graph method, order off the full-rank family.
    std::vector<MyMatrix<Tint>> gens =
        ArithmeticAutomorphismGroup<T, Tint, Tgroup>(eMat, os);
    // read the order from the permutation action on the invariant family
    CanonicVectorFamily<Tint> fam = GetCanonicVectorFamily<T, Tint>(eMat, os);
    MyMatrix<Tint> SHV = fam.get_full();
    int n_row = SHV.rows();
    std::unordered_map<MyVector<Tint>, Tidx> MapRow;
    for (int r = 0; r < n_row; r++) {
      MapRow[GetMatrixRow(SHV, r)] = static_cast<Tidx>(r);
    }
    std::vector<Telt> ListPermGens;
    for (auto &eGen : gens) {
      MyMatrix<Tint> Mtr = eGen.transpose();
      std::vector<Tidx> ePerm(n_row);
      for (int r = 0; r < n_row; r++) {
        MyVector<Tint> v = Mtr * GetMatrixRow(SHV, r);
        ePerm[r] = MapRow.at(v);
      }
      ListPermGens.push_back(Telt(ePerm));
    }
    Tgroup grp(ListPermGens, n_row);
    mpz_class aut = grp.size();
    bool divides = (aut % w == 0);
    os << "class " << (i + 1) << ": roots=" << rs.n_roots
       << " rank=" << rs.rank << " type=" << type << " |W|=" << w
       << " |Aut|=" << aut << " W|Aut=" << (divides ? "yes" : "NO") << "\n";
    if (!divides) {
      std::cerr << "ROOTSYS TEST: |W| does not divide |Aut| on class "
                << (i + 1) << "\n";
      n_error++;
    }
  }
  if (n_error > 0) {
    std::cerr << "ROOTSYS TEST: " << n_error << " errors\n";
    throw TerminalException{1};
  }
  os << "ROOTSYS TEST: all " << ListGram.size() << " root systems valid\n";
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 2) {
      std::cerr << "TEST_RootSystem [FileListGram]\n";
      return -1;
    }
    std::string FileListGram = argv[1];
    using T = mpq_class;
    process<T>(FileListGram, std::cerr);
    std::cerr << "Normal termination of TEST_RootSystem\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TEST_RootSystem\n";
    exit(e.eVal);
  }
  runtime(time);
}
