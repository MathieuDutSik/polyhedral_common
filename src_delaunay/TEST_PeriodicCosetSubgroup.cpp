// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "PeriodicDelaunay.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

/*
  The subgroup of the transformations preserving a periodic point set is
  checked against the direct computation: a group of transformations is
  built by hand, all its elements are enumerated, and the elements that
  preserve the point set (as decided by PeriodicCosetPermutation, itself
  checked against a membership test in TEST_PeriodicCVP) must be exactly
  the elements of the subgroup returned.
 */

using T = mpq_class;
using Tring = mpz_class;
using Tidx = uint32_t;
using Telt = permutalib::SingleSidedPerm<Tidx>;
using Tint_grp = mpz_class;
using Tgroup = permutalib::Group<Telt, Tint_grp>;
using Ttrans = PeriodicAffineTransform<Tring>;

void check(bool test, std::string const &context) {
  if (!test) {
    std::cerr << "TEST_PeriodicCosetSubgroup failed: " << context << "\n";
    throw TerminalException{1};
  }
}

int main() {
  try {
    std::ostream &os = std::cerr;
    int n = 2;
    // The point set Z^2 + {(0,0), (1/2,1/2)}, of denominator 2.
    MyMatrix<T> Cosets(2, n);
    Cosets(0, 0) = 0;
    Cosets(0, 1) = 0;
    Cosets(1, 0) = T(1) / T(2);
    Cosets(1, 1) = T(1) / T(2);
    PeriodicPointSet<Tring> pps = PeriodicPointSetFromRational<Tring>(Cosets);
    check(pps.N == 2, "denominator of the point set");
    //
    // Two transformations of denominator 2, acting on a set of four
    // formal points so that the permutation determines the transformation:
    //   g1: the exchange of the coordinates, which preserves the set,
    //   g2: the translation by (1/2, 0), which does not (it maps the
    //       coset (1/2,1/2) to (0,1/2), which is not a coset).
    MyMatrix<Tring> A_swap(n, n);
    A_swap(0, 0) = 0;
    A_swap(0, 1) = 1;
    A_swap(1, 0) = 1;
    A_swap(1, 1) = 0;
    MyVector<Tring> w_zero = ZeroVector<Tring>(n);
    Ttrans g1{A_swap, w_zero, Tring(2)};
    MyMatrix<Tring> A_id = IdentityMat<Tring>(n);
    MyVector<Tring> w_half(n);
    w_half(0) = 1;
    w_half(1) = 0;
    Ttrans g2{A_id, w_half, Tring(2)};
    check(IsPeriodicPointSetPreserved(pps, g1), "g1 preserves the set");
    check(!IsPeriodicPointSetPreserved(pps, g2), "g2 does not preserve it");
    //
    // The group is realized on four points, the permutation of the points
    // determining which product of g1 and g2 is meant. Both generators are
    // involutions and they commute, so the group is (Z/2)^2 acting on the
    // four pairs of exponents.
    //   point i encodes (a, b) with i = a + 2 b, meaning g1^a g2^b.
    auto get_perm = [&](int i_gen) -> Telt {
      std::vector<Tidx> eList(4);
      for (int i = 0; i < 4; i++) {
        int a = i % 2, b = i / 2;
        if (i_gen == 0) {
          a = 1 - a;
        } else {
          b = 1 - b;
        }
        eList[i] = static_cast<Tidx>(a + 2 * b);
      }
      return Telt(eList);
    };
    std::vector<Telt> LGen{get_perm(0), get_perm(1)};
    Tgroup GRP(LGen, 4);
    check(GRP.size() == 4, "the group has order 4");
    // The transformation realized by a permutation, read off the image of
    // the identity point.
    auto f_trans = [&](Telt const &eElt) -> Ttrans {
      int img = OnPoints(0, eElt);
      int a = img % 2, b = img / 2;
      Ttrans ret = IdentityPeriodicAffineTransform<Tring>(n, Tring(2));
      if (a == 1) {
        ret = ret * g1;
      }
      if (b == 1) {
        ret = ret * g2;
      }
      return ret;
    };
    Tgroup GRPsub =
        PeriodicCosetPreservingSubgroup<Tring, Tgroup, decltype(f_trans)>(
            pps, GRP, f_trans, os);
    // Direct computation: the elements whose transformation preserves the
    // point set. The subgroup acts on the same four points, so the two can
    // be compared element by element.
    size_t n_pres = 0;
    for (auto &eElt : GRP.get_all_element()) {
      if (IsPeriodicPointSetPreserved(pps, f_trans(eElt))) {
        n_pres++;
      }
    }
    std::cerr << "|GRP|=" << GRP.size() << " |GRPsub|=" << GRPsub.size()
              << " n_preserving=" << n_pres << "\n";
    check(UniversalScalarConversion<size_t, Tint_grp>(GRPsub.size()) == n_pres,
          "the subgroup has the number of preserving elements");
    // Every element of the subgroup preserves the point set.
    for (auto &eElt : GRPsub.get_all_element()) {
      check(IsPeriodicPointSetPreserved(pps, f_trans(eElt)),
            "an element of the subgroup preserves the point set");
    }
    check(n_pres == 2, "exactly the two elements without g2");
    //
    // The stabilizer of a Delaunay polytope of the same point set for the
    // identity form. The point set Z^2 + {(0,0), (1/2,1/2)} is the square
    // lattice rotated and scaled, its Delaunay cells are squares, and the
    // square of vertices (0,0), (1,0), (1/2,1/2), (1/2,-1/2) -- of
    // numerators (0,0), (2,0), (1,1), (1,-1) over N = 2 -- has for
    // stabilizer the dihedral group of order 8.
    {
      MyMatrix<T> GramMat = IdentityMat<T>(n);
      // Short vectors generating Z^2, with both signs so that the set is
      // preserved by the symmetries.
      MyMatrix<T> SHV(4, n);
      SHV(0, 0) = 1;  SHV(0, 1) = 0;
      SHV(1, 0) = -1; SHV(1, 1) = 0;
      SHV(2, 0) = 0;  SHV(2, 1) = 1;
      SHV(3, 0) = 0;  SHV(3, 1) = -1;
      MyMatrix<Tring> EXT_scaled(4, n);
      EXT_scaled(0, 0) = 0;  EXT_scaled(0, 1) = 0;
      EXT_scaled(1, 0) = 2;  EXT_scaled(1, 1) = 0;
      EXT_scaled(2, 0) = 1;  EXT_scaled(2, 1) = 1;
      EXT_scaled(3, 0) = 1;  EXT_scaled(3, 1) = -1;
      Tgroup GRPstab =
          PeriodicDelaunay_Stabilizer<T, Tring, Tgroup>(GramMat, pps, SHV,
                                                        EXT_scaled, os);
      std::cerr << "|stabilizer of the square|=" << GRPstab.size() << "\n";
      check(GRPstab.size() == 8, "the square has a stabilizer of order 8");
      //
      // Equivalence. The square translated by the coset (1/2,1/2), of
      // numerators shifted by (1,1), is another Delaunay cell of the same
      // point set, so the two must be equivalent, and the equivalence
      // found must preserve the point set and realize the correspondence.
      MyMatrix<Tring> EXT_shift(4, n);
      for (int i = 0; i < 4; i++) {
        EXT_shift(i, 0) = EXT_scaled(i, 0) + 1;
        EXT_shift(i, 1) = EXT_scaled(i, 1) + 1;
      }
      std::optional<PeriodicAffineTransform<Tring>> opt_equiv =
          PeriodicDelaunay_TestEquivalence<T, Tring, Tgroup>(
              GramMat, pps, SHV, EXT_scaled, EXT_shift, os);
      check(opt_equiv.has_value(), "the translated square is equivalent");
      check(IsPeriodicPointSetPreserved(pps, *opt_equiv),
            "the equivalence preserves the point set");
      // The equivalence maps the first vertex set onto the second.
      {
        MyMatrix<T> M = transform_traits<PeriodicAffineTransform<Tring>>::
            to_field<T>(*opt_equiv);
        std::set<MyVector<Tring>> set2;
        for (int i = 0; i < 4; i++) {
          set2.insert(GetMatrixRow(EXT_shift, i));
        }
        T N_T = UniversalScalarConversion<T, Tring>(pps.N);
        for (int i = 0; i < 4; i++) {
          // The image of the point, back in numerator form.
          MyVector<Tring> img(n);
          for (int j = 0; j < n; j++) {
            T eSum = M(0, j + 1);
            for (int k = 0; k < n; k++) {
              T coord =
                  UniversalScalarConversion<T, Tring>(EXT_scaled(i, k)) / N_T;
              eSum += coord * M(k + 1, j + 1);
            }
            img(j) = UniversalScalarConversion<Tring, T>(eSum * N_T);
          }
          check(set2.count(img) == 1,
                "the equivalence maps a vertex onto a vertex of the second");
        }
      }
      // A square that is NOT a cell of the point set: the same square
      // translated by (1/2, 0), whose numerators are shifted by (1,0). It
      // is isometric to the first but its vertices are not points of the
      // set, so no equivalence may be returned.
      MyMatrix<Tring> EXT_bad(4, n);
      for (int i = 0; i < 4; i++) {
        EXT_bad(i, 0) = EXT_scaled(i, 0) + 1;
        EXT_bad(i, 1) = EXT_scaled(i, 1);
      }
      std::optional<PeriodicAffineTransform<Tring>> opt_bad =
          PeriodicDelaunay_TestEquivalence<T, Tring, Tgroup>(
              GramMat, pps, SHV, EXT_scaled, EXT_bad, os);
      check(!opt_bad.has_value(),
            "no equivalence to a square outside the point set");
    }
    std::cerr << "Normal termination of TEST_PeriodicCosetSubgroup\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of TEST_PeriodicCosetSubgroup\n";
    exit(e.eVal);
  }
}
