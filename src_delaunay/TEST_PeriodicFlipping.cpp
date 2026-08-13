// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "PeriodicDelaunay.h"
#include "IsoDelaunayDomains.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

/*
  The flipping of a periodic Delaunay tessellation across a wall of its
  iso-Delaunay domain, checked against the direct enumeration on the other
  side of the wall.

  The tessellation of a generic form is built, the defining inequalities of
  its iso-Delaunay domain are collected and, for every irredundant wall,
  FlippingLtype produces the tessellation of the adjacent domain. What the
  flip announces is then confronted with what the point set says: at a form
  taken beyond the wall, every cell of the flipped tessellation has to be a
  Delaunay cell of the point set, and the direct enumeration there has to
  find the same number of orbits.
 */

using T = mpq_class;
using Tring = mpz_class;
using Tidx = uint32_t;
using Telt = permutalib::SingleSidedPerm<Tidx>;
using Tint_grp = mpz_class;
using Tgroup = permutalib::Group<Telt, Tint_grp>;

void check(bool test, std::string const &context) {
  if (!test) {
    std::cerr << "TEST_PeriodicFlipping failed: " << context << "\n";
    throw TerminalException{1};
  }
}

/*
  The enumeration of the periodic Delaunay cells of a form, in the shape the
  iso-Delaunay machinery consumes. The short vectors of the two forms used
  below are the four unit vectors, which is why they are hardcoded here.
 */
DelaunayTesselationIneq<T, Tgroup>
GetTesselation(MyMatrix<T> const &GramMat, PeriodicPointSet<Tring> const &pps,
               std::vector<std::vector<T>> const &ListGram, size_t &n_orbit,
               std::ostream &os) {
  int n = GramMat.rows();
  MyMatrix<T> SHV(4, n);
  SHV(0, 0) = 1;  SHV(0, 1) = 0;
  SHV(1, 0) = -1; SHV(1, 1) = 0;
  SHV(2, 0) = 0;  SHV(2, 1) = 1;
  SHV(3, 0) = 0;  SHV(3, 1) = -1;
  MyMatrix<Tring> Graver = UniversalMatrixConversion<Tring, T>(SHV);
  PolyHeuristicSerial<Tint_grp> AllArr =
      AllStandardHeuristicSerial<T, Tint_grp>(n + 1, os);
  PeriodicDataDelaunay<T, Tring, Tgroup> data =
      GetPeriodicDataDelaunay<T, Tring, Tgroup>(GramMat, pps, SHV, Graver,
                                                AllArr, os);
  PeriodicDataDelaunayFunc<T, Tring, Tgroup> data_func{data};
  auto f_incorrect =
      [&]([[maybe_unused]] PeriodicDelaunay_Obj<Tring, Tgroup> const &x)
      -> bool { return false; };
  auto opt = EnumerateAndStore_Serial(data_func, f_incorrect, 0);
  check(opt.has_value(), "the enumeration terminated");
  n_orbit = opt->size();
  return BuildPeriodicDelaunayTesselationIneq<T, Tring, Tgroup>(*opt, pps,
                                                                ListGram, os);
}

int main() {
  HumanTime time;
  try {
    std::ostream &os = std::cerr;
    int n = 2;
    // The point set Z^2 + {(0,0), (1/3,0)}, whose scaled frame is
    // 3 Z^2 + {(0,0), (1,0)}.
    MyMatrix<T> Cosets(2, n);
    Cosets(0, 0) = 0;
    Cosets(0, 1) = 0;
    Cosets(1, 0) = T(1) / T(3);
    Cosets(1, 1) = 0;
    PeriodicPointSet<Tring> pps = PeriodicPointSetFromRational<Tring>(Cosets);
    // The T-space of all the symmetric matrices of order n.
    std::vector<MyMatrix<T>> ListMat;
    std::vector<std::vector<T>> ListGram;
    for (int i = 0; i < n; i++) {
      for (int j = i; j < n; j++) {
        MyMatrix<T> eMat = ZeroMatrix<T>(n, n);
        eMat(i, j) = 1;
        eMat(j, i) = 1;
        ListMat.push_back(eMat);
        ListGram.push_back(GetLineVector(eMat));
      }
    }
    int dimSpace = ListMat.size();
    auto get_matrix = [&](MyVector<T> const &V) -> MyMatrix<T> {
      MyMatrix<T> eMat = ZeroMatrix<T>(n, n);
      for (int u = 0; u < dimSpace; u++) {
        MatAddMul(eMat, V(u), ListMat[u]);
      }
      return eMat;
    };
    // The generic form and its coordinates in the T-space.
    MyVector<T> Gvec(dimSpace);
    Gvec(0) = 1;
    Gvec(1) = T(1) / T(7);
    Gvec(2) = 1;
    MyMatrix<T> GramMat = get_matrix(Gvec);
    size_t n_orbit = 0;
    DelaunayTesselationIneq<T, Tgroup> DTI =
        GetTesselation(GramMat, pps, ListGram, n_orbit, os);
    std::cerr << "|orbits|=" << n_orbit << "\n";
    check(IsPeriodicDelaunayTesselation<T, Tring, Tgroup>(DTI, GramMat, pps, os),
          "the tessellation of the form is its Delaunay tessellation");
    // The defining inequalities of the iso-Delaunay domain, whose
    // irredundant ones are the walls.
    std::vector<FullAdjInfo<T>> ListIneq =
        ComputeListIneqFromTesselationIneq<T, Tgroup>(DTI);
    MyMatrix<T> FAC = GetFACineq(ListIneq);
    std::vector<int> ListIrred = get_non_redundant_indices(FAC, os);
    MyMatrix<T> FACred = SelectRow(FAC, ListIrred);
    std::cerr << "|FAC|=" << FAC.rows() << " |FACred|=" << FACred.rows()
              << "\n";
    check(FACred.rows() > 0, "the domain has walls");
    // The form is in the interior of the domain, so it satisfies every
    // inequality strictly.
    for (int i = 0; i < FAC.rows(); i++) {
      check(GetMatrixRow(FAC, i).dot(Gvec) > 0,
            "the form is in the interior of its iso-Delaunay domain");
    }
    PolyHeuristicSerial<Tint_grp> AllArr =
        AllStandardHeuristicSerial<T, Tint_grp>(n + 1, os);
    RecordDualDescOperation<T, Tgroup> rddo(AllArr, os);
    size_t n_flip = 0;
    for (int i_irred = 0; i_irred < FACred.rows(); i_irred++) {
      // A point of the relative interior of the wall, then a form beyond
      // it: the interior of the adjacent domain is reached by going past
      // the wall from the form.
      MyVector<T> WallPt = GetSpaceInteriorPointFacet(FACred, i_irred, os);
      MyVector<T> Beyond = WallPt + (WallPt - Gvec) / T(100);
      MyMatrix<T> WallMat = get_matrix(WallPt);
      MyMatrix<T> BeyondMat = get_matrix(Beyond);
      if (!IsPositiveDefinite(WallMat, os) ||
          !IsPositiveDefinite(BeyondMat, os)) {
        // The wall is at the boundary of the cone of the positive definite
        // forms, so there is no adjacent domain to flip to.
        std::cerr << "   wall " << i_irred << ": not flippable\n";
        continue;
      }
      FullAdjInfo<T> const &eRecIneq = ListIneq[ListIrred[i_irred]];
      DelaunayTesselationIneq<T, Tgroup> DTIadj = FlippingLtype<T, Tgroup>(
          DTI, GramMat, eRecIneq.ListAdjInfo, ListGram, rddo);
      // The flipped tessellation is the Delaunay tessellation of the point
      // set beyond the wall, and the original one is not: the flip really
      // crossed the wall.
      check(IsPeriodicDelaunayTesselation<T, Tring, Tgroup>(DTIadj, BeyondMat, pps, os),
            "the flipped tessellation is the one beyond the wall");
      check(!IsPeriodicDelaunayTesselation<T, Tring, Tgroup>(DTI, BeyondMat, pps, os),
            "the tessellation before the flip is not the one beyond the wall");
      size_t n_orbit_dir = 0;
      DelaunayTesselationIneq<T, Tgroup> DTIdir =
          GetTesselation(BeyondMat, pps, ListGram, n_orbit_dir, os);
      std::cerr << "   wall " << i_irred << ": |flipped|=" << DTIadj.l_dels.size()
                << " |direct|=" << n_orbit_dir << "\n";
      check(DTIadj.l_dels.size() == n_orbit_dir,
            "the flip finds the orbits of the direct enumeration");
      // The two tessellations of the same form have the same multiset of
      // vertex counts.
      std::vector<int> l_flip, l_dir;
      for (auto &eDel : DTIadj.l_dels) {
        l_flip.push_back(eDel.EXT.rows());
      }
      for (auto &eDel : DTIdir.l_dels) {
        l_dir.push_back(eDel.EXT.rows());
      }
      std::sort(l_flip.begin(), l_flip.end());
      std::sort(l_dir.begin(), l_dir.end());
      check(l_flip == l_dir,
            "the flip and the direct enumeration agree on the cell sizes");
      n_flip++;
    }
    check(n_flip > 0, "at least one wall was flipped");
    //
    // The walk over the iso-Delaunay domains: what a flip produces has to
    // be a valid input for the next one, which a single flip from a freshly
    // enumerated tessellation does not exercise. The walk crosses a wall,
    // computes the walls of the domain it lands in and crosses one of
    // those, checking at every step that the tessellation is the Delaunay
    // tessellation of the point set for the form reached.
    DelaunayTesselationIneq<T, Tgroup> DTIcur = DTI;
    MyVector<T> Gcur = Gvec;
    size_t n_step = 0;
    for (int i_step = 0; i_step < 6; i_step++) {
      MyMatrix<T> GcurMat = get_matrix(Gcur);
      std::vector<FullAdjInfo<T>> ListIneqCur =
          ComputeListIneqFromTesselationIneq<T, Tgroup>(DTIcur);
      MyMatrix<T> FACcur = GetFACineq(ListIneqCur);
      std::vector<int> IrredCur = get_non_redundant_indices(FACcur, os);
      MyMatrix<T> FACredCur = SelectRow(FACcur, IrredCur);
      // The walls are taken in turn, so that the walk does not bounce
      // between the same two domains.
      std::optional<std::pair<MyVector<T>, FullAdjInfo<T>>> opt_wall;
      int n_wall = FACredCur.rows();
      for (int u = 0; u < n_wall; u++) {
        int i_wall = (i_step + u) % n_wall;
        MyVector<T> WallPt = GetSpaceInteriorPointFacet(FACredCur, i_wall, os);
        MyVector<T> Beyond = WallPt + (WallPt - Gcur) / T(100);
        if (IsPositiveDefinite(get_matrix(WallPt), os) &&
            IsPositiveDefinite(get_matrix(Beyond), os)) {
          opt_wall = std::pair<MyVector<T>, FullAdjInfo<T>>{
              Beyond, ListIneqCur[IrredCur[i_wall]]};
          break;
        }
      }
      if (!opt_wall) {
        // Every wall of the domain reached is at the boundary of the cone
        // of the positive definite forms.
        break;
      }
      MyVector<T> const &Beyond = opt_wall->first;
      MyMatrix<T> BeyondMat = get_matrix(Beyond);
      DelaunayTesselationIneq<T, Tgroup> DTInew = FlippingLtype<T, Tgroup>(
          DTIcur, GcurMat, opt_wall->second.ListAdjInfo, ListGram, rddo);
      check(IsPeriodicDelaunayTesselation<T, Tring, Tgroup>(DTInew, BeyondMat, pps, os),
            "the tessellation reached by the walk is a Delaunay tessellation");
      size_t n_orbit_dir = 0;
      GetTesselation(BeyondMat, pps, ListGram, n_orbit_dir, os);
      check(DTInew.l_dels.size() == n_orbit_dir,
            "the walk agrees with the direct enumeration");
      std::cerr << "   step " << i_step << ": |walls|=" << n_wall
                << " |cells|=" << DTInew.l_dels.size() << "\n";
      DTIcur = DTInew;
      Gcur = Beyond;
      n_step++;
    }
    check(n_step > 1, "the walk crossed more than one wall");
    std::cerr << "|steps of the walk|=" << n_step << "\n";
    std::cerr << "|flipped walls|=" << n_flip << "\n";
    std::cerr << "Normal termination of TEST_PeriodicFlipping\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of TEST_PeriodicFlipping\n";
    exit(e.eVal);
  }
  runtime(time);
}
