// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_POLY_DATABASE_ORBITS_COMMON_H_
#define SRC_DUALDESC_POLY_DATABASE_ORBITS_COMMON_H_

// Shared, storage-agnostic logic for the DatabaseOrbits (on-disk persistence)
// and NoSaveDatabaseOrbits (no persistence, e.g. WASM) classes. Everything that
// does not touch the filesystem lives here, so the two derived classes only
// have to provide their persistence layer. This is a CRTP base: the few common
// methods that read from the (possibly absent) database dispatch to the derived
// class via derived(). No file I/O headers are included here, so the base can
// be used on platforms without a filesystem.

// clang-format off
#include "Balinski_basic.h"
#include <string>
// clang-format on

template <typename TbasicBank, typename Derived>
struct DatabaseOrbitsCommon {
public:
  using Tgroup = typename TbasicBank::Tgroup;
  using T = typename TbasicBank::T;
  using Tint = typename TbasicBank::Tint;
  Tint CritSiz;
  TbasicBank &bb;

protected:
  bool NeedToFlush;
  bool AdvancedTerminationCriterion;
  std::ostream &os;
  std::string strPresChar;
  HumanTime time;

  DatabaseOrbitsCommon(TbasicBank &_bb,
                       const bool &_AdvancedTerminationCriterion,
                       std::ostream &_os)
      : CritSiz(_bb.EXT.cols() - 2), bb(_bb), NeedToFlush(true),
        AdvancedTerminationCriterion(_AdvancedTerminationCriterion), os(_os) {
    strPresChar = "|EXT|=" + std::to_string(bb.nbRow) + "/" +
                  std::to_string(bb.nbCol) +
                  " |GRP|=" + std::to_string(bb.GRP.size());
  }
  Derived &derived() { return static_cast<Derived &>(*this); }
  Derived const &derived() const {
    return static_cast<Derived const &>(*this);
  }

public:
  void print_status() const {
#ifdef TRACK_RUN
    os << "RDD: Status : orbit=(" << bb.foc.nbOrbit << "," << bb.foc.nbOrbitDone
       << "," << (bb.foc.nbOrbit - bb.foc.nbOrbitDone) << ") facet=("
       << bb.foc.TotalNumber << "," << (bb.foc.TotalNumber - bb.foc.nbUndone)
       << "," << bb.foc.nbUndone << ")"
       << " " << strPresChar << "\n\n";
#endif
  }
  vectface get_runtime_testcase() const {
    size_t n_orbit = derived().preload_nb_orbit();
    size_t n_target = 100;
    int nbRow = bb.nbRow;
#ifdef DEBUG_RECURSIVE_DUAL_DESC
    os << "RDD: get_runtime_testcase n_orbit=" << n_orbit
       << " n_target=" << n_target << " nbRow=" << nbRow << "\n";
#endif
    if (n_orbit == 0) {
      vectface vf(nbRow);
      for (size_t i = 0; i < n_target; i++) {
        Face f = RandomFace(nbRow);
        vf.push_back(f);
      }
      return vf;
    } else {
      vectface vfo = derived().ReadDatabase(n_target);
      return vectface_reduction(vfo, nbRow);
    }
  }
  DatabaseAction determine_action_database(std::string const &choice) {
    if (choice == "load")
      return DatabaseAction::simple_load;
    if (choice == "guess")
      return DatabaseAction::guess;
    CanonicStrategy choice_i = bb.convert_string_method(choice);
    if (bb.canonic_method == choice_i)
      return DatabaseAction::simple_load;
    return DatabaseAction::recompute_and_shuffle;
  }
  void FuncInsert(Face const &face) { bb.FuncInsert(face); }
  void FuncInsertPair(Face const &face) { bb.FuncInsertPair(face); }
  void FuncPutOrbitAsDone(size_t const &i_orb) {
    bb.FuncPutOrbitAsDone(i_orb);
    print_status();
  }
  Face ComputeIntersectionUndone() const {
    size_t n_row = bb.EXT.rows();
    Face eSetReturn(n_row);

    // don't do full computation if many orbit remaining
    // for some polytopes only the last orbit sets eSetReturn = 0
    // resulting in large slowdowns here
    // alternative fix: enumerate in decending order
    if (bb.foc.nbOrbit - bb.foc.nbOrbitDone > 1000)
      return eSetReturn;

    for (size_t i_row = 0; i_row < n_row; i_row++)
      eSetReturn[i_row] = 1;
    typename TbasicBank::iterator_face iter = bb.begin_face_undone();
    while (iter != bb.end_face_undone()) {
      eSetReturn &= OrbitIntersection(bb.GRP, *iter);
      if (eSetReturn.count() == 0) {
        return eSetReturn;
      }
      iter++;
    }
    return eSetReturn;
  }
  size_t FuncNumberOrbit() const { return bb.foc.nbOrbit; }
  bool IsFinished() const { return bb.foc.nbOrbit == bb.foc.nbOrbitDone; }
  DataFacet<T, Tgroup> FuncGetMinimalUndoneOrbit() {
    DataFacet<T, Tgroup> data = bb.FuncGetMinimalUndoneOrbit();
#ifdef DEBUG_RECURSIVE_DUAL_DESC
    os << "RDD: " << strPresChar << " Considering orbit " << data.SelectedOrbit
       << " |inc|=" << data.eInc.count() << " |stab|=" << data.Stab.size()
       << "\n";
#endif
    return data;
  }
  bool GetTerminationStatus() const {
    auto get_val = [&]() -> bool {
      if (bb.foc.nbOrbitDone > 0) {
        if (bb.foc.nbUndone <= CritSiz) {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
          os << "RDD: Termination by classic Balinski criterion nbUndone="
             << bb.foc.nbUndone << "\n";
#endif
          return true;
        }
        Face eSetUndone = ComputeIntersectionUndone();
        if (eSetUndone.count() > 0) {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
          os << "RDD: Termination by linear programming criterion "
                "|eSetUndone|="
             << eSetUndone.count() << "\n";
#endif
          return true;
        }
      }
      if (AdvancedTerminationCriterion)
        return EvaluationConnectednessCriterion_Serial(bb, os);
      return false;
    };
    if (get_val()) {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
      os << "RDD: End of computation, nbObj=" << bb.foc.TotalNumber
         << " |EXT|=" << bb.nbRow << "/" << bb.nbCol
         << " time=" << time.const_eval() << "\n";
#endif
      return true;
    }
    return false;
  }
  UndoneOrbitInfo<Tint> GetTerminationInfo() const {
    return {bb.foc.nbOrbitDone, bb.foc.nbUndone, ComputeIntersectionUndone()};
  }
};

// clang-format off
#endif  // SRC_DUALDESC_POLY_DATABASE_ORBITS_COMMON_H_
// clang-format on
