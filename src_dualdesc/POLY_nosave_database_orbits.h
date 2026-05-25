// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_POLY_DATABASE_ORBITS_H_
#define SRC_DUALDESC_POLY_DATABASE_ORBITS_H_

// Holds the `NoSaveDatabaseOrbits<TbasicBank>` class extracted from
// POLY_RecursiveDualDesc.h. The class wraps a `TbasicBank` (face-orbit
// store) with no on-disk persistence. This is to be used for context
// without save such as WASM platform. Macros used internally
// (DEBUG_RECURSIVE_DUAL_DESC, TIMINGS_RECURSIVE_DUAL_DESC, TRACK_DATABASE,
// TRACK_RUN) are intentionally shared with POLY_RecursiveDualDesc.h --
// when this file is included from there, all four toggles are already
// configured.

// clang-format off
#include "Balinski_basic.h"
#include <string>
#include <unordered_set>
// clang-format on

template <typename TbasicBank> struct NoSaveDatabaseOrbits {
public:
  using Tgroup = typename TbasicBank::Tgroup;
  using T = typename TbasicBank::T;
  using Telt = typename Tgroup::Telt;
  using Tint = typename TbasicBank::Tint;
  Tint CritSiz;
  TbasicBank &bb;

private:
  /* TRICK 7: Using separate files for faces and status allow us to gain
     locality. The faces are written one by one while the access to status is
     random */
  bool NeedToFlush;
  bool AdvancedTerminationCriterion;
  std::ostream &os;
  size_t delta;
  std::string strPresChar;
  HumanTime time;

public:
  // method encodes the algorithm used for the database and essentially applies
  // only to the canonic.
  // ---If method file is absent then we assume it was computed with the
  // default.
  // ---Otherwise we read it.
  void write_method([[maybe_unused]] std::string const &eFileMethod, [[maybe_unused]] int const &method) const {
  }
  int read_method([[maybe_unused]] std::string const &eFileMethod) const {
    return bb.get_default_strategy();
  }
  bool is_database_present() const {
    return false;
  }
  NoSaveDatabaseOrbits() = delete;
  NoSaveDatabaseOrbits(const NoSaveDatabaseOrbits<TbasicBank> &) = delete;
  NoSaveDatabaseOrbits(NoSaveDatabaseOrbits<TbasicBank> &&) = delete;
  NoSaveDatabaseOrbits &operator=(const NoSaveDatabaseOrbits<TbasicBank> &) = delete;
  void print_status() const {
#ifdef TRACK_RUN
    os << "RDD: Status : orbit=(" << bb.foc.nbOrbit << "," << bb.foc.nbOrbitDone
       << "," << (bb.foc.nbOrbit - bb.foc.nbOrbitDone) << ") facet=("
       << bb.foc.TotalNumber << "," << (bb.foc.TotalNumber - bb.foc.nbUndone)
       << "," << bb.foc.nbUndone << ")"
       << " " << strPresChar << "\n\n";
#endif
  }
  NoSaveDatabaseOrbits(TbasicBank &bb,
                       const bool &_AdvancedTerminationCriterion, std::ostream &os)
      : CritSiz(bb.EXT.cols() - 2), bb(bb),
        AdvancedTerminationCriterion(_AdvancedTerminationCriterion), os(os) {
    strPresChar = "|EXT|=" + std::to_string(bb.nbRow) + "/" +
                  std::to_string(bb.nbCol) +
                  " |GRP|=" + std::to_string(bb.GRP.size());
    delta = bb.delta;
    NeedToFlush = true;
    int val = bb.get_default_strategy();
    bb.the_method = val;
  }
  size_t preload_nb_orbit() const {
    return 0;
  }
  void LoadDatabase() {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
    os << "RDD: No database present\n";
#endif
    print_status();
  }
  vectface ReadDatabase([[maybe_unused]] size_t const &n_read) const {
    vectface vfo(bb.delta + 1);
#ifdef DEBUG_RECURSIVE_DUAL_DESC
    os << "RDD: No database present\n";
#endif
    return vfo;
  }
  vectface get_runtime_testcase() const {
    size_t n_orbit = preload_nb_orbit();
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
      vectface vfo = ReadDatabase(n_target);
      return vectface_reduction(vfo, nbRow);
    }
  }
  int determine_action_database(std::string const &choice) {
    if (choice == "load")
      return DATABASE_ACTION__SIMPLE_LOAD;
    if (choice == "guess")
      return DATABASE_ACTION__GUESS;
    int choice_i = bb.convert_string_method(choice);
    if (bb.the_method == choice_i)
      return DATABASE_ACTION__SIMPLE_LOAD;
    return DATABASE_ACTION__RECOMPUTE_AND_SHUFFLE;
  }
  void set_method(int const &the_method) {
    bb.the_method = the_method;
  }
  void DirectAppendDatabase(vectface &&vf) {
    bb.clear();
    size_t n_orbit = vf.size();
    for (size_t i_orbit = 0; i_orbit < n_orbit; i_orbit++) {
      Face f = vf[i_orbit];
      Face f_red(bb.delta);
      for (size_t u = 0; u < bb.delta; u++) {
        bool val = f[u];
        f_red[u] = val;
      }
      bool status = f[bb.delta];
      std::pair<Face, Tint> eEnt = bb.foc.recConvert.ConvertFace(f_red);
      bb.InsertListOrbitEntry(f_red, i_orbit);
      bb.InsertEntryDatabase(eEnt, status, i_orbit);
    }
    print_status();
  }
  ~NoSaveDatabaseOrbits() {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
    os << "RDD: Clean closing of the NoSaveDatabaseOrbits\n";
#endif
  }
  // FuncListOrbitIncidence() {
  FaceOrbitsizeTableContainer<Tint> GetListFaceOrbitsize() {
    NeedToFlush = false;
    return bb.GetListFaceOrbitsize();
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
#endif  // SRC_DUALDESC_POLY_DATABASE_ORBITS_H_
// clang-format on
