// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_POLY_NOSAVE_DATABASE_ORBITS_H_
#define SRC_DUALDESC_POLY_NOSAVE_DATABASE_ORBITS_H_

// Holds the `NoSaveDatabaseOrbits<TbasicBank>` class extracted from
// POLY_RecursiveDualDesc.h. The class wraps a `TbasicBank` (face-orbit
// store) with no on-disk persistence. This is to be used for context
// without save such as WASM platform. The storage-agnostic logic is shared
// with DatabaseOrbits through the DatabaseOrbitsCommon base. Macros used
// internally (DEBUG_RECURSIVE_DUAL_DESC, TIMINGS_RECURSIVE_DUAL_DESC,
// TRACK_DATABASE, TRACK_RUN) are intentionally shared with
// POLY_RecursiveDualDesc.h -- when this file is included from there, all four
// toggles are already configured.

// clang-format off
#include "POLY_database_orbits_common.h"
#include <string>
// clang-format on

template <typename TbasicBank>
struct NoSaveDatabaseOrbits
    : public DatabaseOrbitsCommon<TbasicBank, NoSaveDatabaseOrbits<TbasicBank>> {
public:
  using Base =
      DatabaseOrbitsCommon<TbasicBank, NoSaveDatabaseOrbits<TbasicBank>>;
  using typename Base::Tint;
  using Base::bb;
  using Base::NeedToFlush;
  using Base::os;
  using Base::print_status;

public:
  // method encodes the algorithm used for the database and essentially applies
  // only to the canonic.
  // ---If method file is absent then we assume it was computed with the
  // default.
  // ---Otherwise we read it.
  void write_method([[maybe_unused]] std::string const &eFileMethod,
                    [[maybe_unused]] CanonicStrategy const &method) const {}
  CanonicStrategy read_method([[maybe_unused]] std::string const &eFileMethod) const {
    return bb.get_default_strategy();
  }
  bool is_database_present() const { return false; }
  NoSaveDatabaseOrbits() = delete;
  NoSaveDatabaseOrbits(const NoSaveDatabaseOrbits<TbasicBank> &) = delete;
  NoSaveDatabaseOrbits(NoSaveDatabaseOrbits<TbasicBank> &&) = delete;
  NoSaveDatabaseOrbits &
  operator=(const NoSaveDatabaseOrbits<TbasicBank> &) = delete;
  NoSaveDatabaseOrbits(TbasicBank &bb,
                       const bool &_AdvancedTerminationCriterion,
                       std::ostream &os)
      : Base(bb, _AdvancedTerminationCriterion, os) {
    CanonicStrategy val = bb.get_default_strategy();
    bb.canonic_method = val;
  }
  size_t preload_nb_orbit() const { return 0; }
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
  void set_method(CanonicStrategy const &canonic_method) { bb.canonic_method = canonic_method; }
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
};

// clang-format off
#endif  // SRC_DUALDESC_POLY_NOSAVE_DATABASE_ORBITS_H_
// clang-format on
