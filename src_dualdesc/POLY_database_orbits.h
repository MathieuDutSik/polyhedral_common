// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_POLY_DATABASE_ORBITS_H_
#define SRC_DUALDESC_POLY_DATABASE_ORBITS_H_

// Holds the `DatabaseOrbits<TbasicBank>` class extracted from
// POLY_RecursiveDualDesc.h. The class wraps a `TbasicBank` (face-orbit
// store) with optional on-disk persistence (eFileEXT / eFileGRP /
// eFileNB / eFileFB / eFileFF / eFileMethod) so that long enumerations
// can be resumed after a crash. Macros used internally
// (DEBUG_RECURSIVE_DUAL_DESC, TIMINGS_RECURSIVE_DUAL_DESC, TRACK_DATABASE,
// TRACK_RUN) are intentionally shared with POLY_RecursiveDualDesc.h --
// when this file is included from there, all four toggles are already
// configured.

// clang-format off
#include "Balinski_basic.h"
#include "Basic_file.h"
#include "basic_datafile.h"
#include <string>
#include <unordered_set>
// clang-format on

template <typename TbasicBank> struct DatabaseOrbits {
public:
  using Tgroup = typename TbasicBank::Tgroup;
  using T = typename TbasicBank::T;
  using Telt = typename Tgroup::Telt;
  using Tint = typename TbasicBank::Tint;
  Tint CritSiz;
  TbasicBank &bb;

private:
  std::string MainPrefix;
  std::string eFileEXT, eFileGRP, eFileNB, eFileFB, eFileFF, eFileMethod;
  /* TRICK 7: Using separate files for faces and status allow us to gain
     locality. The faces are written one by one while the access to status is
     random */
  bool SavingTrigger;
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
  void write_method(std::string const &eFileMethod, int const &method) const {
    std::ofstream os_file(eFileMethod);
    os_file << method;
  }
  int read_method(std::string const &eFileMethod) const {
#ifdef TRACK_DATABASE
    os << "SavingTrigger=" << SavingTrigger << " eFileMethod=" << eFileMethod
       << "\n";
#endif
    if (SavingTrigger) {
#ifdef TRACK_DATABASE
      os << "Running the save system\n";
#endif
      if (!IsExistingFile(eFileMethod)) {
#ifdef TRACK_DATABASE
        os << "The file does not exists\n";
#endif
        int the_method = bb.get_default_strategy();
        write_method(eFileMethod, the_method);
        return the_method;
      } else {
#ifdef TRACK_DATABASE
        os << "The file exists\n";
#endif
        std::ifstream is_file(eFileMethod);
        int method;
        is_file >> method;
        return method;
      }
    } else {
      return bb.get_default_strategy();
    }
  }
  bool is_database_present() const {
    if (IsExistingFile(eFileEXT) == false) {
      return false;
    }
    // verify that EXT file is same as bb.EXT
    MyMatrix<T> EXT = ReadMatrixFile<T>(eFileEXT);
    if (EXT == bb.EXT) {
      return true;
    }
    // else database got changed e.g. due to method change
    // remove it
    // optional future function: check for equivalence and convert
    if (SavingTrigger) {
#ifdef TRACK_DATABASE
      os << "Database got changed, removing old one\n";
#endif
      RemoveFileIfExist(eFileNB);
      RemoveFileIfExist(eFileFB);
      RemoveFileIfExist(eFileFF);
      RemoveFileIfExist(eFileEXT);
      RemoveFileIfExist(eFileGRP);
      RemoveFileIfExist(eFileMethod);
    }
    return false;
  }
  DatabaseOrbits() = delete;
  DatabaseOrbits(const DatabaseOrbits<TbasicBank> &) = delete;
  DatabaseOrbits(DatabaseOrbits<TbasicBank> &&) = delete;
  DatabaseOrbits &operator=(const DatabaseOrbits<TbasicBank> &) = delete;
  void print_status() const {
#ifdef TRACK_RUN
    os << "RDD: Status : orbit=(" << bb.foc.nbOrbit << "," << bb.foc.nbOrbitDone
       << "," << (bb.foc.nbOrbit - bb.foc.nbOrbitDone) << ") facet=("
       << bb.foc.TotalNumber << "," << (bb.foc.TotalNumber - bb.foc.nbUndone)
       << "," << bb.foc.nbUndone << ")"
       << " " << strPresChar << "\n\n";
#endif
  }
  DatabaseOrbits(TbasicBank &bb, const std::string &MainPrefix,
                 const bool &_SavingTrigger,
                 const bool &_AdvancedTerminationCriterion, std::ostream &os)
      : CritSiz(bb.EXT.cols() - 2), bb(bb), SavingTrigger(_SavingTrigger),
        AdvancedTerminationCriterion(_AdvancedTerminationCriterion), os(os) {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
    os << "RDD: MainPrefix=" << MainPrefix << "\n";
#endif
    eFileEXT = MainPrefix + ".ext";
    eFileGRP = MainPrefix + ".grp";
    eFileNB = MainPrefix + ".nb";
    eFileFB = MainPrefix + ".fb";
    eFileFF = MainPrefix + ".ff";
    eFileMethod = MainPrefix + ".method";
    strPresChar = "|EXT|=" + std::to_string(bb.nbRow) + "/" +
                  std::to_string(bb.nbCol) +
                  " |GRP|=" + std::to_string(bb.GRP.size());
    delta = bb.delta;
    NeedToFlush = true;
    int val = read_method(eFileMethod);
#ifdef DEBUG_RECURSIVE_DUAL_DESC
    os << "RDD: read_method val=" << val << "\n";
#endif
    bb.the_method = val;
    if (SavingTrigger && !is_database_present()) {
      if (!FILE_IsFileMakeable(eFileEXT)) {
        std::cerr << "Error in DatabaseOrbits: File eFileEXT=" << eFileEXT
                  << " is not makeable\n";
        throw TerminalException{1};
      }
      initial_writes();
    }
  }
  void initial_writes() {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
    os << "RDD: Creating the initial files (NB, FB, FF) with zero state\n";
#endif
    FileNumber fn(eFileNB, true);
    FileBool fb(eFileFB);
    FileFace ff(eFileFF, bb.delta);
    std::vector<uint8_t> V_empty; // empty write, maybe useless.
    fn.setval(0);
    ff.direct_write(V_empty);
    fb.direct_write(V_empty);
    write_method(eFileMethod, bb.the_method);
    std::ofstream os_grp(eFileGRP);
    os_grp << bb.GRP;
    WriteMatrixFile(eFileEXT, bb.EXT);
  }
  size_t preload_nb_orbit() const {
    if (SavingTrigger && is_database_present()) {
      FileNumber fn(eFileNB, false);
      return fn.getval();
    }
    return 0;
  }
  void LoadDatabase() {
    if (is_database_present()) {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
      os << "RDD: Opening existing files (NB, FB, FF)\n";
#endif
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      MicrosecondTime time;
#endif
      FileNumber fn(eFileNB, false);
      size_t n_orbit = fn.getval();
#ifdef DEBUG_RECURSIVE_DUAL_DESC
      os << "RDD: Loading database with n_orbit=" << n_orbit << "\n";
#endif
      FileBool fb(eFileFB, n_orbit);
      FileFace ff(eFileFF, bb.delta, n_orbit);
      for (size_t i_orbit = 0; i_orbit < n_orbit; i_orbit++) {
        Face f = ff.getface(i_orbit);
        std::pair<Face, Tint> eEnt = bb.foc.recConvert.ConvertFace(f);
        bool status = fb.getbit(i_orbit);
        bb.InsertListOrbitEntry(f, i_orbit);
        bb.InsertEntryDatabase(eEnt, status, i_orbit);
      }
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      os << "|RDD: Loading Database|=" << time << "\n";
#endif
    } else {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
      os << "RDD: No database present\n";
#endif
    }
    print_status();
  }
  vectface ReadDatabase(size_t const &n_read) const {
    vectface vfo(bb.delta + 1);
    if (is_database_present()) {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
      os << "RDD: Opening existing files (NB, FB, FF)\n";
#endif
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      MicrosecondTime time;
#endif
      FileNumber fn(eFileNB, false);
      size_t n_orbit = fn.getval();
#ifdef DEBUG_RECURSIVE_DUAL_DESC
      os << "RDD: Reading database with n_orbit=" << n_orbit << "\n";
#endif
      FileBool fb(eFileFB, n_orbit);
      FileFace ff(eFileFF, bb.delta, n_orbit);
      size_t n_read_eff = std::min(n_read, n_orbit);
      for (size_t i_orbit = 0; i_orbit < n_read_eff; i_orbit++) {
        Face f = ff.getface(i_orbit);
        bool status = fb.getbit(i_orbit);
        Face f_insert(bb.delta + 1);
        for (size_t u = 0; u < bb.delta; u++) {
          f_insert[u] = f[u];
        }
        f_insert[bb.delta] = status;
        vfo.push_back(f_insert);
      }
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
      os << "|RDD: Reading Database|=" << time << "\n";
#endif
    } else {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
      os << "RDD: No database present\n";
#endif
    }
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
    if (SavingTrigger) {
      write_method(eFileMethod, the_method);
    }
  }
  void DirectAppendDatabase(vectface &&vf) {
    bb.clear();
    size_t n_orbit = vf.size();
    size_t len_ff = 0;
    size_t len_fb = 0;
    if (SavingTrigger) {
      len_ff = (n_orbit * bb.delta + 7) / 8;
      len_fb = (n_orbit + 7) / 8;
    }
    std::vector<uint8_t> ListOrbit_ff(len_ff);
    std::vector<uint8_t> V_status(len_fb);
    size_t pos_ff = 0;
    for (size_t i_orbit = 0; i_orbit < n_orbit; i_orbit++) {
      Face f = vf[i_orbit];
      Face f_red(bb.delta);
      for (size_t u = 0; u < bb.delta; u++) {
        bool val = f[u];
        f_red[u] = val;
        if (SavingTrigger) {
          setbit_vector(ListOrbit_ff, pos_ff, val);
          pos_ff++;
        }
      }
      bool status = f[bb.delta];
      if (SavingTrigger) {
        setbit_vector(V_status, i_orbit, status);
      }
      std::pair<Face, Tint> eEnt = bb.foc.recConvert.ConvertFace(f_red);
      bb.InsertListOrbitEntry(f_red, i_orbit);
      bb.InsertEntryDatabase(eEnt, status, i_orbit);
    }
    if (SavingTrigger) {
      FileNumber fn(eFileNB, true);
      FileBool fb(eFileFB);
      FileFace ff(eFileFF, bb.delta);
      fn.setval(n_orbit);
      ff.direct_write(ListOrbit_ff);
      fb.direct_write(V_status);
    }
    print_status();
  }
  ~DatabaseOrbits() {
    /* TRICK 5: The destructor does NOT destroy the database! This is because it
       can be used in another call. Note that the returning of the list of orbit
       does destroy the database and this gives a small window in which bad
       stuff can happen.
     */
    if (SavingTrigger && NeedToFlush) {
      flush();
    }
#ifdef DEBUG_RECURSIVE_DUAL_DESC
    os << "RDD: Clean closing of the DatabaseOrbits\n";
#endif
  }
  void flush() const {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
    os << "RDD: Doing the flushing operation\n";
#endif
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
    MicrosecondTime time;
#endif
    FileNumber fn(eFileNB, true);
    FileBool fb(eFileFB);
    FileFace ff(eFileFF, bb.delta);
    ff.direct_write(bb.foc.ListOrbit);
    size_t nbOrbit = bb.foc.nbOrbit;
    fn.setval(nbOrbit);
    size_t len = (nbOrbit + 7) / 8;
    std::vector<uint8_t> V_status(len, 255);
    auto iter = bb.begin_index_undone();
    while (iter != bb.end_index_undone()) {
      size_t pos = *iter;
      setbit_vector(V_status, pos, false);
      iter++;
    }
    fb.direct_write(V_status);
#ifdef TIMINGS_RECURSIVE_DUAL_DESC
    os << "|RDD: flush|=" << time << "\n";
#endif
  }
  // FuncListOrbitIncidence() {
  FaceOrbitsizeTableContainer<Tint> GetListFaceOrbitsize() {
    NeedToFlush = false;
    if (SavingTrigger) {
      RemoveFileIfExist(eFileNB);
      RemoveFileIfExist(eFileFB);
      RemoveFileIfExist(eFileFF);
      RemoveFileIfExist(eFileEXT);
      RemoveFileIfExist(eFileGRP);
      RemoveFileIfExist(eFileMethod);
    }
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
