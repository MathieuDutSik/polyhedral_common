// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_POLY_DATABASE_ORBITS_H_
#define SRC_DUALDESC_POLY_DATABASE_ORBITS_H_

// Holds the `DatabaseOrbits<TbasicBank>` class extracted from
// POLY_RecursiveDualDesc.h. The class wraps a `TbasicBank` (face-orbit
// store) with optional on-disk persistence (eFileEXT / eFileGRP /
// eFileNB / eFileFB / eFileFF / eFileMethod) so that long enumerations
// can be resumed after a crash. The storage-agnostic logic is shared with
// NoSaveDatabaseOrbits through the DatabaseOrbitsCommon base. Macros used
// internally (DEBUG_RECURSIVE_DUAL_DESC, TIMINGS_RECURSIVE_DUAL_DESC,
// TRACK_DATABASE, TRACK_RUN) are intentionally shared with
// POLY_RecursiveDualDesc.h -- when this file is included from there, all four
// toggles are already configured.

// clang-format off
#include "POLY_database_orbits_common.h"
#include "Basic_file.h"
#include "basic_datafile.h"
#include <string>
// clang-format on

template <typename TbasicBank>
struct DatabaseOrbits
    : public DatabaseOrbitsCommon<TbasicBank, DatabaseOrbits<TbasicBank>> {
public:
  using Base = DatabaseOrbitsCommon<TbasicBank, DatabaseOrbits<TbasicBank>>;
  using typename Base::T;
  using typename Base::Tint;
  using Base::bb;
  using Base::NeedToFlush;
  using Base::os;
  using Base::print_status;

private:
  std::string eFileEXT, eFileGRP, eFileNB, eFileFB, eFileFF, eFileMethod;
  bool SavingTrigger;

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
  void remove_database_files() const {
    RemoveFileIfExist(eFileNB);
    RemoveFileIfExist(eFileFB);
    RemoveFileIfExist(eFileFF);
    RemoveFileIfExist(eFileEXT);
    RemoveFileIfExist(eFileGRP);
    RemoveFileIfExist(eFileMethod);
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
      remove_database_files();
    }
    return false;
  }
  DatabaseOrbits() = delete;
  DatabaseOrbits(const DatabaseOrbits<TbasicBank> &) = delete;
  DatabaseOrbits(DatabaseOrbits<TbasicBank> &&) = delete;
  DatabaseOrbits &operator=(const DatabaseOrbits<TbasicBank> &) = delete;
  DatabaseOrbits(TbasicBank &bb, const std::string &MainPrefix,
                 const bool &_SavingTrigger,
                 const bool &_AdvancedTerminationCriterion, std::ostream &os)
      : Base(bb, _AdvancedTerminationCriterion, os),
        SavingTrigger(_SavingTrigger) {
#ifdef DEBUG_RECURSIVE_DUAL_DESC
    os << "RDD: MainPrefix=" << MainPrefix << "\n";
#endif
    eFileEXT = MainPrefix + ".ext";
    eFileGRP = MainPrefix + ".grp";
    eFileNB = MainPrefix + ".nb";
    eFileFB = MainPrefix + ".fb";
    eFileFF = MainPrefix + ".ff";
    eFileMethod = MainPrefix + ".method";
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
      remove_database_files();
    }
    return bb.GetListFaceOrbitsize();
  }
};

// clang-format off
#endif  // SRC_DUALDESC_POLY_DATABASE_ORBITS_H_
// clang-format on
