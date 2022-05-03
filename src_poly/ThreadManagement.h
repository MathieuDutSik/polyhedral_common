#ifndef SRC_POLY_THREADMANAGEMENT_H_
#define SRC_POLY_THREADMANAGEMENT_H_

#include "BasicThreadInclude.h"
#include "facet.h"
//#include "tbb/tbb.h"
#include "tbb/concurrent_vector.h"
#include <string>

/*
  The complicacy of this function is that we have to deal with the
  fact that some of the thread can be stuck and so numberSpannedThread
  can be larger than MaxNumberThread.
  The number of running thread are changed by decNRT / incNRT */
struct MainProcessor {
public:
  // no accidental construction, i.e. temporaries and the like
  MainProcessor() = delete;

  // no copy
  MainProcessor(const MainProcessor &) = delete;

  // no assign
  MainProcessor &operator=(const MainProcessor &) = delete;

  // no move
  MainProcessor(MainProcessor &&) = delete;

  // lock is of course needed but established somewhere else
  void SpanThreadNoLock() {
    std::ostringstream convert;
    convert << numberSpannedThread;
    std::string str3 = "log" + convert.str();
    //    os.imbue(LogFacet);
    std::ofstream *os_ptr = nullptr;
    int pos = LogErr.size();
    LogErr.push_back(os_ptr);
    LogErr[pos] = new std::ofstream(str3);
    *(LogErr[pos]) << "First message in the logs\n";
    *(LogErr[pos]) << std::unitbuf;
    ListStatus.push_back(0);
    numberSpannedThread++;
  }
  MainProcessor(int inpMaxNumThread) {
    numberSpannedThread = 0;
    MaxNumThread = inpMaxNumThread;
    numberRunningThread = 0;
  }
  ~MainProcessor() {
    for (int iThr = 0; iThr < numberSpannedThread; iThr++) {
      *(LogErr[iThr]) << "Last message in the logs\n";
      LogErr[iThr]->close();
      delete LogErr[iThr];
    }
  }
  int MPU_NumberFree() { return MaxNumThread - numberRunningThread; }
  void incNRT(int const &TheId) {
    std::lock_guard<std::mutex> lk(mulNRT);
    numberRunningThread++;
    *(LogErr[TheId]) << "Call to incNRT: Now numberRunningThread="
                     << numberRunningThread << "\n";
    *(LogErr[TheId]) << "MaxNumThread=" << MaxNumThread
                     << " numberSpannedThread=" << numberSpannedThread << "\n";
  }
  void decNRT(int const &TheId) {
    std::lock_guard<std::mutex> lk(mulNRT);
    numberRunningThread--;
    *(LogErr[TheId]) << "Call to decNRT: Now numberRunningThread="
                     << numberRunningThread << "\n";
    *(LogErr[TheId]) << "MaxNumThread=" << MaxNumThread
                     << " numberSpannedThread=" << numberSpannedThread << "\n";
  }
  int MPU_GetId() {
    std::lock_guard<std::mutex> lk(mul);
    for (int iProc = 0; iProc < numberSpannedThread; iProc++)
      if (ListStatus[iProc] == 0) {
        *(LogErr[iProc]) << "Reusing existing stream with iProc=" << iProc
                         << "\n";
        ListStatus[iProc] = 1;
        numberRunningThread++;
        return iProc;
      }
    if (numberRunningThread < MaxNumThread) {
      int iProc = numberSpannedThread;
      SpanThreadNoLock();
      *(LogErr[iProc]) << "Using a newly spanned stream with iProc=" << iProc
                       << "\n";
      ListStatus[iProc] = 1;
      numberRunningThread++;
      return iProc;
    }
    return -1;
  }
  void MPU_Terminate(int eId) {
    std::lock_guard<std::mutex> lk(mul);
    *(LogErr[eId]) << "Terminating eId=" << eId << "\n";
    ListStatus[eId] = 0;
    *(LogErr[eId]) << "numberSpannedThread=" << numberSpannedThread << "\n";
    numberRunningThread--;
  }
  std::ostream &GetO(int eId) { return *(LogErr[eId]); }

private:
  std::mutex mul;
  std::mutex mulNRT;
  int numberSpannedThread;
  int MaxNumThread;
  int numberRunningThread;
  tbb::concurrent_vector<std::ofstream *> LogErr;
  tbb::concurrent_vector<int> ListStatus;
};

struct LimitNumberProcessor {
public:
  // no accidental construction, i.e. temporaries and the like
  LimitNumberProcessor() = delete;

  // no copy
  LimitNumberProcessor(const LimitNumberProcessor &) = delete;

  // no assign
  LimitNumberProcessor &operator=(const LimitNumberProcessor &) = delete;

  // no move
  LimitNumberProcessor(LimitNumberProcessor &&) = delete;

  LimitNumberProcessor(int inpMaxNumThread) {
    numberThread = 0;
    MaxNumThread = inpMaxNumThread;
  }
  ~LimitNumberProcessor() = default;
  void MPU_Acquire() {
    std::unique_lock<std::mutex> lk(mul);
    data_cond.wait(lk, [&] { return numberThread < MaxNumThread - 1; });
    numberThread++;
  }
  void MPU_Terminate() {
    std::lock_guard<std::mutex> lk(mul);
    numberThread--;
  }

private:
  std::mutex mul;
  std::condition_variable data_cond;
  int numberThread;
  int MaxNumThread;
};

// clang-format off
#endif  // SRC_POLY_THREADMANAGEMENT_H_
// clang-format on
