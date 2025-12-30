// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_PERFECTFORMSMPI_H_
#define SRC_PERFECT_PERFECTFORMSMPI_H_

// clang-format off
#include "LatticeStabEquiCan.h"
#include "Namelist.h"
#include "PerfectMPI_types.h"
#include "PerfectForm.h"
#include "rational.h"
#include "Positivity.h"
#include "POLY_lrslib.h"
#include <unordered_map>
#include "hash_functions.h"
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/optional.hpp>
// clang-format on

template <typename T, typename Tint>
std::vector<TypePerfectExch<Tint>>
GetAdjacentObjects(TypePerfectExch<Tint> const &eObjIn, std::ostream &os) {
  MyMatrix<T> eMat_T = UniversalMatrixConversion<T, Tint>(eObjIn.eMat);
  Tshortest<T, Tint> eRec = T_ShortestVector<T, Tint>(eMat_T, os);
  int n = eRec.SHV.cols();
  int nbShort = eRec.SHV.rows() / 2;
  int dimSymm = n * (n + 1) / 2;
  MyMatrix<Tint> SHVred(nbShort, n);
  for (int iShort = 0; iShort < nbShort; iShort++)
    for (int i = 0; i < n; i++)
      SHVred(iShort, i) = eRec.SHV(2 * iShort, i);
  MyMatrix<T> ConeClassical = GetNakedPerfectConeClassical<T, Tint>(SHVred);
  vectface ListIncd = lrs::DualDescription_incd(ConeClassical);
  MyVector<T> Wvect = GetSymmetricMatrixWeightVector<T>(n);
  std::vector<TypePerfectExch<Tint>> ListAdjMat;
  for (auto &eIncd : ListIncd) {
    MyVector<T> eFacet = FindFacetInequality(ConeClassical, eIncd);
    MyVector<T> Vexpand(dimSymm);
    for (int i = 0; i < dimSymm; i++)
      Vexpand(i) = eFacet(i) / Wvect(i);
    MyMatrix<T> eMatDir = VectorToSymmetricMatrix(Vexpand, n);
    std::pair<MyMatrix<T>, Tshortest<T, Tint>> ePairAdj =
        Flipping_Perfect<T, Tint>(eMat_T, eMatDir, os);
    int incd = ePairAdj.second.SHV.rows() / 2;
    MyMatrix<T> eMat2 =
        ComputeCanonicalForm<T, Tint>(ePairAdj.first, std::cerr).Mat;
    MyMatrix<T> eMat3 = RemoveFractionMatrix(eMat2);
    MyMatrix<Tint> eMat4 = UniversalMatrixConversion<Tint, T>(eMat3);
    TypePerfectExch<Tint> RecMat{incd, eMat4};
    ListAdjMat.emplace_back(RecMat);
  }
  return ListAdjMat;
}

template <typename T, typename Tint>
void EnumeratePerfectCones_MPI(boost::mpi::communicator &comm,
                               FullNamelist const &eFull) {
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  int n = BlockDATA.get_int("n");
  int MaxNumberFlyingMessage = BlockDATA.get_int("MaxNumberFlyingMessage");
  int MaxIncidenceTreating = BlockDATA.get_int("MaxIncidenceTreating");
  int MaxStoredUnsentMatrices = BlockDATA.get_int("MaxStoredUnsentMatrices");

  int irank = comm.rank();
  int size = comm.size();
  std::string eFileO = "LOG_" + std::to_string(irank);
  std::ofstream os(eFileO);
  os << "PERF_MPI_EnumeratePerfectCones: Initial log entry for rank " << irank
     << std::endl;

  // Initialize the MPI communication structures
  static int tag_new_form = 37;
  std::vector<boost::mpi::request> ListRequest(MaxNumberFlyingMessage);
  std::vector<int> RequestStatus(MaxNumberFlyingMessage, 0);

  auto GetFreeIndex = [&]() -> int {
    for (int u = 0; u < MaxNumberFlyingMessage; u++) {
      if (RequestStatus[u] == 0)
        return u;
      boost::optional<boost::mpi::status> stat = ListRequest[u].test();
      if (stat) {
        if (stat->error() != 0) {
          std::cerr << "MPI communication error" << std::endl;
          throw TerminalException{1};
        }
        RequestStatus[u] = 0;
        return u;
      }
    }
    return -1;
  };

  // Data structures for managing perfect forms
  std::unordered_map<TypePerfectExch<Tint>, int> ListCasesDone;
  std::unordered_map<TypePerfectExch<Tint>, int> ListCasesNotDone;
  std::vector<PairExch<Tint>> ListMatrixUnsent;
  int idxMatrixCurrent = 0;

  auto fInsert = [&](PairExch<Tint> const &ePair) -> void {
    TypePerfectExch<Tint> ePerfect = ePair.ePerfect;
    auto it1 = ListCasesDone.find(ePerfect);
    if (it1 != ListCasesDone.end()) {
      os << "Already processed entry=" << ePair.eIndex << " END" << std::endl;
      return;
    }
    auto it2 = ListCasesNotDone.find(ePerfect);
    if (it2 != ListCasesNotDone.end()) {
      os << "Already in queue entry=" << ePair.eIndex << " END" << std::endl;
      return;
    }
    ListCasesNotDone[ePerfect] = idxMatrixCurrent;
    os << "Inserting new perfect form " << ePerfect
       << " idxMatrix=" << idxMatrixCurrent << " from " << ePair.eIndex
       << " END" << std::endl;
    idxMatrixCurrent++;
  };

  auto fSendMatrix = [&](PairExch<Tint> const &ePair, int const &u) -> void {
    int res = IntegerDiscriminantInvariant(ePair.ePerfect.eMat, size);
    ListRequest[u] = comm.isend(res, tag_new_form, ePair);
    RequestStatus[u] = 1;
  };

  auto ClearUnsentAsPossible = [&]() -> void {
    int pos = ListMatrixUnsent.size() - 1;
    while (true) {
      if (pos == -1)
        break;
      int idx = GetFreeIndex();
      if (idx == -1)
        break;
      fSendMatrix(ListMatrixUnsent[pos], idx);
      ListMatrixUnsent.pop_back();
      pos--;
    }
  };

  auto fInsertUnsent = [&](PairExch<Tint> const &ePair) -> void {
    int res = IntegerDiscriminantInvariant(ePair.ePerfect.eMat, size);
    if (res == irank) {
      fInsert(ePair);
    } else {
      ListMatrixUnsent.push_back(ePair);
      if (static_cast<int>(ListMatrixUnsent.size()) < MaxStoredUnsentMatrices) {
        ClearUnsentAsPossible();
      }
    }
  };

  auto GetNextUndone = [&]() -> boost::optional<TypePerfectExch<Tint>> {
    auto it = ListCasesNotDone.begin();
    if (it != ListCasesNotDone.end()) {
      return boost::optional<TypePerfectExch<Tint>>(it->first);
    }
    return {};
  };

  auto SetMatrixAsDone = [&](TypePerfectExch<Tint> const &TheMat) -> void {
    int idxMatrix = ListCasesNotDone.at(TheMat);
    ListCasesNotDone.erase(TheMat);
    ListCasesDone[TheMat] = idxMatrix;
  };

  // Initialize with identity matrix for dimension n
  if (irank == 0) {
    MyMatrix<Tint> IdentMat = IdentityMat<Tint>(n);
    int incd = 2 * n; // Number of shortest vectors for identity matrix
    TypePerfectExch<Tint> eRecMat{incd, IdentMat};
    TypeIndex eIndex{irank, idxMatrixCurrent, 0};
    PairExch<Tint> ePair{eRecMat, eIndex};
    fInsert(ePair);
  }

  os << "Starting main enumeration loop" << std::endl;

  // Main enumeration loop
  while (true) {
    // Check for incoming messages
    boost::optional<boost::mpi::status> prob = comm.iprobe();
    if (prob) {
      if (prob->tag() == tag_new_form) {
        PairExch<Tint> ePair;
        comm.recv(prob->source(), prob->tag(), ePair);
        fInsert(ePair);
      }
    } else {
      // Process next undone matrix if queue isn't full
      int nMatrixUnsent = ListMatrixUnsent.size();
      if (nMatrixUnsent < MaxStoredUnsentMatrices) {
        boost::optional<TypePerfectExch<Tint>> eReq = GetNextUndone();
        if (eReq) {
          SetMatrixAsDone(*eReq);

          os << "Processing perfect form with " << eReq->incd
             << " shortest vectors" << std::endl;

          // Compute adjacent perfect forms
          std::vector<TypePerfectExch<Tint>> ListAdjacentObject =
              GetAdjacentObjects<T, Tint>(*eReq, os);

          int nbAdjacent = ListAdjacentObject.size();
          os << "Found " << nbAdjacent << " adjacent perfect forms END"
             << std::endl;

          int iAdj = 0;
          for (auto &eObj : ListAdjacentObject) {
            if (eObj.incd <= MaxIncidenceTreating) {
              TypeIndex eIndex{irank, ListCasesDone[*eReq], iAdj};
              PairExch<Tint> ePair{eObj, eIndex};
              fInsertUnsent(ePair);
            }
            iAdj++;
          }
        }
      }
    }

    ClearUnsentAsPossible();

    // Check termination condition
    if (ListCasesNotDone.empty() && ListMatrixUnsent.empty()) {
      // Signal completion to other processes
      comm.barrier();
      break;
    }

    // Periodic progress report
    if (idxMatrixCurrent % 100 == 0) {
      std::cerr << "Rank " << irank << ": processed " << ListCasesDone.size()
                << " forms, " << ListCasesNotDone.size() << " pending"
                << std::endl;
    }
  }

  os << "Completed enumeration. Total forms found: " << ListCasesDone.size()
     << std::endl;
  std::cerr << "Rank " << irank << " completed with " << ListCasesDone.size()
            << " perfect forms" << std::endl;
}

// clang-format off
#endif  // SRC_PERFECT_PERFECTFORMSMPI_H_
// clang-format on
