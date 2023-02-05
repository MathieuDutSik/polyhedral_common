// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_TEMP_MPI_H_
#define SRC_POLY_TEMP_MPI_H_

#include "GRP_GroupFct.h"
#include "MPI_functions.h"
#include "MPQ_Matrix.h"
#include <vector>

struct MPIworkingstructure {
  int Nproc;
  std::vector<int> StatusProc;
};

template <typename T> struct EquivalenceSetFunc {
  typedef typename equiv_info<T>::equiv_type Tequiv;
  typedef typename invariant_info<T>::invariant_type Tinv;
  std::function<EquivInfo<Tequiv>(T const &, T const &)> fEquiv;
  std::function<Tinv(T const &)> fInv;
  std::function<int(Tinv const &)> fHash;
};

template <typename T>
void InfiniteLoopWrite(std::vector<int> const &ListNode) {}

template <typename T> void BankingLoop(std::vector<int> const &ListNodes) {
  while (1) {
  }
}

template <>
void MPI_SEND_vector<mpq_class>(vector<mpq_class> a, int dest, int tag,
                                MPI_Comm comm) {
  MPI_SEND_mpq_vector(a, dest, tag, comm);
}

template <>
void MPI_RECV_vector(vector<mpq_class> &a, int src, int tag, MPI_Comm comm) {
  void MPI_RECV_mpq_vector(a, src, tag, comm);
}

template <>
void MPI_SEND_vector(vector<double> &a, int dest, int tag, MPI_Comm comm) {
  std::vector<int> eVectSendC(1, a.size());
  (int)MPI_Send(eVectSendC.data(), 1, MPI_INT, dest, tag, comm);
  (int)MPI_Send(a.data(), eSize, MPI_DOUBLE_PRECISION, dest, 1002, comm);
}

template <>
void MPI_RECV_vector<double>(vector<double> &a, int src, int tag,
                             MPI_Comm comm) {
  int *eVectRecvC;
  double *eVectRecvDouble;
  if ((eVectRecvC = (int *)malloc(1 * sizeof(int))) == 0) {
    throw TerminalException{1};
  }
  MPI_Recv(eVectRecvC, 1, MPI_INT, src, tag, comm, &status);
  eSize = eVectRecvC[0];
  free(eVectRecvC);

  if ((eVectRecvDouble = (double *)malloc(eSize * sizeof(double))) == 0) {
    throw TerminalException{1};
  }
  MPI_Recv(eVectRecvDouble, eSize, MPI_DOUBLE_PRECISION, src, 1002, comm,
           &status);
  a.clear();
  for (i = 0; i < eSize; i++)
    a.push_back(eVectRecvDouble[i]);
  free(eVectRecvDouble);
}

template <typename T>
void MPI_SEND_MyMatrix(MyMatrix<T> *eMat, int dest, int tag, MPI_Comm comm) {
  vector<T> eVect;
  int *eVectSendC;
  int nbRow, nbCol, ret, nbEnt, i;
  if ((eVectSendC = (int *)malloc(2 * sizeof(int))) == 0) {
    throw TerminalException{1};
  }
  nbRow = eMat->nbRow;
  nbCol = eMat->nbCol;
  eVectSendC[0] = nbRow;
  eVectSendC[1] = nbCol;
  ret = MPI_Send(eVectSendC, 2, MPI_INT, dest, 1292, comm);
  free(eVectSendC);

  nbEnt = nbRow * nbCol;
  for (i = 0; i < nbEnt; i++)
    eVect.push_back(eMat->ListElt[i]);
  MPI_SEND_vector(eVect, dest, tag, comm);
}

template <typename T>
void MPI_RECV_MyMatrix(MyMatrix<T> *eMat, int src, int tag, MPI_Comm comm) {
  VectGMP eVect;
  int *eVectRecvC;
  MPI_Status status;
  int nbRow, nbCol, nbEnt, i;
  if ((eVectRecvC = (int *)malloc(2 * sizeof(int))) == 0) {
    throw TerminalException{1};
  }
  MPI_Recv(eVectRecvC, 2, MPI_INT, src, 1292, comm, &status);
  nbRow = eVectRecvC[0];
  nbCol = eVectRecvC[1];
  free(eVectRecvC);

  MATRIX_Allocate(eMat, nbRow, nbCol);
  MPI_RECV_vector(eVect, src, tag, comm);
  nbEnt = nbRow * nbCol;
  for (i = 0; i < nbEnt; i++)
    eMat->ListElt[i] = eVect[i];
}

// clang-format off
#endif  // SRC_POLY_TEMP_MPI_H_
// clang-format on
