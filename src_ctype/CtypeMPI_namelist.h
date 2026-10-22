// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_CTYPE_CTYPEMPI_NAMELIST_H_
#define SRC_CTYPE_CTYPEMPI_NAMELIST_H_

// clang-format off
#include <map>
#include <string>
// clang-format on

// The keys common to CTYP_MPI_Enumeration_c and CTYP_MPI_EnumerationAdjacencies
// for the dimension and the MPI exchanges between the processes
void Insert_CtypeMPI_ExchangeKeys(
    std::map<std::string, std::string> &ListIntValues_doc) {
  ListIntValues_doc["n"] = "The dimension n of the lattices.\n\
It must be the dimension of the C-types stored in the netCDF files";
  ListIntValues_doc["MaxNumberFlyingMessage"] = "Default: 100\n\
The maximum number of MPI messages being sent at the same time";
  ListIntValues_doc["MaxStoredUnsentMatrices"] = "Default: 1000\n\
When more entries than this wait to be sent, no new adjacency\n\
computation is started until the buffers are cleared";
  ListIntValues_doc["MaxRunTimeSecond"] = "Default: -1\n\
The maximum runtime in seconds after which no new adjacency computation\n\
is started. The state is in the netCDF files, so the run can be resumed.\n\
If negative then there is no maximum runtime";
  ListIntValues_doc["TimeForDeclaringItOver"] = "Default: 120\n\
The number of seconds without any operation after which a process\n\
declares itself finished to the other processes";
  ListIntValues_doc["MpiBufferSize"] = "Default: 20\n\
The maximum number of entries sent in a single MPI message";
}

// clang-format off
#endif  // SRC_CTYPE_CTYPEMPI_NAMELIST_H_
// clang-format on
