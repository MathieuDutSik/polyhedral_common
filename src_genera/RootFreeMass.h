// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_ROOTFREEMASS_H_
#define SRC_GENERA_ROOTFREEMASS_H_

// clang-format off
#include "GenusMass.h"
#include "LocalDensityOverlattice.h"
#include "RootEmbeddings.h"
#include <algorithm>
#include <utility>
#include <vector>
// clang-format on

#ifdef SANITY_CHECK
#define SANITY_CHECK_ROOT_FREE_MASS
#endif

/*
  The masses x_T = sum_{L in G, rho(L) = T} 1 / |Aut L| of the lattices of a
  genus G by the type T of their root system, x_0 being the mass of the
  root-free lattices, from local data only.

  For every ADE type R, M_R = sum_{L in G} r(R, L) / |Aut L| =
  sum_T r(R, T) x_T, and the matrix r(R, T) is triangular for "R embeds in
  T" with diagonal r(R, R) = |Aut R|. Decreasing number of roots is a
  linear extension of that order (a proper root sublattice has fewer
  roots; rank alone is not, A_1^4 lies in D_4), so

      x_R = (M_R - sum_{T above R} r(R, T) x_T) / |Aut R|

  in that order, and x_0 = mass(G) - sum_{T != 0} x_T. The M_R come from
  Siegel's formula below the full rank (MassRatioFromDensitiesFull times
  GenusMass) and from the count of overlattices at the full rank
  (FullRankRepresentationMass). The system r(R, T) does not depend on the
  genus and is built once.
 */

template <typename Tint> struct RootMassSystem {
  int max_rank;
  // the types, by decreasing number of roots, the zero type last
  std::vector<AdeType> types;
  // above[i]: the j < i with r(types[i], types[j]) != 0, and that number
  std::vector<std::vector<std::pair<size_t, Tint>>> above;
  std::vector<Tint> diagonal;
};

template <typename Tint>
RootMassSystem<Tint> BuildRootMassSystem(int max_rank,
                                         RootEmbeddingCache<Tint> &cache) {
  RootMassSystem<Tint> sys;
  sys.max_rank = max_rank;
  sys.types = AdeTypesUpToRank(max_rank);
  sys.types.push_back(AdeType());
  std::stable_sort(sys.types.begin(), sys.types.end(),
                   [](AdeType const &a, AdeType const &b) -> bool {
                     return AdeNumberRoots(a) > AdeNumberRoots(b);
                   });
  size_t n_type = sys.types.size();
  sys.above.resize(n_type);
  for (size_t i = 0; i < n_type; i++) {
    sys.diagonal.push_back(
        RootEmbeddingNumber<Tint>(sys.types[i], sys.types[i], cache));
    for (size_t j = 0; j < i; j++) {
      Tint r = RootEmbeddingNumber<Tint>(sys.types[i], sys.types[j], cache);
      if (r != 0) {
        sys.above[i].push_back({j, r});
      }
    }
  }
#ifdef SANITY_CHECK_ROOT_FREE_MASS
  for (size_t i = 0; i < n_type; i++) {
    for (size_t j = i + 1; j < n_type; j++) {
      Tint r = RootEmbeddingNumber<Tint>(sys.types[i], sys.types[j], cache);
      if (r != 0) {
        std::cerr << "RFM: " << AdeTypeName(sys.types[i]) << " embeds in "
                  << AdeTypeName(sys.types[j])
                  << " which has no more roots\n";
        throw TerminalException{1};
      }
    }
  }
#endif
  return sys;
}

template <typename T> struct RootFreeMassResult {
  T mass;
  T x0;
  // x_T for the types of the system, in its order
  std::vector<T> x;
};

template <typename T> struct RootFreeMassCaches {
  LocalDensityCache<T> density;
  FullRankOverlatticeCache<T> full_rank;
};

template <typename T, typename Tint>
RootFreeMassResult<T> ComputeRootFreeMass(MyMatrix<T> const &Tgram,
                                          RootMassSystem<Tint> const &sys,
                                          RootFreeMassCaches<T> &caches) {
  int n = Tgram.rows();
  if (n > sys.max_rank) {
    std::cerr << "RFM: the system is built for rank " << sys.max_rank
              << " and the genus has rank " << n << "\n";
    throw TerminalException{1};
  }
  RootFreeMassResult<T> res;
  res.mass = GenusMass<T, Tint>(Tgram);
  GenusDensityDataFull<T> data =
      ComputeGenusDensityDataFull<T, Tint>(Tgram, caches.density.siegel);
  size_t n_type = sys.types.size();
  res.x.resize(n_type);
  for (size_t i = 0; i < n_type; i++) {
    AdeType const &R = sys.types[i];
    int m = AdeRank(R);
    T M_R;
    if (m > n) {
      M_R = T(0);
    } else if (m == 0) {
      M_R = res.mass;
    } else if (m == n) {
      M_R = FullRankRepresentationMass<T, Tint>(AdeGramMatrix<T>(R), Tgram,
                                                caches.full_rank);
    } else {
      M_R = res.mass * MassRatioFromDensitiesFull<T, Tint>(
                           AdeGramMatrix<T>(R), data, caches.density);
    }
    for (auto const &[j, r] : sys.above[i]) {
      M_R -= UniversalScalarConversion<T, Tint>(r) * res.x[j];
    }
    res.x[i] = M_R / UniversalScalarConversion<T, Tint>(sys.diagonal[i]);
#ifdef SANITY_CHECK_ROOT_FREE_MASS
    if (res.x[i] < 0) {
      std::cerr << "RFM: negative mass " << res.x[i] << " for the type "
                << AdeTypeName(R) << "\n";
      throw TerminalException{1};
    }
#endif
  }
  res.x0 = res.x[n_type - 1];
  return res;
}

// clang-format off
#endif  // SRC_GENERA_ROOTFREEMASS_H_
// clang-format on
