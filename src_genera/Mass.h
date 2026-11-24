// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_MASS_H_
#define SRC_GENERA_MASS_H_

// clang-format off
#include "MassIsometry.h"
#include "MassGroupHomCount.h"
#include "MassRootFree.h"
#include "MassHermitian.h"
// clang-format on

/*
  The masses, the header to include. The layers, each including the one
  above it:

  - the masses of the components, one header per kind of centralizer:
    MassGenus (+-1, End_G = Q), MassCyclotomicHermitian (over Q(zeta_d)),
    MassAbelianHermitian (a CM abelian center), MassAbelianQuadraticLocal,
    MassQuadraticHeckeL and MassAbelianQuadratic (a totally real abelian
    center), MassQuaternionic, MassAlbertTypeTwo, MassUnitaryInnerForm;
    and MassHermitian (Hermitian lattices over an imaginary quadratic
    field) and MassRootFree (the masses by root system);
  - MassComponents: the decomposition of L under an isometry or an
    abelian group;
  - MassGlue: the discriminant modules of the components and the orbits
    of the glue;
  - MassMultiplicityLocalFactor: the local factors at the primes dividing
    |G|;
  - MassGroupComponents: the decomposition under a non-abelian group;
  - MassIsometry: the mass of (L, G);
  - MassGroupHomCount: the total mass of the genera of pairs.
 */

// clang-format off
#endif  // SRC_GENERA_MASS_H_
// clang-format on
