// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
#include "QuantizationDeformation.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

// Given a positive definite form Q, enumerate the orbit representatives of the
// integer vectors v with |v_i| <= bound under the symmetry group of Q
// (v ~ w iff w = +/- U^T v for U in Aut(Q), the equivalence that makes the
// deformations Q + t v v^T isometric). For the first K orbits, sorted by the
// invariant v^T Q^{-1} v, compute the second derivative of the normalized
// quantizer G(Q + t v v^T) at t = 0. Different orbits should give different
// values.
template <typename T, typename Tint, typename Tgroup>
void process(std::string const &Qfile, int bound, int K, std::ostream &os) {
  MyMatrix<T> Q = ReadMatrixFile<T>(Qfile);
  int n = Q.rows();
  MyMatrix<T> Qinv = Inverse(Q);
  std::vector<MyMatrix<Tint>> autom =
      ArithmeticAutomorphismGroup<T, Tint, Tgroup>(Q, os);
  os << "QORB: n=" << n << " |Aut(Q) generators|=" << autom.size() << "\n";
  //
  // Enumerate integer vectors with |v_i| <= bound (antipodally reduced).
  //
  std::vector<MyVector<Tint>> candidates;
  std::unordered_set<MyVector<Tint>> seen_cand;
  MyVector<Tint> v = ZeroVector<Tint>(n);
  int range = 2 * bound + 1;
  long total = 1;
  for (int i = 0; i < n; i++) {
    total *= range;
  }
  for (long code = 0; code < total; code++) {
    long c = code;
    bool is_zero = true;
    for (int i = 0; i < n; i++) {
      int digit = c % range;
      c /= range;
      v(i) = digit - bound;
      if (v(i) != 0) {
        is_zero = false;
      }
    }
    if (is_zero) {
      continue;
    }
    MyVector<Tint> cv = sign_canonicalize_vector(v);
    if (seen_cand.insert(cv).second) {
      candidates.push_back(cv);
    }
  }
  //
  // Classify the candidates into orbits.
  //
  std::unordered_set<MyVector<Tint>> assigned;
  struct OrbitRec {
    MyVector<Tint> rep;
    long orbit_size;
    T invariant;
  };
  std::vector<OrbitRec> orbits;
  for (auto &cand : candidates) {
    if (assigned.count(cand) > 0) {
      continue;
    }
    std::unordered_set<MyVector<Tint>> orb =
        orbit_vector_deformation(autom, cand);
    for (auto &w : orb) {
      assigned.insert(w);
    }
    MyVector<T> cand_T = UniversalVectorConversion<T, Tint>(cand);
    T invariant = cand_T.dot(Qinv * cand_T);
    orbits.push_back(OrbitRec{cand, static_cast<long>(orb.size()), invariant});
  }
  std::sort(orbits.begin(), orbits.end(),
            [](OrbitRec const &a, OrbitRec const &b) -> bool {
              return a.invariant < b.invariant;
            });
  os << "QORB: number of vector orbits (|v_i|<=" << bound
     << ")=" << orbits.size() << "\n";
  //
  // Compute the second derivative for the first K orbits.
  //
  int nb = std::min<int>(K, orbits.size());
  os << "QORB: computing G''(0) for the first " << nb << " orbits\n";
  for (int i = 0; i < nb; i++) {
    OrbitRec const &orb = orbits[i];
    MyVector<T> v_T = UniversalVectorConversion<T, Tint>(orb.rep);
    MyMatrix<T> H = v_T * v_T.transpose();
    DeformationDerivatives<T> der =
        compute_deformation_derivatives<T, Tint, Tgroup>(Q, H, os);
    os << "QORB: orbit " << i << " v=" << StringVectorGAP(orb.rep)
       << " |orbit|=" << orb.orbit_size << " vTQinvV=" << orb.invariant << "\n";
    os << "QORB:   SecMoment''(0)=" << der.S2 << " G''(0)=" << der.G2 << "\n";
  }
}

template <typename T, typename Tint>
void process_B(std::string const &Qfile, int bound, int K) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tint_grp = mpz_class;
  using Tgroup = permutalib::Group<Telt, Tint_grp>;
  return process<T, Tint, Tgroup>(Qfile, bound, K, std::cerr);
}

void process_C(std::string const &arith, std::string const &Qfile, int bound,
               int K) {
  if (arith == "gmp") {
    using T = mpq_class;
    using Tint = mpz_class;
    return process_B<T, Tint>(Qfile, bound, K);
  }
  if (arith == "safe") {
    using T = Rational<SafeInt64>;
    using Tint = SafeInt64;
    return process_B<T, Tint>(Qfile, bound, K);
  }
  std::cerr << "LATT_QuantizationDeformationOrbits: no match for arith=" << arith
            << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr
          << "LATT_QuantizationDeformationOrbits [arith] [Qfile] [bound] [K]\n";
      std::cerr << "     arith   gmp, safe\n";
      std::cerr << "     Qfile   the positive definite form Q\n";
      std::cerr << "     bound   enumerate integer vectors with |v_i| <= bound\n";
      std::cerr << "     K       number of orbits (by invariant) to evaluate\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string Qfile = argv[2];
    int bound = ParseScalar<int>(argv[3]);
    int K = ParseScalar<int>(argv[4]);
    process_C(arith, Qfile, bound, K);
    std::cerr << "Normal termination of LATT_QuantizationDeformationOrbits\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_QuantizationDeformationOrbits\n";
    exit(e.eVal);
  }
  runtime(time);
}
