// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#include "LatticeDelaunay.h"
#include "IsoDelaunayDomains.h"
#include "Rigidity.h"
#include "QuantizationIntegral.h"
#include "QuantizationDeformation.h"
#include "FreeVectors.h"
#include "Permutation.h"
#include "Group.h"
#include <boost/archive/text_oarchive.hpp>
// clang-format on

template <typename T> struct HessianPlusResult {
  int n;
  int N;               // dim Sym^n = n(n+1)/2
  bool found;          // a rank-one direction with strictly positive curvature
  bool spanned;        // the sampled rank-one directions span Sym^n
  int nbEval;          // moment-derivative evaluations used
  MyVector<T> witness; // the direction v with beta_vv > 0 (found=true)
  T value;             // beta_vv on the witness (found=true)
};

// Weaker, cheaper companion of compute_hessian_signature: decide whether the
// quantizer Hessian at Q has ANY positive direction, without assembling the whole
// Hessian or its signature. It scans rank-one directions v v^T (exactly as the
// full Hessian enumerates them -- shell by shell of increasing dual norm
// v^T Q^{-1} v, one moment-derivative evaluation per Aut(Q)-orbit, since the
// Hessian value along the direction,
//   beta_vv = <B, DM[B]> - (2/n) qv mv + (S/n) qv^2 + (S/n^2) qv^2,
// with B = v v^T, qv = v^T Q^{-1} v, mv = v^T M0 v, is constant on an orbit) and
// returns as soon as one is strictly positive. A '+' is a sound, EXACT witness
// that the Hessian is not negative semidefinite (Q is not a local maximum of the
// normalized quantizer). This is the intended use on the hard, low-symmetry forms
// where the full signature needs too many orbit evaluations: if a positive
// direction shows up early we stop after a handful. If the sampled directions span
// Sym^n with none positive (spanned=true) no rank-one '+' exists; the full Hessian
// could still have a '+' on a higher-rank direction, so this test is deliberately
// weaker.
template <typename T, typename Tint, typename Tgroup>
HessianPlusResult<T> HasHessianPlus(MyMatrix<T> const &Q, std::ostream &os) {
  int n = Q.rows();
  int N = n * (n + 1) / 2;
  HessianPlusResult<T> res;
  res.n = n;
  res.N = N;
  res.found = false;
  res.spanned = false;
  res.nbEval = 0;
  res.value = T(0);
  MyMatrix<T> Qinv = Inverse(Q);
  QuantizationResult<T> q0 = quant_at_gram<T, Tint, Tgroup>(Q, os);
  MyMatrix<T> M0 = q0.SecMomentMat;
  T S = q0.SecMoment;
  std::vector<MyMatrix<Tint>> autom =
      ArithmeticAutomorphismGroup<T, Tint, Tgroup>(Q, os);
  T Tn(n);
  T detQ = DeterminantMat(Q);
  MyMatrix<T> Qadj = detQ * Qinv; // adjugate: integral, positive definite
  CVPSolver<T, Tint> solver(Qadj, os);
  T incr = GetSmallestIncrement(Qadj);
  T norm = incr;
  std::unordered_set<MyVector<Tint>> assigned;
  std::vector<MyVector<T>> span_rows; // rows spanning the reached subspace of Sym^n
  while (!res.spanned) {
    std::vector<MyVector<Tint>> ListVect = solver.fixed_norm_vectors(norm);
    for (auto &v0 : ListVect) {
      MyVector<Tint> start = SignCanonicalizeVector(v0);
      if (!assigned.insert(start).second) {
        continue;
      }
      // The whole Aut(Q)-orbit of "start" shares the same beta_vv; enumerate it
      // (marks the members and lets the spanning test use them for free).
      std::vector<MyVector<Tint>> orbit{start};
      size_t head = 0;
      while (head < orbit.size()) {
        MyVector<Tint> v = orbit[head];
        head++;
        for (auto &U : autom) {
          MyVector<Tint> w = U * v;
          MyVector<Tint> cw = SignCanonicalizeVector(w);
          if (assigned.insert(cw).second) {
            orbit.push_back(cw);
          }
        }
      }
      // One moment-derivative evaluation for the orbit; check the Hessian diagonal.
      MyVector<T> vRep = UniversalVectorConversion<T, Tint>(start);
      MyMatrix<T> B = vRep * vRep.transpose();
      MyMatrix<T> DM = compute_moment_derivative_jet<T, Tint, Tgroup>(Q, B, os);
      res.nbEval++;
      T qv = vRep.dot(Qinv * vRep);
      T mv = vRep.dot(M0 * vRep);
      T dm = vRep.dot(DM * vRep);
      T beta_vv = dm - (T(2) / Tn) * qv * mv + (S / Tn) * qv * qv +
                  (S / (Tn * Tn)) * qv * qv;
      if (beta_vv > T(0)) {
        res.found = true;
        res.witness = vRep;
        res.value = beta_vv;
        return res;
      }
      // No '+': this orbit only advances the Sym^n spanning bound.
      for (auto &m : orbit) {
        MyVector<T> vm = UniversalVectorConversion<T, Tint>(m);
        MyMatrix<T> Bm = vm * vm.transpose();
        MyVector<T> row = SymmetricMatrixToVector(Bm);
        int cur = static_cast<int>(span_rows.size());
        MyMatrix<T> Test(cur + 1, N);
        for (int r = 0; r < cur; r++) {
          Test.row(r) = span_rows[r].transpose();
        }
        Test.row(cur) = row.transpose();
        if (RankMat(Test) == cur + 1) {
          span_rows.push_back(row);
        }
        if (static_cast<int>(span_rows.size()) == N) {
          res.spanned = true;
          break;
        }
      }
      if (res.spanned) {
        break;
      }
    }
    norm += incr;
  }
  return res;
}

template <typename T>
void WriteHessianPlusGAP(std::ostream &os_out, HessianPlusResult<T> const &res) {
  os_out << "return rec(n:=" << res.n << ", dimSpace:=" << res.N
         << ", hasHessianPlus:=" << (res.found ? "true" : "false")
         << ", spanned:=" << (res.spanned ? "true" : "false")
         << ", nbEval:=" << res.nbEval;
  if (res.found) {
    os_out << ", witness:=" << StringVectorGAP(res.witness)
           << ", value:=" << res.value;
  }
  os_out << ");\n";
}

template <typename T, typename Tint, typename Tgroup>
void process_A(FullNamelist const &eFull, std::ostream &os) {
  using TintGroup = typename Tgroup::Tint;
  SingleBlock const &BlockSYSTEM = eFull.get_block("SYSTEM");
  SingleBlock const &BlockDATA = eFull.get_block("DATA");
  SingleBlock const &BlockQUERIES = eFull.get_block("QUERIES");
  //
  std::string GRAMfile = BlockDATA.get_string("GRAMfile");
  MyMatrix<T> GramMat = ReadMatrixFile<T>(GRAMfile);
  int dimEXT = GramMat.rows() + 1;
  //
  std::string FileDualDesc = BlockDATA.get_string("FileDualDescription");
  PolyHeuristicSerial<TintGroup> AllArr =
      Read_AllStandardHeuristicSerial_File<T, TintGroup>(FileDualDesc, dimEXT,
                                                         os);
  DataLattice<T, Tint, Tgroup> data =
      get_data_lattice<T, Tint, Tgroup>(eFull, AllArr, os);
  //
  // Compute (or load from cache) the Delaunay tesselation.
  //
  std::string CacheFile = BlockDATA.get_string("CacheFile");
  int max_runtime_second = BlockSYSTEM.get_int("max_runtime_second");
  DelaunayTesselation<Tint, Tgroup> DT =
      get_delaunay_tessellation_serial<T, Tint, Tgroup>(
          data, CacheFile, max_runtime_second, os);
  // Field-type view of the tessellation for the queries doing rational
  // geometry on it (quantization, isotropy, free vectors, rigidity).
  DelaunayTesselation<T, Tgroup> DT_T =
      ConvertTesselationScalar<T, Tint, Tgroup>(DT);
  //
  // Write the tesselation in the requested format.
  //
  std::string OutFormat = BlockSYSTEM.get_string("OutFormat");
  std::string OutFile = BlockSYSTEM.get_string("OutFile");
  auto f = [&](std::ostream &os_out) -> void {
    WriteDelaunayTesselation(OutFormat, os_out, GramMat, DT);
  };
  FILE_PrintStderrStdoutFile(OutFile, f);
  //
  // Additional computations on the tesselation (the QUERIES block). Each one is
  // skipped when its file entry is "null".
  //
  std::string FileQuantization = BlockQUERIES.get_string("FileQuantization");
  if (FileQuantization != "null") {
    QuantizationResult<T> qres =
        ComputeQuantizationIntegral<T, Tint, Tgroup>(data, DT_T, os);
    std::ofstream os_out(FileQuantization);
    os_out << "return ";
    WriteQuantizationGAP(os_out, qres);
    os_out << ";\n";
  }
  // Isotropy (extremal / "white" quantizer) test: GramMat times the Voronoi
  // cell second-moment matrix is a scalar multiple of the identity. The defect
  // matrix and the boolean are written to FileIsotropy.
  std::string FileIsotropy = BlockQUERIES.get_string("FileIsotropy");
  if (FileIsotropy != "null") {
    IsotropyResult<T> ires =
        ComputeIsotropy<T, Tint, Tgroup>(data, DT_T, GramMat, os);
    std::ofstream os_out(FileIsotropy);
    WriteIsotropyGAP(os_out, ires);
  }
  std::string FileFreeVectors = BlockQUERIES.get_string("FileFreeVectors");
  if (FileFreeVectors != "null") {
    FreeVectorsResult<Tint> fres =
        compute_free_vectors<T, Tint, Tgroup>(GramMat, DT_T, os);
    std::ofstream os_out(FileFreeVectors);
    WriteFreeVectorsGAP(os_out, fres);
  }
  std::string FileRigidityDegree = BlockQUERIES.get_string("FileRigidityDegree");
  if (FileRigidityDegree != "null") {
    int rigidity =
        ComputeRigidityDegreeLattice<T, Tint, Tgroup>(GramMat, DT_T, os);
    std::ofstream os_out(FileRigidityDegree);
    os_out << "return " << rigidity << ";\n";
  }
  std::string FileIsoDelaunayDomain =
      BlockQUERIES.get_string("FileIsoDelaunayDomain");
  if (FileIsoDelaunayDomain != "null") {
    // Only meaningful when the iso-Delaunay domain has full L-type
    // dimension, i.e. the L-space (= lineality of the L-type cone) has
    // dim n*(n+1)/2 — equivalently, every iso-Delaunay-domain inequality
    // is non-degenerate and the polyhedron is full-dim in canonical Sym^n.
    // Lower rigidity means the lattice sits on a wall and the iso-Delaunay
    // representation would be ill-defined.
    int n = GramMat.rows();
    int sym_dim = (n * (n + 1)) / 2;
    int rigidity =
        ComputeRigidityDegreeLattice<T, Tint, Tgroup>(GramMat, DT_T, os);
    if (rigidity != sym_dim) {
      std::cerr << "LATT_SerialComputeDelaunay: FileIsoDelaunayDomain "
                   "requires rigidity == n*(n+1)/2 = "
                << sym_dim << " but rigidity=" << rigidity << "\n";
      throw TerminalException{1};
    }
    // Compute SHV_T (full-rank invariant family), then bundle DT + GramMat
    // + SHV_T as an IsoDelaunayDomain and write it via boost::serialization
    // (text_oarchive — same format used by SerializeMatrix and friends).
    MyMatrix<Tint> SHV =
        ExtractInvariantVectorFamilyZbasis<T, Tint>(GramMat, os);
    MyMatrix<T> SHV_T = UniversalMatrixConversion<T, Tint>(SHV);
    // The rigidity check above guarantees the tesselation is generic over the
    // full canonical space, so its Voronoi inequalities are well defined.
    LinSpaceMatrix<T> LinSpa = ComputeCanonicalSpace<T>(n);
    std::vector<std::vector<Tint>> ListGramRing =
        GetListGramRing(LinSpa.ListLineMat);
    DelaunayTesselationIneq<Tint, Tgroup> DTI =
        BuildDelaunayTesselationIneq(DT, ListGramRing, os);
    IsoDelaunayDomain<T, Tint, Tgroup> x_iso{DTI, GramMat, SHV_T};
    std::ofstream ofs(FileIsoDelaunayDomain);
    boost::archive::text_oarchive oa(ofs);
    oa << x_iso;
  }
  // Second differential of G along Q + t H for a symmetric direction H read
  // from FileDeformation; the record is written to FileDeformation + ".output".
  std::string FileDeformation = BlockQUERIES.get_string("FileDeformation");
  if (FileDeformation != "null") {
    MyMatrix<T> H = ReadMatrixFile<T>(FileDeformation);
    DeformationDerivatives<T> der =
        compute_deformation_derivatives<T, Tint, Tgroup>(GramMat, H, os);
    std::ofstream os_out(FileDeformation + ".output");
    WriteDeformationGAP(os_out, der);
  }
  // Orbit scan of the rank-one directions v v^T: G''(0) for the first
  // DeformationNumberOrbit orbits of integer vectors under Aut(Q), ordered by the
  // invariant v^T Q^{-1} v. We get the orbit representatives shell by shell with
  // get_k_short_orbit_vectors, enumerating by the (integer) dual norm
  // v^T adj(Q) v = det(Q) v^T Q^{-1} v, then evaluate each one.
  std::string FileDeformationOrbits =
      BlockQUERIES.get_string("FileDeformationOrbits");
  if (FileDeformationOrbits != "null") {
    size_t Korb = BlockQUERIES.get_int("DeformationNumberOrbit");
    MyMatrix<T> Qinv = Inverse(GramMat);
    MyMatrix<T> Qadj = DeterminantMat(GramMat) * Qinv; // adj(Q): integral, PD
    std::vector<MyMatrix<Tint>> autom =
        ArithmeticAutomorphismGroup<T, Tint, Tgroup>(GramMat, os);
    std::vector<std::pair<MyVector<Tint>, size_t>> reps =
        get_k_short_orbit_vectors<T, Tint>(Qadj, autom, Korb, os);
    std::ofstream os_out(FileDeformationOrbits);
    os_out << "return rec(nbEvaluated:=" << reps.size() << ",\nListOrbit:=[";
    bool IsFirst = true;
    for (auto &[v, orbit_size] : reps) {
      MyVector<T> v_T = UniversalVectorConversion<T, Tint>(v);
      MyMatrix<T> H = v_T * v_T.transpose();
      DeformationDerivatives<T> der =
          compute_deformation_derivatives<T, Tint, Tgroup>(GramMat, H, os);
      T invariant = v_T.dot(Qinv * v_T);
      if (!IsFirst) {
        os_out << ",\n";
      }
      IsFirst = false;
      os_out << "rec(v:=" << StringVectorGAP(v) << ", OrbitSize:=" << orbit_size
             << ", vTQinvV:=" << invariant
             << ", SecMoment0:=" << jet_deriv(der.secmoment, 0)
             << ", SecMoment2:=" << jet_deriv(der.secmoment, 2)
             << ", det0:=" << jet_deriv(der.det, 0)
             << ", det1:=" << jet_deriv(der.det, 1)
             << ", det2:=" << jet_deriv(der.det, 2)
             << ", Gpp:=" << jet_deriv(der.G, 2) << ")";
    }
    os_out << "]);\n";
  }
  // Hessian of the normalized quantizer constant G at Q and its signature, via
  // the moment-derivative method on a rank-one basis v v^T built shell by shell
  // of increasing dual norm v^T Q^{-1} v.
  std::string FileHessian = BlockQUERIES.get_string("FileHessian");
  if (FileHessian != "null") {
    HessianResult<T> hres =
        compute_hessian_signature<T, Tint, Tgroup>(GramMat, os);
    std::ofstream os_out(FileHessian);
    WriteHessianGAP(os_out, hres);
  }
  // Weaker but cheaper: does the Hessian have any positive direction? Scans
  // rank-one directions and stops at the first one with positive curvature.
  std::string FileHessianPlus = BlockQUERIES.get_string("FileHessianPlus");
  if (FileHessianPlus != "null") {
    HessianPlusResult<T> hpres = HasHessianPlus<T, Tint, Tgroup>(GramMat, os);
    std::ofstream os_out(FileHessianPlus);
    WriteHessianPlusGAP(os_out, hpres);
  }
}

template <typename T, typename Tint> void process_B(FullNamelist const &eFull) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tint_grp = mpz_class;
  using Tgroup = permutalib::Group<Telt, Tint_grp>;
  return process_A<T, Tint, Tgroup>(eFull, std::cerr);
}

void process_C(FullNamelist const &eFull) {
  std::string arithmetic = GetNamelistStringEntry(eFull, "DATA", "arithmetic");
  if (arithmetic == "gmp") {
    using T = mpq_class;
    using Tint = mpz_class;
    return process_B<T, Tint>(eFull);
  }
#ifdef ENABLE_ALL_NUMERICAL_TYPES
  if (arithmetic == "gmp_boost") {
    using T = boost::multiprecision::mpq_rational;
    using Tint = boost::multiprecision::mpz_int;
    return process_B<T, Tint>(eFull);
  }
  if (arithmetic == "multi_boost") {
    using T = boost::multiprecision::cpp_rational;
    using Tint = boost::multiprecision::cpp_int;
    return process_B<T, Tint>(eFull);
  }
#endif
  std::cerr << "LATT_SerialComputeDelaunay: Failed to find a matching type for "
               "arithmetic="
            << arithmetic << "\n";
  std::cerr << "Available types: gmp, gmp_boost, multi_boost\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    FullNamelist eFull = NAMELIST_GetStandard_SERIAL_COMPUTE_DELAUNAY();
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_SerialComputeDelaunay [file.nml]\n";
      std::cerr << "With file.nml a namelist file\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      return -1;
    }
    std::string eFileName = argv[1];
    NAMELIST_ReadNamelistFile(eFileName, eFull);
    process_C(eFull);
    std::cerr << "Normal termination of LATT_SerialComputeDelaunay\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_SerialComputeDelaunay\n";
    exit(e.eVal);
  }
  runtime(time);
}
