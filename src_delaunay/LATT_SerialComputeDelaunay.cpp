// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheorySafeInt.h"
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
  DelaunayTesselation<T, Tgroup> DT =
      get_delaunay_tessellation_serial<T, Tint, Tgroup>(
          data, CacheFile, max_runtime_second, os);
  //
  // Write the tesselation in the requested format.
  //
  std::string OutFormat = BlockSYSTEM.get_string("OutFormat");
  std::string OutFile = BlockSYSTEM.get_string("OutFile");
  auto f = [&](std::ostream &os_out) -> void {
    WriteDelaunayTesselation(OutFormat, os_out, GramMat, DT);
  };
  print_stderr_stdout_file(OutFile, f);
  //
  // Additional computations on the tesselation (the QUERIES block). Each one is
  // skipped when its file entry is "null".
  //
  std::string FileQuantization = BlockQUERIES.get_string("FileQuantization");
  if (FileQuantization != "null") {
    QuantizationResult<T> qres =
        ComputeQuantizationIntegral<T, Tint, Tgroup>(data, DT, os);
    std::ofstream os_out(FileQuantization);
    os_out << "return ";
    WriteQuantizationGAP(os_out, qres);
    os_out << ";\n";
  }
  std::string FileFreeVectors = BlockQUERIES.get_string("FileFreeVectors");
  if (FileFreeVectors != "null") {
    FreeVectorsResult<Tint> fres =
        compute_free_vectors<T, Tint, Tgroup>(GramMat, DT, os);
    std::ofstream os_out(FileFreeVectors);
    WriteFreeVectorsGAP(os_out, fres);
  }
  std::string FileRigidityDegree = BlockQUERIES.get_string("FileRigidityDegree");
  if (FileRigidityDegree != "null") {
    int rigidity =
        ComputeRigidityDegreeLattice<T, Tint, Tgroup>(GramMat, DT, os);
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
        ComputeRigidityDegreeLattice<T, Tint, Tgroup>(GramMat, DT, os);
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
    IsoDelaunayDomain<T, Tint, Tgroup> x_iso{DT, GramMat, SHV_T};
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
    std::vector<MyVector<Tint>> reps =
        get_k_short_orbit_vectors<T, Tint>(Qadj, autom, Korb, os);
    std::ofstream os_out(FileDeformationOrbits);
    os_out << "return rec(nbEvaluated:=" << reps.size() << ",\nListOrbit:=[";
    bool IsFirst = true;
    for (auto &v : reps) {
      MyVector<T> v_T = UniversalVectorConversion<T, Tint>(v);
      MyMatrix<T> H = v_T * v_T.transpose();
      DeformationDerivatives<T> der =
          compute_deformation_derivatives<T, Tint, Tgroup>(GramMat, H, os);
      T invariant = v_T.dot(Qinv * v_T);
      long orbit_size = orbit_elements(autom, v).size();
      if (!IsFirst) {
        os_out << ",\n";
      }
      IsFirst = false;
      os_out << "rec(v:=" << StringVectorGAP(v) << ", OrbitSize:=" << orbit_size
             << ", vTQinvV:=" << invariant
             << ", SecMoment0:=" << jet_derivative(der.secmoment, 0)
             << ", SecMoment2:=" << jet_derivative(der.secmoment, 2)
             << ", det0:=" << der.det.coeffs[0]
             << ", det1:=" << der.det.coeffs[1]
             << ", det2:=" << der.det.coeffs[2]
             << ", Gpp:=" << jet_derivative(der.G, 2) << ")";
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
  if (arithmetic == "safe") {
    using T = Rational<SafeInt64>;
    using Tint = SafeInt64;
    return process_B<T, Tint>(eFull);
  }
  std::cerr << "LATT_SerialComputeDelaunay: Failed to find a matching type for "
               "arithmetic="
            << arithmetic << "\n";
  std::cerr << "Available types: gmp, gmp_boost, multi_boost, safe\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
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
