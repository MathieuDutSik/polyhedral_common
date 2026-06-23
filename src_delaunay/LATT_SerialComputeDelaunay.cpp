// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
#include "LatticeDelaunay.h"
#include "QuantizationIntegral.h"
#include "FreeVectors.h"
#include "Permutation.h"
#include "Group.h"
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
