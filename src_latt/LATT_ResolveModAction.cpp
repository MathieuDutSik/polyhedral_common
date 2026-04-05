// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheorySafeInt.h"
#include "Zp_action.h"
// clang-format on

template <typename T>
void compute_orbit_partition(std::string const &list_matrix_file,
                             std::string const& mod_val_string,
                             std::string const &OutFormat,
                             std::ostream &os_out) {
  std::vector<MyMatrix<T>> l_gens = ReadListMatrixFile<T>(list_matrix_file);
  int dim = l_gens[0].rows();
  T mod_val = ParseScalar<T>(mod_val_string);
  ResultModEnumeration rme = get_partition(dim, l_gens, mod_val, std::cerr);
  if (OutFormat == "GAP") {
    os_out << "return [";
    for (size_t i = 0; i < rme.orbit_sizes.size(); i++) {
      if (i > 0)
        os_out << ",\n";
      os_out << rme.orbit_sizes[i];
    }
    os_out << "];\n";
    return;
  }
  std::cerr << "Failed to find a matching entry for OutFormat\n";
  throw TerminalException{1};
}

void process(std::string const &arith,
             std::string const &list_matrix_file,
             std::string const &mod_val_string,
             std::string const &OutFormat, std::ostream &os_out) {
  if (arith == "mpz_class") {
    using T = mpz_class;
    return compute_orbit_partition<T>(list_matrix_file, mod_val_string, OutFormat, os_out);
  }
  std::cerr << "Failed to find a matching entry for arith\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 4 && argc != 6) {
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_ResolveModAction [arith] [list_matrix_file] [mod_val] "
                   "[OutFormat] [OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_ResolveModAction [arith] [list_matrix_file] [mod_val]\n";
      std::cerr << "\n";
      std::cerr << "    where\n";
      std::cerr << "arith: mpz_class\n";
      std::cerr << "list_matrix_file: The input matrix file\n";
      std::cerr << "mod_val: The modulo considered\n";
      std::cerr << "OutFormat: GAP\n";
      std::cerr << "OutFile: Output file or stderr/stdout (optional, "
                   "default: stderr)\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string list_matrix_file = argv[2];
    std::string mod_val_string = argv[3];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      OutFile = argv[5];
    }
    auto f=[&](std::ostream& os_out) -> void {
      return process(arith, list_matrix_file, mod_val_string, OutFormat, os_out);
    };
    print_stderr_stdout_file(OutFile, f);
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something wrong happened in the computation\n";
    exit(e.eVal);
  }
  runtime(time);
}
