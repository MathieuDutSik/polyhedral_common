// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "basic_datafile.h"
// clang-format on

int main(int argc, char *argv[]) {
  try {
    if (argc != 3) {
      std::cerr << "Read_MPI_sizes [Prefix] [n_proc]\n";
      return -1;
    }
    using T_unused = int;
    std::string prefix = argv[1];
    std::string n_proc_str = argv[2];
    size_t n_proc = ParseScalar<int>(n_proc_str);
    for (size_t i_proc=0; i_proc<n_proc; i_proc++) {
      std::string file_name = prefix + "_nproc" + std::to_string(n_proc) + "_rank" + std::to_string(i_proc) + ".nb";
      bool overwrite = false;
      FileData<T_unused> fdata(file_name, overwrite);
      std::vector<size_t> l_sizes = fdata.read_all_sizes();
      std::cerr << "i_proc=" << i_proc << " l_sizes =";
      for (auto & val : l_sizes) {
        std::cerr << " " << val;
      }
      std::cerr << "\n";
    }
    std::cerr << "Normal end of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something went wrong\n";
    exit(e.eVal);
  }
}
