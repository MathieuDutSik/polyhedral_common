#include "POLY_netcdf_file.h"

void CreateFile(std::string const &eFile, Face f) {
  FileBool fb(eFile);
  for (size_t i = 0; i < f.size(); i++) {
    bool val = f[i];
    fb.setbit(i, val);
  }
}

Face ReadFile(std::string const &eFile, int m) {
  FileBool fb(eFile, m);
  Face f(m);
  for (size_t i = 0; i < f.size(); i++) {
    bool val = fb.getbit(i);
    f[i] = val;
  }
  return f;
}

void test_specific_size(int const &m) {
  std::string TestFile =
      "/tmp/testbool_m" + std::to_string(m) + "_" + random_string(20);
  std::cerr << "TestFile = " << TestFile << "\n";
  //
  Face f(m);
  for (int i = 0; i < m; i++)
    f[i] = random() % 2;
  //
  CreateFile(TestFile, f);
  Face f_read = ReadFile(TestFile, m);
  if (f != f_read) {
    std::cerr << "Inconsistency in f != f_read\n";
    throw TerminalException{1};
  }
}

void test_specific_size_randaccess(int const &m) {
  std::string TestFile =
      "/tmp/testbool_m" + std::to_string(m) + "_" + random_string(20);
  std::cerr << "TestFile = " << TestFile << "\n";
  FileBool fb(TestFile);
  //
  Face f(m);
  std::vector<int> Status(m, 0);
  for (size_t iter = 0; iter < 1000; iter++) {
    size_t pos = random() % m;
    bool val = random() % 2;
    f[pos] = val;
    Status[pos] = 1;
    fb.setbit(pos, val);
    //
    for (int i = 0; i < m; i++) {
      if (Status[i] == 1) {
        bool val1 = fb.getbit(i);
        bool val2 = f[i];
        if (val1 != val2) {
          std::cerr << "Inconsistency in val1 and val2\n";
          throw TerminalException{1};
        }
      }
    }
  }
}

int main() {
  //
  for (size_t iter = 0; iter < 100; iter++) {
    int m = 100 + random() % 200;
    test_specific_size(m);
    test_specific_size_randaccess(m);
  }
  std::cerr << "Normal termination of the program\n";
}
