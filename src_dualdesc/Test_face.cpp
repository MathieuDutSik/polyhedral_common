// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "basic_datafile.h"
// clang-format on

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

void test_specific_size(int const &m, int const &n) {
  std::vector<Face> ListFace;
  for (int i = 0; i < n; i++) {
    Face f(m);
    for (int i = 0; i < m; i++) {
      bool rnd = random() % 2;
      f[i] = rnd;
    }
    ListFace.push_back(f);
    std::cerr << "i=" << i << " f=";
    for (int i = 0; i < m; i++) {
      std::cerr << f[i];
    }
    std::cerr << "\n";
  }
  //
  std::string TestFile = "/tmp/testbool_n" + std::to_string(n) + "_m" +
                         std::to_string(m) + "_" + random_string(20);
  RemoveFileIfExist(TestFile);
  //  std::cerr << "TestFile=" << TestFile << "\n";
  {
    FileFace ff(TestFile, m);
    for (int i = 0; i < n; i++)
      ff.setface(i, ListFace[i]);
  }
  //
  FileFace ff(TestFile, m, n);
  for (int i = 0; i < n; i++) {
    Face f = ff.getface(i);
    if (f != ListFace[i]) {
      std::cerr << "Inconsistency at i=" << i << "\n";
      throw TerminalException{1};
    }
  }
}

std::string get_string(Face f) {
  std::string estr;
  for (size_t i = 0; i < f.size(); i++)
    estr += std::to_string(f[i]);
  return estr;
}

void test_specific_size_randaccess(int const &m, int const &n) {
  std::vector<Face> ListFace(n);
  std::vector<int> Status(n, 0);
  //
  std::string TestFile = "/tmp/testbool_n" + std::to_string(n) + "_m" +
                         std::to_string(m) + "_" + random_string(20);
  RemoveFileIfExist(TestFile);
  FileFace ff(TestFile, m);
  for (size_t iter = 0; iter < 1000; iter++) {
    size_t pos = random() % n;
    Face f(m);
    for (int i = 0; i < m; i++)
      f[i] = random() % 2;
    ListFace[pos] = f;
    Status[pos] = 1;
    std::cerr << "iter=" << iter << " pos=" << pos << " f=" << get_string(f)
              << "\n";
    ff.setface(pos, f);
    //
    for (int i = 0; i < n; i++) {
      if (Status[i] == 1) {
        Face f = ff.getface(i);
        if (f != ListFace[i]) {
          std::cerr << "Inconsistency at i=" << i << " iter=" << iter << "\n";
          throw TerminalException{1};
        }
      }
    }
  }
}

int main() {
  for (size_t i = 0; i < 100; i++) {
    int n = 10 + random() % 20;
    int m = 20 + random() % 20;
    std::cerr << "n=" << n << " m=" << m << "\n";
    //    test_specific_size(m, n);
    test_specific_size_randaccess(m, n);
  }
  std::cerr << "Normal end of the program\n";
}
