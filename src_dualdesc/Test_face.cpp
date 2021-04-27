#include "POLY_netcdf_file.h"




void CreateFile(std::string const& eFile, Face f)
{
  FileBool fb(eFile);
  for (size_t i=0; i<f.size(); i++) {
    bool val = f[i];
    fb.setbit(i, val);
  }
}


Face ReadFile(std::string const& eFile, int m)
{
  FileBool fb(eFile, m);
  Face f(m);
  for (size_t i=0; i<f.size(); i++) {
    bool val = fb.getbit(i);
    f[i] = val;
  }
  return f;
}




void test_specific_size(int const& m, int const& n)
{
  std::vector<Face> ListFace;
  for (int i=0; i<n; i++) {
    Face f(m);
    for (int i=0; i<m; i++) {
      bool rnd = rand() % 2;
      f[i] = rnd;
    }
    ListFace.push_back(f);
    std::cerr << "i=" << i << " f=";
    for (int i=0; i<m; i++) {
      std::cerr << f[i];
    }
    std::cerr << "\n";
  }
  //
  std::string TestFile = "/tmp/testbool_n" + std::to_string(n) + "_m" + std::to_string(m) + "_" + random_string(20);
  RemoveFileIfExist(TestFile);
  //  std::cerr << "TestFile=" << TestFile << "\n";
  {
    FileFace ff(TestFile, m);
    for (int i=0; i<n; i++)
      ff.setface(i, ListFace[i]);
  }
  //
  FileFace ff(TestFile, m, n);
  for (int i=0; i<n; i++) {
    Face f = ff.getface(i);
    if (f != ListFace[i]) {
      std::cerr << "Inconsistency at i=" << i << "\n";
      throw TerminalException{1};
    }
  }

}











int main()
{
  for (int i=0; i<1000; i++) {
    std::timespec ts;
    std::timespec_get(&ts, TIME_UTC);
    unsigned val = ts.tv_nsec;
    std::cerr << "i=" << i << " val=" << val << "\n";
  }

  
  for (size_t i=0; i<100; i++) {
    int n = 10 + rand()% 20;
    int m = 20 + rand()% 20;
    std::cerr << "n=" << n << " m=" << m << "\n";
    test_specific_size(m, n);
  }
  std::cerr << "Normal end of the program\n";
}
