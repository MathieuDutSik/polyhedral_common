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




int main()
{
  std::string TestFile = "/tmp/testbool_" + random_string(20);
  std::cerr << "TestFile = " << TestFile << "\n";
  //
  int m = 24;
  Face f(m);
  for (int i=0; i<m; i++) {
    bool rnd = rand() % 2;
    std::cerr << "i=" << i << " rnd=" << rnd << "\n";
    f[i] = rnd;
  }
  //
  CreateFile(TestFile, f);
  Face f_read = ReadFile(TestFile, m);
  if (f != f_read) {
    std::cerr << "Reading OR writing failed\n";
  } else {
    std::cerr << "Correct reading/writing\n";
  }
}
