// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NonLatticePacking.h"
#include <iostream>
// clang-format on

/*
  Search tool for m-periodic sphere packings in double precision, after
  Andreanov--Kallus (see NonLatticePacking.h and the accompanying
  NonLatticePacking.tex in Paper_non_lattice). Configuration files use the
  same format as the covering search: a header "n m", the n rows of Q, the
  m rows of the cosets (the first being zero).

  The best lattice packing densities, for reference (they are also the best
  m-periodic densities for n = 3, 4, 5 and m <= 2 by Andreanov--Kallus):
    n = 3 : fcc     phi = pi/sqrt(18)  = 0.7404804897
    n = 4 : D4      phi = pi^2/16      = 0.6168502751
    n = 5 : D5      phi = pi^2 sqrt(2)/30 = 0.4652576133
 */

using namespace non_lattice_packing;

int main(int argc, char *argv[]) {
  try {
    if (argc < 2) {
      std::cerr << "PERIODIC_NonLatticePacking [mode] ...\n";
      std::cerr << "\n";
      std::cerr << "modes:\n";
      std::cerr << "  evaluate [config]\n";
      std::cerr << "      packing density phi, lambda^2 and the number of "
                   "shortest pair vectors\n";
      std::cerr << "  descend [config] [out_config] [rounds]\n";
      std::cerr << "      local maximization of the packing density\n";
      std::cerr << "  certify [config]\n";
      std::cerr << "      first-order extremeness test (the chart form of "
                   "Andreanov-Kallus algebraic extremeness) and the count "
                   "of flat (floating) directions\n";
      std::cerr << "  multistart [n] [m] [count] [out_config] [seed] "
                   "[rounds]\n";
      std::cerr << "      repeated descents from random configurations, "
                   "keeping the best\n";
      return -1;
    }
    std::string mode = argv[1];
    if (mode == "evaluate" && argc == 3) {
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      PackingResult pr = PackingDensity(conf);
      std::cout << "n=" << conf.n << " m=" << conf.m << "\n";
      printf("lambda2=%.13f n_short=%d\n", pr.lambda2, pr.n_short);
      printf("phi=%.13f\n", pr.phi);
      return 0;
    }
    if (mode == "descend" && (argc == 4 || argc == 5)) {
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      int rounds = argc == 5 ? atoi(argv[4]) : 15;
      DescendResult res = DescendPacking(conf, rounds, std::cerr, true);
      if (!res.success) {
        std::cerr << "the descent could not evaluate any configuration\n";
        return 1;
      }
      printf("phi=%.13f\n", res.theta);
      WriteConfigFile(argv[3], res.conf);
      return 0;
    }
    if (mode == "certify" && argc == 3) {
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      PackingResult pr = PackingDensity(conf);
      PackingCertificate ct = PackingCertify(conf);
      printf("phi=%.13f lambda2=%.13f\n", pr.phi, pr.lambda2);
      printf("active=%d rank=%d/%d rate=%.6e extreme=%d floating=%d\n",
             ct.n_active, ct.grad_rank, ct.dim, ct.rate, ct.extreme ? 1 : 0,
             ct.n_floating);
      return 0;
    }
    if (mode == "multistart" && (argc >= 6 && argc <= 8)) {
      int n = atoi(argv[2]);
      int m = atoi(argv[3]);
      int count = atoi(argv[4]);
      unsigned seed = argc >= 7 ? unsigned(atol(argv[6])) : 1u;
      int rounds = argc == 8 ? atoi(argv[7]) : 12;
      std::mt19937_64 gen(seed);
      std::normal_distribution<double> gauss(0.0, 1.0);
      std::uniform_real_distribution<double> unif(0.0, 1.0);
      double best = -1;
      for (int st = 0; st < count; st++) {
        PeriodicConfig conf;
        conf.n = n;
        conf.m = m;
        Eigen::MatrixXd R = Eigen::MatrixXd::Identity(n, n);
        for (int i = 0; i < n; i++)
          for (int j = 0; j < n; j++)
            R(i, j) += 0.25 * gauss(gen);
        conf.Q = R.transpose() * R;
        conf.C.resize(m, n);
        conf.C.row(0).setZero();
        for (int t = 1; t < m; t++)
          for (int j = 0; j < n; j++)
            conf.C(t, j) = unif(gen);
        DescendResult res = DescendPacking(conf, rounds, std::cerr, false);
        if (res.success) {
          printf("start %d: phi=%.13f%s\n", st, res.theta,
                 res.theta > best ? "  *" : "");
          fflush(stdout);
          if (res.theta > best) {
            best = res.theta;
            WriteConfigFile(argv[5], res.conf);
          }
        } else {
          printf("start %d: failed\n", st);
          fflush(stdout);
        }
      }
      printf("best=%.13f\n", best);
      return 0;
    }
    std::cerr << "unrecognized arguments; run without arguments for usage\n";
    return -1;
  } catch (std::exception const &e) {
    std::cerr << "Error in PERIODIC_NonLatticePacking: " << e.what() << "\n";
    return 1;
  }
}
