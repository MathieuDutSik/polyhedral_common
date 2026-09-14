// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "JointCoveringDouble.h"
#include <iostream>
// clang-format on

/*
  Search tool for periodic sphere coverings in double precision: covering
  density evaluation, joint local descent on (Q, C), and multistart, for
  any dimension and any number of cosets. See JointCoveringDouble.h for
  the method; candidates found here are to be re-verified by the exact
  machinery.

  Configuration files: a header "n m", then the n rows of Q, then the m
  rows of the cosets (the first being zero).
 */

using namespace joint_covering_double;

int main(int argc, char *argv[]) {
  try {
    if (argc < 2) {
      std::cerr << "PERIODIC_JointCoveringSearch [mode] ...\n";
      std::cerr << "\n";
      std::cerr << "modes:\n";
      std::cerr << "  evaluate [config]\n";
      std::cerr << "      covering density, covering radius and cell-class "
                   "count of the configuration\n";
      std::cerr << "  descend [config] [out_config] [rounds]\n";
      std::cerr << "      joint local descent from the configuration; the "
                   "best found is written to out_config\n";
      std::cerr << "  multistart [n] [m] [count] [out_config] [seed] [rounds]\n";
      std::cerr << "  multistart-alt [n] [m] [count] [out_config] [seed]\n";
      std::cerr << "  certify [config]\n";
      std::cerr << "      LP rigidity certificate at the configuration: is it "
                   "jammed (no first-order descent direction)?\n";
      std::cerr << "  descend-lp [config] [out_config] [rounds]\n";
      std::cerr << "      joint (Q,C) descent driven by the LP direction step; "
                   "reports whether the limit point is rigid (jammed)\n";
      std::cerr << "  polish [config] [out_config] [rounds]\n";
      std::cerr << "      L-BFGS then LP finish: the recommended optimizer; "
                   "reports rigidity of the limit\n";
      std::cerr << "  multistart-lp [n] [m] [count] [out_config] [seed] "
                   "[max_seconds]\n";
      std::cerr << "      LP-direction descents from random seeds, reporting "
                   "the rigid configurations found\n";
      std::cerr << "      alternating SDP/minimax descent from well-rounded "
                   "seeds (fast)\n";
      std::cerr << "      repeated descents from random configurations, "
                   "keeping the best in out_config\n";
      return -1;
    }
    std::string mode = argv[1];
    if (mode == "evaluate" && argc == 3) {
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      DensityResult dr = CoveringDensity(conf);
      std::cout << "n=" << conf.n << " m=" << conf.m
                << " classes=" << dr.cells.size() << "\n";
      std::cout << "mu2=" << dr.mu2 << "\n";
      printf("theta=%.13f\n", dr.theta);
      return 0;
    }
    if (mode == "qstepcheck" && argc == 3) {
      // gradient check of the Q-step barrier, then a full Q-step: verify it
      // improves (or holds) the covering density at fixed cosets.
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      DensityResult dr = CoveringDensity(conf);
      QStepData qd = BuildQStepData(dr.cells, conf.C);
      VectorXd x(qd.dim);
      for (int u = 0; u < qd.dim; u++) {
        x(u) = conf.Q(qd.basis[u].first, qd.basis[u].second);
      }
      double val;
      VectorXd g;
      MatrixXd H;
      QStepBarrier(qd, x, 100.0, val, &g, &H);
      double eps = 1e-6, worst = 0;
      for (int u = 0; u < qd.dim; u++) {
        VectorXd xp = x, xm = x;
        xp(u) += eps;
        xm(u) -= eps;
        double vp, vm;
        QStepBarrier(qd, xp, 100.0, vp, nullptr, nullptr);
        QStepBarrier(qd, xm, 100.0, vm, nullptr, nullptr);
        double fd = (vp - vm) / (2 * eps);
        worst = std::max(worst, std::abs(fd - g(u)) /
                                    std::max(1.0, std::abs(fd)));
      }
      printf("worst relative gradient error = %.3e\n", worst);
      MatrixXd Qopt = QStep(dr.cells, conf.C, conf.Q);
      PeriodicConfig c2 = conf;
      c2.Q = Qopt;
      DensityResult d2 = CoveringDensity(c2);
      printf("density before Q-step = %.12f\n", dr.theta);
      printf("density after  Q-step = %.12f\n", d2.theta);
      return 0;
    }
    if (mode == "gradcheck" && argc == 3) {
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      DensityResult dr = CoveringDensity(conf);
      Packing pk{conf.n, conf.m, conf.n * (conf.n + 1) / 2};
      Eigen::LLT<MatrixXd> llt(conf.Q);
      MatrixXd Lo = llt.matrixL();
      VectorXd x = Pack(pk, Lo, conf.C);
      VectorXd g;
      double beta = 300.0;
      double F0 = ObjectiveGradient(pk, dr.cells, beta, x, g);
      double eps = 1e-6;
      double worst = 0;
      for (int k = 0; k < x.size(); k++) {
        VectorXd xp = x, xm = x, gd;
        xp(k) += eps;
        xm(k) -= eps;
        double fd = (ObjectiveGradient(pk, dr.cells, beta, xp, gd) -
                     ObjectiveGradient(pk, dr.cells, beta, xm, gd)) /
                    (2 * eps);
        double err = std::abs(fd - g(k)) / std::max(1.0, std::abs(fd));
        worst = std::max(worst, err);
      }
      printf("F=%.12f dim=%d worst relative gradient error=%.3e\n", F0,
             int(x.size()), worst);
      return 0;
    }
    if (mode == "joint" && (argc == 4 || argc == 5)) {
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      int rounds = argc == 5 ? atoi(argv[4]) : 60;
      DescendResult res = JointDescent(conf, rounds, std::cerr, true);
      if (!res.success) { std::cerr << "no config evaluated\n"; return 1; }
      printf("theta=%.13f\n", res.theta);
      WriteConfigFile(argv[3], res.conf);
      return 0;
    }
    if (mode == "descend-alt2" && (argc == 4 || argc == 5)) {
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      int rounds = argc == 5 ? atoi(argv[4]) : 30;
      DescendResult res = DescendAlt2(conf, rounds, std::cerr, true);
      if (!res.success) {
        std::cerr << "the descent could not evaluate any configuration\n";
        return 1;
      }
      printf("theta=%.13f\n", res.theta);
      WriteConfigFile(argv[3], res.conf);
      return 0;
    }
    if (mode == "descend-alt" && (argc == 4 || argc == 5)) {
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      int rounds = argc == 5 ? atoi(argv[4]) : 20;
      DescendResult res = DescendAlt(conf, rounds, std::cerr, true);
      if (!res.success) {
        std::cerr << "the descent could not evaluate any configuration\n";
        return 1;
      }
      printf("theta=%.13f\n", res.theta);
      WriteConfigFile(argv[3], res.conf);
      return 0;
    }
    if (mode == "descend" && (argc == 4 || argc == 5)) {
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      int rounds = argc == 5 ? atoi(argv[4]) : 15;
      DescendResult res = Descend(conf, rounds, std::cerr, true);
      if (!res.success) {
        std::cerr << "the descent could not evaluate any configuration\n";
        return 1;
      }
      printf("theta=%.13f\n", res.theta);
      WriteConfigFile(argv[3], res.conf);
      return 0;
    }
    if (mode == "multistart-alt" && (argc == 6 || argc == 7)) {
      int n = atoi(argv[2]);
      int m = atoi(argv[3]);
      int count = atoi(argv[4]);
      unsigned seed = argc == 7 ? unsigned(atol(argv[6])) : 1u;
      std::mt19937_64 gen(seed);
      std::normal_distribution<double> gauss(0.0, 1.0);
      std::uniform_real_distribution<double> unif(0.0, 1.0);
      double best = 1e30;
      for (int st = 0; st < count; st++) {
        PeriodicConfig conf;
        conf.n = n;
        conf.m = m;
        // well-rounded seed: R = I + small gaussian, Q = R^T R near the
        // identity, so the tessellation stays cheap and the descent starts
        // in the region where covering-optimal (well-rounded) forms live
        Eigen::MatrixXd R = Eigen::MatrixXd::Identity(n, n);
        for (int i = 0; i < n; i++) {
          for (int j = 0; j < n; j++) {
            R(i, j) += 0.25 * gauss(gen);
          }
        }
        conf.Q = R.transpose() * R;
        conf.C.resize(m, n);
        conf.C.row(0).setZero();
        for (int t = 1; t < m; t++) {
          for (int j = 0; j < n; j++) {
            conf.C(t, j) = unif(gen);
          }
        }
        DescendResult res = DescendAlt(conf, 25, std::cerr, false);
        if (res.success) {
          printf("start %d: theta=%.13f%s\n", st, res.theta,
                 res.theta < best ? "  *" : "");
          fflush(stdout);
          if (res.theta < best) {
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
    if (mode == "multistart" && (argc >= 6 && argc <= 8)) {
      int n = atoi(argv[2]);
      int m = atoi(argv[3]);
      int count = atoi(argv[4]);
      unsigned seed = argc >= 7 ? unsigned(atol(argv[6])) : 1u;
      int rounds = argc == 8 ? atoi(argv[7]) : 15;
      std::mt19937_64 gen(seed);
      std::normal_distribution<double> gauss(0.0, 1.0);
      std::uniform_real_distribution<double> unif(0.0, 1.0);
      double best = 1e30;
      for (int st = 0; st < count; st++) {
        PeriodicConfig conf;
        conf.n = n;
        conf.m = m;
        Eigen::MatrixXd A(n, n);
        for (int i = 0; i < n; i++) {
          for (int j = 0; j < n; j++) {
            A(i, j) = gauss(gen);
          }
        }
        conf.Q = A * A.transpose() +
                 0.05 * Eigen::MatrixXd::Identity(n, n);
        conf.C.resize(m, n);
        conf.C.row(0).setZero();
        for (int t = 1; t < m; t++) {
          for (int j = 0; j < n; j++) {
            conf.C(t, j) = unif(gen);
          }
        }
        DescendResult res = Descend(conf, rounds, std::cerr, false);
        if (res.success) {
          printf("start %d: theta=%.13f%s\n", st, res.theta,
                 res.theta < best ? "  *" : "");
          fflush(stdout);
          if (res.theta < best) {
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
    if (mode == "certify" && argc == 3) {
      // The pure LP rigidity certificate at a fixed configuration: no descent.
      // Reports the scale-free steepest-descent rate of log Theta over the
      // joint (Q, C) directions. rate ~ 0 means no direction lowers the
      // density to first order -- the configuration is rigid (jammed).
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      DensityResult dr = CoveringDensity(conf);
      auto basis = SymBasis(conf.n);
      int dQ = basis.size();
      int dim = dQ + conf.n * (conf.m - 1);
      std::vector<double> R2;
      std::vector<Eigen::VectorXd> g;
      JointGrad(dr.cells, conf.Q, conf.C, basis, R2, g);
      double mu2 = *std::max_element(R2.begin(), R2.end());
      Eigen::VectorXd gdet = Eigen::VectorXd::Zero(dim);
      Eigen::MatrixXd Qinv = conf.Q.inverse();
      for (int u = 0; u < dQ; u++) {
        int k = basis[u].first, l = basis[u].second;
        gdet(u) = (k == l) ? Qinv(k, k) : 2.0 * Qinv(k, l);
      }
      // A small box so the linearization stays valid (a large box lets the
      // linear model "see" descent that would in fact activate other radii):
      // the reported rate is then the genuine one-sided directional derivative.
      double rho_cert = 1e-3;
      DirLPResult lp = DirectionLP(R2, g, gdet, mu2, conf.n, rho_cert);
      double rate = -lp.predicted / rho_cert;
      printf("theta=%.13f\n", dr.theta);
      printf("orbits=%d active=%d rate=%.6e rigid=%d\n",
             (int)dr.cells.size(), lp.n_active, rate,
             rate < 1e-7 ? 1 : 0);
      return 0;
    }
    if (mode == "descend-lp" && (argc == 4 || argc == 5)) {
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      int rounds = argc == 5 ? atoi(argv[4]) : 200;
      DescendResult res = JointDescentLP(conf, rounds, std::cerr, true);
      if (!res.success) {
        std::cerr << "the descent could not evaluate any configuration\n";
        return 1;
      }
      printf("theta=%.13f\n", res.theta);
      printf("rigid=%d stationarity=%.3e active=%d\n", res.rigid ? 1 : 0,
             res.stationarity, res.n_active);
      WriteConfigFile(argv[3], res.conf);
      return 0;
    }
    if (mode == "polish" && (argc == 4 || argc == 5)) {
      // The recommended optimizer: L-BFGS descends fast to the neighbourhood
      // of a local minimum, then the LP direction step finishes the descent
      // along the non-smooth ridge and certifies whether the limit is rigid.
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      int rounds = argc == 5 ? atoi(argv[4]) : 20;
      DescendResult a = Descend(conf, rounds, std::cerr, true);
      if (!a.success) {
        std::cerr << "the descent could not evaluate any configuration\n";
        return 1;
      }
      std::cerr << "L-BFGS theta=" << a.theta << "; LP finishing...\n";
      DescendResult b = JointDescentLP(a.conf, 300, std::cerr, true);
      DescendResult r = (b.success && b.theta <= a.theta) ? b : a;
      printf("theta=%.13f\n", r.theta);
      printf("rigid=%d stationarity=%.3e active=%d\n", b.rigid ? 1 : 0,
             b.stationarity, b.n_active);
      WriteConfigFile(argv[3], r.conf);
      return 0;
    }
    if (mode == "multistart-lp" && (argc >= 6 && argc <= 8)) {
      int n = atoi(argv[2]);
      int m = atoi(argv[3]);
      int count = atoi(argv[4]);
      unsigned seed = argc >= 7 ? unsigned(atol(argv[6])) : 1u;
      // Optional wall-clock budget (seconds): when given, run starts until the
      // budget is exhausted rather than a fixed count.
      double max_seconds = argc == 8 ? atof(argv[7]) : 0.0;
      // The lattice record to beat in this dimension (best lattice covering).
      double record = 0.0;
      if (n == 3) record = 1.4635030689668180;
      if (n == 4) record = 1.7655285081493524;
      if (n == 5) record = 2.1242859089916246;
      std::mt19937_64 gen(seed);
      std::normal_distribution<double> gauss(0.0, 1.0);
      std::uniform_real_distribution<double> unif(0.0, 1.0);
      double best = 1e30;
      int n_rigid = 0;
      int n_record = 0;
      auto t_start = std::chrono::steady_clock::now();
      for (int st = 0;; st++) {
        double elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(
                             std::chrono::steady_clock::now() - t_start).count();
        if (max_seconds > 0.0) {
          if (elapsed >= max_seconds) break;
        } else {
          if (st >= count) break;
        }
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
        DescendResult a = Descend(conf, 15, std::cerr, false);
        if (!a.success) {
          printf("start %d: failed\n", st);
          fflush(stdout);
          continue;
        }
        // finish on the non-smooth ridge and certify rigidity
        DescendResult res = JointDescentLP(a.conf, 200, std::cerr, false, 90);
        if (!res.success || res.theta > a.theta) res = a;
        if (res.rigid) {
          n_rigid++;
          char fn[4096];
          snprintf(fn, sizeof(fn), "%s.rigid.%d", argv[5], n_rigid);
          WriteConfigFile(fn, res.conf);
        }
        bool beats = record > 0.0 && res.theta < record - 1e-9;
        if (beats) {
          n_record++;
          char fn[4096];
          snprintf(fn, sizeof(fn), "%s.record.%d", argv[5], n_record);
          WriteConfigFile(fn, res.conf);
        }
        printf("start %d: theta=%.13f rigid=%d active=%d rate=%.2e%s%s\n", st,
               res.theta, res.rigid ? 1 : 0, res.n_active, res.stationarity,
               res.theta < best ? "  *" : "",
               beats ? "  <<< BEATS RECORD" : "");
        fflush(stdout);
        if (res.theta < best) {
          best = res.theta;
          WriteConfigFile(argv[5], res.conf);
        }
      }
      printf("best=%.13f rigid_count=%d record_beats=%d (record=%.13f)\n", best,
             n_rigid, n_record, record);
      return 0;
    }
    std::cerr << "unrecognized arguments; run without arguments for usage\n";
    return -1;
  } catch (std::exception const &e) {
    std::cerr << "Error in PERIODIC_JointCoveringSearch: " << e.what() << "\n";
    return 1;
  }
}
