#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {
 public:
  MPC();

  virtual ~MPC();

  double delta_prev {0};
  double a_prev {0.1};
  int latency = 2;

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<vector<double>> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

#endif /* MPC_H */
