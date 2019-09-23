#include <vector>
#include "Eigen-3.3/Eigen/QR"
#include "helpers.h"
#include "matplotlibcpp.h"
#include "custom_MPC.h"

namespace plt = matplotlibcpp;

using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

int main() {
  CustomMPC mpc;
  int iters = 500;

  double x = 288.0;
  double y = -178.0;
  double psi = 2.0;
  double v = 39;
  /**
   * calculate the cross track error
   * since the cross track error is the perpendicular dist from point to ref
   * and the stright line is a horizontal one:
   */
  double cte = sqrt( pow(x - 287.39, 2) + pow(y - (-178.82), 2) );
  /**
   * calculate the orientation error
   * due to the sign starting at 0, the orientation error is -f'(x)
   * derivative of ( coeffs[1] * x + coeffs[0] ) -> coeffs[1]
   */
  double epsi = psi - 1.9603;

  VectorXd state(6);
  state << x, y, psi, v, cte, epsi;

  vector<double> x_vals = {state[0]};
  vector<double> y_vals = {state[1]};
  vector<double> psi_vals = {state[2]};
  vector<double> v_vals = {state[3]};
  vector<double> cte_vals = {state[4]};
  vector<double> epsi_vals = {state[5]};
  vector<double> delta_vals = {};
  vector<double> a_vals = {};

  for (size_t i = 0; i < iters; ++i) {
    cout << "Iteration " << i << endl;

    auto vars = mpc.Solve(state);

    x_vals.push_back(vars[0]);
    y_vals.push_back(vars[1]);
    psi_vals.push_back(vars[2]);
    v_vals.push_back(vars[3]);
    cte_vals.push_back(vars[4]);
    epsi_vals.push_back(vars[5]);

    delta_vals.push_back(vars[6]);
    a_vals.push_back(vars[7]);

    state << vars[0], vars[1], vars[2], vars[3], vars[4], vars[5];
    cout << "x = " << vars[0] << endl;
    cout << "y = " << vars[1] << endl;
    cout << "psi = " << vars[2] << endl;
    cout << "v = " << vars[3] << endl;
    cout << "cte = " << vars[4] << endl;
    cout << "epsi = " << vars[5] << endl;
    cout << "delta = " << vars[6] << endl;
    cout << "a = " << vars[7] << endl;
    cout << endl;
  }

  // Plot values
  // NOTE: feel free to play around with this.
  // It's useful for debugging!
  plt::subplot(3, 2, 1);
  plt::title("CTE");
  plt::plot(cte_vals);
  plt::subplot(3, 2, 2);
  plt::title("Delta (Radians)");
  plt::plot(delta_vals);
  plt::subplot(3, 2, 3);
  plt::title("Speed [m/s]");
  plt::plot(v_vals);
  plt::subplot(3, 2, 4);
  plt::title("pos x [m]");
  plt::plot(x_vals);
  plt::subplot(3, 2, 5);
  plt::title("pos y [m]");
  plt::plot(y_vals);
  plt::subplot(3, 2, 6);
  plt::title("acceleration [m/s^2]");
  plt::plot(a_vals);

  plt::show();
}