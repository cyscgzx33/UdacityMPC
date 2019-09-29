#include "Ipopt_MPC.h"
#include <math.h>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "IpIpoptApplication.hpp"

using Eigen::VectorXd;

/**
 * Set N and dt
 */
size_t N = 50;
double dt = 0.02;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

// NOTE: feel free to play around with this or do something completely different
double ref_v = 20;

// The solver takes all the state variables and actuator
// variables in a singular vector. Thus, we should to establish
// when one variable starts and another ends to make our lifes easier.
size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
size_t v_start = psi_start + N;
size_t cte_start = v_start + N;
size_t epsi_start = cte_start + N;
size_t delta_start = epsi_start + N;
size_t a_start = delta_start + N - 1;

class FG_eval {
 public:
  VectorXd coeffs;
  // CppAD::vector<AD<double>> cl_x_;
  // CppAD::vector<AD<double>> cl_y_;
  // CppAD::vector<AD<double>> cl_phi_;
  std::vector<double> cl_x_;
  std::vector<double> cl_y_;
  std::vector<double> cl_phi_;

  // Coefficients of the fitted polynomial.
  FG_eval(VectorXd coeffs) { this->coeffs = coeffs; }
  FG_eval(std::vector<double>& cl_x, std::vector<double>& cl_y, std::vector<double>& cl_phi) 
  { 
    // for (int i = 0; i < cl_phi.size; i++)
    // {
    //   cl_x_.push_back( cl_x[i] );
    //   cl_y_.push_back( cl_y[i] );
    //   cl_phi_.push_back( cl_phi[i] );
    // }
    cl_x_ = cl_x;
    cl_y_ = cl_y;
    cl_phi_ = cl_phi;
  }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  // `fg` is a vector containing the cost and constraints.
  // `vars` is a vector containing the variable values (state & actuators).
  void operator()(ADvector& fg, const ADvector& vars) {
    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    fg[0] = 0;

    // Reference State Cost
    /**
     * Define the cost related the reference and the state in 3 parts:
     */

    // Part I: reference part
    for (int t = 0; t < N; t++)
    {
      fg[0] += CppAD::pow( vars[cte_start + t], 2 );  // cte has a (const) reference value equals to 0
      fg[0] += CppAD::pow( vars[epsi_start + t], 2 ); // epsi has a (const) reference value equals to 0
      fg[0] += CppAD::pow( vars[v_start + t] - ref_v , 2 ); // v has a (const) reference value equals to ref_v
    }

    // Part II: actuator part
    for (int t = 0; t < N - 1; t++)
    {
      fg[0] += CppAD::pow( vars[delta_start + t], 2 );  
      fg[0] += CppAD::pow( vars[a_start + t], 2 );
    }

    // Part III: actuator changing rate part
    for (int t = 0; t < N - 2; t++)
    {
      fg[0] += CppAD::pow( vars[delta_start + t + 1] - vars[delta_start + t], 2 );
      fg[0] += CppAD::pow( vars[a_start + t + 1] - vars[a_start + t], 2 );
    }

    //
    // Setup Constraints
    //
    // NOTE: In this section you'll setup the model constraints.

    // Initial constraints
    //
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`.
    // This bumps up the position of all the other values.
    fg[1 + x_start]     =  vars[x_start];
    fg[1 + y_start]     =  vars[y_start];
    fg[1 + psi_start]   =  vars[psi_start];
    fg[1 + v_start]     =  vars[v_start];
    fg[1 + cte_start]   =  vars[cte_start];
    fg[1 + epsi_start]  =  vars[epsi_start];

    // The rest of the constraints
    for (int t = 1; t < N; ++t) {
      /**
       *   Grab the rest of the states at t+1 and t.
       *   We have given you parts of these states below.
       */

      // consider states at time t and t + 1
      AD<double> x1        =  vars[x_start + t];
      AD<double> y1        =  vars[y_start + t];
      AD<double> psi1      =  vars[psi_start + t];
      AD<double> v1        =  vars[v_start + t];
      AD<double> cte1      =  vars[cte_start + t];
      AD<double> epsi1     =  vars[epsi_start + t];

      AD<double> x0        =  vars[x_start + t - 1];
      AD<double> y0        =  vars[y_start + t];
      AD<double> psi0      =  vars[psi_start + t - 1];
      AD<double> v0        =  vars[v_start + t - 1];
      AD<double> cte0      =  vars[cte_start + t - 1];
      AD<double> epsi0     =  vars[epsi_start + t - 1];

      // only consider the actuators at time t
      AD<double> delta0    =  vars[delta_start + t -1];
      AD<double> a0        =  vars[a_start + t - 1];

      // generate some assistant variables
      // AD<double> f_x0      =  coeffs[1] * x0 + coeffs[0];
      // AD<double> psi_des0  =  CppAD::atan( coeffs[1] );

      // Here's `x` to get you started.
      // The idea here is to constraint this value to be 0.
      
      /** 
       *  [Description]
       *  Recall the equations for the kinematics model:
       *  
       *  x_[t+1]   =  x[t] + v[t] * cos(psi[t]) * dt
       *  y_[t+1]   =  y[t] + v[t] * sin(psi[t]) * dt
       *  psi_[t+1] =  psi[t] + v[t] / Lf * delta[t] * dt
       *  v_[t+1]   =  v[t] + a[t] * dt
       * 
       * */
      
      // NOTE: The use of `AD<double>` and use of `CppAD`!
      // CppAD can compute derivatives and pass these to the solver.

      /**
       *   Setup the rest of the model constraints
       */

      // calculate cte & epsi based on (x, y, psi) and (cl_x, cl_y, cl_phi)
      std::vector<AD<double>> dist(cl_phi_.size(), 0);
      std::vector<int> ind( cl_phi_.size() );
      for (int i = 0; i < cl_phi_.size(); i++)
      {
        dist[i] = (x1 - cl_x_[i]) * (x1 - cl_x_[i]) + (y1 - cl_y_[i]) * (y1 - cl_y_[i]);
        ind[i] = i;
      }
      CppAD::index_sort(dist, ind);
      int closest_ind = ind[0];

      fg[1 + x_start + t]     =  x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start + t]     =  y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_start + t]   =  psi1 - (psi0 + v0 / Lf * delta0 * dt);
      fg[1 + v_start + t]     =  v1 - (v0 + a0 * dt);
      fg[1 + cte_start + t]   =  cte1 - CppAD::sqrt( CppAD::pow(x1 - cl_x_[closest_ind], 2) + CppAD::pow(y1 - cl_y_[closest_ind], 2) );
      fg[1 + epsi_start + t]  =  epsi1 - psi1 + cl_phi_[closest_ind];
    }
  }
};

//
// IpoptMPC class definition
//

IpoptMPC::IpoptMPC()
{
  readRoadmapFromCSV("/home/honda/git/UdacityMPC/mpc_to_line/roadmap.csv");

  // parse waypoints_ to center line (x, y, phi)
  for (auto& wp : waypoints_)
  {
    cl_x_.push_back( wp[4] );
    cl_y_.push_back( wp[5] );
    cl_phi_.push_back( atan(wp[6]) ); // slope to direction angle in [rad]
  }
}
IpoptMPC::~IpoptMPC() {}

// std::vector<double> CustomMPC::Solve(const VectorXd& x0, const VectorXd& coeffs) {
std::vector<double> IpoptMPC::Solve(const VectorXd& x0) {
  typedef CPPAD_TESTVECTOR(double) Dvector;

  double x    =  x0[0];
  double y    =  x0[1];
  double psi  =  x0[2];
  double v    =  x0[3];
  double cte  =  x0[4];
  double epsi =  x0[5];

  // number of independent variables
  // N timesteps == N - 1 actuations
  size_t n_vars = N * 6 + (N - 1) * 2;
  // Number of constraints
  size_t n_constraints = N * 6;

  // Initial value of the independent variables.
  // Should be 0 except for the initial values.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; ++i) {
    vars[i] = 0.0;
  }
  // Set the initial variable values
  vars[x_start]     =  x;
  vars[y_start]     =  y;
  vars[psi_start]   =  psi;
  vars[v_start]     =  v;
  vars[cte_start]   =  cte;
  vars[epsi_start]  =  epsi;

  // Lower and upper limits for x
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (int i = 0; i < delta_start; ++i) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  // NOTE: Feel free to change this to something else.
  for (int i = delta_start; i < a_start; ++i) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
  for (int i = a_start; i < n_vars; ++i) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for constraints
  // All of these should be 0 except the initial
  // state indices.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; ++i) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  // Object that computes objective and constraints
  FG_eval fg_eval(cl_x_, cl_y_, cl_phi_);

  // options
  std::string options;
  options += "Integer print_level  0\n";
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  //
  // Check some of the solution values
  //
  bool ok = true;
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;
  return {solution.x[x_start + 1],   solution.x[y_start + 1],
          solution.x[psi_start + 1], solution.x[v_start + 1],
          solution.x[cte_start + 1], solution.x[epsi_start + 1],
          solution.x[delta_start],   solution.x[a_start]};
}

// constructor
IpoptMPC::IpoptMPC() {}

// destructor
IpoptMPC::~IpoptMPC() {}

// return the size of the problem
bool IpoptMPC::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
                            Ipopt::Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style)
{
  // The problem described in Ipopt_MPC.hpp has ( 4 * N + 2 * (N - 1) ) ??? variables
  n = 6 * N - 2;

  // The problem described in Ipopt_MPC.hpp has ( 4 * (N - 1) + 2 * (N - 1) )???  constraints 
  m = 6 * N - 6;

  // non-zeros in the jacobian example
  nnz_jac_g = N; // ???

  // related to hessian
  nnz_h_lag = N; // ???

  // use the C style indexing (0-based)
  index_style = TNLP::C_STYLE;

  return true;
}

bool IpoptMPC::get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
                               Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u)
{
  // here we need to assert to make sure the variable numbers and constraint numbers
  // are what we think they are
  assert(n == 6 * N - 2); // ???
  assert(m == 6 * N - 6); // ???

  // lower bound & upper bound
  for (Ipopt::Index i = 0; i < n; i++)
  {
    // x_l[i] = ???
    // x_u[i] = ???
  }

  // the first constraint g1 has a lower bound of 25
  g_l[0] = 25; // ???
  // the first constraint g1 has NO upper bound, here we set it to 2e19.
  // Ipopt interprets any number greater than nlp_upper_bound_inf as
  // infinity. The default value of nlp_upper_bound_inf and nlp_lower_bound_inf
  // is 1e19 and can be changed through ipopt options.
  g_u[0] = 2e19; // ???

  // the second constraint g2 is an equality constraint, so we set the
  // upper and lower bound to the same value
  g_l[1] = g_u[1] = 40.0; // ???

  return true;
}

// return the initial point for the problem
bool IpoptMPC::get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
                                  bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
                                  Ipopt::Index m, bool init_lambda,
                                  Ipopt::Number* lambda)
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the dual variables
  // if you wish
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  // initialize to the given starting point
  /**
   * [Simple Reference]
   * x[0] = 1.0;
   * x[1] = 5.0;
   * x[2] = 5.0;
   * x[3] = 1.0;
   */

  return true;
}

// returns the value of the objective function
bool IpoptMPC::eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value)
{
  assert(n == 6 * N - 2); // ???

  /**
   * [Simple Reference]
   * obj_value = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
   */

  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool IpoptMPC::eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f)
{ 
  // assert variable numbers
  assert(n == 6 * N - 2); // ???

  /**
   * [Simple Reference]
   * grad_f[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
   * grad_f[1] = x[0] * x[3];
   * grad_f[2] = x[0] * x[3] + 1;
   * grad_f[3] = x[0] * (x[0] + x[1] + x[2]);
   */

  return true;
}

// return the value of the constraints: g(x)
bool IpoptMPC::eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g)
{
  assert(n == 6 * N - 2); // ???
  assert(m == 6 * N - 6); // ???

  /**
   *  [Simple Reference]
   *  g[0] = x[0] * x[1] * x[2] * x[3];
   *  g[1] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3];
   */

  return true;
}

// return the sturcture or values of the jacobian
bool IpoptMPC::eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                          Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index *jCol,
                          Ipopt::Number* values)
{
  // [Simple Reference]
  /* if (values == NULL) {
    // return the structure of the jacobian

    // this particular jacobian is dense
    iRow[0] = 0;
    jCol[0] = 0;
    iRow[1] = 0;
    jCol[1] = 1;
    iRow[2] = 0;
    jCol[2] = 2;
    iRow[3] = 0;
    jCol[3] = 3;
    iRow[4] = 1;
    jCol[4] = 0;
    iRow[5] = 1;
    jCol[5] = 1;
    iRow[6] = 1;
    jCol[6] = 2;
    iRow[7] = 1;
    jCol[7] = 3;
  }
  else {
    // return the values of the jacobian of the constraints

    values[0] = x[1]*x[2]*x[3]; // 0,0
    values[1] = x[0]*x[2]*x[3]; // 0,1
    values[2] = x[0]*x[1]*x[3]; // 0,2
    values[3] = x[0]*x[1]*x[2]; // 0,3

    values[4] = 2*x[0]; // 1,0
    values[5] = 2*x[1]; // 1,1
    values[6] = 2*x[2]; // 1,2
    values[7] = 2*x[3]; // 1,3
  } */

  return true;
}

//return the structure or values of the hessian
bool IpoptMPC::eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                      Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
                      bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
                      Ipopt::Index* jCol, Ipopt::Number* values)
{
  // [Simple Reference]
  /* if (values == NULL) {
    // return the structure. This is a symmetric matrix, fill the lower left
    // triangle only.

    // the hessian for this problem is actually dense
    Ipopt::Index idx=0;
    for (Ipopt::Index row = 0; row < 4; row++) {
      for (Ipopt::Index col = 0; col <= row; col++) {
        iRow[idx] = row;
        jCol[idx] = col;
        idx++;
      }
    }

    assert(idx == nele_hess);
  }
  else {
    // return the values. This is a symmetric matrix, fill the lower left
    // triangle only

    // fill the objective portion
    values[0] = obj_factor * (2*x[3]); // 0,0

    values[1] = obj_factor * (x[3]);   // 1,0
    values[2] = 0.;                    // 1,1

    values[3] = obj_factor * (x[3]);   // 2,0
    values[4] = 0.;                    // 2,1
    values[5] = 0.;                    // 2,2

    values[6] = obj_factor * (2*x[0] + x[1] + x[2]); // 3,0
    values[7] = obj_factor * (x[0]);                 // 3,1
    values[8] = obj_factor * (x[0]);                 // 3,2
    values[9] = 0.;                                  // 3,3


    // add the portion for the first constraint
    values[1] += lambda[0] * (x[2] * x[3]); // 1,0

    values[3] += lambda[0] * (x[1] * x[3]); // 2,0
    values[4] += lambda[0] * (x[0] * x[3]); // 2,1

    values[6] += lambda[0] * (x[1] * x[2]); // 3,0
    values[7] += lambda[0] * (x[0] * x[2]); // 3,1
    values[8] += lambda[0] * (x[0] * x[1]); // 3,2

    // add the portion for the second constraint
    values[0] += lambda[1] * 2; // 0,0

    values[2] += lambda[1] * 2; // 1,1

    values[5] += lambda[1] * 2; // 2,2

    values[9] += lambda[1] * 2; // 3,3
  } */

  return true;
}

void IpoptMPC::finalize_solution(Ipopt::SolverReturn status,
                                 Ipopt::Index n, const Ipopt::Number* x, const Ipopt::Number* z_L, const Ipopt::Number* z_U,
                                 Ipopt::Index m, const Ipopt::Number* g, const Ipopt::Number* lambda,
                                 Ipopt::Number obj_value,
                                 const Ipopt::IpoptData* ip_data,
                                 Ipopt::IpoptCalculatedQuantities* ip_cq)
{
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution.

  // [Simple reference]
  // For this example, we write the solution to the console
  /* std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
  for (Ipopt::Index i=0; i<n; i++) {
     std::cout << "x[" << i << "] = " << x[i] << std::endl;
  }

  std::cout << std::endl << std::endl << "Solution of the bound multipliers, z_L and z_U" << std::endl;
  for (Ipopt::Index i=0; i<n; i++) {
    std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
  }
  for (Ipopt::Index i=0; i<n; i++) {
    std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
  }

  std::cout << std::endl << std::endl << "Objective value" << std::endl;
  std::cout << "f(x*) = " << obj_value << std::endl;

  std::cout << std::endl << "Final value of the constraints:" << std::endl;
  for (Ipopt::Index i=0; i<m ;i++) {
    std::cout << "g(" << i << ") = " << g[i] << std::endl;
  } */
}