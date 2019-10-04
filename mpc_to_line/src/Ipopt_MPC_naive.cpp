#include "Ipopt_MPC_naive.h"
#include "IpIpoptApplication.hpp"
#include <cassert>
#include <iostream>
#include <math.h>
#include <vector>

using namespace Ipopt;

/** A bunch of global variables **/
//@{

/**
 * Set N and dt
 */
Index N = 2;
// double dt = 0.02;
double dt = 0.5;

// related to map
// int map_sz = 3;
// std::vector<std::vector<double>> waypoints = { {279.02,-182.26,294.22,-176.02,287.39,-178.82,1.9603},
//                                                {276.62,-172.5,291.29,-167.3,284.2,-169.82,1.9115},
//                                                {272.91,-160.58,287.16,-155.43,279.94,-158.03,1.9174} };
// std::vector<double> x0     =  { 288.0, -178.0, 2.0, 19.0 };
// std::vector<double> cl_x   =  { 287.39, 284.2, 279.97 };
// std::vector<double> cl_y   =  { -178.82, -169.82, -158.03 };
// std::vector<double> cl_phi =  { 1.9603, 1.9115, 1.9174 };


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
double v_ref = 20;

// The solver takes all the state variables and actuator
// variables in a singular vector. Thus, we should to establish
// when one variable starts and another ends to make our lifes easier.
Index x_start = 0;
Index y_start = x_start + N;
Index psi_start = y_start + N;
Index v_start = psi_start + N;
Index delta_start = v_start + N;
Index a_start = delta_start + N - 1;

//@}

IpoptMPCNaive::IpoptMPCNaive() {}

// destructor
IpoptMPCNaive::~IpoptMPCNaive() {}

// return the size of the problem
bool IpoptMPCNaive::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // The problem described in Ipopt_MPC.hpp has ( 4 * N + 2 * (N - 1) ) variables
  n = 6 * N - 2;

  // The problem described in Ipopt_MPC.hpp has ( 4 * (N - 1) ) EUQALITY (state function) constraints 
  // and ( 2 * (N - 1) + 6 * N ) INEQUALITY (state/actuator) boundaries
  m = 4 * (N - 1);

  // non-zeros in the jacobian example
  nnz_jac_g = 15 * (N - 1); // 4 + 4 + 4 + 3

  // related to hessian
  nnz_h_lag = 6 * N - 2 + 3 * (N - 1); // crucial change: still need proving

  // use the C style indexing (0-based)
  index_style = C_STYLE;

  // verbosely test
  std::cout << "get_nlp_info() seems to be initialzied correct" << std::endl;

  return true;
}

bool IpoptMPCNaive::get_bounds_info(Index n, Number* x_l, Number* x_u,
                               Index m, Number* g_l, Number* g_u)
{
  // here we need to assert to make sure the variable numbers and constraint numbers
  // are what we think they are
  assert(n == 6 * N - 2);
  assert(m == 4 * N - 4);

  /* lower bound & upper bound */
  // Part I:
  // set all non-actuators upper and lower limits
  // to the max negative and positive values
  for (Index i = x_start; i < delta_start; i++)
  {
    x_l[i] = -1.0e19;
    x_u[i] = 1.0e19;
  }
  // Part II:
  // set all actuators upper and lower limits
  // according to their physics properties
  for (Index i = delta_start; i < a_start; i++)
  {
    x_l[i] = -0.436332;
    x_u[i] = 0.436332;
  }
  for (Index i = a_start; i < n; i++) // n: # of variables
  {
    x_l[i] = -1.0;
    x_u[i] = 1.0;
  }

  // the first constraint g1 has a lower bound of 25
  // g_l[0] = 25;
  // the first constraint g1 has NO upper bound, here we set it to 2e19.
  // Ipopt interprets any number greater than nlp_upper_bound_inf as
  // infinity. The default value of nlp_upper_bound_inf and nlp_lower_bound_inf
  // is 1e19 and can be changed through ipopt options.
  // g_u[0] = 2e19;

  // the second constraint g2 is an equality constraint, so we set the
  // upper and lower bound to the same value
  for (Index i = 0; i < 4 * N - 4; i++)
  {
    g_l[i] = g_u[i] = 0.0; // every constraint is an equality constraint, w/ lower & upper bound equal to 0 
    // g_l[i] = -400000.0;  // dummy test
    // g_u[i] = 400000.0;
  }

  // verbosely test
  std::cout << "get_bounds_info() seems to be initialized correct" << std::endl;

  return true;
}

// return the initial point for the problem
bool IpoptMPCNaive::get_starting_point(Index n, bool init_x, Number* x,
                                  bool init_z, Number* z_L, Number* z_U,
                                  Index m, bool init_lambda,
                                  Number* lambda)
{
  std::vector<double> x0     =  { 288.0, -178.0, 2.0, 19.0 };
  // verbosely test
  std::cout << "get_starting_point() starts to work" << std::endl;

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
  // Initial value of the idependent variables
  // Should be 0 except for the initial values
  for (int i = 0; i < n; i++)
  {
    x[i] = 0.0;
  }
  // Set the initial variable values
  x[x_start]   = x0[0];
  x[y_start]   = x0[1];
  x[psi_start] = x0[2];
  x[v_start]   = x0[3];

  x[x_start + 1]   =  x[x_start] + x[v_start] * cos( x[psi_start] ) * dt; 
  x[y_start + 1]   =  x[y_start] + x[v_start] * sin( x[psi_start] ) * dt; 
  x[psi_start + 1] =  x[psi_start] + x[v_start] / Lf * x[delta_start] * dt;
  x[v_start + 1]   =  x[v_start] + x[a_start] * dt;

  // verbosely test
  std::cout << "get_starting_point() seems to be initialized correct" << std::endl;

  return true;
}

// returns the value of the objective function
bool IpoptMPCNaive::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  // verbosely test
  std::cout << "eval_f() starts to work" << std::endl;

  assert(n == 6 * N - 2);

  /**
   * [Simple Reference]
   * obj_value = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
   */

  obj_value = 0.0;
  obj_value += pow( x[x_start] - 287.39, 2 ) + pow( x[y_start] + 178.82, 2 );
  obj_value += pow( x[x_start + 1] - 284.2, 2 ) + pow( x[y_start + 1] + 169.82, 2 );
  obj_value += pow( x[psi_start] - 1.9603, 2 ) + pow(x[psi_start + 1] - 1.9115, 2 );
  obj_value += pow( x[v_start] - v_ref, 2 ) + pow(x[v_start + 1] - v_ref, 2 );
  obj_value += pow( x[delta_start], 2 ) + pow( x[a_start], 2 );
  
  std::cout << "obj_value = " << obj_value << std::endl;
  // verbosely test
  std::cout << "eval_f() seems to be initialzied correct" << std::endl;

  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool IpoptMPCNaive::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  // verbosely test
  std::cout << "eval_grad_f() starts to work" << std::endl;

  // assert variable numbers
  assert(n == 6 * N - 2);

  /**
   * [Simple Reference]
   * grad_f[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
   * grad_f[1] = x[0] * x[3];
   * grad_f[2] = x[0] * x[3] + 1;
   * grad_f[3] = x[0] * (x[0] + x[1] + x[2]);
   */
  grad_f[x_start]       =   2 * ( x[x_start] - 287.39 );
  grad_f[x_start + 1]   =   2 * ( x[x_start + 1] - 284.2 );
  grad_f[y_start]       =   2 * ( x[y_start] + 178.82 );
  grad_f[y_start + 1]   =   2 * ( x[y_start + 1] + 169.82 );
  grad_f[psi_start]     =   2 * ( x[psi_start] - 1.9603 );
  grad_f[psi_start + 1] =   2 * ( x[psi_start + 1] - 1.9115 );
  grad_f[v_start]       =   2 * ( x[v_start] - v_ref );
  grad_f[v_start + 1]   =   2 * ( x[v_start + 1] - v_ref );
  grad_f[delta_start]   =   2 * ( x[delta_start] );
  grad_f[a_start]       =   2 * ( x[a_start] );

  // verbosely test
  std::cout << "eval_grad_f() seems to be initialzied correct" << std::endl;

  return true;
}

// return the value of the constraints: g(x)
bool IpoptMPCNaive::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  // verbosely test
  std::cout << "eval_g_f() starts to work " << std::endl;

  assert(n == 6 * N - 2);
  assert(m == 4 * N - 4);

  /**
   *  [Simple Reference]
   *  g[0] = x[0] * x[1] * x[2] * x[3];
   *  g[1] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3];
   */

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

  for (Index i = 0; i < N - 1; i++)
  {
    g[4 * i]       =  x[x_start + i]   + x[v_start + i] * cos( x[psi_start + i] ) * dt  - x[x_start + i + 1]   ;
    g[4 * i + 1]   =  x[y_start + i]   + x[v_start + i] * sin( x[psi_start + i] ) * dt  - x[x_start + i + 1]   ;
    g[4 * i + 2]   =  x[psi_start + i] + x[v_start + i] / Lf * x[delta_start + i] * dt  - x[psi_start + i + 1] ;
    g[4 * i + 3]   =  x[v_start + i]   + x[a_start + i] * dt                            - x[v_start + i + 1]   ;
  }

  // verbosely test
  std::cout << "eval_g_f() seems to be initialized correct " << std::endl;

  return true;
}

// return the sturcture or values of the jacobian
bool IpoptMPCNaive::eval_jac_g(Index n, const Number* x, bool new_x,
                          Index m, Index nele_jac, Index* iRow, Index *jCol,
                          Number* values)
{
  // verbosely test
  std::cout << "eval_jac_g() starts to work" << std::endl;

  if (values == NULL)
  {
    for (Index i = 0; i < N - 1; i++)
    { 
      /* @iRow = x_start + i, it contains 4 elements */ 
      // x(t)
      iRow[i * 15] = x_start + i;
      jCol[i * 15] = x_start + i; 
      // v(t)
      iRow[i * 15 + 1] = x_start + i;
      jCol[i * 15 + 1] = v_start + i;
      // psi(t)
      iRow[i * 15 + 2] = x_start + i;
      jCol[i * 15 + 2] = psi_start + i;
      // x(t+1)
      iRow[i * 15 + 3] = x_start + i;
      jCol[i * 15 + 3] = x_start + i + 1;

      /* @iRow = y_start + i, it contains 4 elements */ 
      // y(t)
      iRow[i * 15 + 4] = y_start + i;
      jCol[i * 15 + 4] = y_start + i; 
      // v(t)
      iRow[i * 15 + 5] = y_start + i;
      jCol[i * 15 + 5] = v_start + i;
      // psi(t)
      iRow[i * 15 + 6] = y_start + i;
      jCol[i * 15 + 6] = psi_start + i;
      // y(t+1)
      iRow[i * 15 + 7] = y_start + i;
      jCol[i * 15 + 7] = y_start + i + 1; 

      /* @iRow = psi_start + i, it contains 4 elements */ 
      // psi(t)
      iRow[i * 15 + 8] = psi_start + i;
      jCol[i * 15 + 8] = psi_start + i; 
      // v(t)
      iRow[i * 15 + 9] = psi_start + i;
      jCol[i * 15 + 9] = v_start + i;
      // delta(t)
      iRow[i * 15 + 10] = psi_start + i;
      jCol[i * 15 + 10] = delta_start + i;
      // psi(t+1)
      iRow[i * 15 + 11] = psi_start + i;
      jCol[i * 15 + 11] = psi_start + i + 1;

      /* @iRow = v_start + i, it contains 3 elements */ 
      // v(t)
      iRow[i * 15 + 12] = v_start + i;
      jCol[i * 15 + 12] = v_start + i; 
      // a(t)
      iRow[i * 15 + 13] = v_start + i;
      jCol[i * 15 + 13] = a_start + i;
      // v(t+1)
      iRow[i * 15 + 14] = v_start + i;
      jCol[i * 15 + 14] = v_start + i + 1;
    }

    // verbosely test
    std::cout << "eval_jac_g() structure part seems to be initialized correct " << std::endl;
  }
  else
  { 
    // verbosely test
    std::cout << "eval_jac_g() starts to work on values part " << std::endl;
    for (Index i = 0; i < N - 1; i++)
    { 
      /* @iRow = x_start + i, it contains 4 elements */ 
      // x(t)
      values[i * 15] = 1.0;
      // v(t)
      values[i * 15 + 1] = cos( x[psi_start + i] ) * dt; // cos( psi(t) ) * dt
      // psi(t)
      values[i * 15 + 2] = -sin( x[psi_start + i] ) * x[v_start + i] * dt; // -v(t) * sin( psi(t) ) * dt
      // x(t+1)
      values[i * 15 + 3] = -1.0;

      /* @iRow = y_start + i, it contains 4 elements */ 
      // y(t)
      values[i * 15 + 4] = 1.0;
      // v(t)
      values[i * 15 + 5] = sin( x[psi_start + i] ) * dt; // sin( psi(t) ) * dt
      // psi(t)
      values[i * 15 + 6] = cos( x[psi_start + i] ) * x[v_start + i] * dt; // v(t) * cos( psi(t) ) * dt
      // y(t+1)
      values[i * 15 + 7] = -1.0;

      /* @iRow = psi_start + i, it contains 4 elements */ 
      // psi(t)
      values[i * 15 + 8] = 1.0;
      // v(t)
      values[i * 15 + 9] = x[delta_start + i] * dt / Lf; // delta(t) * dt / Lf
      // delta(t)
      values[i * 15 + 10] = x[v_start + i] * dt / Lf; // v(t) * dt / Lf
      // psi(t+1)
      values[i * 15 + 11] = -1.0;

      /* @iRow = v_start + i, it contains 3 elements */ 
      // v(t)
      values[i * 15 + 12] = 1.0;
      // a(t)
      values[i * 15 + 13] = dt;
      // v(t+1)
      values[i * 15 + 14] = -1.0;
    }

    // verbosely test
    std::cout << "eval_jac_g() values part seems to be initialized correct " << std::endl;
  }
  
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
bool IpoptMPCNaive::eval_h(Index n, const Number* x, bool new_x,
                             Number obj_factor, Index m, const Number* lambda,
                             bool new_lambda, Index nele_hess, Index* iRow,
                             Index* jCol, Number* values)
{ 
  // verbosely test
  std::cout << "eval_h() starts to work" << std::endl;

  // [Simple Reference]
  /* if (values == NULL) {
    // return the structure. This is a symmetric matrix, fill the lower left
    // triangle only.

    // the hessian for this problem is actually dense
    Index idx=0;
    for (Index row = 0; row < 4; row++) {
      for (Index col = 0; col <= row; col++) {
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

void IpoptMPCNaive::finalize_solution(SolverReturn status,
                                 Index n, const Number* x, const Number* z_L, const Number* z_U,
                                 Index m, const Number* g, const Number* lambda,
                                 Number obj_value,
                                 const IpoptData* ip_data,
                                 IpoptCalculatedQuantities* ip_cq)
{
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution.

  // verbosely test
  std::cout << "finalize_solution() starts to work" << std::endl;

  // parse the solution to each state vector & actuator vector

  // [Simple reference]
  // For this example, we write the solution to the console
  std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
  for (Index i=0; i<n; i++) {
     std::cout << "x[" << i << "] = " << x[i] << std::endl;
  }

  std::cout << std::endl << std::endl << "Solution of the bound multipliers, z_L and z_U" << std::endl;
  for (Index i=0; i<n; i++) {
    std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
  }
  for (Index i=0; i<n; i++) {
    std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
  }

  std::cout << std::endl << std::endl << "Objective value" << std::endl;
  std::cout << "f(x*) = " << obj_value << std::endl;

  std::cout << std::endl << "Final value of the constraints:" << std::endl;
  for (Index i=0; i<m ;i++) {
    std::cout << "g(" << i << ") = " << g[i] << std::endl;
  }
}