// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: TDNaive.cpp 2005 2011-06-06 12:55:16Z stefan $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#include "apply_ipopt_kinematic_model.hpp"
#include "IpIpoptApplication.hpp"
#include <cassert>
#include <iostream>
#include <algorithm>

using namespace Ipopt;

/* global variable */
const double dt = 0.02;
const double v_ref = 20.0;
const double Lf = 2.67;


// constructor
TDNaive::TDNaive()
{}

//destructor
TDNaive::~TDNaive()
{}

// returns the size of the problem
bool TDNaive::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                             Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // The problem described in TDNaive.hpp has 4 variables, x[0] through x[3]
  n = 10;

  // one equality constraint and one inequality constraint
  m = 4;

  // in this example the jacobian is dense and contains 8 nonzeros
  nnz_jac_g = 15;

  // the hessian is also dense and has 16 total nonzeros, but we
  // only need the lower left corner (since it is symmetric)
  nnz_h_lag = 12;

  // use the C style indexing (0-based)
  index_style = TNLP::C_STYLE;

  return true;
}

// returns the variable bounds
bool TDNaive::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{
  // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.
  assert(n == 10);
  assert(m == 4);

  // the variables have lower bounds of 1
  for (Index i=0; i<10; i++) {
    x_l[i] = -500;
  }

  // the variables have upper bounds of 5
  for (Index i=0; i<10; i++) {
    x_u[i] = 500;
  }

  // the first constraint g1 has a lower bound of 25
  g_l[0] = 0.0;
  // the first constraint g1 has NO upper bound, here we set it to 2e19.
  // Ipopt interprets any number greater than nlp_upper_bound_inf as
  // infinity. The default value of nlp_upper_bound_inf and nlp_lower_bound_inf
  // is 1e19 and can be changed through ipopt options.
  g_u[0] = 0.0;

  // the second constraint g2 is an equality constraint, so we set the
  // upper and lower bound to the same value
  g_l[1] = g_u[1] = 0.0;
  g_l[2] = g_u[2] = 0.0;
  g_l[3] = g_u[3] = 0.0;

  return true;
}

// returns the initial point for the problem
bool TDNaive::get_starting_point(Index n, bool init_x, Number* x,
                                   bool init_z, Number* z_L, Number* z_U,
                                   Index m, bool init_lambda,
                                   Number* lambda)
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the dual variables
  // if you wish
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  // initialize to the given starting point
  x[0] = 287.0;
  x[1] = 282.1;
  x[2] = -178.0;
  x[3] = -168.2;
  x[4] = 1.95;
  x[5] = 1.90;
  x[6] = 18.0;
  x[7] = 18.4;
  x[8] = 0.5;
  x[9] = 2.0;
  return true;
}

// returns the value of the objective function
bool TDNaive::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  assert(n == 10);

  // map info
  std::vector<double> cl_x   =  { 287.39, 284.2, 279.97 };
  std::vector<double> cl_y   =  { -178.82, -169.82, -158.03 };
  std::vector<double> cl_phi =  { 1.9603, 1.9115, 1.9174 };

  // find closest point
  std::vector<double> dist_x0(3, 0.0);
  std::vector<double> dist_x1(3, 0.0);
  for (int i = 0; i < 3; i++)
  {
    dist_x0[i] = pow(cl_x[i] - x[0], 2) + pow(cl_y[i] - x[2], 2);
    dist_x1[i] = pow(cl_x[i] - x[1], 2) + pow(cl_y[i] - x[3], 2);
  }
  auto closest_x0 = std::min_element( dist_x0.begin(), dist_x0.end() );
  auto closest_x1 = std::min_element( dist_x1.begin(), dist_x1.end() );
  int cl_idx_x0 = closest_x0 - dist_x0.begin();
  int cl_idx_x1 = closest_x1 - dist_x1.begin(); 

  obj_value = 0.0;

  // position 0
  obj_value += pow(x[0] - cl_x[cl_idx_x0], 2) + pow(x[2] - cl_y[cl_idx_x0], 2);
  obj_value += pow(x[4] - cl_phi[cl_idx_x0], 2);
  obj_value += pow(x[6] - v_ref, 2);
  // position 1
  obj_value += pow(x[1] - cl_x[cl_idx_x1], 2) + pow(x[3] - cl_y[cl_idx_x1], 2);
  obj_value += pow(x[5] - cl_phi[cl_idx_x1], 2);
  obj_value += pow(x[7] - v_ref, 2);
  // actuators
  obj_value += x[8] * x[8];
  obj_value += x[9] * x[9];

  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool TDNaive::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  assert(n == 10);

  // map info
  std::vector<double> cl_x   =  { 287.39, 284.2, 279.97 };
  std::vector<double> cl_y   =  { -178.82, -169.82, -158.03 };
  std::vector<double> cl_phi =  { 1.9603, 1.9115, 1.9174 };

  // find closest point
  std::vector<double> dist_x0(3, 0.0);
  std::vector<double> dist_x1(3, 0.0);
  for (int i = 0; i < 3; i++)
  {
    dist_x0[i] = pow(cl_x[i] - x[0], 2) + pow(cl_y[i] - x[2], 2);
    dist_x1[i] = pow(cl_x[i] - x[1], 2) + pow(cl_y[i] - x[3], 2);
  }
  auto closest_x0 = std::min_element( dist_x0.begin(), dist_x0.end() );
  auto closest_x1 = std::min_element( dist_x1.begin(), dist_x1.end() );
  int cl_idx_x0 = closest_x0 - dist_x0.begin();
  int cl_idx_x1 = closest_x1 - dist_x1.begin(); 

  grad_f[0] = 2 * ( x[0] - cl_x[cl_idx_x0] );
  grad_f[1] = 2 * ( x[1] - cl_x[cl_idx_x1] );
  grad_f[2] = 2 * ( x[2] - cl_y[cl_idx_x0] );
  grad_f[3] = 2 * ( x[3] - cl_y[cl_idx_x1] );
  grad_f[4] = 2 * ( x[4] - cl_phi[cl_idx_x0] );
  grad_f[5] = 2 * ( x[5] - cl_phi[cl_idx_x1] );
  grad_f[6] = 2 * ( x[6] - v_ref );
  grad_f[7] = 2 * ( x[7] - v_ref );
  grad_f[8] = 2 * x[8];
  grad_f[9] = 2 * x[9];

  return true;
}

// return the value of the constraints: g(x)
bool TDNaive::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  assert(n == 10);
  assert(m == 4);

  g[0] = x[0] + x[6] * cos( x[4] ) * dt - x[1];
  g[1] = x[2] + x[6] * sin( x[4] ) * dt - x[3];
  g[2] = x[4] + x[6] * x[8] / Lf * dt - x[5];
  g[3] = x[6] + x[9] * dt - x[7];

  return true;
}

// return the structure or values of the jacobian
bool TDNaive::eval_jac_g(Index n, const Number* x, bool new_x,
                         Index m, Index nele_jac, Index* iRow, Index *jCol,
                         Number* values)
{
  if (values == NULL) {
    // return the structure of the jacobian
    // this particular jacobian is dense
    iRow[0] = 0;
    jCol[0] = 0;
    iRow[1] = 0;
    jCol[1] = 1;
    iRow[2] = 0;
    jCol[2] = 4;
    iRow[3] = 0;
    jCol[3] = 5;
    iRow[4] = 1;
    jCol[4] = 2;
    iRow[5] = 1;
    jCol[5] = 3;
    iRow[6] = 1;
    jCol[6] = 4;
    iRow[7] = 1;
    jCol[7] = 6;
    iRow[8] = 2;
    jCol[8] = 4;
    iRow[9] = 2;
    jCol[9] = 5;
    iRow[10] = 2;
    jCol[10] = 6;
    iRow[11] = 2;
    jCol[11] = 8;
    iRow[12] = 3;
    jCol[12] = 6;
    iRow[13] = 3;
    jCol[13] = 7;
    iRow[14] = 3;
    jCol[14] = 9;
  }
  else {
    // return the values of the jacobian of the constraints

    values[0] = 1.0; // 0,0
    values[1] = -1.0; // 0,1
    values[2] = -sin( x[4] ) * x[6] * dt; // 0,4
    values[3] = cos( x[4] ) * dt; // 0,6

    values[4] = 1.0; // 1,2
    values[5] = -1.0; // 1,3
    values[6] = cos( x[4] ) * x[6] * dt; // 1,4
    values[7] = sin( x[4] ) * dt; // 1,6

    values[8] = 1.0; // 2,4
    values[9] = -1.0; // 2,5
    values[10] = x[8] * dt / Lf; // 2,6
    values[11] = x[6] * dt / Lf; // 2,8

    values[12] = 1.0; // 3,6
    values[13] = -1.0; // 3,7
    values[14] = dt; // 3,9
  }

  return true;
}

//return the structure or values of the hessian
bool TDNaive::eval_h(Index n, const Number* x, bool new_x,
                       Number obj_factor, Index m, const Number* lambda,
                       bool new_lambda, Index nele_hess, Index* iRow,
                       Index* jCol, Number* values)
{

  if (values == NULL) {
    // return the structure. This is a symmetric matrix, fill the lower left
    // triangle only.

    // the hessian for this problem is actually dense
    Index idx=0;
    for (Index row = 0; row < 10; row++) {
        iRow[idx] = row;
        jCol[idx] = row;
        idx++;
    }

    // element 10
    iRow[idx] = 6;
    jCol[idx] = 4;
    idx++;

    // element 11
    iRow[idx] = 8;
    jCol[idx] = 6;
    idx++;

    assert(idx == nele_hess);
  }
  else {
    // return the values. This is a symmetric matrix, fill the lower left
    // triangle only

    // fill the objective portion
    values[0] += obj_factor * 2; // 0,0
    values[1] += obj_factor * 2; // 1,1
    values[2] += obj_factor * 2; // 2,2
    values[3] += obj_factor * 2; // 3,3
    values[4] += obj_factor * 2; // 4,4
    values[5] += obj_factor * 2; // 0,0
    values[6] += obj_factor * 2; // 1,1
    values[7] += obj_factor * 2; // 2,2
    values[8] += obj_factor * 2; // 3,3
    values[9] += obj_factor * 2; // 4,4

    values[4] += lambda[0] * x[6] * cos( x[4] ) * dt; // 4,4
    values[4] += lambda[1] * x[6] * ( -sin( x[4] ) ) * dt; // 4,4

    values[10] += lambda[0] * ( -sin( x[4] ) ) * dt; // 6,4
    values[10] += lambda[1] * cos( x[4] ) * dt; // 6,4

    values[11] += dt / Lf;
  }
  return true;
}

void TDNaive::finalize_solution(SolverReturn status,
                                  Index n, const Number* x, const Number* z_L, const Number* z_U,
                                  Index m, const Number* g, const Number* lambda,
                                  Number obj_value,
                                  const IpoptData* ip_data,
                                  IpoptCalculatedQuantities* ip_cq)
{
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution.

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


int main(int argv, char* argc[])
{
  // Create a new instance of your nlp
  //  (use a SmartPtr, not raw)
  SmartPtr<TNLP> mynlp = new TDNaive();

  // Create a new instance of IpoptApplication
  //  (use a SmartPtr, not raw)
  // We are using the factory, since this allows us to compile this
  // example with an Ipopt Windows DLL
  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

  // Change some options
  // Note: The following choices are only examples, they might not be
  //       suitable for your optimization problem.
  app->Options()->SetNumericValue("tol", 1e-7);
  app->Options()->SetStringValue("mu_strategy", "adaptive");
  app->Options()->SetStringValue("output_file", "ipopt.out");
  // The following overwrites the default name (ipopt.opt) of the
  // options file
  // app->Options()->SetStringValue("option_file_name", "hs071.opt");

  // Initialize the IpoptApplication and process the options
  ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Solve_Succeeded) {
    std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
    return (int) status;
  }

  // Ask Ipopt to solve the problem
  status = app->OptimizeTNLP(mynlp);

  if (status == Solve_Succeeded) {
    std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
  }
  else {
    std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
  }

  // As the SmartPtrs go out of scope, the reference count
  // will be decremented and the objects will automatically
  // be deleted.

  return (int) status;
}