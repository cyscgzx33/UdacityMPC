// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: ODNaive.cpp 2005 2011-06-06 12:55:16Z stefan $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#include "apply_ipopt_example.hpp"
#include "IpIpoptApplication.hpp"
#include <cassert>
#include <iostream>
#include <algorithm>

using namespace Ipopt;

/* global variable */
const double dt = 1.0;

// constructor
ODNaive::ODNaive()
{}

//destructor
ODNaive::~ODNaive()
{}

// returns the size of the problem
bool ODNaive::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                             Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // The problem described in ODNaive.hpp has 4 variables, x[0] through x[3]
  n = 5;

  // one equality constraint and one inequality constraint
  m = 2;

  // in this example the jacobian is dense and contains 8 nonzeros
  nnz_jac_g = 6;

  // the hessian is also dense and has 16 total nonzeros, but we
  // only need the lower left corner (since it is symmetric)
  nnz_h_lag = 5;

  // use the C style indexing (0-based)
  index_style = TNLP::C_STYLE;

  return true;
}

// returns the variable bounds
bool ODNaive::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{
  // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.
  assert(n == 5);
  assert(m == 2);

  // the variables have lower bounds of 1
  for (Index i=0; i<5; i++) {
    x_l[i] = -2.0;
  }

  // the variables have upper bounds of 5
  for (Index i=0; i<5; i++) {
    x_u[i] = 4.0;
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

  return true;
}

// returns the initial point for the problem
bool ODNaive::get_starting_point(Index n, bool init_x, Number* x,
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
  x[0] = 0.0;
  x[1] = 2.0;
  x[2] = 0.5;
  x[3] = 1.0;
  x[4] = 0.0;

  return true;
}

// returns the value of the objective function
bool ODNaive::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  assert(n == 5);

  std::vector<double> map = {1.0, 2.2, 3.0};
  std::vector<double> dist(3, 0.0);

  for (int i = 0; i < 3; i++)
    dist[i] = pow(map[i] - x[1], 2);
  auto closest = std::min_element( dist.begin(), dist.end() );
  int cl_idx = closest - dist.begin();

  obj_value = 0.0;

  obj_value += x[0] * x[0];
  // obj_value += (x[1] - 2.0) * (x[1] - 2.0);
  obj_value += (x[1] - map[cl_idx]) * (x[1] - map[cl_idx]);
  obj_value += (x[2] - 1.0) * (x[2] - 1.0);
  obj_value += (x[3] - 1.0) * (x[3] - 1.0);
  obj_value += x[4] * x[4];

  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool ODNaive::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  assert(n == 5);

  std::vector<double> map = {1.0, 2.2, 3.0};
  std::vector<double> dist(3, 0.0);

  for (int i = 0; i < 3; i++)
    dist[i] = pow(map[i] - x[1], 2);
  auto closest = std::min_element( dist.begin(), dist.end() );
  int cl_idx = closest - dist.begin();

  grad_f[0] = 2 * x[0];
  // grad_f[1] = 2 * (x[1] - 2.0);
  grad_f[1] = 2 * (x[1] - map[cl_idx]);
  grad_f[2] = 2 * (x[2] - 1.0);
  grad_f[3] = 2 * (x[3] - 1.0);
  grad_f[4] = 2 * x[4];

  return true;
}

// return the value of the constraints: g(x)
bool ODNaive::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  assert(n == 5);
  assert(m == 2);

  g[0] = x[0] + x[2] * dt - x[1];
  g[1] = x[2] + x[4] * dt - x[3];

  return true;
}

// return the structure or values of the jacobian
bool ODNaive::eval_jac_g(Index n, const Number* x, bool new_x,
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
    jCol[2] = 2;
    iRow[3] = 1;
    jCol[3] = 2;
    iRow[4] = 1;
    jCol[4] = 3;
    iRow[5] = 1;
    jCol[5] = 4;
  }
  else {
    // return the values of the jacobian of the constraints

    values[0] = 1.0; // 0,0
    values[1] = -1.0; // 0,1
    values[2] = dt; // 0,2

    values[3] = 1.0; // 0,3
    values[4] = -1.0; // 1,0
    values[5] = dt; // 1,1
  }

  return true;
}

//return the structure or values of the hessian
bool ODNaive::eval_h(Index n, const Number* x, bool new_x,
                       Number obj_factor, Index m, const Number* lambda,
                       bool new_lambda, Index nele_hess, Index* iRow,
                       Index* jCol, Number* values)
{

  if (values == NULL) {
    // return the structure. This is a symmetric matrix, fill the lower left
    // triangle only.

    // the hessian for this problem is actually dense
    Index idx=0;
    for (Index row = 0; row < 5; row++) {
        iRow[idx] = row;
        jCol[idx] = row;
        idx++;
    }

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
  }
  return true;
}

void ODNaive::finalize_solution(SolverReturn status,
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
  SmartPtr<TNLP> mynlp = new ODNaive();

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