// Author: Yusen Chen   2019-11-08

#include "apply_ipopt_n_steps_kinematics_model.hpp"
#include "IpIpoptApplication.hpp"
#include <cassert>
#include <iostream>
#include <algorithm>

using namespace Ipopt;

/* global variable */

// physical properties
const double v_ref = 20.0;
const double Lf = 2.67;

// iteration steps
const double dt    =  5.0;         // longer time: 0.02, 0.2, 2, 2.5, 4, 
const int N        =  2;           // steps
const int x_st     =  0;           // x var start idx
const int y_st     =  N;           // y var start idx
const int phi_st   =  2 * N;       // phi var start idx
const int v_st     =  3 * N;       // v var start idx
const int delta_st =  4 * N;       // delta var start idx
const int a_st     =  5 * N - 1;   // a var start idx


/* Definition of NStepsKM class methods & attributes */

// constructor
NStepsKM::NStepsKM()
{}

//destructor
NStepsKM::~NStepsKM()
{}

// Status: done
// returns the size of the problem
bool NStepsKM::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                             Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // The problem described in NStepsKM.hpp has 4 variables, x[0] through x[3]
  n = 6 * N - 2; // 4 * N + 2 * (N - 1) = 6 * N - 2

  // equality constraints 
  // (it can be equality or inequality ones, here we only have eqaulity ones)
  m = 4 * N - 4;

  // TODO(done): modify it
  // in this example the jacobian is dense and contains 8 nonzeros
  nnz_jac_g = 15 * N - 15; // 15 * (N - 1)

  // TODO: modify it
  // the hessian is also dense and has 16 total nonzeros, but we
  // only need the lower left corner (since it is symmetric)
  nnz_h_lag = 7 * N - 4; // (6 * N - 2) + ( 2 * (N - 1) ) 

  // use the C style indexing (0-based)
  index_style = TNLP::C_STYLE;

  return true;
}

// Status: done
// returns the variable bounds
bool NStepsKM::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{
  // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.
  assert( n == (6 * N - 2) );
  assert( m == (4 * N - 4) );

  // the variables have lower bounds & upper bounds
  for (Index i = 0; i < n; i++) 
  {
    x_l[i] = -500;
    x_u[i] = 500;
  }

  // the constraints are equality ones, i.e., all of them has val equals to 0
  for (Index i = 0; i < m; i++)
  {
      g_l[i] = 0.0;
      g_u[i] = 0.0;
  }

  return true;
}

// Status: pending
// returns the initial point for the problem
bool NStepsKM::get_starting_point(Index n, bool init_x, Number* x,
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

// Status: pending
// returns the value of the objective function
bool NStepsKM::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  assert( n == (6 * N - 2) );

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
  // pay more attention on the second position tracking performance
  obj_value += 10 * pow(x[1] - cl_x[cl_idx_x1], 2) + 10 * pow(x[3] - cl_y[cl_idx_x1], 2);
  obj_value += 10 * pow(x[5] - cl_phi[cl_idx_x1], 2);
  obj_value += pow(x[7] - v_ref, 2);
  // actuators
  obj_value += x[8] * x[8];
  obj_value += x[9] * x[9];

  return true;
}

// Status: pending
// return the gradient of the objective function grad_{x} f(x)
bool NStepsKM::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
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
  grad_f[1] = 2 * ( x[1] - cl_x[cl_idx_x1] ) * 10;
  grad_f[2] = 2 * ( x[2] - cl_y[cl_idx_x0] );
  grad_f[3] = 2 * ( x[3] - cl_y[cl_idx_x1] ) * 10;
  grad_f[4] = 2 * ( x[4] - cl_phi[cl_idx_x0] );
  grad_f[5] = 2 * ( x[5] - cl_phi[cl_idx_x1] ) * 10;
  grad_f[6] = 2 * ( x[6] - v_ref );
  grad_f[7] = 2 * ( x[7] - v_ref );
  grad_f[8] = 2 * x[8];
  grad_f[9] = 2 * x[9];

  return true;
}

// Status: done
// return the value of the constraints: g(x)
bool NStepsKM::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  assert( n == (6 * N - 2) );
  assert( m == (4 * N - 4) );

  for (Index i = 0; i < N - 1; i++)
  {
      g[4 * i]      =  x[x_st + i] + x[v_st + i] * cos( x[phi_st + i] ) * dt - x[x_st + i + 1];
      g[4 * i + 1]  =  x[y_st + i] + x[v_st + i] * sin( x[phi_st + i] ) * dt - x[y_st + i + 1];
      g[4 * i + 2]  =  x[phi_st + i] + x[v_st + i] * x[delta_st + i] * dt / Lf - x[phi_st + i + 1];
      g[4 * i + 3]  =  x[v_st + i] + x[a_st + i] * dt - x[v_st + i + 1];
  }

  return true;
}

// Status: done
// return the structure or values of the jacobian
bool NStepsKM::eval_jac_g(Index n, const Number* x, bool new_x,
                         Index m, Index nele_jac, Index* iRow, Index *jCol,
                         Number* values)
{
  if (values == NULL) {
    // return the structure of the jacobian
    // this particular jacobian is dense
    for (Index i = 0; i < N - 1; i++)
    {
        iRow[i * 15] = 4 * i;
        jCol[i * 15] = x_st + i;
        iRow[i * 15 + 1] = 4 * i;
        jCol[i * 15 + 1] = x_st + i + 1;
        iRow[i * 15 + 2] = 4 * i;
        jCol[i * 15 + 2] = phi_st + i;
        iRow[i * 15 + 3] = 4 * i;
        jCol[i * 15 + 3] = v_st + i;
        iRow[i * 15 + 4] = 4 * i + 1;
        jCol[i * 15 + 4] = y_st + i;
        iRow[i * 15 + 5] = 4 * i + 1;
        jCol[i * 15 + 5] = y_st + i + 1;
        iRow[i * 15 + 6] = 4 * i + 1;
        jCol[i * 15 + 6] = phi_st + i;
        iRow[i * 15 + 7] = 4 * i + 1;
        jCol[i * 15 + 7] = v_st + i;
        iRow[i * 15 + 8] = 4 * i + 2;
        jCol[i * 15 + 8] = phi_st + i;
        iRow[i * 15 + 9] = 4 * i + 2;
        jCol[i * 15 + 9] = phi_st + i + 1;
        iRow[i * 15 + 10] = 4 * i + 2;
        jCol[i * 15 + 10] = v_st + i;
        iRow[i * 15 + 11] = 4 * i + 2;
        jCol[i * 15 + 11] = delta_st + i;
        iRow[i * 15 + 12] = 4 * i + 3;
        jCol[i * 15 + 12] = v_st + i;
        iRow[i * 15 + 13] = 4 * i + 3;
        jCol[i * 15 + 13] = v_st + i + 1;
        iRow[i * 15 + 14] = 4 * i + 3;
        jCol[i * 15 + 14] = a_st + i;
    }
  }
  else {
    // return the values of the jacobian of the constraints
    for (Index i = 0; i < N - 1; i++)
    {
        values[i * 15] = 1.0;
        values[i * 15 + 1] = -1.0;
        values[i * 15 + 2] = -sin( x[phi_st + i] ) * x[v_st + i] * dt;
        values[i * 15 + 3] = cos( x[phi_st + i] ) * dt;

        values[i * 15 + 4] = 1.0;
        values[i * 15 + 5] = -1.0;
        values[i * 15 + 6] = cos( x[phi_st + i] ) * x[v_st + i] * dt;
        values[i * 15 + 7] = sin( x[phi_st + i] ) * dt;

        values[i * 15 + 8] = 1.0;
        values[i * 15 + 9] = -1.0;
        values[i * 15 + 10] = x[delta_st + i] * dt / Lf;
        values[i * 15 + 11] = x[v_st + i] * dt / Lf;

        values[i * 15 + 12] = 1.0;
        values[i * 15 + 13] = -1.0;
        values[i * 15 + 14] = dt;
    }
  }

  return true;
}

// Status: done
//return the structure or values of the hessian
bool NStepsKM::eval_h(Index n, const Number* x, bool new_x,
                       Number obj_factor, Index m, const Number* lambda,
                       bool new_lambda, Index nele_hess, Index* iRow,
                       Index* jCol, Number* values)
{
  if (values == NULL) {
    // return the structure. This is a symmetric matrix, fill the lower left
    // triangle only.

    // the hessian for this problem is actually dense
    Index idx = 0;

    // for the grad_grad_f part:
    // formulate an diagonal matrix first
    for (Index row = 0; row < 6 * N - 2; row++)
    {
      iRow[idx] = row;
      jCol[idx] = row;
      idx++;
    }

    // for the grad_grad_g part:
    // formulate each non-zero element
    for (Index i = 0; i <  N - 1; i++)
    {
      // for grad_grad_g_(i) & grad_grad_g_(i+1)
      iRow[idx] = v_st + i;
      jCol[idx] = phi_st + i;
      idx++;

      // for grad_grad_g(i+2)
      iRow[idx] = delta_st + i;
      jCol[idx] = v_st + i;
      idx++;
    }

    assert(idx == nele_hess);
  }
  else {
    // return the values. This is a symmetric matrix, fill the lower left
    // triangle only

    Index idx = 0;
    // for the grad_grad_f part:
    // assign values affected by grad_grad_f
    // Note (important): adjust the params after changing the weights of cost function 
    for (Index row = 0; row < 6 * N - 2; row++)
    {
      values[idx] += obj_factor * 2;
      idx++;
    }

    // for the grad_grad_g part:
    // assign values affected by grad_grad_g
    for (Index i = 0; i <  N - 1; i++)
    {
      // for diagnal parts
      // Note: the "idx" of each values happens to be equal to phi_st + i
      values[phi_st + i] += lambda[4 * i] * x[v_st + i] * cos( x[phi_st + i] ) * dt;
      values[phi_st + i] += lambda[4 * i  + 1] * x[v_st + i] * ( -sin( x[phi_st + i] ) ) * dt;

      // for non-diagnal parts
      // Note: here we can only rely on "idx" to count where to fill in values
      // grad_grad_g_(i) & grad_grad_g_(i+1)
      values[idx] += lambda[4 * i] * ( -sin( x[phi_st + i] ) ) * dt;
      values[idx] += lambda[4 * i + 1] * cos( x[phi_st + i] ) * dt;
      idx++;

      // for grad_grad_g(i+2)
      values[idx] += lambda[4 * i + 2] * dt / Lf;
      idx++;
    }

    assert(idx == nele_hess);
  }

  return true;
}

void NStepsKM::finalize_solution(SolverReturn status,
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
  SmartPtr<TNLP> mynlp = new NStepsKM();

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