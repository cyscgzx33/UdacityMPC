// Author: Yusen Chen   2019-11-08

#include "apply_ipopt_n_steps_kinematics_model.hpp"
#include "IpIpoptApplication.hpp"
#include <cassert>
#include <iostream>
#include <algorithm>

using namespace Ipopt;

/* global variable */

// physical properties
const double v_ref = 20.0; // [20.0]
const double Lf = 2.67;

// iteration steps
const double dt    =  1.0;         // longer time: 0.02, 0.2, 2, 2.5, 4, 5
const int N        =  3;           // steps
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
  nnz_h_lag = 8 * N - 4; // (6 * N - 2) + ( 2 * (N - 1) ) 

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
  // Example: N = 2
  /*   
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
  */

  // Example: N = 3
  x[0] = 287.0;
  x[1] = 282.1;
  x[2] = 279.0;
  x[3] = -178.0;
  x[4] = -168.2;
  x[5] = -159.2;
  x[6] = 1.95;
  x[7] = 1.90;
  x[8] = 1.91;
  x[9] = 10.0;
  x[10] = 11.0;
  x[11] = 12.0;
  x[12] = 0.5;
  x[13] = 0.5;
  x[14] = 2.0;
  x[15] = 2.0;

  return true;
}

// Status: done
// returns the value of the objective function
bool NStepsKM::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  assert( n == (6 * N - 2) );

  // init objective value f
  obj_value = 0.0;

  // map info
  std::vector<double> cl_x   =  { 287.39,  284.2,   279.97,  276.07,  271.63,  265.64,  257.88,  251.96,  245.37, 238.04,  236.1 };
  std::vector<double> cl_y   =  { -178.82, -169.82, -158.03, -146.81, -133.26, -117.35, -93.156, -76.336, -56.4,  -35.048, -26.611 };
  std::vector<double> cl_phi =  { 1.9603,  1.9115,  1.9174,  1.9027,  1.8878,  1.9306,  1.8812,  1.9093,  1.89,   1.9015,  1.7966 };

  // for finding closest point
  int sz = cl_x.size();
  std::vector<double> dist(sz, 0.0);

  for (int i = 0; i < N; i++) // totally N steps
  {
      for (int j = 0; j < sz; j++) // the waypoints size equals to sz
      {
          // the distances from point ( pos_X[i], pos_Y[i] ) to point sequence ( cl_x[j], cl_y[j] )
          dist[j] = pow(cl_x[j] - x[x_st + i], 2) + pow(cl_y[j] - x[y_st + i], 2);
      }

      // find the index belongs to the smallest element
      auto closest_pt = std::min_element( dist.begin(), dist.end() );
      int closest_cl_idx = closest_pt - dist.begin();

      /* assign the objective value f */
      // states penalization
      // (1) pos_X and pos_Y
      obj_value += pow(x[x_st + i] - cl_x[closest_cl_idx], 2) + pow(x[y_st + i] - cl_y[closest_cl_idx], 2);
      // (2) heading angle
      obj_value += pow(x[phi_st + i] - cl_phi[closest_cl_idx], 2);
      // (3) reference velocity
      obj_value += pow(x[v_st + i] - v_ref, 2);
      
      // actuation penalization (ignore this step when i == N - 1)
      if (i < N - 1)
      { 
        // (1) steering angle delta
        obj_value += x[delta_st + i] * x[delta_st + i];
        // (2) acceleration
        obj_value += x[a_st + i] * x[a_st + i];
      }
  }
  
  return true;
}

// Status: done
// return the gradient of the objective function grad_{x} f(x)
bool NStepsKM::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  assert( n == (6 * N - 2) );

  // map info
  std::vector<double> cl_x   =  { 287.39,  284.2,   279.97,  276.07,  271.63,  265.64,  257.88,  251.96,  245.37, 238.04,  236.1 };
  std::vector<double> cl_y   =  { -178.82, -169.82, -158.03, -146.81, -133.26, -117.35, -93.156, -76.336, -56.4,  -35.048, -26.611 };
  std::vector<double> cl_phi =  { 1.9603,  1.9115,  1.9174,  1.9027,  1.8878,  1.9306,  1.8812,  1.9093,  1.89,   1.9015,  1.7966 };

  // for finding closest point
  int sz = cl_x.size();
  std::vector<double> dist(sz, 0.0);

  for (int i = 0; i < N; i++) // totally N steps
  {
      for (int j = 0; j < sz; j++) // the waypoints size equals to sz
      {
          // the distances from point ( pos_X[i], pos_Y[i] ) to point sequence ( cl_x[j], cl_y[j] )
          dist[j] = pow(cl_x[j] - x[x_st + i], 2) + pow(cl_y[j] - x[y_st + i], 2);
      }

      // find the index belongs to the smallest element
      auto closest_pt = std::min_element( dist.begin(), dist.end() );
      int closest_cl_idx = closest_pt - dist.begin();

      /* assign the objective value f */
      // states penalization
      // (1) pos_X and pos_Y
      grad_f[x_st + i]   =  2 * ( x[x_st + i] - cl_x[closest_cl_idx] );
      grad_f[y_st + i]   =  2 * ( x[y_st + i] - cl_y[closest_cl_idx] );
      // (2) heading angle
      grad_f[phi_st + i] =  2 * ( x[phi_st + i] - cl_phi[closest_cl_idx] );
      // (3) reference velocity
      grad_f[v_st + i]   =  2 * ( x[v_st + i] - v_ref );

      // actuation penalization (ignore this step when i == N - 1)
      if (i < N - 1)
      {
        // (1) steering angle delta
        grad_f[delta_st + i] =  2 * x[delta_st + i];
        // (2) acceleration
        grad_f[a_st + i]     =  2 * x[a_st + i];
      }
  }
  
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