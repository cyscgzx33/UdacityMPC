// Copyright (C) 2005, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: hs071_nlp.hpp 1861 2010-12-21 21:34:47Z stefan $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-09

#ifndef __TD_NAIVE_HPP__
#define __TD_NAIVE_HPP__

#include "IpTNLP.hpp"

using namespace Ipopt;
// std::vector<double> cl_x   =  { 287.39, 284.2, 279.97 };
// std::vector<double> cl_y   =  { -178.82, -169.82, -158.03 };
// std::vector<double> cl_phi =  { 1.9603, 1.9115, 1.9174 };

/**
 *
 * Problem ODNaive looks like this
 *     x[0]: x(0)
 *     x[1]: x(1)
 *     x[2]: y(0)
 *     x[3]: y(1)
 *     x[4]: psi(0)
 *     x[5]: psi(1)
 *     x[6]: v(0)
 *     x[7]: v(1)
 *     x[8]: delta(0)
 *     x[9]: a(0)    
 * 
 * 
 *     min   (x0 - x0_cl)^2 + (x1 - x1_cl)^2 + (x2 - x2_cl)^2 + (x3 - x3_cl)^2 +
 *           (x4 - x4_cl)^2 + (x5 - x5_cl)^2 + (x6 - v_ref)^2 + (x7 - v_ref)^2 +
 *           x8^2 + x9^2
 *     s.t.  x0 + x6 * cos(x4) * dt - x1       =  0
 *           x2 + x6 * sin(x4) * dt - x3       =  0
 *           x4 + x6 / Lf * x8 * dt - x5  =  0
 *           x6 + x9 * dt - x7            =  0
 *
 *     Starting point:
 *        x = (287.0, 282.1, -178.0, -168.2, 1.95, 1.90, 18, 18.4, 0.5, 2.0)
 *
 *     Optimal solution:
 *        x = ?
 *
 *
 */

class TDNaive : public TNLP
{
public:
  /** default constructor */
  TDNaive();

  /** default destructor */
  virtual ~TDNaive();

  /**@name Overloaded from TNLP */
  //@{
  /** Method to return some info about the nlp */
  virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style);

  /** Method to return the bounds for my problem */
  virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                               Index m, Number* g_l, Number* g_u);

  /** Method to return the starting point for the algorithm */
  virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                  bool init_z, Number* z_L, Number* z_U,
                                  Index m, bool init_lambda,
                                  Number* lambda);

  /** Method to return the objective value */
  virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

  /** Method to return the gradient of the objective */
  virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

  /** Method to return the constraint residuals */
  virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

  /** Method to return:
   *   1) The structure of the jacobian (if "values" is NULL)
   *   2) The values of the jacobian (if "values" is not NULL)
   */
  virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                          Index m, Index nele_jac, Index* iRow, Index *jCol,
                          Number* values);

  /** Method to return:
   *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
   *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
   */
  virtual bool eval_h(Index n, const Number* x, bool new_x,
                      Number obj_factor, Index m, const Number* lambda,
                      bool new_lambda, Index nele_hess, Index* iRow,
                      Index* jCol, Number* values);

  //@}

  /** @name Solution Methods */
  //@{
  /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
  virtual void finalize_solution(SolverReturn status,
                                 Index n, const Number* x, const Number* z_L, const Number* z_U,
                                 Index m, const Number* g, const Number* lambda,
                                 Number obj_value,
                                 const IpoptData* ip_data,
                                 IpoptCalculatedQuantities* ip_cq);
  //@}

private:
  /**@name Methods to block default compiler methods.
   * The compiler automatically generates the following three methods.
   *  Since the default compiler implementation is generally not what
   *  you want (for all but the most simple classes), we usually
   *  put the declarations of these methods in the private section
   *  and never implement them. This prevents the compiler from
   *  implementing an incorrect "default" behavior without us
   *  knowing. (See Scott Meyers book, "Effective C++")
   * 
   */
  //@{
  //  TDNaive();
  TDNaive(const TDNaive&);
  TDNaive& operator=(const TDNaive&);
  //@}
};


#endif /* __TD_NAIVE_HPP__ */