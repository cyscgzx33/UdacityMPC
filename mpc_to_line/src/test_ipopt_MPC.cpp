#include <iostream>
#include <IpTNLP.hpp>
#include <memory>
#include "Ipopt_MPC.h"

int main(int argv, char* argc[])
{

  /** MPC problem formulation **/
  // iteration times
  int iters = 50;
  // initial state
  std::vector<double> x_0{ 288.0, -178.0, 2.0, 39.0 };

  // state iterating
  for (size_t i = 0; i < iters; i++)
  {
      std::cout << "Iteration #" << i << ": " << std::endl;
    // Create a new instance of your nlp
    //  (use a SmartPtr, not raw)
    std::shared_ptr<IpoptMPC> kinematic_mpc_ipopt( new IpoptMPC(x_0) );
    Ipopt::SmartPtr<Ipopt::TNLP> kinematic_mpc = &(*kinematic_mpc_ipopt);
    // Ipopt::SmartPtr<Ipopt::TNLP> kinematic_mpc = new IpoptMPC(x_0);

    // Create a new instance of IpoptApplication
    //  (use a SmartPtr, not raw)
    // We are using the factory, since this allows us to compile this
    // example with an Ipopt Windows DLL
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();

    // Change some options
    // Note: The following choices are only examples, they might not be
    //       suitable for your optimization problem.
    app->Options()->SetNumericValue("tol", 1e-7);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("output_file", "ipopt.out");
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");

    // The following overwrites the default name (ipopt.opt) of the
    // options file
    // app->Options()->SetStringValue("option_file_name", "hs071.opt");

    // Initialize the IpoptApplication and process the options
    Ipopt::ApplicationReturnStatus status;
    status = app->Initialize();
    if (status != Ipopt::Solve_Succeeded) {
      std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
      return (int) status;
    }

    // Ask Ipopt to solve the problem
    status = app->OptimizeTNLP(kinematic_mpc);

    x_0 = kinematic_mpc_ipopt->getSolutionState();

    if (status == Ipopt::Solve_Succeeded) {
      std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
    }
    else {
      std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
    }

    // As the SmartPtrs go out of scope, the reference count
    // will be decremented and the objects will automatically
    // be deleted.
  }

  return 0;
}