#include <iostream>
#include <IpTNLP.hpp>
#include <memory>
#include "Ipopt_MPC_one_step.h"

int main(int argv, char* argc[])
{

    /** MPC single step problem formulation **/
    // initial state
    // std::vector<double> x_0{ 288.0, -178.0, 2.0, 19.0 };

    // Create a new instance of your nlp
    // (use a SmartPtr, not raw)
    // std::vector<std::vector<double>> wps { {279.02,-182.26,294.22,-176.02,287.39,-178.82,1.9603},
    //                                        {276.62,-172.5,291.29,-167.3,284.2,-169.82,1.9115},
    //                                        {272.91,-160.58,287.16,-155.43,279.94,-158.03,1.9174} };
    // Ipopt::SmartPtr<Ipopt::TNLP> kinematic_mpc = new IpoptMPCOneStep(x_0, wps);
    Ipopt::SmartPtr<Ipopt::TNLP> kinematic_mpc = new IpoptMPCOneStep();

    // Create a new instance of IpoptApplication
    //  (use a SmartPtr, not raw)
    // We are using the factory, since this allows us to compile this
    // example with an Ipopt Windows DLL
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();

    // Change some options
    // Note: The following choices are only examples, they might not be
    //       suitable for your optimization problem.
    app->Options()->SetNumericValue("tol", 1e-7); // default: 1e-7
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("output_file", "ipopt.out");
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");

    // The following overwrites the default name (ipopt.opt) of the
    // options file
    // app->Options()->SetStringValue("option_file_name", "hs071.opt");

    // Initialize the IpoptApplication and process the options
    Ipopt::ApplicationReturnStatus status;
    status = app->Initialize();

    // verbosely test
    std::cout << "status = " << status << std::endl;

    if (status != Ipopt::Solve_Succeeded) {
        std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
        return (int) status;
    }

    // verbosely test
    std::cout << "Seems like initialization is correctly done!" << std::endl;

    // Ask Ipopt to solve the problem
    status = app->OptimizeTNLP(kinematic_mpc);

    if (status == Ipopt::Solve_Succeeded) {
        std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
    }
    else {
        std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
    }

    // As the SmartPtrs go out of scope, the reference count
    // will be decremented and the objects will automatically
    // be deleted.

    return 0;
}