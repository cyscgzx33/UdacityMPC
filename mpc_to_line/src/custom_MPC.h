#ifndef COSTOM_MPC_H
#define COSTOM_MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

class CustomMPC {
  public:
    CustomMPC();

    virtual ~CustomMPC();

    // Solve the model given an initial state.
    // Return the next state and actuations as a vector.
    std::vector<double> Solve(const Eigen::VectorXd &x0, 
                              const Eigen::VectorXd &coeffs);
    std::vector<double> Solve(const Eigen::VectorXd& x0);

  private:
    std::vector<std::vector<double>> waypoints_;
    std::vector<double> cl_x_; // center line x
    std::vector<double> cl_y_; // center line y
    std::vector<double> cl_phi_; // center line direction phi

    void readRoadmapFromCSV(const std::string& roadmap_file_name)
    {
        // Parse the roadmap file and add waypoints
        std::istringstream rm_stream(roadmap_file_name);
        std::string line;
        while (std::getline(rm_stream, line))
            parseRoadMapLine(line);
    }

    void parseRoadMapLine(const std::string& line)
    {
        // Parse each line in roadmap file and add waypoint vector
        std::istringstream line_stream(line);
        std::string number;
        std::vector<double> wp;
        while (std::getline(line_stream, number, ','))
            wp.push_back(std::stof(number));
        waypoints_.push_back(wp);
    }

};

#endif  // COSTOM_MPC_H