#define _USE_MATH_DEFINES

#include <vector>
#include <fstream>
#include <ios>
#include <limits>
#include <string>
#include <tuple>
#include <utility>

#include "./headers/inputoutput.h"
#include "./headers/utils.h"
#include "./headers/numerical.h"


int main () {
    // Initialize the variables to be used
    vector<double> dimensions;
    vector<int> nodes;
    bool isTestCase;
    vector<double> BC_temperatures;
    int coordinate_sys = 1; // fixed for now to the cartesian coordinate system

    // Output intro messages to the console
    initiateProgram();
    
    // Start getting input parameters 
    std::tie(dimensions, nodes, isTestCase, BC_temperatures) = requestInputs();
    
    // Create mesh accordingly
    vector<vector<double>> mesh = generate_mesh(dimensions[0], dimensions[1], nodes[0], nodes[1], coordinate_sys);
    
    // Use the 4-point finite difference method to provide numerical approximation to the problem
    vector<double> t_sim = finite_diff_4_pts(dimensions[0], dimensions[1], nodes[0], nodes[1], coordinate_sys, BC_temperatures);
    
    // Implement unit test case (if requested) and output the results
    if(isTestCase) { 
        // Calculate the analytical solution
        vector<double> t_analytical = analytical_sol(mesh, dimensions[0], dimensions[1], BC_temperatures[0], BC_temperatures[3]);

        // Calculate the error between analytical and numerical values
        vector<vector<double>> error = test_fdm_error(t_analytical, t_sim);

        // Save the total results (with analytical and error vectors)
        outputResults(mesh, t_sim, t_analytical, error[1]);
    }
    else {
        // Save the total results
        outputResults(mesh, t_sim);
    }

    // Plot a graph of the results for the simulated temperature values
    plotResults(dimensions, nodes, t_sim);
}